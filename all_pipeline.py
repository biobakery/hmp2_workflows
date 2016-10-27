import os
import re
import json
import glob
import gzip
import operator
from bz2 import BZ2File
from tarfile import TarFile
from datetime import datetime
from itertools import imap, izip
from collections import Counter
from StringIO import StringIO
from os.path import dirname, basename

from bunch import Bunch
from toolz import concat, groupby, get, concatv
import biom
import biom.table
from biom.table import DenseOTUTable
from biom.parse import (
    OBS_META_TYPES,
    parse_biom_table,
    parse_classic_table_to_rich_table
    )

from anadama.sge import SGEPowerup
from anadama.runcontext import RunContext
from anadama.deps import FileDependency
from anadama.deps import ExecutableDependency
from anadama.deps import GlobDependency
from anadama.deps import HugeFileDependency as HFD
from anadama.deps import DirectoryDependency
from anadama.helpers import system
from anadama.helpers import rm
from anadama.helpers import rm_r
from anadama.reporters import LoggerReporter
from anadama.util import fname

snd = operator.itemgetter(1)

picard_dir = "/seq/picard_aggregation"
ignores_file = "/seq/ibdmdb/centos6/hmp2_v2/samples_ignore"

metadata_json = FileDependency("/seq/ibdmdb/data_deposition/HMP2"
                               "/Metadata/json/joined_subcoll.json")
output_dirs = Bunch(
    metagenomics = Bunch(
        cleanseq = "/seq/ibdmdb/processing/HMP2/WGS/1508/kbayer/anadama_products",
        taxprof = "/seq/ibdmdb/public/HMP2/WGS/1508/kbayer/anadama_products",
        h2 = "/seq/ibdmdb/public/HMP2/WGS/1508/kbayer/anadama_products",
        ),
    metatranscriptomics = Bunch(
        cleanseq = "/seq/ibdmdb/processing/HMP2_MT/WGS/1524/rschwager/anadama_products",
        taxprof = "/seq/ibdmdb/public/HMP2_MT/WGS/1524/rschwager/anadama_products",
        h2 = "/seq/ibdmdb/public/HMP2_MT/WGS/1524/rschwager/anadama_products",
        ),
    amplicon = Bunch(
        cleanseq = "/seq/ibdmdb/processing/HMP2/16S/1547/rschwager/anadama_products",
        taxprof = "/seq/ibdmdb/public/HMP2/16S/1547/rschwager/anadama_products",
        h2 = "/seq/ibdmdb/public/HMP2/16S/1547/rschwager/anadama_products",
        )
    )

conf = Bunch(otudb="/seq/ibdmdb/centos6/var/lib/gg_13_5_otus/rep_set/97_otus.udb",
             chimeradb="/seq/ibdmdb/centos6/var/lib/uchime/gold.fa",
             taxonomy="/seq/ibdmdb/centos6/var/lib/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt",
             copyno="/seq/ibdmdb/centos6/picrust-1.0.0/picrust/data/16S_13_5_precalculated.tab.gz")

bin = Bunch(fqjoin="fastq-join",
            usearch="usearch8",
            picrust1="normalize_by_copy_number.py",
            picrust2="predict_metagenomes.py",
            samtools="samtools",
            trimmomatic="/seq/ibdmdb/centos6/var/lib/Trimmomatic-0.33/trimmomatic-0.33.jar",
            knead="kneaddata",
            metaphlan2="metaphlan2.py",
            humann2="humann2")

for k in bin:
    bin[k] = ExecutableDependency(bin[k])



kneaddb = Bunch(dna="/seq/ibdmdb/centos6/var/lib/Homo_sapiens_Bowtie2_v0.1",
                rna="/seq/ibdmdb/centos6/var/lib/mrna")
kneaddb_deps = Bunch([ (k, GlobDependency(v+"/*")) for k, v in kneaddb.items() ])

h2ntdb = DirectoryDependency("/seq/ibdmdb/centos6/var/lib/humann2/chocophlan")
h2pdb  = DirectoryDependency("/seq/ibdmdb/centos6/var/lib/humann2/uniref")


def findfile(row):
    if not row[0].startswith("G"):
        return
    projdir = os.path.join(picard_dir, row[0])
    sampdirs = [ os.path.join(projdir, x) for x in os.listdir(projdir) ]
    sampdirs = filter(os.path.isdir, sampdirs)
    if len(sampdirs) > 1:
        d = dict( zip(map(os.path.basename, sampdirs), sampdirs) )
        sampdir = d.get(row[1])
        if not sampdir:
            import pdb; pdb.set_trace()
            raise Exception("more than one sample per project")
    else:
        sampdir = sampdirs[0]
    bams = glob.glob(os.path.join(sampdir, "current", "*.bam"))
    if len(bams) > 1:
        import pdb; pdb.set_trace()
        raise Exception("more than one bamfile per sample")
    return bams[0]


def apply_ignores(rows):
    if not ignores_file or not os.path.exists(ignores_file):
        return rows
    with open(ignores_file, 'r') as f:
        ignores = set(iter(f))
    return [ row for row in rows
             if row[1] not in ignores ]


def fasta_sequences(seqs_f):
    id = None
    lines = iter(seqs_f)
    line = next(lines).strip()
    while True:
        if line.startswith(">"):
            id = line.split(None, 1)[0]
            id = id.replace(">", "")
            line = next(lines).strip()
        else:
            seq = str()
            while not line.startswith(">"):
                seq += line
                try:
                    line = next(lines).strip()
                except StopIteration as e:
                    yield (id, seq)
                    raise e
            yield (id, seq)
            id, seq = str(), str()


def fields(fname):
    with open(fname, 'r') as f:
        for line in f:
            yield line.strip().split("\t")


def cumsum(it):
    sm = 0
    for item in it:
        sm += item
        yield sm


def find_cutoff(fname, ratio=0.95):
    with open(fname, 'r') as f:
        reads = izip(*[iter(f)]*4)
        cntr = Counter( imap(len, imap(snd, reads)) )
    hist = sorted(cntr.iteritems(), reverse=True)
    read_length_counts = map(snd, hist)
    read_depth = sum(read_length_counts)
    for i, sm in enumerate(cumsum(read_length_counts)):
        if float(sm)/read_depth > ratio:
            return hist[i][0]


def split(fq, r1, r2):
    def _actually_split(task):
        in_f, outr1, outr2 = open(fq, 'r'), open(r1, 'w'), open(r2, 'w')
        with in_f, outr1, outr2:
            reads = izip(*[iter(in_f)]*4)
            for read in reads:
                if read[0].endswith('/1\n'):
                    outr1.write("".join(read))
                else:
                    outr2.write("".join(read))

    return _actually_split


def truncate(joinfq, truncated, logfile):
    def _actually_truncate(task):
        cutoff = find_cutoff(joinfq)
        system(["usearch8", "-fastx_truncate", joinfq, "-trunclen", 
                cutoff, "-fastaout", truncated], stdout=logfile)(None)
    return _actually_truncate

get_otu_id = lambda s: re.search(r'OTU_(\d+)', s).group(1)

def count_otus(taxonomy_txt, nochimera_fa, refout_uc, mapout_uc):
    taxmap = dict(fields(taxonomy_txt))
    with open(nochimera_fa, 'r') as f:
        seqmap = dict(fasta_sequences(f))
    otu_tax = {}
    for row in fields(refout_uc):
        if row[0] != "H":
            continue
        query, target = row[8], row[9]
        otu_id = get_otu_id(query)
        otu_tax[otu_id] = [taxmap[target], seqmap[query], 0, target]
    for row in fields(mapout_uc):
        if row[0] != "H":
            continue
        otu_id = get_otu_id(row[9])
        if otu_id in otu_tax:
            otu_tax[otu_id][2] += 1
    ret = {}
    for otu_id, (tax, seq, cnt, tax_id) in otu_tax.iteritems():
        if tax_id not in ret:
            ret[tax_id] = [tax, [otu_id], [seq], cnt]
        else:
            ret[tax_id][1].append(otu_id)
            ret[tax_id][2].append(seq)
            ret[tax_id][3] += cnt
    return ret


def write_otu_txt(tax_seqs, otu_txt):
    with open(otu_txt, 'w') as f:
        print >> f, "\t".join(("Taxonomy", "Count"))
        for tax_id, (t, _, _, count) in tax_seqs.iteritems():
            print >> f, "\t".join((t+":"+str(tax_id), str(count)))


def write_otu_biom(tax_seqs, otu_biom):
    sample_ids = [fname.rmext(otu_biom)]
    obs_metadata = list()
    data = list()
    tax_ids = list()
    for tax_id, (t, _, seqs, count) in tax_seqs.iteritems():
        data.append([count])
        obs_metadata.append({"taxonomy": t.split("; "), "otusequences": seqs})
        tax_ids.append(tax_id)
    biom_table = biom.table.table_factory(
        data, sample_ids, tax_ids, 
        sample_metadata=None, observation_metadata=obs_metadata,
        table_id="ibdmdb_org_otu_table", constructor=biom.table.DenseOTUTable
        )
    with open(otu_biom, 'w') as f:
        json.dump(biom_table.getBiomFormatObject("ibdmdb_org_otu_table"), f)


def tabulate_otu_table(taxonomy_txt, nochimera_fa, mapout_uc, 
                       refout_uc, otu_txt, otu_biom):
    def _acutally_tabulate_otu_table(task):
        otu_seqs = count_otus(taxonomy_txt, nochimera_fa, refout_uc, mapout_uc)
        write_otu_txt(otu_seqs, otu_txt)
        write_otu_biom(otu_seqs, otu_biom)
    return _acutally_tabulate_otu_table


def convert(ctx, rawbam, rawfq):
    ctx.add_task(
        [ system(["samtools", "bam2fq", str(rawbam),], stdout_clobber=rawfq) ],
        depends=[rawbam, bin.samtools],
        targets=HFD(rawfq),
        name="Convert "+str(rawbam)
        )
 

def otupick(ctx, rawfq, otu_biom, otu_txt, log, truncqual=10, pct_id=0.97):
    to_del = {}
    r1     = to_del[0]  = fname.mangle(rawfq, tag="r1", ext='fq')
    r2     = to_del[1]  = fname.mangle(rawfq, tag="r1", ext='fq')
    joinfq = to_del[3]  = fname.mangle(rawfq, tag="join", ext='fq')
    trimfq = to_del[4]  = fname.mangle(rawfq, tag="trim", ext='fq')
    truncd = to_del[5]  = fname.mangle(rawfq, tag="truncated", ext="fa")
    derep  = to_del[6]  = fname.mangle(rawfq, tag="derep", ext="fa")
    sortd  = to_del[7]  = fname.mangle(rawfq, tag="sort", ext="fa")
    rawotu = to_del[8]  = fname.mangle(rawfq, tag="raw_otus", ext="fa")
    nochim = to_del[9]  = fname.mangle(rawfq, tag="nonchimeric", ext="fa")
    mapout = to_del[10]  = fname.mangle(rawfq, tag="map", ext="uc")
    refout = to_del[11]  = fname.mangle(rawfq, tag="refout", ext="uc")
    
    ctx.add_task(
        [ split(rawfq, r1, r2),
          system(["fastq-join", r1, r2, "-o", "/dev/null", 
                                        "-o", "/dev/null", 
                                        "-o", joinfq], 
                 stdout=log),
          system(["usearch8", "-fastq_filter", joinfq, "-fastqout", trimfq, 
                  "-fastq_truncqual", truncqual], stderr=log),
          truncate(trimfq, truncd, log),
          system(["usearch8", "-derep_fulllength", truncd, "-fastaout", derep,
                  "-sizeout"], stderr=log),
          system(["usearch8", "-sortbysize", derep, "-fastaout", sortd,
                  "-minsize", 2], stderr=log),
          system(["usearch8", "-cluster_otus", sortd, "-otus", rawotu, 
                  "-uparseout", "/dev/null", "-relabel", "OTU_", "-sizein", 
                  "-sizeout"], stderr=log),
          system(["usearch8", "-uchime_ref", rawotu, "-nonchimeras", nochim, 
                  "-db", conf.chimeradb, "-strand", "plus"], stderr=log),
          system(["usearch8", "-usearch_global", joinfq, "-strand", "plus",
                  "-db", nochim, "-uc", mapout, "-id", pct_id], stderr=log),
          system(["usearch8", "-usearch_global", nochim, "-uc", refout, 
                  "-db", conf.otudb, "-top_hit_only", 
                  "-strand", "both", "-id", pct_id], stderr=log),
          tabulate_otu_table(conf.taxonomy, nochim, mapout, refout, 
                             otu_txt, otu_biom),
          rm(to_del.values())
          ],
        depends=[bin.fqjoin, bin.usearch, HFD(rawfq), 
                 HFD(conf.chimeradb), HFD(conf.otudb)],
        targets=[log, otu_biom, otu_txt],
        name="OTU Pick "+str(rawfq)
        )


def _drop_unknown_otu_ids(otu_biom, filtered_biom, copy_number_fname):
    def _actually_drop(task):
        idx = set([ row.strip().split('\t')[0]
                    for row in gzip.open(copy_number_fname) ])
        filter_func = lambda a, otu_id, c: str(otu_id) in idx
        with open(otu_biom) as f, open(filtered_biom, 'w') as f_out:
            try:
                table = parse_biom_table(f)
            except Exception:
                table = parse_classic_table_to_rich_table(
                    f, None, None, OBS_META_TYPES['taxonomy'], DenseOTUTable)
            table = table.filterObservations(filter_func)
            json.dump( table.getBiomFormatObject("AnADAMA"), f_out )
    return _actually_drop


def picrust(ctx, otu_biom, tarbz):
    tmpdir = fname.rmext(str(tarbz), all=True)
    fltrd = fname.mangle(otu_biom, dir=tmpdir, tag="filtered", ext="biom")
    normed = fname.mangle(otu_biom, dir=tmpdir, tag="normalized", ext="biom")
    pred = fname.mangle(otu_biom, dir=tmpdir, tag="picrust", ext="biom")

    ctx.add_task(
        [ _drop_unknown_otu_ids(otu_biom, fltrd, conf.copyno),
          system(["normalize_by_copy_number.py", "-i", fltrd, "-o", normed]),
          system(["predict_metagenomes.py", "-i", normed, "-o", pred]),
          system(["tar", "-cjf", tarbz, tmpdir]),
          rm_r(tmpdir) ],
        depends=[otu_biom, bin.picrust1, bin.picrust2, HFD(conf.copyno)],
        targets=tarbz,
        name="Predict function"+str(otu_biom)
        )


def knead(ctx, bamfile, output_clean, output_stats, threads=2, 
          dbs=[kneaddb.dna], db_deps=[kneaddb_deps.dna]):
    dbs = list(concat([ ['-db', d] for d in dbs ]))
    tmpdir = fname.rmext(output_clean, all=True)
    rawfq = fname.mangle(str(bamfile), dir=tmpdir, ext="fastq")
    kneadout = os.path.join(tmpdir, basename(output_clean))
    kneadlog = fname.mangle(kneadout, ext="log")
    ctx.grid_add_task(
        [ system(["samtools", "bam2fq", str(bamfile)], stdout_clobber=rawfq),
          system(["kneaddata", '-i', rawfq, '-o', tmpdir]+dbs+[
                  "--trimmomatic", dirname(str(bin.trimmomatic)), 
                  "--threads", threads,
                  "--output-prefix", fname.rmext(basename(output_clean))]),
          system(["mv", '-f', kneadout, output_clean ]),
          system(["mv", '-f', kneadlog, output_stats ]),
          rm_r(tmpdir) ],
        depends=[bamfile, bin.samtools, bin.knead]+db_deps,
        targets=[HFD(output_clean), output_stats],
        name="Decontaminate "+str(bamfile),
        mem=4000,
        time=90,
        cores=threads
        )


def metaphlan2(ctx, fastq, output_biom, output_txt, output_strainbam):
    sam = fname.mangle(output_strainbam, ext="sam")
    ctx.grid_add_task(
        [ system(["metaphlan2.py", "--no_map", "--input_type", "fastq",
                  "-s", sam, "--biom", output_biom, fastq, output_txt]),
          system(["samtools", "view", "-b", sam, '-o', output_strainbam]),
          rm(sam) ],
        depends=[HFD(fastq), bin.metaphlan2],
        targets=[output_biom, output_txt, output_strainbam],
        name="Profile taxonomy "+str(fastq),
        mem=2500,
        time=90
        )



def humann2(ctx, fastq, tax_tsv, output_tarbz, threads=2):
    tmpdir = fname.rmext(str(fastq), all=True)
    base = fname.rmext(basename(str(fastq)), all=True)
    logfile = fname.mangle("sample.log", dir=tmpdir)
    def _h2(task):
        try:
            system(["humann2", "--output-basename", base, "--log-level", "DEBUG",
                    "--remove-temp-output", "--threads", str(threads),
                    "--taxonomic-profile", tax_tsv, "--input", fastq, 
                    "--output", tmpdir, "--o-log", logfile])(task)
            system(["tar", "-cjf", output_tarbz, tmpdir])(task)
        finally:
            rm_r(tmpdir)(task)

    ctx.grid_add_task(
        _h2,
        depends=[HFD(fastq), tax_tsv, bin.humann2, h2ntdb, h2pdb],
        targets=output_tarbz,
        name="Profile function "+str(fastq),
        mem=12000,
        time=8*60,
        cores=int(threads)
        )


def _load_metaphlan_biom(fn):
    sn = re.sub(r'.biom$', '', basename(fn))
    with open(fn) as f:
        ret = StringIO(f.read().replace("Metaphlan2_Analysis", sn))
    return parse_biom_table(ret)


def _merge_metaphlan_bioms(in_bioms, out_merged_biom):
    bs = [ _load_metaphlan_biom(fn) for fn in in_bioms 
           if os.stat(fn).st_size > 0 ]
    for b in bs[1:]:
        bs[0] = bs[0].merge(b)
    with open(out_merged_biom, 'w') as f:
        bs[0].getBiomFormatJsonString("ibdmdb.org", f)


def merge_metaphlan_bioms(ctx, in_bioms, out_merged_biom):
    ctx.add_task(lambda t: _merge_metaphlan_bioms(in_bioms, out_merged_biom),
                 depends=in_bioms, targets=out_merged_biom,
                 name="Merge metphlan bioms")

    
def stack_chart(ctx, merged_biom, outdir):
    tmpdir = fname.mangle("stacktmp", dir=dirname(merged_biom))
    chartdir = fname.mangle("charts", dir=tmpdir)

    def _run(task):
        try:
            system(["summarize_taxa_through_plots.py", "-i", merged_biom, 
                    "-o", chartdir])(task)
            system(["mkdir", outdir])(task)
            system(["find", chartdir, "-name", '*_legend.pdf', 
                    '-exec', 'pdftoppm', '-singlefile', '-png', 
                    '{}', '{}', ';'])(task)
            system(["find", chartdir, "-name", "*.png",
                    "-exec", "mv", "{}", outdir, ";"])(task)
        finally:
            rm_r(tmpdir)

    ctx.add_task(_run, depends=merged_biom,
                 name="Plot taxonomy as stacked bars" )


def _load_table(f):
    for line in f:
        name, value = line.strip().split("\t",1)
        try:
            yield name, float(value)
        except:
            continue


def _geth2tables(tarbz):
    with TarFile(fileobj=BZ2File(tarbz)) as f:
        infos = f.getmembers()
        ret = dict([(i.name, dict(_load_table(f.extractfile(i))))
                    for i in infos if i.name.endswith(".tsv")])
    return ret


def _merge_write_h2tab(tables, outfname):
    names = set()
    data = dict()
    for colname, tab in tables:
        names.update(tab.keys())
        data[colname] = tab
        colnames = list(sorted(data.keys()))
    with open(outfname, 'w') as outf:
        print >> outf, "\t".join(["#pathway_or_genefamily"]+colnames)
        for name in names:
            values = [str(data[colname].get(name, 0)) for colname in colnames]
            print >> outf, "\t".join([name]+values)


def _merge_humann2_tables(tables, out_genefam, out_pathabund, out_pathcov):
    alltables = map(_geth2tables, tables)
    gf, pa, pc = list(), list(), list()
    for t in alltables:
        for k in t:
            if k.endswith("_genefamilies.tsv"):
                gf.append((basename(k).replace("_genefamilies.tsv", ""), t[k]))
            elif k.endswith("_pathabundance.tsv"):
                pa.append((basename(k).replace("_pathabundance.tsv", ""), t[k]))
            elif k.endswith("_pathcoverage.tsv"):
                pc.append((basename(k).replace("_pathcoverage.tsv", ""), t[k]))
    _merge_write_h2tab(gf, out_genefam)
    _merge_write_h2tab(pa, out_pathabund)
    _merge_write_h2tab(pc, out_pathcov)

def merge_humann2_tables(ctx, tables, out_genefam, out_pathabund, out_pathcov):
    ctx.add_task(lambda t: _merge_humann2_tables(tables,),
                 depends=tables, targets=[gf, pa, pc], 
                 name="Merge Humann2 tables")
                                                    

def load_metaphlan_table(fname):
    with open(fname, 'r') as f:
        ret = _load_table(f)
    return ret

def _merge_metaphlan_tables(in_fnames, out_fname):
    cnames = set()
    data = dict()
    with open(out_fname, 'w') as outf:
        samplename = lambda x: os.path.basename(x).split("_", 1)[0]
        for fn in in_fnames:
            l = list(load_metaphlan_table(fn))
            if not l:
                continue
            data[fn] = dict(l)
            cnames.update(data[fn].keys())
        all_fnames = list(sorted(data.keys()))
        print >> outf, "\t".join(["#Clade_name"]+map(samplename, all_fnames))
        for cname in cnames:
            values = [str(data[fn].get(cname, 0)) for fn in all_fnames]
            print >> outf, "\t".join([cname]+values)

            
def merge_metaphlan_tables(ctx, in_fnames, out_fname):
    ctx.add_task(lambda t: _merge_metaphlan_tables(in_fnames, out_fname),
                 depends=in_fnames, targets=out_fname,
                 name="Merge metaphlan taxonomic profiles")

    
def _last_meta_name(fn):
    prev_line = str()
    with open(fn) as f:
        for line in f:
            if re.search(r'[Bb]acteria|[Aa]rchaea.*\s+\d', line):
                return prev_line.split('\t')[0]
            prev_line = line
        return prev_line.split('\t')[0]
            

def _sample_id(fn):
    id_ = str()
    with open(fn) as f:
        for line in f:
            if line.startswith("#"):
                id_ = line.split('\t')[0]
                continue
            else:
                return id_ or line.split('\t')[0]

    
def pcoa_chart(ctx, merged_tax, chart_fname):
    tmpdir = fname.mangle("pcoatmp", dir=os.path.dirname(merged_tax))
    
    def _pcoa(task):
        try:
            system(["scriptPcoa.py", "--meta", _last_meta_name(merged_tax),
                    "--id", _sample_id(merged_tax), "--noShape", 
                    "--outputFile", fname.mangle("chart.png", dir=tmpdir),
                    merged_tax])(task)
            system(["find", tmpdir, '-name', "*.png", 
                    "-exec", "mv", "{}", chart_fname, ";"])(task)
        finally:
            rm_r(tmpdir)

    ctx.add_task(_pcoa,
                 targets=chart_fname, depends=merged_tax,
                 name="Plot taxonomy as ordination")


def compress(ctx, fastq):
    ctx.add_task(system(["pbzip2", fastq]),
                 depends=[fastq],
                 targets=str(fastq)+".bz2",
                 name="Compress "+str(fastq))


def _tabsave(headers, row, output):
    def _actually_tabsave(task):
        with open(output, 'w') as f:
            print >> f, "\t".join( map(str, headers) )
            print >> f, "\t".join( map(str, row) )
            
    return _actually_tabsave


def save_metadata(ctx, metadata_fname, headers, data_row, output_tsv):
    ctx.add_task( _tabsave(headers, data_row, output_tsv),
                  depends=metadata_fname,
                  targets=output_tsv,
                  name="Save metadata slice "+output_tsv )


def process_amplicon(ctx, data, metadata_json, headers):
    for row in data:
        ctx.already_exists(HFD(row.file))
        dirconf = output_dirs.amplicon

        metadata_slice = fname.mangle(row.file, dir=dirconf.taxprof, ext="tsv")
        cleanfq = fname.mangle(row.file, dir=dirconf.cleanseq, ext="fastq")
        cleanlog = fname.mangle(row.file, dir=dirconf.cleanseq, tag="clean", ext="log")
        biom = fname.mangle(row.file, dir=dirconf.taxprof, ext="biom")
        taxtxt = fname.mangle(row.file, dir=dirconf.taxprof, tag="tax", ext="txt")
        tarbz = fname.mangle(row.file, dir=dirconf.h2, tag="picrust", ext="tar.bz2")

        ctx.already_exists(HFD(row.file))
        save_metadata(ctx, metadata_json, headers, row.data, metadata_slice)
        convert(ctx, HFD(row.file), cleanfq)
        otupick(ctx, cleanfq, biom, taxtxt, cleanlog)
        picrust(ctx, biom, tarbz)

        
def process_wgs_mtx(ctx, wgs_data, mtx_data, metadata_json, headers, datatype_index):
    if len(wgs_data) > 1:
        raise Exception("too many wgs data for sample "+wgs_data[0].data[3])
    elif len(wgs_data) < 1:
        return
    for row in concatv(wgs_data, mtx_data):
        dirconf = output_dirs[row.data[datatype_index]]
        metadata_slice = fname.mangle(row.file, dir=dirconf.taxprof, ext="tsv")
        save_metadata(ctx, metadata_json, headers, row.data, metadata_slice)
    cleandir, tpdir, h2dir = get(["cleanseq", "taxprof", "h2"], output_dirs.metagenomics)
    ctx.already_exists(HFD(wgs_data[0].file))
    cleanfq = fname.mangle(wgs_data[0].file, dir=cleandir, ext="fastq")
    knead(ctx, HFD(wgs_data[0].file), cleanfq,
          fname.mangle(wgs_data[0].file, dir=cleandir, tag="clean", ext="log"))
    wgs_taxtxt = fname.mangle(wgs_data[0].file, dir=tpdir, tag="tax", ext="txt")
    metaphlan2(ctx, cleanfq, 
               fname.mangle(wgs_data[0].file, dir=tpdir, ext="biom"),
               wgs_taxtxt, 
               fname.mangle(wgs_data[0].file, dir=tpdir, tag="strain", ext="bam"))
    humann2(ctx, cleanfq, wgs_taxtxt, 
            fname.mangle(wgs_data[0].file, dir=h2dir, tag="humann2", ext="tar.bz2"))

    for mtx_row in mtx_data:
        dirconf = output_dirs.metatranscriptomics
        cleanfq = fname.mangle(mtx_row.file, dir=dirconf.cleanseq, ext="fastq")
        cleanlog = fname.mangle(mtx_row.file, dir=dirconf.cleanseq, tag="clean", ext="log")
        biom = fname.mangle(mtx_row.file, dir=dirconf.taxprof, ext="biom")
        strainbam = fname.mangle(mtx_row.file, dir=dirconf.taxprof, tag="strain", ext="bam")
        taxtxt = fname.mangle(mtx_row.file, dir=dirconf.taxprof, tag="tax", ext="txt")
        tarbz = fname.mangle(mtx_row.file, dir=dirconf.h2, tag="humann2", ext="tar.bz2")

        ctx.already_exists(HFD(mtx_row.file))
        knead(ctx, HFD(mtx_row.file), cleanfq, cleanlog, 
              dbs=kneaddb.values(), db_deps=kneaddb_deps.values())
        metaphlan2(ctx, cleanfq, biom, taxtxt, strainbam)
        humann2(ctx, cleanfq, wgs_taxtxt, tarbz)


def generate_figures(ctx):
    today = datetime.now().strftime("%F")
    for d in output_dirs.values():
        in_tax = [fn for fn in glob.glob(d.taxprof+"/*_tax.txt")
                  if os.stat(fn).st_size > 0]
        in_biom = glob.glob(d.taxprof+"/*.biom")
        in_tarbz = glob.glob(d.h2+"/*.tar.bz2")
        for fn in in_tax+in_biom+in_tarbz:
            ctx.already_exists(fn)
        merged_tax = fname.mangle("merged-metaphlan-table", tag=today, 
                                  dir=d.taxprof, ext="tsv")
        merged_biom = fname.mangle(merged_tax, ext="biom")
        gf = fname.mangle("merged_genefamilies", dir=d.h2, tag=today, ext="tsv")
        pa = fname.mangle("merged_pathabundance", dir=d.h2, tag=today, ext="tsv")
        pc = fname.mangle("merged_pathcoverage", dir=d.h2, tag=today, ext="tsv")
        figdir = fname.mangle("figures", dir=d.taxprof)
        merge_metaphlan_tables(ctx, in_tax, merged_tax)
        merge_metaphlan_bioms(ctx, in_biom, merged_biom)
        merge_humann2_tables(ctx, in_tarbz, gf, pa, pc)
        stack_chart(ctx, merged_biom, fname.mangle("stacked", dir=figdir))
        pcoa_chart(ctx, merged_tax, fname.mangle("taxonomy_pcoa.png", dir=figdir))
        pcoa_chart(ctx, gf, fname.mangle("genefamily_pcoa.png", dir=figdir))
        pcoa_chart(ctx, pa, fname.mangle("pathabundance_pcoa.png", dir=figdir))
        
    
def main():
    with open(str(metadata_json), 'r') as f:
        d = Bunch(json.load(f))
    d.data = apply_ignores(d.data)
    datatype_idx = d.headers.index("data_type")
    fnames = filter(bool, map(findfile, d.data))
    d.data = [ Bunch(data=r, file=fn) for r, fn in zip(d.data, fnames) ]
    subjs = groupby(lambda b: b.data[3], d.data)
    
    ctx = RunContext(grid_powerup=SGEPowerup(queue="long", tmpdir="/seq/ibdmdb/tmp/anasge"))
    reporter = LoggerReporter(ctx, loglevel_str='debug', 
                              logfile="/seq/ibdmdb/anadama_run.log")
    ctx.already_exists(metadata_json)
    ctx.already_exists(*bin.values())
    ctx.already_exists(h2ntdb, h2pdb)
    ctx.already_exists(*kneaddb_deps.values())
    ctx.already_exists(*map(HFD, conf.values()))

    for subj in subjs.values():
        bytype = groupby(lambda s: s.data[4], subj)
        if "amplicon" in bytype:
            process_amplicon(ctx, bytype['amplicon'], metadata_json, d.headers)
        process_wgs_mtx(ctx, bytype.get("metagenomics", []), 
                        bytype.get("metatranscriptomics", []),
                        metadata_json, d.headers, datatype_idx)

    try:
        ctx.go(reporter=reporter, n_grid_parallel=35, n_parallel=2)
    except anadama.runcontext.RunFailed:
        pass

    generate_figures(ctx)
    ctx.go(reporter=reporter, n_grid_parallel=0, n_parallel=2)
    # ctx.cli()


if __name__ == '__main__':
    main()
