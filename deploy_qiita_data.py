import os
import sys
import json
import shutil
from os.path import basename
from glob import glob
from collections import OrderedDict

from bunch import Bunch
from toolz import groupby
from toolz import count as iter_len

from all_pipeline import output_dirs, metadata_json
from anadama.util.fname import rmext

with open(str(metadata_json)) as f:
    meta = Bunch(json.load(f))

grp = Bunch(groupby(meta.headers.index("data_type"), meta.data))

qiitadir = "/seq/ibdmdb/centos6/qiita-data/HMP2"
qiitatypes = {
    "amplicon": "16S",
    "metagenomics": "WGS",
    "metatranscriptomics": "WTS"
    }

default_prep_kvs = lambda : OrderedDict([
            ("sample_name", ""),
            ("BARCODE",""),
            ("EXPERIMENT_CENTER","Broad Institute"),
            ("EXPERIMENT_DESIGN_DESCRIPTION","Fecal samples in time from CD, UC, and healthy controls"),
            ("EXPERIMENT_TITLE","The_Inflammatory_Bowel_Disease_Multi_omics_Database"),
            ("KEY_SEQ",""),
            ("LIBRARY_CONSTRUCTION_PROTOCOL","30 cycle PCR"),
            ("LINKER",""),
            ("PLATFORM","Illumina MiSeq"),
            ("PRIMER",""),
            ("REGION","V4"),
            ("RUN_CENTER","Broad Institute"),
            ("RUN_DATE","unknown"),
            ("RUN_PREFIX",""),
            ("SAMPLE_CENTER","Broad Institute"),
            ("STUDY_CENTER","Broad Institute"),
            ("samp_size",".1,g"),
            ("sequencing_meth","Illumina Dye Technology"),
            ("target_gene","16S rRNA"),
            ("target_subfragment","V4"),
            ("pcr_primers",""),
            ("emp_status","NOT_EMP"),
            ("center_name","Broad Institute"),
            ("center_project_name","ibdmdb"),
            ("linkerprimersequence", "GATC"),
            ("barcodesequence", "GATC")
        ])

default_samp_kvs = lambda : OrderedDict([
            ("sample_name",""),
            ("collection_timestamp",""),
            ("physical_location","Broad Insitute 415 Main Street Cambridge, MA 02142"),
            ("taxon_id","9606"),
            ("description",""),
            ("scientific_name","Homo sapiens"),
            ("sample_type","Stool"),
            ("physical_specimen_remaining","TRUE"),
            ("dna_extracted","TRUE"),
            ("latitude",""),
            ("longitude",""),
            ("host_subject_id",""),
            ("preprocessed_name", ""),
            ("has_extracted_data", "TRUE"),
            ("required_sample_info_status", "completed"),
            ("has_physical_specimen", "TRUE"),
        ])

latlon = {"C": (34.07533, -118.39253), #c.sinai
          "H": (39.13787, -84.50358),  #cin. children's
          "M": (42.36279, -71.06807),  #mass. general
          "E": (33.79322, -84.32333)}  #emory


def output(out_dicts, outtsv_fname, headers):
    with open(outtsv_fname, 'w') as f:
        print >> f, "\t".join(headers)
        for d in out_dicts:
            print >> f, "\t".join(d.itervalues())


def _count_fastq(fname):
    with open(fname, 'r') as f:
        l = iter_len(iter(f))
    return l/4


def count_fastq(fname):
    rcfile = fname+".readcount"
    rc = None
    if os.path.exists(rcfile):
        try:
            with open(rcfile) as f:
                rc = int(f.read().strip())
        except Exception as e:
            print "Exception parsing {}: {}".format(rcfile, e)

    if not rc:
        rc = _count_fastq(fname)
        with open(rcfile, 'w') as f:
            print >> f, rc
    return rc


def count_seqs(datatype):
    fqdir = output_dirs[datatype].cleanseq
    biomdir = output_dirs[datatype].taxprof
    fqsizes = dict([ (basename(fname), count_fastq(fname)) 
                     for fname in glob(fqdir+"/*.fastq") ])
    biomsizes = dict([ (basename(fname), os.stat(fname).st_size)
                       for fname in glob(biomdir+"/*.biom") ])
    for biom_fn in biomsizes.keys():
        if not biomsizes[biom_fn] > 0:
            fqsizes.pop(biom_fn.replace("biom", "fastq"))
        if "merged" in biom_fn:
            biomsizes.pop(biom_fn)
    for name in fqsizes.keys():
        if name.replace("fastq","biom") not in biomsizes:
            fqsizes.pop(name)
    return fqsizes, biomsizes


def prep_create(datatype, sizes):
    md = groupby(1, grp[datatype])
    qt = qiitatypes[datatype]
    outtsv_fname = os.path.join(qiitadir, qt, "qiita_prep_"+qt.lower()+".tsv")
    
    def _out_dicts():
        for name, size in sizes.iteritems():
            md_row = md[rmext(name, True)][0]
            d = default_prep_kvs()
            d['sample_name'] = md_row[1]
            d['RUN_PREFIX'] = name
            d['samp_size'] = str(size)
            yield d
                
    return output(_out_dicts(), outtsv_fname, default_prep_kvs().keys())


def samp_create(datatype, sizes):
    md = groupby(1, grp[datatype])
    qt = qiitatypes[datatype]
    outtsv_fname = os.path.join(qiitadir, qt, "qiita_samples_"+qt.lower()+".tsv")
    
    def _out_dicts():
        for name, size in sizes.iteritems():
            md_row = md[rmext(name, True)][0]
            d = default_prep_kvs()
            d['sample_name'] = md_row[1]
            d['latitude'], d['longitude'] = map(
                str, latlon.get(md_row[3][0], ("", "")) )
            d['preprocessed_name'] = name
            d['collection_timestamp'] = "2000-01-01 12:00:00"
            d['description'] = md_row[1]
            d['host_subject_id'] = md_row[2]
            yield d
                
    return output(_out_dicts(), outtsv_fname, default_prep_kvs().keys())
    

def process(datatype):
    qdir = os.path.join(qiitadir, qiitatypes[datatype])
    fqsizes, biomsizes = count_seqs(datatype)
    for f in os.listdir(qdir):
        os.remove(os.path.join(qdir, f))
    samp_create(datatype, fqsizes)
    prep_create(datatype, fqsizes)
    biomdir = output_dirs[datatype].taxprof
    for name in biomsizes:
        shutil.copy(os.path.join(biomdir, name), qdir)

        
def main():
    process("amplicon")
    process("metagenomics")
    process("metatranscriptomics")

if __name__ == '__main__':
    sys.exit(main())
