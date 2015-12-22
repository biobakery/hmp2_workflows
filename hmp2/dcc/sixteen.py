import os
import sys
import hashlib
from operator import attrgetter
import multiprocessing

import cutlass
import cutlass.mimarks
from toolz import groupby, compose, first
import anadama.util

from .subject import fields
from .visit import SaveError

def default_mimarks_dict():
    return dict([ (k, v()) for k, v in
                  cutlass.mimarks.MIMARKS._fields.iteritems() ])

class settings:
    library_method_text = \
        """ Samples are diluted and run in a 30 cycle PCR that specifically
        amplifies the V4 region of the 16S ribosomal subunit, with the
        primers making use of the highly conserved regions on either
        side. PCR primers are barcoded by well, in a dual plate format,
        for a total of 192 available barcodes. Samples are quantified
        using the caliper, and then pooled together. Finally, the pooled
        sample is cleaned up and selected for size by using a combination
        of purification columns and a Pippen Prep SAGE protocol, and
        handed off for MiSeq. """
    prep_id = "1"
    

def default_prep(sample):
    prep = cutlass.SixteenSDnaPrep()
    prep.comment = "Broad IBDMDB default 16S dna prep"
    prep.ncbi_taxon_id = "408170" # human gut metagenomes
    prep.lib_layout = "fragment"
    prep.lib_selection = "Random"
    prep.prep_id = settings.prep_id
    prep.mimarks = default_mimarks_dict()
    prep.mimarks['target_gene'] = "16S"
    prep.mimarks['lib_const_meth'] = settings.library_method_text
    prep.mimarks['target_subfragment'] = "V4"
    prep.sequencing_center = "Broad Institute"
    prep.sequencing_contact = "tpoon@broadintitute.org"
    prep.storage_duration = 365 # 1 year in days
    prep.links['prepared_from'] = [sample.id]
    return prep

seqtype_docs = {
    "raw": ("Raw 16S ", "Raw, untrimmed "),
    "trimmed": ("Trimmed 16S","Trimmed ")
}

def parse_record(r, sample, fname, md5sum, seqtype="raw"):
    prep = default_prep(sample)
    com, doc = seqtype_docs[seqtype]
    seq = cutlass.SixteenSRawSeqSet()
    seq.comment = com+str(sample)
    seq.checksums = {"md5": md5sum}
    seq.exp_length = 175
    seq.format = anadama.util.guess_seq_filetype(fname)
    seq.format_doc = doc+seq.format
    seq.local_file = fname
    seq.seq_model = "Illumina MiSeq"
    seq.size = os.stat(fname).st_size
    seq.study = "ibd"
    return (prep, seq)
    

md5_buffer_len = 4096
def md5(fname):
    m = hashlib.md5()
    with open(fname, 'rb') as f:
        while True:
            buf = f.read(md5_buffer_len)
            if not buf:
                break
            m.update(buf)
    return m.hexdigest()

def rm_common_prefix(list_fnames):
    p = os.path.commonprefix(list_fnames)
    return p, [ fname.lstrip(p) for fname in list_fnames ]

def find_fname(substr, fnames, commonprefix):
    matches = [fname for fname in fnames if substr in fname]
    if not matches:
        return
    return commonprefix+matches[0]
    
def from_file(fname, samples, seq_fnames, sample_id_idx=3,
              fname_idx=0, nest=True, n_procs=8, parse_record=parse_record):
    samples = groupby(compose(first, attrgetter("tags")), samples)
    pool = multiprocessing.Pool(n_procs)
    md5s = pool.map(md5, seq_fnames)
    md5s = dict([(f,m) for f,m in zip(seq_fnames, md5s)])
    commonprefix, fnames = rm_common_prefix(seq_fnames)
    for record in fields(fname):
        sampleid = record[sample_id_idx]
        if sampleid not in samples:
            continue
        seq_fname = find_fname(record[fname_idx], fnames, commonprefix)
        if not seq_fname:
            continue
        ret = parse_record(record, samples.get(sampleid)[0],
                           seq_fname, md5s[seq_fname])
        if nest:
            yield ret
        else:
            for item in ret:
                yield item


def save_all(ps_groups):
    valids = [ p.is_valid() for p, _ in ps_groups ]
    if not all(valids):
        valid_errors = [ None if valid else ps[0].validate()
                         for valid, ps in zip(valids, ps_groups) ]
        raise SaveError(valid_errors)
    for i, (p, seq) in enumerate(ps_groups):
        p.save()
        seq.links['sequenced_from'] = [p.id]
        valids[i] = seq.is_valid()
    if not all(valids):
        valid_errors = [ None if valid else ps[1].validate()
                         for valid, ps in zip(valids, ps_groups) ]
        raise SaveError(valid_errors)
    for _, seq in ps_groups:
        seq.save()
        


def main():
    pass

if __name__ == '__main__':
    sys.exit(main())
