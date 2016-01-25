import sys
import csv
from collections import OrderedDict
from os.path import basename

from anadama.util import guess_seq_filetype
from toolz import count as iter_len

from . import fname_search

default_kvs = lambda : OrderedDict([
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

def output(out_dicts, outtsv_fname, rmkeys):
    with open(outtsv_fname, 'w') as f:
        d = default_kvs()
        for k in rmkeys:
            d.pop(k, None)
        print >> f, "\t".join(d.iterkeys())
        for d in out_dicts:
            print >> f, "\t".join(d.itervalues())
    

def count_fastq(fname):
    with open(fname, 'r') as f:
        l = iter_len(iter(f))
    return l/4


def count_fasta(fname):
    with open(fname, 'r') as f:
        l = iter_len(iter(f))
    return l/2


def readcount(fname):
    ftype = guess_seq_filetype(fname)
    if ftype == "fastq":
        return count_fastq(fname)
    elif ftype == "fasta":
        return count_fasta(fname)
    else:
        raise ValueError("Unrecognized sequence type: "+ftype)


def create(fna_fnames, incsv_fname, outtsv_fname, skip_header=True,
           rmkeys=(), sample_name_idx=18, fname_prefix_idx=16):
    sizes = dict([ (fname, readcount(fname)) for fname in fna_fnames ])
    def _out_dicts():
        with open(incsv_fname, 'r') as in_f:
            reader = csv.reader(in_f)
            if skip_header:
                next(reader)
            for line in reader:
                d = default_kvs()
                for k in rmkeys:
                    d.pop(k, None)
                d['sample_name'] = line[sample_name_idx]
                seq_fname = fname_search(fna_fnames, line[fname_prefix_idx])
                d['RUN_PREFIX'] = basename(seq_fname)
                if seq_fname not in sizes:
                    msg = ("Unable to find sequence file "
                           "matching "+line[fname_prefix_idx])
                    print >> sys.stderr, msg
                    continue
                d['samp_size'] = str(sizes[seq_fname])
                yield d

    return output(_out_dicts(), outtsv_fname, rmkeys)
