# -*- coding: utf-8 -*-
import os

import cutlass
import cutlass.mims

import anadama.util

from .sixteen import seq_len_mode


def default_mims_dict():
    return dict([ (k, v()) for k, v in
                  cutlass.mims.MIMS._fields.iteritems() ])


class settings:
    library_method_text = \
        """Metagenomic DNA samples were quantified by Quant-iT PicoGreen dsDNA
        Assay (Life Technologies) and normalized to a concentration of
        50 pg μL-1. Illumina sequencing libraries were prepared from
        100-250 pg DNA using the Nextera XT DNA Library Preparation
        kit (Illumina) according to the manufacturer’s recommended
        protocol, with reaction volumes scaled accordingly. Batches of
        24, 48, or 96 libraries were pooled by transferring equal
        volumes of each library using a Labcyte Echo 550 liquid
        handler. Insert sizes and concentrations for each pooled
        library were determined using an Agilent Bioanalyzer DNA 1000
        kit (Agilent Technologies)."""
    prep_id = "2"
    sequencing_contact = "tpoon@broadintitute.org"
    sequencing_center = "Broad Institute"
    ncbi_taxon_id = "408170" # human gut metagenomes
    lib_layout = "fragment"
    lib_selection = "Random"
    prep_comment = "Broad IBDMDB default WGS dna prep"
    seq_comment = "Raw WGS Sequence set "
    storage_duration = 365 # 1 year in days
    sequencer_model = "Illumina HiSeq"
    

def default_prep(sample):
    maybe_extant_preps = list(sample.wgsDnaPreps())
    if maybe_extant_preps:
        prep = maybe_extant_preps[0]
    else:
        prep = cutlass.WGSDnaPrep()
    prep.comment = settings.prep_comment
    prep.ncbi_taxon_id = settings.ncbi_taxon_id
    prep.lib_layout = settings.lib_layout
    prep.lib_selection = settings.lib_selection
    prep.prep_id = settings.prep_id
    prep.mims = default_mims_dict()
    prep.mims['lib_const_meth'] = settings.library_method_text
    prep.sequencing_center = settings.sequencing_center
    prep.sequencing_contact = settings.sequencing_contact
    prep.storage_duration = settings.storage_duration
    prep.links['prepared_from'] = [sample.id]
    return prep



def parse_record(r, sample, fname, md5sum, seqtype="raw"):
    prep = default_prep(sample)
    maybe_seqs = list(prep.raw_seq_sets())
    if maybe_seqs:
        seq = maybe_seqs[0]
    else:
        seq = cutlass.WgsRawSeqSet()
    seq.comment = settings.seq_comment + str(sample)
    seq.checksums = {"md5": md5sum}
    seq.exp_length = seq_len_mode(fname)
    seq.format = anadama.util.guess_seq_filetype(fname)
    seq.format_doc = "Raw WGS "+seq.format
    seq.local_file = fname
    seq.seq_model = settings.sequencer_model
    seq.size = os.stat(fname).st_size
    seq.study = "ibd"
    return (prep, seq)
