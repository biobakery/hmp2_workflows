
import cutlass
import cutlass.mims

from .sixteen import from_file


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
    

def default_prep(sample):
    prep = cutlass.SixteenSDnaPrep()
    prep.comment = "Broad IBDMDB default WGS dna prep"
    prep.ncbi_taxon_id = "408170" # human gut metagenomes
    prep.lib_layout = "fragment"
    prep.lib_selection = "Random"
    prep.prep_id = settings.prep_id
    prep.mims = default_mims_dict()
    prep.mims['lib_const_meth'] = settings.library_method_text
    prep.sequencing_center = "Broad Institute"
    prep.sequencing_contact = "tpoon@broadintitute.org"
    prep.storage_duration = 365 # 1 year in days
    prep.links['prepared_from'] = [sample.id]
    return prep


def parse_record(r, sample, fname, md5sum, seqtype="raw"):
    pass
