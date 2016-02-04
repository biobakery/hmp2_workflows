import re
import os

from anadama.pipelines import Pipeline
from anadama_workflows import settings as aw_settings

from .dcc import submit
from .qiita import prep 
from .qiita import sample
from . import participant

class ExternalMetadataPipeline(Pipeline):
    name = "ExternalMetadata"
    products = {
        "hospital_metadata":   list(),
        "wgs_sample_metadata": list(),
        "wts_sample_metadata": list(),
        "six_sample_metadata": list(),

        "six_raw_seqs": list(),
        "wgs_raw_seqs": list(),
        "wts_raw_seqs": list(),

        "six_otu_tables": list(),
        "wgs_otu_tables": list(),
        "wts_otu_tables": list(),

        "participant_metadata": list(),
        "qiita_preps": list(),
        "qiita_samples": list(),
    }
    default_options = { # first field's index is 1 (not zero-based)
        "participant_join": {
        },
        "qiita_join_16s": {
            "fname_prefix_idx": 1,
            "sample_name_idx": 1,
            "subject_idx": 18,
            "collection_time_idx": -2,
            "description_idx": 1
        },
        "qiita_join_wgs": {
            "fname_prefix_idx": 2,
            "sample_name_idx": 2,
            "subject_idx": 21,
            "collection_time_idx": 27,
            "description_idx": 2
        },
        "qiita_join_wts": {
            "fname_prefix_idx": 2,
            "sample_name_idx": 1,
            "subject_idx": 9,
            "collection_time_idx": 15,
            "description_idx": 1
        },
        "dcc_upload": {
        }
    }

    def __init__(self,
                 hospital_metadata=list(),
                 wgs_sample_metadata=list(),
                 wts_sample_metadata=list(),
                 six_sample_metadata=list(),

                 six_raw_seqs=list(),
                 wgs_raw_seqs=list(),
                 wts_raw_seqs=list(),
                 
                 six_otu_tables=list(),
                 wgs_otu_tables=list(),
                 wts_otu_tables=list(),

                 participant_metadata=list(),
                 qiita_preps=list(),
                 qiita_samples=list(),
                 products_dir=None,
                 workflow_options=dict(),
                 *args, **kwargs):

        super(ExternalMetadataPipeline, self).__init__(*args, **kwargs)

        self.add_products(
            hospital_metadata=hospital_metadata,
            wgs_sample_metadata=wgs_sample_metadata,
            wts_sample_metadata=wts_sample_metadata,
            six_sample_metadata=six_sample_metadata,
            
            six_raw_seqs=six_raw_seqs,
            wgs_raw_seqs=wgs_raw_seqs,
            wts_raw_seqs=wts_raw_seqs,
                 
            six_otu_tables=six_otu_tables,
            wgs_otu_tables=wgs_otu_tables,
            wts_otu_tables=wts_otu_tables,
            
            participant_metadata=list(),
            qiita_preps=list(),
            qiita_samples=list(),
        )

        if not products_dir:
            products_dir = aw_settings.workflows.product_directory
        self.products_dir = os.path.abspath(products_dir)

        self.options = self.default_options.copy()
        self.options.update(workflow_options)

        for w_k, v in self.options.iteritems():
            for k, i in v.iteritems():
                if k.endswith("idx"):
                    self.options[w_k][k] = i-1


    def _configure(self):
        if not os.path.isdir(self.products_dir):
            os.mkdir(self.products_dir)

        yield {
            "name": "submit_to_dcc",
            "uptodate": [lambda *a, **k: False],
            "targets": [None],
            "actions": [lambda: submit(self.hospital_metadata[0], 
                                       self.six_sample_metadata[0],
                                       self.six_raw_seqs,
                                       self.wgs_sample_metadata[0],
                                       self.wgs_raw_seqs,
                                       self.wts_sample_metadata[0],
                                       self.wts_raw_seqs,
                                       self.options['dcc_upload']['dcc_user'],
                                       self.options['dcc_upload']['dcc_pw'])],
        }

        pjoin = lambda s: os.path.join(self.products_dir, s)
        yield qiita_task(pjoin("qiita_samples_16s.tsv"),
                         pjoin("qiita_prep_16s.tsv"),
                         self.six_raw_seqs, self.six_sample_metadata[0],
                         self.options['qiita_join_16s'])
        yield qiita_task(pjoin("qiita_samples_wgs.tsv"),
                         pjoin("qiita_prep_wgs.tsv"),
                         self.wgs_raw_seqs, self.wgs_sample_metadata[0],
                         self.options['qiita_join_wgs'],
                         rmkeys=("PRIMER", "REGION", "target_gene",
                                  "LINKER", "PRIMER", "KEY_SEQ", "BARCODE",))
        yield qiita_task(pjoin("qiita_samples_wts.tsv"),
                         pjoin("qiita_prep_wts.tsv"),
                         self.wts_raw_seqs, self.wts_sample_metadata[0],
                         self.options['qiita_join_wts'],
                         rmkeys=("PRIMER", "REGION", "target_gene",
                                  "LINKER", "PRIMER", "KEY_SEQ", "BARCODE",))

        six_args = (self.six_sample_metadata[0], self.six_otu_tables,
                    self.options['qiita_join_16s']['subject_idx'], # subj_idx
                    lambda r: r[0]+".bam.merged.fasta_mangled_otus_otu_table.biom",
                    visit_id_gen(21),)
        wgs_args = (self.wgs_sample_metadata[0], self.wgs_otu_tables,
                    self.options['qiita_join_wgs']['subject_idx'],
                    lambda r: r[1]+".fastq.metaphlan2.biom",
                    visit_id_gen(24),)
        yield {
            "name":"participant_metadata",
            "targets": [pjoin("participant_metadata.tsv")],
            "actions": [lambda: participant.create(pjoin("participant_metadata.tsv"),
                                                   six_args, wgs_args)]
        }

def visit_id_gen(idx):
    def get(record):
        visit_str = record[idx]
        if "Baseline" in visit_str:
            return 1
        else:
            return re.sub(r'\D', '', visit_str)

    return get


def qiita_task(sample_out, prep_out, raw_seqfiles, sample_metadata,
               options_dict, rmkeys=() ):
    yield {
        "name": "qiita_join:sample: "+sample_out,
        "targets": [sample_out],
        "file_dep": [sample_metadata],
        "actions": [
            lambda *a, **kw: sample.create(
                raw_seqfiles,
                sample_metadata,
                sample_out, True,
                options_dict['sample_name_idx'],
                options_dict['subject_idx'],
                options_dict['fname_prefix_idx'],
                options_dict['collection_time_idx'],
                options_dict['description_idx']
            )
        ]
    }

    yield {
        "name": "qiita_join:prep: "+prep_out,
        "targets": [prep_out],
        "file_dep": [sample_metadata],
        "actions": [
            lambda *a, **kw: prep.create(
                raw_seqfiles,
                sample_metadata,
                prep_out, True, rmkeys,
                options_dict['sample_name_idx'],
                options_dict['fname_prefix_idx'],
            )
        ]
    }
