import os

from anadama.pipelines import Pipeline
from anadama_workflows import settings as aw_settings

from .qiita import prep 
from .qiita import sample

class ExternalMetadataPipeline(Pipeline):
    name = "ExternalMetadata"
    products = {
        "wgs_sample_metadata": list(),
        "wms_sample_metadata": list(),
        "six_sample_metadata": list(),

        "six_raw_seqs": list(),
        "wgs_raw_seqs": list(),
        "wms_raw_seqs": list(),

        "six_otu_tables": list(),
        "wgs_otu_tables": list(),
        "wms_otu_tables": list(),

        "participant_metadata": list(),
        "qiita_preps": list(),
        "qiita_samples": list(),
    }
    default_options = { # first field's index is 1 (not zero-based)
        "participant_join": {
        },
        "qiita_join": {
            "fname_prefix_idx": 1,
            "sample_name_idx": 1,
            "subject_idx": 2,
            "collection_time_idx": -2,
            "description_idx": 1
        },
        "dcc_upload": {
        }
    }

    def __init__(self,
                 wgs_sample_metadata=list(),
                 wms_sample_metadata=list(),
                 six_sample_metadata=list(),

                 six_raw_seqs=list(),
                 wgs_raw_seqs=list(),
                 wms_raw_seqs=list(),
                 
                 six_otu_tables=list(),
                 wgs_otu_tables=list(),
                 wms_otu_tables=list(),

                 participant_metadata=list(),
                 qiita_preps=list(),
                 qiita_samples=list(),
                 products_dir=None,
                 workflow_options=dict(),
                 *args, **kwargs):

        super(ExternalMetadataPipeline, self).__init__(*args, **kwargs)

        self.add_products(
            wgs_sample_metadata=wgs_sample_metadata,
            wms_sample_metadata=wms_sample_metadata,
            six_sample_metadata=six_sample_metadata,
            
            six_raw_seqs=six_raw_seqs,
            wgs_raw_seqs=wgs_raw_seqs,
            wms_raw_seqs=wms_raw_seqs,
                 
            six_otu_tables=six_otu_tables,
            wgs_otu_tables=wgs_otu_tables,
            wms_otu_tables=wms_otu_tables,
            
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

        pjoin = lambda s: os.path.join(self.products_dir, s)
        yield qiita_task(pjoin("qiita_samples_16s.tsv"),
                         pjoin("qiita_prep_16s.tsv"),
                         self.six_raw_seqs, self.six_sample_metadata[0],
                         self.options['qiita_join'])
        yield qiita_task(pjoin("qiita_samples_wgs.tsv"),
                         pjoin("qiita_prep_wgs.tsv"),
                         self.wgs_raw_seqs, self.wgs_sample_metadata[0],
                         self.options['qiita_join'],
                         rmkeys=("PRIMER", "REGION", "target_gene",
                                  "LINKER", "PRIMER", "KEY_SEQ", "BARCODE",))
        yield qiita_task(pjoin("qiita_samples_wms.tsv"),
                         pjoin("qiita_prep_wms.tsv"),
                         self.wms_raw_seqs, self.wms_sample_metadata[0],
                         self.options['qiita_join'],
                         rmkeys=("PRIMER", "REGION", "target_gene",
                                  "LINKER", "PRIMER", "KEY_SEQ", "BARCODE",))



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
