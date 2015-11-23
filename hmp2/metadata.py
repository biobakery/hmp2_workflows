import re
import os
import operator

from anadama.pipelines import Pipeline

from anadama_workflows import settings as aw_settings

from . import join_metadata

class BroadMetadataPipeline(Pipeline):
    name = "BroadMetadata"
    products = {
        "hospital_sample_csv": list(),
        "hospital_subject_csv": list(),
        "broad_16s_csv": list(),
        "broad_wgs_csv": list(),
        "broad_wts_csv": list(),

        "wgs_sample_metadata": list(),
        "wms_sample_metadata": list(),
        "six_sample_metadata": list(),

        "otu_tables_wgs": list(),
        "otu_tables_wms": list(),
        "otu_tables_16s": list(),
    }
    default_options = { # first field's index is 1 (not zero-based)
        "hospital_join": {
            "subject_subject_idx": 1,
            "sample_subject_idx": 2,
            "sample_stool_idx": 9,
        },
        "16s_join": {
            "sample_idx": 1
        },
        "wgs_join": {
            "sample_idx": 2
        },
        "wts_join": {
            "sample_idx": 1
        }
    }

    def __init__(self,
                 hospital_subject_csv=list(),
                 hospital_sample_csv=list(),
                 broad_16s_csv=list(),
                 broad_wgs_csv=list(),
                 broad_wts_csv=list(),
                 otu_tables_wgs=list(),
                 otu_tables_wms=list(),
                 otu_tables_16s=list(),
                 products_dir=None,
                 workflow_options=dict(),
                 *args, **kwargs):

        super(BroadMetadataPipeline, self).__init__(*args, **kwargs)

        self.add_products(
            hospital_subject_csv = hospital_subject_csv,
            hospital_sample_csv  = hospital_sample_csv,
            broad_16s_csv        = broad_16s_csv,
            broad_wgs_csv        = broad_wgs_csv,
            broad_wts_csv        = broad_wts_csv,
            otu_tables_wgs       = otu_tables_wgs,
            otu_tables_wms       = otu_tables_wms,
            otu_tables_16s       = otu_tables_16s,

            wgs_sample_metadata  = list(),
            wms_sample_metadata  = list(),
            six_sample_metadata  = list(),
        )

        if not products_dir:
            products_dir = aw_settings.workflows.product_directory
        self.products_dir = os.path.abspath(products_dir)

        self.options = self.default_options.copy()
        self.options.update(workflow_options)

        for w_k, v in self.options.iteritems():
            for k, i in v.iteritems():
                self.options[w_k][k] = i-1

    def _configure(self):
        if not os.path.isdir(self.products_dir):
            os.mkdir(self.products_dir)

        with open(self.hospital_subject_csv[0], 'r') as f:
            hospital_sample_width = iter(f).next().count(',')

        a, b = self.hospital_subject_csv[0], self.hospital_sample_csv[0]
        t = hospital_joined = os.path.join(self.products_dir,
                                          "hospital_joined.csv")
        k_a, k_b = map(self.options['hospital_join'].get,
                       ("subject_subject_idx", "sample_subject_idx"))
        k_a, k_b = map(operator.itemgetter, map(int, (k_a, k_b)))
        yield {
            "name": "join_hospital_metadata",
            "file_dep": [a, b],
            "targets": [t],
            "actions": [
                lambda *args, **kwargs: join_metadata.main(a, b, k_a, k_b, t)
            ]
        }

        deps_16s = (hospital_joined, self.broad_16s_csv[0])
        self.six_sample_metadata = [os.path.join(self.products_dir,
                                                 "16s_joined.csv")]
        other_idx = self.options[
            'hospital_join']['sample_stool_idx'] + hospital_sample_width + 1

        yield {
            "name": "16s_join",
            "file_dep": deps_16s,
            "targets": self.six_sample_metadata,
            "actions": [
                lambda *ar, **kw: join_metadata.join_16s(
                    deps_16s[1], deps_16s[0], self.six_sample_metadata[0],
                    six_idx=self.options['16s_join']['sample_idx'],
                    other_idx=other_idx)
            ]
        }
        
        deps_wgs = (hospital_joined, self.broad_wgs_csv[0])
        self.wgs_sample_metadata = [os.path.join(self.products_dir,
                                                 "wgs_joined.csv")]
        yield {
            "name": "wgs_join",
            "file_dep": deps_wgs,
            "targets": self.wgs_sample_metadata,
            "actions": [
                lambda *ar, **kw: join_metadata.join_wgs(
                    deps_wgs[1], deps_wgs[0], self.wgs_sample_metadata[0],
                    six_idx=self.options['wgs_join']['sample_idx'],
                    other_idx=other_idx)
            ]
        }

        
        deps_wts = (hospital_joined, self.broad_wts_csv[0])
        self.wts_sample_metadata = [os.path.join(self.products_dir,
                                                 "wts_joined.csv")]
        yield {
            "name": "wts_join",
            "file_dep": deps_wts,
            "targets": self.wts_sample_metadata,
            "actions": [
                lambda *ar, **kw: join_metadata.join_wgs(
                    deps_wts[1], deps_wts[0], self.wts_sample_metadata[0],
                    six_idx=self.options['wts_join']['sample_idx'],
                    other_idx=other_idx)
            ]
        }

        
