# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.metagenome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles HMP2 WGS metagenome data and metadata files.

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import itertools
import os
import tempfile

from itertools import chain

from anadama2 import Workflow

from biobakery_workflows.tasks.shotgun import (quality_control, 
                                               taxonomic_profile,
                                               functional_profile)

from biobakery_workflows.utilities import (find_files,
                                           sample_names as get_sample_names,
                                           create_folders,
                                           name_files)
from hmp2_workflows.tasks.common import (verify_files, 
                                         stage_files,
                                         tar_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.file_conv import (bam_to_fastq,
                                            batch_convert_tsv_to_biom)
from hmp2_workflows.tasks.metadata import (generate_sample_metadata, 
                                           add_metadata_to_tsv)
from hmp2_workflows.utils.misc import (parse_cfg_file, 
                                       create_merged_md5sum_file)
from hmp2_workflows.utils.files import create_project_dirs
                                      

def parse_cli_arguments():
    """Parses any command-line arguments passed into the workflow.

    Args:
        None
    Requires:
        None
    Returns:
        anadama2.Workflow: The workflow object for this pipeline
        anadama2.cli.Configuration: Arguments passed into this workflow.
    """
    workflow = Workflow(version='1.0', description='A workflow to handle HMP2 '
                        'WGS data.', remove_options=['input', 'output'])
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('metadata-file', desc='Accompanying metadata file '
                           'for the provided data files.', default=None)
    workflow.add_argument('threads', desc='number of threads/cores for each '
                          'task to use', default=1)

    return workflow


def main(workflow):
    args = workflow.parse_args()

    conf = parse_cfg_file(args.config_file, section='mgx')
    manifest = parse_cfg_file(args.manifest_file)

    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')
    contaminate_db = conf.get('databases').get('knead_dna')

    if data_files and data_files.get('mgx'):
        input_files = data_files.get('mgx').get('input')
        sample_names = get_sample_names(input_files)

        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                            project,
                                            creation_date,
                                            'WGS')
        (deposition_dir, processing_dir, public_dir) = project_dirs
        base_depo_dir = os.path.abspath(os.path.join(deposition_dir, '..'))

        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_depo_dir)

        deposited_files = stage_files(workflow,
                                      input_files,
                                      deposition_dir,
                                      symlink=True)

        fastq_files = bam_to_fastq(workflow, deposited_files, 
                                   processing_dir, args.threads)

        (cleaned_fastqs, read_counts) = quality_control(workflow, 
                                                        fastq_files, 
                                                        processing_dir,
                                                        args.threads,
                                                        contaminate_db)
        
        ## Generate taxonomic profile output. Output are stored in a list 
        ## and are the following:   
        ##
        ##      * Merged taxonomic profile 
        ##      * Individual taxonomic files
        ##      * metaphlan2 SAM files 
        #tax_profile_outputs = taxonomic_profile(workflow,
        #                                        cleaned_fastqs,
        #                                        processing_dir,
        #                                        args.threads)
        tax_profile_outputs = ('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/taxonomic_profiles.tsv', [], [])
        tax_profile_outputs[1].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDK_taxonomic_profile.tsv')
        tax_profile_outputs[1].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDM_taxonomic_profile.tsv')
        tax_profile_outputs[1].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDI_taxonomic_profile.tsv')
        tax_profile_outputs[1].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDU_taxonomic_profile.tsv')
        tax_profile_outputs[1].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDQ_taxonomic_profile.tsv')
    
        tax_profile_outputs[2].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDM_bowtie2.sam')
        tax_profile_outputs[2].append('/seq/ibdmdb/carze_test/processing/HMP2/2017-05-15/WGS/metaphlan2/MSM5LLDM_bowtie2.sam')

        ## Generate functional profile output using humann2. Outputs are the 
        ## the following:
        ## 
        ##      * Merged gene families
        ##      * Merged relative abundances
        ##      * Merged pathway abundances
        func_profile_outputs = functional_profile(workflow,
                                                  cleaned_fastqs,
                                                  processing_dir,
                                                  args.threads,
                                                  tax_profile_outputs[1])

        biom_files = batch_convert_tsv_to_biom(workflow, tax_profile_outputs[1])
        tax_biom_files = stage_files(workflow,
                                     biom_files,
                                     processing_dir)

        kneaddata_log_files = name_files(sample_names,
                                         processing_dir,
                                         subfolder = 'kneaddata',
                                         extension = 'log')

        pub_raw_dir = os.path.join(public_dir, 'raw')
        pub_tax_profile_dir = os.path.join(public_dir, 'tax_profile')
        pub_func_profile_dir = os.path.join(public_dir, 'func_profile')
        map(create_folders, [pub_raw_dir, pub_tax_profile_dir, 
                             pub_func_profile_dir])
    
        tax_profile_pcl = add_metadata_to_tsv(workflow,
                                              [tax_profile_outputs[0]],
                                              args.metadata_file,
                                              conf.get('metadata_id_col'),
                                              ['_taxonomic_profile', '_functional_profile'],
                                              conf.get('target_metadata_cols'))

        pub_files = [stage_files(workflow, files, target_dir) for (files, target_dir) 
                     in [(cleaned_fastqs, pub_raw_dir), 
                         ([tax_profile_outputs[0]], pub_tax_profile_dir),
                         (tax_biom_files, pub_tax_profile_dir),
                         (tax_profile_pcl, pub_tax_profile_dir),
                         (func_profile_outputs, pub_func_profile_dir),
                         (kneaddata_log_files, pub_raw_dir)]]

        pub_metadata_files = generate_sample_metadata(workflow,
                                                      'metagenomics',
                                                      sample_names,
                                                      args.metadata_file,
                                                      public_dir)

        norm_genefamilies = name_files(sample_names, 
                                       processing_dir, 
                                       subfolder = 'genes',
                                       tag = 'genefamilies_relab',
                                       extension = 'tsv')
        norm_ecs_files = name_files(sample_names,
                                    processing_dir,
                                    subfolder = 'ecs',
                                    tag = 'genefamilies_ecs_relab',
                                    extension = 'tsv')
        norm_path_files = name_files(sample_names,
                                     processing_dir,
                                     subfolder = 'pathways',
                                     tag = 'pathabundance_relab',
                                     extension = 'tsv')

        func_tar_files = []
        for (sample, gene_file, ecs_file, path_file) in zip(sample_names,
                                                            norm_genefamilies,
                                                            norm_ecs_files,
                                                            norm_path_files):
            tar_path = os.path.join(pub_func_profile_dir, 
                                    "%s_humann2.tgz" % sample)
            func_tar_file = tar_files(workflow, 
                                      [gene_file, ecs_file, path_file],
                                      tar_path)
            func_tar_files.append(func_tar_file)


        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
