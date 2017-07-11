# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.wmgx
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    workflow.add_argument('threads-kneaddata', desc='OPTIONAL. A specific '
                          'number of threads/cores to use just for the '
                          'kneaddata task.', default=None)
    workflow.add_argument('threads-metaphlan', desc='OPTIONAL. A specific '
                          'number of threads/cores to use just for the '
                          'metaphlan2 task.', default=None)
    workflow.add_argument('threads-humann', desc='OPTIONAL. A specific '
                          'number of threads/cores to use just for the humann2 '
                          'task.', default=None)
    return workflow


def main(workflow):
    args = workflow.parse_args()

    conf = parse_cfg_file(args.config_file, section='MGX')
    manifest = parse_cfg_file(args.manifest_file)

    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')
    contaminate_db = conf.get('databases').get('knead_dna')

    if data_files and data_files.get('MGX'):
        input_files = data_files.get('MGX').get('input')
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

        qc_threads = args.threads_kneaddata if args.threads_kneaddata else args.threads
        (cleaned_fastqs, read_counts) = quality_control(workflow, 
                                                        deposited_files, 
                                                        processing_dir,
                                                        qc_threads,
                                                        contaminate_db,
                                                        remove_intermediate_output=True)
        
        ## Generate taxonomic profile output. Output are stored in a list 
        ## and are the following:   
        ##
        ##      * Merged taxonomic profile 
        ##      * Individual taxonomic files
        ##      * metaphlan2 SAM files 
        tax_threads = args.threads_metaphlan if args.threads_metaphlan else args.threads
        tax_profile_outputs = taxonomic_profile(workflow,
                                                cleaned_fastqs,
                                                processing_dir,
                                                tax_threads,
                                                input_extension="*.fastq")

        ## Generate functional profile output using humann2. Outputs are the 
        ## the following:
        ## 
        ##      * Merged gene families
        ##      * Merged relative abundances
        ##      * Merged pathway abundances
        func_threads = args.threads_humann if args.threads_humann else args.threads
        func_profile_outputs = functional_profile(workflow,
                                                  cleaned_fastqs,
                                                  processing_dir,
                                                  func_threads,
                                                  tax_profile_outputs[1],
                                                  remove_intermediate_output=True)

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
 
        knead_read_counts = os.path.join(processing_dir, 
                                         'counts', 
                                         'kneaddata_read_count_table.tsv') 
        tax_profile_pcl = add_metadata_to_tsv(workflow,
                                              [tax_profile_outputs[0]],
                                              args.metadata_file,
                                              'metagenomics',
                                              conf.get('metadata_id_col'),
                                              conf.get('analysis_col_patterns'),
                                              conf.get('target_metadata_cols'),
                                              aux_files=[knead_read_counts])
        func_profile_pcl = add_metadata_to_tsv(workflow,
                                               [func_profile_outputs[0]],
                                               args.metadata_file,
                                               'metagenomics',
                                               conf.get('metadata_id_col'),
                                               conf.get('analysis_col_patterns'),
                                               conf.get('target_metadata_cols'),
                                               aux_files=[knead_read_counts])

        pub_files = [stage_files(workflow, files, target_dir) for (files, target_dir) 
                     in [(cleaned_fastqs, pub_raw_dir), 
                         ([tax_profile_outputs[0]], pub_tax_profile_dir),
                         (tax_profile_outputs[1], pub_tax_profile_dir),
                         (tax_biom_files, pub_tax_profile_dir),
                         (tax_profile_pcl, pub_tax_profile_dir),
                         (func_profile_outputs, pub_func_profile_dir),
                         (func_profile_pcl, pub_func_profile_dir),
                         (kneaddata_log_files, pub_raw_dir)]]

        norm_genefamilies = name_files(sample_names, 
                                      processing_dir, 
                                       subfolder = 'genes',
                                       tag = 'genefamilies_relab',
                                       extension = 'tsv')
        norm_ecs_files = name_files(sample_names,
                                    processing_dir,
                                    subfolder = 'ecs',
                                    tag = 'ecs_relab',
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
                                      tar_path,
                                      depends=func_profile_outputs)
            func_tar_files.append(func_tar_file)


        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
