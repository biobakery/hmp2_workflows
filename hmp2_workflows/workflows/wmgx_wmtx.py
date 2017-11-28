# -*- coding: utf-8 -*-

"""
hmp2_workflows.wmgx_wmtx
~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles analysis and dissemination of HMP2 
metagenome and metatranscriptome data.

Copyright (c) 2016 Harvard School of Public Health

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

import os

from anadama2 import Workflow

from biobakery_workflows.utilities import (create_folders,
                                           sample_names,
                                           paired_files,
                                           name_files)
from biobakery_workflows.tasks.shotgun import (quality_control, 
                                               taxonomic_profile,
                                               functional_profile)

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                         tar_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.file_conv import bam_to_fastq
from hmp2_workflows.utils.files import (find_files, match_tax_profiles, 
                                        create_project_dirs)
from hmp2_workflows.utils.misc import parse_cfg_file
from hmp2_workflows.tasks.metadata import add_metadata_to_tsv


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
    workflow = Workflow(version='0.1', description='A workflow to handle '
                        'analysis of metatranscriptomic data.',
                        remove_options=['input', 'output'])
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

    conf_mtx = parse_cfg_file(args.config_file, section='MTX')
    conf_mgx = parse_cfg_file(args.config_file, section='MGX')
    manifest = parse_cfg_file(args.manifest_file)

    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')

    contaminate_db = conf_mtx.get('databases').get('knead_dna')
    mtx_db = conf_mtx.get('databases').get('knead_mtx')
    rrna_db = conf_mtx.get('databases').get('knead_rrna')
    adapter_sequences = conf_mtx.get('adapter_sequences')

    qc_threads = args.threads_kneaddata if args.threads_kneaddata else args.threads
    tax_threads = args.threads_metaphlan if args.threads_metaphlan else args.threads
    func_threads = args.threads_humann if args.threads_humann else args.threads

    if data_files and data_files.get('MTX', {}).get('input'):
        input_files_mtx = data_files.get('MTX').get('input')
        file_extension_mtx = data_files.get('MTX').get('file_ext')
        pair_identifier_mtx = data_files.get('MTX').get('pair_identifier')
        input_tax_profiles = []

        project_dirs_mtx = create_project_dirs([conf_mtx.get('deposition_dir'),
                                                conf_mtx.get('processing_dir'),
                                                conf_mtx.get('public_dir')],
                                               project,
                                               creation_date,
                                               'MTX')
        public_dir_mtx = project_dirs_mtx[-1]

        deposited_files_mtx = stage_files(workflow,
                                          input_files_mtx,
                                          project_dirs_mtx[0],
                                          symlink=True)

        if file_extension_mtx == ".bam":
            ## Need to sort our BAM files to be sure here...
            paired_end_seqs = bam_to_fastq(workflow, 
                                            deposited_files_mtx, 
                                            project_dirs_mtx[1],
                                            paired_end=True,
                                            compress=False,
                                            threads=args.threads)
            pair_identifier_mtx = "_R1"                                            
        else:
            paired_end_seqs = paired_files(deposited_files_mtx, file_extension_mtx, pair_identifier_mtx)   

        adapter_trim_opts = (" --trimmomatic-options \"ILLUMINACLIP:%s:2:30:10:8:TRUE "
                             "SLIDINGWINDOW:4:20 MINLEN:50\"" % adapter_sequences)
        (cleaned_fastqs_mtx, read_counts_mtx) = quality_control(workflow,
                                                                paired_end_seqs,
                                                                file_extension_mtx,
                                                                project_dirs_mtx[1],
                                                                qc_threads,
                                                                [contaminate_db,
                                                                 rrna_db,
                                                                 mtx_db],
                                                                pair_identifier=pair_identifier_mtx,
                                                                additional_options=adapter_trim_opts,
                                                                remove_intermediate_output=True)

        sample_names_mtx = sample_names(cleaned_fastqs_mtx)                                                                

        ##########################################
        #          MGX FILE PROCESSING           #
        ##########################################
        # Ideally we would be passed in a set of corresponding metagenome
        # sequence(s) to go with our metatranscriptomic files but we also
        # have two other scenarios:
        #
        #       1.) No accompanying metagenomic sequences exist; in this
        #           case we will proceed just using the metatranscriptomic
        #           data.
        #       2.) Taxonomic profiles are passed directly in in our MANIFEST
        #           file; here we remove these from our input files and
        #           prevent them from running through the kneaddata ->
        #           metaphlan2 portions of our pipeline
        if data_files.get('MGX', {}).get('input'):
            input_files_wgs = data_files.get('MGX').get('input')
            file_extension_mgx = data_files.get('MGX').get('file_ext')
            pair_identifier_mgx = data_files.get('MGX').get('pair_identifier')
            input_tax_profiles = [in_file for in_file in input_files_wgs
                                  if 'taxonomic_profile.tsv' in in_file]
            input_files_wgs = set(input_files_wgs) - set(input_tax_profiles)

            if input_files_wgs:
                sample_names_wgs = sample_names(input_files_wgs)

                project_dirs_wgs = create_project_dirs([conf_mgx.get('deposition_dir'),
                                                        conf_mgx.get('processing_dir'),
                                                        conf_mgx.get('public_dir')],
                                                    project,
                                                    creation_date,
                                                    'WGS')
                public_dir_wgs = project_dirs_wgs[-1]

                deposited_files_wgs = stage_files(workflow,
                                                  input_files_wgs,
                                                  project_dirs_wgs[0],
                                                  symlink=True)

                if file_extension_mgx == ".bam":
                    ## Need to sort our BAM files to be sure here...
                    paired_end_seqs = bam_to_fastq(workflow, 
                                                    deposited_files_mgx, 
                                                    project_dirs_mgx[1],
                                                    paired_end=True,
                                                    compress=False,
                                                    threads=args.threads)
                    pair_identifier_mgx = "_R1"                                            
                else:
                    paired_end_seqs_mgx = paired_files(deposited_files_mgx, pair_identifier_mgx)  

                (cleaned_fastqs_wgs, read_counts_wgs) = quality_control(workflow,
                                                                        paired_end_seqs_mgx,
                                                                        project_dirs_wgs[1],
                                                                        qc_threads,
                                                                        [contaminate_db,
                                                                        rrna_db],
                                                                        remove_intermediate_output=True)

                tax_profile_outputs_wgs = taxonomic_profile(workflow,
                                                            cleaned_fastqs_wgs,
                                                            project_dirs_wgs[1],
                                                            tax_threads,
                                                            '*.fastq')

                func_profile_outputs_wgs = functional_profile(workflow,
                                                            cleaned_fastqs_wgs,
                                                            project_dirs_wgs[1],
                                                            func_threads,
                                                            tax_profile_outputs_wgs[1],
                                                            remove_intermediate_output=True)
                input_tax_profiles.extend(tax_profile_outputs_wgs[1])

                pub_wgs_raw_dir = os.path.join(public_dir_wgs, 'raw')
                pub_wgs_tax_profile_dir = os.path.join(public_dir_wgs, 'tax_profile')
                pub_wgs_func_profile_dir = os.path.join(public_dir_wgs, 'func_profile')
                map(create_folders, [pub_wgs_raw_dir, pub_wgs_tax_profile_dir,
                                    pub_wgs_func_profile_dir])

                norm_genefamilies_wgs = name_files(sample_names,
                                                project_dirs_wgs[1],
                                                subfolder='genes',
                                                tag='genefamilies_relab',
                                                extension='tsv')
                norm_ecs_files_wgs = name_files(sample_names,
                                                project_dirs_wgs[1],
                                                subfolder='ecs',
                                                tag='genefamilies_ecs_relab',
                                                extension='tsv')
                norm_path_files_wgs = name_files(sample_names,
                                                project_dirs_wgs[1],
                                                subfolder='pathways',
                                                tag='pathabundance_relab',
                                                extension='tsv')

                pcl_files = add_metadata_to_tsv(workflow,
                                                [tax_profile_outputs_wgs[1]] 
                                                + func_profile_outputs_wgs,
                                                'metagenomics',
                                                conf_mgx.get('metadata_id_col'),
                                                conf_mgx.get('analysis_col_patterns'),
                                                conf_mgx.get('target_metadata_cols'))
                                      
                func_tar_files_wgs = []
                for (sample, gene_file, ecs_file, path_file) in zip(sample_names_wgs,
                                                                    norm_genefamilies_wgs,
                                                                    norm_ecs_files_wgs,
                                                                    norm_path_files_wgs):
                    tar_path = os.path.join(pub_wgs_func_profile_dir, 
                                            "%s_humann2.tgz" % sample)
                    func_tar_file = tar_files(workflow,
                                            [gene_file, ecs_file, path_file],
                                            tar_path,
                                            depends=func_profile_outputs_wgs)
                    func_tar_files_wgs.append(func_tar_file)

        ##########################################
        #          MTX FILE PROCESSING           #
        ##########################################
        # Here we want to see if we can create a set of matching cleaned
        # MTX files to corresponding MGX taxonomic profiles. If these exist
        # we want to run functional profiling wit hthe corresponding MGX
        # taxonomic profile otherwise we will run a taxonomic profiling
        # on the MTX sequences and run functional profiling with the produced
        # taxonomic profile.
        func_outs_match_mtx = []
        if input_tax_profiles:
            (matched_fqs, matched_tax_profiles) = match_tax_profiles(cleaned_fastqs_mtx,
                                                                     conf_mtx.get('metadat_id_col'),
                                                                     input_tax_profiles,
                                                                     conf_mtx.get('tax_profile_id'),
                                                                     args.metadata_file)

            func_outs_match_mtx = functional_profile(workflow,
                                                     matched_fqs,
                                                     project_dirs_mtx[1],
                                                     func_threads,
                                                     matched_tax_profiles,
                                                     remove_intermediate_output=True)

            # Reset the remaining MTX files left over here so that we can run them through
            # the metaphlan2 -> humann2 pipeline.
            cleaned_fastqs_mtx = set(cleaned_fastqs_mtx) - set(matched_fqs)

        if cleaned_fastqs_mtx:
            tax_outs_mtx = taxonomic_profile(workflow,
                                             cleaned_fastqs_mtx,
                                             project_dirs_mtx[1],
                                             tax_threads,
                                             '*.fastq')
            func_outs_mtx = functional_profile(workflow,
                                               cleaned_fastqs_mtx,
                                               project_dirs_mtx[1],
                                               func_threads,
                                               tax_outs_mtx[1],
                                               remove_intermediate_output=True)
            func_outs_mtx = list(func_outs_mtx).extend(func_outs_match_mtx)
        else:
            func_outs_mtx = func_outs_match_mtx

        pub_mtx_raw_dir = os.path.join(public_dir_mtx, 'raw')
        pub_mtx_tax_profile_dir = os.path.join(public_dir_mtx, 'tax_profile')
        pub_mtx_func_profile_dir = os.path.join(public_dir_mtx, 'func_profile')
        map(create_folders, [pub_mtx_raw_dir, pub_mtx_tax_profile_dir,
                             pub_mtx_func_profile_dir])

        norm_genefamilies_mtx = name_files(sample_names_mtx,
                                           project_dirs_mtx[1],
                                           subfolder='genes',
                                           tag='genefamilies_relab',
                                           extension='tsv')
        norm_ecs_files_mtx = name_files(sample_names_mtx,
                                        project_dirs_mtx[1],
                                        subfolder='ecs',
                                        tag='genefamilies_ecs_relab',
                                        extension='tsv')
        norm_path_files_mtx = name_files(sample_names_mtx,
                                         project_dirs_mtx[1],
                                         subfolder='pathways',
                                         tag='pathabundance_relab',
                                         extension='tsv')

        func_tar_files_mtx = []
        for (sample, gene_file, ecs_file, path_file) in zip(sample_names_mtx,
                                                            norm_genefamilies_mtx,
                                                            norm_ecs_files_mtx,
                                                            norm_path_files_mtx):
            tar_path = os.path.join(pub_mtx_func_profile_dir,
                                    "%s_humann2.tgz" % sample)
            func_tar_file = tar_files(workflow,
                                      [gene_file, ecs_file, path_file],
                                      tar_path,
                                      depends=func_outs_mtx)
            func_tar_files_mtx.append(func_tar_file)
    
        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
