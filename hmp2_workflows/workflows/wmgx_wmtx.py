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

from biobakery_workflows import utilities as bb_utils
from biobakery_workflows.tasks.shotgun import (quality_control, 
                                               taxonomic_profile,
                                               functional_profile)

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.file_conv import bam_to_fastq
from hmp2_workflows.utils.files import find_files
from hmp2_workflows.utils import (parse_cfg_file, 
                                  find_files,
                                  create_merged_md5checksum_file)
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
    workflow = Workflow(version='0.1', description='A workflow to handle '
                        'analysis of metatranscriptomic data.', 
                        remove_options=['input', 'output'])
    workflow.add_argument('--manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('--config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('metadata-file', desc='Accompanying metadata file '
                           'for the provided data files.', default=None)
    workflow.add_argument('match_wgs_out', desc='Attempt to search our output '
                          'space for matching WGS output data that corresponds '
                          'to our MTX files.')                           
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


def create_seq_mapping_file(mtx_seqs, wgs_seqs, output_dir):
    """Creates a mapping file mapping metatranscriptome samples to 
    metagenome samples for use in biobakery-workflows tasks.

    Mapping file is in format:

        # wts   wms
        wts_1   wms_1
        wts_2   wms_2
 
    Args:
        mtx_seqs (list):  A list of metatranscriptome sample files
        wgs_seqs (list): A list of metagenome sample files
        output_dir (string): The directory to write the sequence mapping
            file too.

    Requires:
        None

    Returns:
        string: The path to the sequence mapping file.       
    """
    mapping_file = os.path.join(output_dir, "mapping_file.tsv")
    
    mapping_fh = open(mapping_file, 'w')
    mapping_fh.write('# wts\t#wms\n')

    mtx_seq_dir = os.path.dirname(mtx_seqs[0])
    wgs_seq_dir = os.path.dirname(wgs_seqs[0])

    mtx_seqs_base = [os.path.splitext(os.path.basename(mtx_seq))[0] 
                     for mtx_seq in mtx_seqs]
    wgs_seqs_base = [os.path.splitext(os.path.basename(wgs_seq))[0] 
                     for wgs_seq in wgs_seqs]
    
    matching_seqs = set(mtx_seqs_base).intersection(set(wgs_seqs_base))
    
    for seq_base in matching_seqs:
        mapping_fh.write('%s.fastq\t%s.fastq\n' % (os.path.join(mtx_seq_dir, seq_base),
                                                   os.path.join(wgs_seq_dir, seq_base)))
 
    mapping_fh.close()

    return mapping_file


def main(workflow):
    args = workflow.parse_args()    

    conf = parse_cfg_file(args.config_file, section='mtx')
    manifest = parse_cfg_file(args.manifest_file)

    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')

    contaminate_db = conf.get('database').get('knead_dna')
    mtx_db = conf.get('database').get('knead_mtx')
    rrna_db = conf.get('database').get('knead_rrna')

    if data_files and data_files.get('MTX'):
        input_files_mtx = data_files.get('MTX')
        sample_names_mtx = bb_utils.sample_names(input_files_mtx)

        project_dirs_mtx = create_project_dirs([conf.get('deposition_dir'),
                                                conf.get('processing_dir'),
                                                conf.get('public_dir')],
                                                project,
                                                creation_date,
                                                'MTX')

        deposited_files_mtx = stage_files(workflow,
                                          input_files_mtx,
                                          project_dirs_mtx[0],
                                          symlink=True)

        (cleaned_fastqs_mtx, read_counts_mtx) = quality_control(workflow, 
                                                                input_files_mtx, 
                                                                project_dirs_mtx[1],
                                                                args.threads,
                                                                [contaminate_db, 
                                                                 rrna_db,
                                                                 mtx_db])

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
        if data_files.get('WGS'):
            input_files_wgs = data_files.get('WGS')
            input_tax_profiles = [in_file for in_file in input_files_wgs
                                  if 'taxonomic_profile.tsv' in in_file]
            input_files_wgs = set(input_files_wgs) - set(input_tax_profiles)

            sample_names_wgs = bb_utils.sample_names(input_files_wgs)

            project_dirs_wgs = create_project_dirs([conf.get('deposition_dir'),
                                                    conf.get('processing_dir'),
                                                    conf.get('public_dir')],
                                                    project,
                                                    creation_date,
                                                    'WGS')

            deposited_files_wgs = stage_files(workflow,
                                            input_files_wgs,
                                            project_dirs_wgs[0],
                                            symlink=True)

            (cleaned_fastqs_wgs, read_counts_wgs) = quality_control(workflow, 
                                                                    input_files_wgs,
                                                                    project_dirs_wgs[1],
                                                                    args.threads,
                                                                    [contaminate_db,
                                                                    rrna_db])

            tax_profile_outputs_wgs = taxonomic_profile(workflow,
                                                        cleaned_fastqs_wgs,
                                                        project_dirs_wgs[1],
                                                        args.threads)
            tax_profile_outputs_wgs[1].extend(input_tax_profiles)                                                        

            func_profile_outputs_wgs = functional_profile(workflow,
                                                        cleaned_fastqs_wgs,
                                                        project_dirs_wgs[1],
                                                        args.threads,
                                                        tax_profile_outputs_wgs[1])

        # Here we want to see if we can create a set of matching cleaned 
        # MTX files to corresponding MGX taxonomic profiles. If these exist 
        # we want to run functional profiling wit hthe corresponding MGX 
        # taxonomic profile otherwise we will run a taxonomic profiling 
        # on the MTX sequences and run functional profiling with the produced 
        # taxonomic profile.
        if tax_profile_outputs:
            (matched_mtx_fq, matched_tax_profile) = match_files(cleaned_fastqs_mtx,
                                                                tax_profile_outputs_wgs[1],
                                                                args.metadata_file)

            func_profile_output_mtx = functional_profile(workflow,
                                                         matched_mtx_fq,
                                                         project_dirs_mtx[1],
                                                         args.threads,
                                                         matched_tax_profile)

            no_match_mtx_fq = set(cleaned_fastqs_mtx) - set(matched_mtx_fq)                                                            
            if no_match_mtx_fq:
                tax_outs_nomatch_mtx = taxonomic_profile(workflow,
                                                         no_match_mtx_fq,
                                                         project_dirs_mtx[1],
                                                         args.threads)
                func_outs_nomatch_mtx = functional_profile(workflow,
                                                           no_match_mtx_fq,
                                                           project_dirs_mtx[1],
                                                           args.threads,
                                                           tax_outs_nomatch_mtx[1])
        else:
            tax_outs_mtx = taxonomic_profile(workflow,
                                             cleaned_fastqs_mtx,
                                             project_dirs_mtx[1],
                                             args.threads)
            func_profile_output_mtx = functional_profile(workflow,
                                                         cleaned_fastqs_mtx,
                                                         project_dirs_mtx[1],
                                                         args.threads,
                                                         tax_outs_mtx[1])

        ## TODO: Handle files that do not have a matching WGS taxonomic profile
        norm_ratio_outputs = norm_ratio(workflow,
                                        func_profile_outputs_wgs[3],
                                        func_profile_outputs_wgs[4],
                                        func_profile_outputs_mtx[3],
                                        func_profile_outputs_mtx[4],
                                        func_profile_outputs_mtx[-1],
                                        project_dirs_mtx[-1],
                                        mapping_file)

        pub_mtx_raw_dir = os.path.join(public_dir, 'raw')
        pub_mtx_tax_profile_dir = os.path.join(public_dir, 'tax_profile')
        pub_mtx_func_profile_dir = os.path.join(public_dir, 'func_profile')
        map(create_folders, [pub_mtx_raw_dir, pub_mtx_tax_profile_dir, 
                             pub_mtx_func_profile_dir])

        pub_wgs_raw_dir = os.path.join(public_dir, 'raw')
        pub_wgs_tax_profile_dir = os.path.join(public_dir, 'tax_profile')
        pub_wgs_func_profile_dir = os.path.join(public_dir, 'func_profile')
        map(create_folders, [pub_wgs_raw_dir, pub_wgs_tax_profile_dir, 
                             pub_wgs_func_profile_dir])

        norm_genefamilies_mtx = name_files(sample_names, 
                                           project_dirs_mtx[1], 
                                           subfolder = 'genes',
                                           tag = 'genefamilies_relab',
                                           extension = 'tsv')
        norm_ecs_files_mtx = name_files(sample_names,
                                        project_dirs_mtx[1],
                                        subfolder = 'ecs',
                                        tag = 'genefamilies_ecs_relab',
                                        extension = 'tsv')
        norm_path_files_mtx = name_files(sample_names,
                                         project_dirs_mtx[1],
                                         subfolder = 'pathways',
                                         tag = 'pathabundance_relab',
                                         extension = 'tsv')

        func_tar_files_mtx = []
        for (sample, gene_file, ecs_file, path_file) in zip(sample_names_mtx,
                                                            norm_genefamilies_mtx,
                                                            norm_ecs_files_mtx,
                                                            norm_path_files_mtx):
            tar_path = os.path.join(pub_mtx_func_profile_dir, 
                                    "%s_humann2.tgz" % sample)
            func_tar_file = tar_files(workflow, 
                                      [gene_file, ecs_file, path_file],
                                      tar_path)
            func_tar_files_mtx.append(func_tar_file)

        norm_genefamilies_wgs = name_files(sample_names, 
                                           project_dirs_wgs[1], 
                                           subfolder = 'genes',
                                           tag = 'genefamilies_relab',
                                           extension = 'tsv')
        norm_ecs_files_wgs = name_files(sample_names,
                                        project_dirs_wgs[1],
                                        subfolder = 'ecs',
                                        tag = 'genefamilies_ecs_relab',
                                        extension = 'tsv')
        norm_path_files_wgs = name_files(sample_names,
                                         project_dirs_wgs[1],
                                         subfolder = 'pathways',
                                         tag = 'pathabundance_relab',
                                         extension = 'tsv')

        func_tar_files_wgs = []
        for (sample, gene_file, ecs_file, path_file) in zip(sample_names_mtx,
                                                            norm_genefamilies_mtx,
                                                            norm_ecs_files_mtx,
                                                            norm_path_files_mtx):
            tar_path = os.path.join(pub_wgs_func_profile_dir, 
                                    "%s_humann2.tgz" % sample)
            func_tar_file = tar_files(workflow, 
                                      [gene_file, ecs_file, path_file],
                                      tar_path)
            func_tar_files_wgs.append(func_tar_file)

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
