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
from hmp2_workflows.tasks.metadata import generate_sample_metadata
from hmp2_workflows.utils.misc import (parse_cfg_file, 
                                       create_project_dirs,
                                       #get_template,
                                       create_merged_md5sum_file)


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

    contaminate_db = conf.get('databases').get('knead_dna')

    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')

    if data_files and data_files.get('mgx'):
        input_files = data_files.get('mgx').get('input')

        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                            project,
                                            creation_date,
                                            'WGS')
        (deposition_dir, processing_dir, public_dir) = project_dirs
        base_depo_dir = os.path.abspath(os.path.join(deposition_dir, '..'))

        ## Our md5sums are going to be in the directories with their 
        ## respective files so we'll want to collate them somehow.
        # input_dirs = [group[0] for group in itertools.groupby(input_files, 
        #                                                      os.path.dirname)]
        #md5sum_files = chain.from_iterable([find_files(input_dir, '.md5') for 
        #                                    input_dir in input_dirs])
        #merged_md5sum_file = tempfile.mktemp('.md5')
        #checksums_file = create_merged_md5sum_file(md5sum_files, 
        #                                           merged_md5sum_file)

        ## Validate files to ensure file integrity via md5 checksums
        #validated_files = verify_files(workflow, input_files, 
        #                              checksums_file)

        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_depo_dir)

        ## WGS files are generated by the Broad and will live under 
        ## the /seq/picard_aggregation directory so we'll want to symlink
        ## them over to our data deposition directory.
        deposited_files = stage_files(workflow,
                                      input_files,
                                      deposition_dir,
                                      symlink=True)

        fastq_files = bam_to_fastq(workflow, deposited_files, 
                                   processing_dir, args.threads)
        fastq_files = deposited_files

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
        tax_profile_outputs = taxonomic_profile(workflow,
                                                cleaned_fastqs,
                                                processing_dir,
                                                args.threads)

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

        ## BIOM files are also distributed on the IBDMDB website so we need 
        ## to make sure we round those all up and move them to the public 
        ## directory
        biom_files = batch_convert_tsv_to_biom(workflow, tax_profile_outputs[1])
        tax_biom_files = stage_files(workflow,
                                     biom_files,
                                     processing_dir)

        ## A bunch of analysis files are going to be in the processing
        ## directory so let's take care of moving the ones we want available 
        ## on the IBDMDB website to the 'public' directory
        pub_raw_dir = os.path.join(public_dir, 'raw')
        create_folders(pub_raw_dir)

        #pub_metadata_dir = os.path.join(public_dir, 'metadata')
        #create_folders(pub_metadata_dir)

        pub_tax_profile_dir = os.path.join(public_dir, 'tax_profile')
        create_folders(pub_tax_profile_dir)

        pub_func_profile_dir = os.path.join(public_dir, 'func_profile')
        create_folders(pub_func_profile_dir)

        pub_cleaned_fastqs = stage_files(workflow,
                                         cleaned_fastqs,
                                         pub_raw_dir)
        pub_merged_tax_profile = stage_files(workflow,
                                             [tax_profile_outputs[0]],
                                             pub_tax_profile_dir)
        pub_tax_profiles = stage_files(workflow,
                                       tax_profile_outputs[1],
                                       pub_tax_profile_dir)
        pub_tax_biom = stage_files(workflow,
                                   tax_biom_files,
                                   pub_tax_profile_dir)
        pub_func_files = stage_files(workflow, 
                                     func_profile_outputs,
                                     pub_func_profile_dir)

        
        sample_names = get_sample_names(input_files)
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

        kneaddata_log_files = name_files(sample_names,
                                         processing_dir,
                                         subfolder = 'kneaddata',
                                         extension = 'log')
        pub_log_files = stage_files(workflow,
                                    kneaddata_log_files,
                                    pub_raw_dir)

        pub_metadata_files = generate_sample_metadata(workflow,
                                                      'metagenomics',
                                                      sample_names,
                                                      args.metadata_file,
                                                      public_dir)

        ## Handle generating all the templates we need here
        #template_files = name_files(sample_names,
        #                            public_dir,
        #                            subfolder='templates',
        #                            extension='html')
        
        #for elts in zip(sample_names, template_file, pub_metadata_files,
        #                pub_tax_profiles, norm_genefamiles, norm_path_files):
        #    sample_name = elts[0]
        #    template_file = elts[1]
        #    metadata_file = elts[2]
        #    tax_profile = elts[3]
        #    features = elts[4]
        #    path_abund = elts[5]

        #    templates = [get_template('header'), 
        #                 get_template('quality_control_single_dna'),
        #                 get_template('taxonomy'),
        #                 get_template('functional_dna')]

        #    doc_task = workflow.add_document(
        #        templates=templates,
        #        depends=[read_counts, metadata_file, tax_profile, features, path_abund],
        #        targets=[template_file],
        #        vars={"sample_name": sample_name,
        #              "dna_read_counts": sample_read_counts,
        #              "taxonomic_profile": tax_profile,
        #              "path_abund": path_abund,
        #              "feature_counts": features,
        #              "contaminate_database": contaminate_db,
        #        }
        #    )

        make_files_web_visible(workflow, 
                               [pub_cleaned_fastqs,
                                pub_log_files,
                                pub_metadata_files,
                                pub_merged_tax_profile,
                                pub_tax_profiles,
                                pub_tax_biom,
                                pub_func_files,
                                func_tar_files])

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
