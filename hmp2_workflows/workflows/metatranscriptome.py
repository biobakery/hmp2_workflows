# -*- coding: utf-8 -*-

"""
hmp2_workflows.wgs
~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles analysis and dissemination of HMP2 
metatranscriptome data.

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

import itertools
import tempfile

from anadama2 import Workflow

from biobakery_workflows import utilities as bb_utils
from biobakery_workflows.tasks.shotgun import (quality_control, 
                                               taxonomic_profile,
                                               functional_profile)

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.file_conv import bam_to_fastq
from hmp2_workflows.utils import (parse_cfg_file, 
                                  create_merged_md5checksum_file)
from hmp2_workflows.utils.files import find_files


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
    workflow.add_argument('threads', desc='number of threads/cores for each '
                          'task to use', default=1)

    return (workflow, workflow.parse_args())


def main(workflow, args):
    conf = parse_cfg_file(args.config_file, section='mtx')

    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')

    contaminate_db = conf.get('database').get('knead_dna')
    metatranscriptome_db = conf.get('database').get('knead_mtx')
    sixs_db = conf.get('database').get('knead_six')

    if data_files and data_files.get('MTX'):
        input_files_mtx = data_files.get('MTX')
        sample_names_mtx = sample_names(input_files_mtx)

        if 'WGS' in data_files:
            wgs_files = data_files.get('WGS')



        deposited_files = stage_data_files(workflow,
                                           input_files_mtx,
                                           conf.deposition_dir,
                                           symlink=True)


        fastq_files = bam_to_fastq(workflow, deposited_files, 
                                   conf.processing_dir, args.threads)



        (cleaned_fastqs, read_counts) = quality_control(workflow, 
                                                        fastq_files, 
                                                        conf.processing_dir,
                                                        args.threads,
                                                        [contaminate_db, 
                                                         sixs_db,
                                                         metatranscriptome_db])


        
        ## Generate taxonomic profile output. Output are stored in a list 
        ## and are the following:   
        ##
        ##      * Merged taxonomic profile 
        ##      * Individual taxonomic files
        ##      * metaphlan2 SAM files 
        tax_profile_outputs = taxonomic_profile(workflow,
                                                cleaned_fastqs,
                                                conf.public_dir,
                                                args.threads)
        
        ## Generate functional profile output using humann2. Outputs are the 
        ## the following:
        ## 
        ##      * Merged gene families
        ##      * Merged relative abundances
        ##      * Merged pathway abundances
        ##
        ## TODO: We probably want to merge WGS and MTX analysis together here 
        ## because we use the WGS taxonomic profiles.
        func_profile_outputs = functional_profile(workflow,
                                                  cleaned_fastqs,
                                                  conf.public_dir,
                                                  args.threads,
                                                  tax_profile_outputs[1])

        ## We get the merged files back from the functional profiling step
        ## but we also need each individual file.
        ## TODO: Figure out the best way to gather up all the individual 
        ## files so we can push them to the website
        ##
        ## Actually we may not have to find each individual file because 
        ## they should still exist in the intermediate dirs somewhere.

        ## TODO: Add static webpage visualization generation in here
        make_files_web_visible(workflow, cleaned_fastqs,
                               tax_profile_outputs[0:1],
                               func_profile_outputs)

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
