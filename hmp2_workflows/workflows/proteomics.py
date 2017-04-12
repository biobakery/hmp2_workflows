# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.proteomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles HMP2 Proteomics data and metadata files.

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

from itertools import chain

from anadama2 import Workflow

from biobakery_workflows.utilities import find_files

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                        make_files_web_visible)
from hmp2_workflows.utils import parse_cfg_file


def parse_cli_arguments():
    """Parses any command-line arguments passed into the workflow.
    
    Args: 
        None
    Requires:
        None
    Returns:
        AnaDAMA2.Workflow: The workflow object for this pipeline.
        AnaDAMA2.cli.Configuration: Arguments passed into this workflow.
    """
    workflow = Workflow(version='0.1', description='A workflow to handle HMP2 '
                        'Proteomics data.', remove_options=['input', 'output'])
    workflow.add_argument('--manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('--config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('--md5-checksums', desc='MD5 checksums for files '
                          'found in the supplied input directory.')

    return (workflow, workflow.parse_args())


def main(workflow, args):
    conf = parse_cfg_file(args.config_file, section='proteomics')

    ## Parse the manifest file containing all data files from this submission
    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')
    
    if data_files and data_files.get('proteomics'):
        input_files = data_files.get('proteomics')

        ## Step #1 - Verify MD5sums of all input data provided to IBDMDB
        validated_files = verify_files(workflow, input_files, 
                                       args.md5_checksums)

        ## Step 2 - Move files over to our deposition directory
        deposited_files = stage_files(workflow,
                                      validated_files,
                                      conf.deposition_dir,
                                      delete=True)
        
        ## Step #3 - Stage files to processing directory
        ##
        ## For the Proteomics data it is ok to symlink these files over from the 
        ## data deposition folder because these files aren't actually processed
        ## but we need them to be in place here to show up on the website.
        files_to_process = stage_files(workflow, 
                                       deposited_files,
                                       conf.processing_dir,
                                       symlink=True)

        output_files = (data_files.get('output_files') if 
            'output_files' in data_files else [])
         
        ## Step #4 - Stage output files to public folder
        public_files = stage_files(workflow, output_files, 
                                   conf.public_dir)
        
        ## Step #5 - Make files web-visible by creating the complete.html file
        ## in each of our output directories.
        make_file_web_visible(workflow, files_to_process, files_to_process, 
                              public_files)


if __name__ == "__main__":
    main(parse_cli_arguments())
