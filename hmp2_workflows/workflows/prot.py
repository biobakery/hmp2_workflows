# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.prot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

import datetime
import os

from itertools import chain

from anadama2 import Workflow

from biobakery_workflows.utilities import find_files, create_folders

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                        make_files_web_visible)
from hmp2_workflows.tasks.metadata import add_metadata_to_tsv
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
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('checksums-file', desc='MD5 checksums for files '
                          'found in the supplied input directory.')
    workflow.add_argument('data_specific_metadata', desc='A collection of '
                          'dataset specific metadata that should be integrated '
                          'with any analysis output (creating a PCL file).')

    return workflow


def main(workflow):
    args = workflow.parse_args()
    conf = parse_cfg_file(args.config_file, section='proteomics')

    ## Parse the manifest file containing all data files from this submission
    manifest = parse_cfg_file(args.manifest_file)
    project = manifest.get('project')
    data_files = manifest.get('submitted_files')

    if data_files and data_files.get('proteomics'):
        (input_files, output_files) = data_files.get('proteomics').values()

        ## Step #1 - Verify MD5sums of all input data provided to IBDMDB
        ## 
        ## Since our proteomics files will be coming from the PNNL our
        ## files won't be in the same location as the Broad files so 
        ## we'll need to get MD5's manually supplied.
        validated_files = verify_files(workflow, input_files, 
                                       args.checksums_file)

        ## Setup the directories where we will be depositing our files
        date_stamp = str(datetime.date.today())
        base_deposition_dir = os.path.join(conf.get('deposition_dir'),
                                           project,
                                           date_stamp)
        deposition_dir = os.path.join(base_deposition_dir, 'proteomics')
        create_folders(deposition_dir)

        processing_dir = os.path.join(conf.get('processing_dir'),
                                      project,
                                      date_stamp,
                                      'proteomics')
        create_folders(processing_dir)

        public_dir = os.path.join(conf.get('public_dir'),
                                  project,
                                  date_stamp,
                                  'proteomics')
        create_folders(public_dir)

        ## Move the manifest file over so we have information about this 
        ## batch of data in the deposition directory
        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_deposition_dir)

        ## Step 2 - Move files over to our deposition directory
        deposited_files = stage_files(workflow,
                                      validated_files,
                                      deposition_dir)

        ## Step #3 - Stage files to processing directory
        ##
        ## For the Proteomics data it is ok to symlink these files over from the 
        ## data deposition folder because these files aren't actually processed
        ## but we need them to be in place here to show up on the website.
        files_to_process = stage_files(workflow, 
                                       deposited_files,
                                       processing_dir,
                                       symlink=True)

        output_files = output_files if output_files else []         

        ## We have a dataset specific metadata file that we can incorporate
        ## into the analysis output.
        if output_files and args.data_specific_metadata:
            output_files = add_metadata_to_tsv(workflow, 
                                               output_files, 
                                               args.data_specific_metadata, 
                                               conf.get('metadata_id_col'),
                                               conf.get('target_metadata_cols', []))

        ## Step #4 - Stage output files to public folder
        public_files = stage_files(workflow, output_files, 
                                   public_dir)
        
        ## TODO: We need to generate metadata files for the output files that
        ## are included with this dataset. Need to talk to George about 
        ## getting the ID-mapped version of these files since they will be 
        ## needed here.

        ## Step #5 - Make files web-visible by creating the complete.html file
        ## in each of our output directories.
        make_files_web_visible(workflow, [files_to_process, public_files])

        ## Step #6 - Once all the files have been staged we can go ahead and
        ## delete the raw files from their original directory as well as the 
        ## MANIFEST file.

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
