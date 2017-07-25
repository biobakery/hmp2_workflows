# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles refreshing HMP2 metadata files and producing
a validated and human-readable copy of all metadata to be disseminated on the 
IBDMDB website.

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

import os
import tempfile

import pandas as pd 

from anadama2 import Workflow

from hmp2_workflows.tasks.metadata import (validate_metadata_file, 
                                           generate_metadata_file)
#                                           make_metadata_human_readable)
from hmp2_workflows.tasks.common import stage_files

from hmp2_workflows.utils.misc import (parse_cfg_file)
                                       #merge_metadata_file)


def parse_cli_arguments():
    """Parses any command-line arguments passed into the workflow.
    
    Args: 
        None
    Requires:
        None
    Returns:
        AnaDAMA2.Workflow: The workflow object for this pipeline.
    """
    workflow = Workflow(version='0.1', description='A workflow to handle '
                        'refreshing and disseminating HMP2 metadata.',
                        remove_options=['input', 'output'])
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('metadata-file', desc='If an existing metadata '
                          'file exists it can be supplied here. This metadata '
                          'file will be appended to instead of a whole new '
                          'metadata file being generated.')
    workflow.add_argument('studytrax-metadata-file', desc='Accompanying '
                          'StudyTrax data all corresponding samples in the '
                          'HMP2 project.')
    workflow.add_argument('broad-sample-tracking-file', desc='Broad Institute '
                           'sample tracking spreadsheet containing status of '
                           'sequence products generated.')
    workflow.add_argument('auxillary-metadata', action='append', default=[],
                          desc='Any auxillary metadata to be appeneded '
                          'to the final metadata table.')                           

    return workflow


def main(workflow):
    args = workflow.parse_args()

    config = parse_cfg_file(args.config_file)
    manifest = parse_cfg_file(args.manifest_file)

    if manifest.get('submitted_files'):
        data_files = manifest.get('submitted_files')

        ## Validate our metadata files
        validate_metadata_file(workflow, args.studytrax_metadata_file,
                               config.get('validators').get('studytrax'))
        validate_metadata_file(workflow, args.broad_sample_tracking_file,
                               config.get('validators').get('broad_sample_status'))
        
        metadata_file = None
        metadata_dir = os.path.join(config.get('public_dir'), 
                                    config.get('project'), 
                                    'metadata')

        manifest_metadata_df = generate_metadata_file(workflow,
                                                      config,
                                                      data_files,
                                                      args.studytrax_metadata_file,
                                                      args.broad_sample_tracking_file,
                                                      args.auxillary_metadata)

        if args.metadata_file:
            metadata_df = pd.read_csv(args.metadata_file)
            #metadata_file = merge_metadata_file(manifest_metadata_df, 
            #                                    metadata_df)

        #metadata_file = make_metadata_human_readable(metadata_file)

        #pub_metadata_file = stage_files(workflow, metadata_file,
        #                                metadata_dir)
        workflow.go()    

if __name__ == "__main__":
    main(parse_cli_arguments())
