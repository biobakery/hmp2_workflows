# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.metabolomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles analysis and dissemination of HMP2 
metabolomics data.

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

from anadama2 import Workflow

from biobakery_workflows.utilities import (find_files,
                                           sample_names as get_sample_names,
                                           create_folders,
                                           name_files)
from hmp2_workflows.tasks.common import (stage_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.metadata import add_metadata_to_tsv
from hmp2_workflows.tasks.file_conv import (excel_to_csv)                                           
from hmp2_workflows.utils.misc import parse_cfg_file
from hmp2_workflows.utils.files import create_project_dirs


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
                        'Metabolomic data.', remove_options=['input', 'output'])
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('metadata-file', desc='Accompanying metadata file '
                           'for the provided data files.', default=None)
    workflow.add_argument('aux_metadata', desc='Any additional metadata '
                          'files that can supply metadata for our ouptut '
                          'PCL files.')                           

    return workflow


def main(workflow):
    args = workflow.parse_args()

    conf = parse_cfg_file(args.config_file, section='mbx')
    manifest = parse_cfg_file(args.manifest_file)

    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')
    dataset_cfg = manifest.get('config')

    if data_files and data_files.get('mbx'):
        input_files = data_files.get('mbx').get('input')
        sample_names = get_sample_names(input_files)

        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                            project,
                                            creation_date,
                                            'Metabolomics')

        (deposition_dir, processing_dir, public_dir) = project_dirs
        base_depo_dir = os.path.abspath(os.path.join(deposition_dir, '..'))

        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_depo_dir)

        deposited_files = stage_files(workflow,
                                      input_files,
                                      deposition_dir)

        # Our metabolite data are just a series of spreadsheets that we are 
        # going to want to append some metadata too. If they are Excel files 
        # we will want to process them and convert to CSV.
        processed_files = excel_to_csv(workflow, 
                                       deposited_files,
                                       processing_dir)

        pcl_files = add_metadata_to_tsv(workflow,
                                        processed_files,
                                        args.metadata_file,
                                        'metabolomics',
                                        conf.get('metadata_id_col'),
                                        metadata_rows=dataset_cfg.get('metadata_rows'),
                                        col_offset=dataset_cfg.get('col_offset'),
                                        target_cols=conf.get('target_metadata_cols', None))

        public_files = stage_files(workflow, pcl_files, public_dir)

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())