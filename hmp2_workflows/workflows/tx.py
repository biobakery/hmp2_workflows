# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.tx
~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles HMP2 host transcriptome data.

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

from biobakery_workflows.utilities import (find_files, create_folders,
                                           paired_files, sample_names)

from hmp2_workflows.tasks.common import verify_files, stage_files, tar_files
from hmp2_workflows.tasks.metadata import add_metadata_to_tsv
from hmp2_workflows.tasks.file_conv import bam_to_fastq

from hmp2_workflows.utils.files import create_project_dirs
from hmp2_workflows.utils.misc import parse_cfg_file


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
                        'host transcriptome data.',
                        remove_options=['input', 'output'])
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('threads', desc='Number of threads to use in '
                          'workflow processing', default=1)

    return workflow


def main(workflow):
    args = workflow.parse_args()
    conf = parse_cfg_file(args.config_file, section='TX')

    ## Parse the manifest file containing all data files from this submission
    manifest = parse_cfg_file(args.manifest_file)
    project = manifest.get('project')
    data_files = manifest.get('submitted_files')
    submission_date = manifest.get('submission_date')

    if data_files and data_files.get('TX', {}).get('input'):
        input_files = data_files.get('TX').get('input')

        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                            project,
                                            submission_date,
                                            'TX')

        deposited_files = stage_files(workflow,
                                      input_files,
                                      project_dirs[0],
                                      symlink=True)

        paired_fastq_files = bam_to_fastq(workflow,
                                          input_files,
                                          project_dirs[1],
                                          paired_end=True,
                                          threads=args.threads)

        paired_fastq_tars = []
        for (mate_1, mate_2) in paired_fastq_files:
            sample_name = sample_names(mate_1, pair_identifier="_R1")
            tar_path = os.path.join(project_dirs[-1], "%s.tar" % sample_name)
            paired_fastq_tar = tar_files(workflow, 
                                         [mate_1, mate_2],
                                         tar_path,
                                         depends=[mate_1, mate_2])
            paired_fastq_tars.append(paired_fastq_tar)



        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())                          

