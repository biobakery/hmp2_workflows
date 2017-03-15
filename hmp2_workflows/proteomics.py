# -*- coding: utf-8 -*-

"""
proteomics.py
~~~~~~~~~~~~~

An AnADAMA2 workflow that handles HMP2 Proteomics data and metadata files.

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

from anadama2 import Workflow

from hmp2_workflows.utils import verify_files, stage_files



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
                        'Proteomics data.')

    workflow.add_argument('--config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('--input-directory', desc='Input directory '
                          'containing Proteomics files to be processed.')
    workflow.add_argument('--md5-checksums', desc='MD5 checksums for files '
                          'found in the supplied input directory.')
    workflow.add_argument('-m' '--metadata-file', desc='Metadata file '
                          'assosciated with provided Proteomics files')

    return (workflow, workflow.parse_args())


def main(workflow, args):
   valid_files = verify_files(workflow, input_directory, 


aa __name__ == "__main__":
    main(parse_cli_arguments())
