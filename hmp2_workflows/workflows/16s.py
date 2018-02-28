# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.16s
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles analysis and dissemination of HMP2 
16S amplicon data.

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
import os
import tempfile

from itertools import chain

from anadama2 import Workflow

from biobakery_workflows.tasks.sixteen_s import (quality_control,
                                                 demultiplex,
                                                 taxonomic_profile,
                                                 merge_samples_and_rename,
                                                 functional_profile)
from biobakery_workflows.utilities import (find_files,
                                           create_folders,
                                           sample_names as get_sample_names)

from hmp2_workflows.tasks.common import (stage_files,
                                         make_files_web_visible)
from hmp2_workflows.utils.misc import (parse_cfg_file, 
                                       create_merged_md5sum_file)
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
                        'analysis of 16S amplicon data.', 
                        remove_options=['input', 'output'])
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
    conf = parse_cfg_file(args.config_file, section='16S')

    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')
    project = manifest.get('project')
    creation_date = manifest.get('submission_date')

    gg_tax = conf.get('databases').get('gg_taxonomy')
    gg_usearch = conf.get('databases').get('gg_usearch')
    gg_fasta = conf.get('databases').get('gg_fasta')

    if data_files and data_files.get('16S', {}).get('input'):
        input_files = data_files.get('16S').get('input')
        barcode_file = data_files.get('16S').get('barcode_file')
        pair_identifier = data_files.get('16S').get('pair_identifier')
        index_identifier = data_files.get('16S').get('index_identifier')

        index_files = [in_file for in_file in input_files if
                       index_identifier in in_file]
        input_files = set(input_files) - set(index_files)                       

        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                           project,
                                           creation_date,
                                           '16S')

        base_depo_dir = os.path.abspath(os.path.join(project_dirs[0], '..'))
        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_depo_dir)
        sequence_files = stage_files(workflow,
                                     input_files,
                                     project_dirs[0],
                                     symlink=True)
 
        # An entry point into this pipeline is any analysis conducted by Baylor/CMMR
        # which will require a slight branching of the pipeline utilized.

        # We are making a very fraught assumption here that if only one sequence file 
        # is passed in alongside a centroid FASTA file we are dealing with any files
        # generated by CMMR
        otu_table = data_files.get('16S').get('otu_table')
        centroid_fasta = data_files.get('16S').get('centroid_fasta')
        if otu_table and centroid_fasta and len(sequence_files) < 1:
            merged_fastq = sequence_files[0]
        else:
            if barcode_file:
                sequence_files = demultiplex(workflow,
                                            input_files,
                                            project_dirs[1],
                                            barcode_file,
                                            index_files,
                                            conf.get('min_pred_qc_score'),
                                            pair_identifier)

            merged_fastq = merge_samples_and_rename(workflow,
                                                    sequence_files,
                                                    project_dirs[1],
                                                    pair_identifier,
                                                    args.threads)

            qc_fasta_outs = quality_control(workflow,
                                            merged_fastq,
                                            project_dirs[1],
                                            args.threads,
                                            conf.get('maxee'),
                                            conf.get('min_trunc_len_max'))

            if not otu_table:
                closed_ref_tsv = taxonomic_profile(workflow,
                                                   qc_fasta_outs[0],
                                                   qc_fasta_outs[1],
                                                   qc_fasta_outs[2],
                                                   project_dirs[1],
                                                   args.threads,
                                                   conf.get('percent_identity'),
                                                   gg_usearch,
                                                   gg_fasta,
                                                   gg_tax,
                                                   conf.get('min_size'))
                
                predict_metagenomes_tsv = functional_profile(workflow, 
                                                             closed_ref_tsv, 
                                                             project_dirs[1])

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
