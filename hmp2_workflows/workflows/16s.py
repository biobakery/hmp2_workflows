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
                                                 taxonomic_profile,
                                                 functional_profile)
from biobakery_workflows.utilities import find_files

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.file_conv import bam_to_fastq
from hmp2_workflows.utils import (parse_cfg_file, 
                                  create_project_dirs,
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
    workflow = Workflow(version='0.1', description='A workflow to handle '
                        'analysis of 16S amplicon data.', 
                        remove_options=['input', 'output'])
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('threads', desc='number of threads/cores for each '
                          'task to use', default=1)

    return workflow


def main(workflow):
    args = workflow.parse_args()
    conf = parse_cfg_file(args.config_file, section='16S')

    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')
    project = manifest.get('project')

    contaminate_db = conf.get('databases').get('knead_dna')
    metatranscriptome_db = conf.get('databases').get('knead_mtx')
    sixs_db = conf.get('databases').get('knead_six')
    gg_tax = conf.get('databases').get('gg_taxonomy')
    gg_usearch = conf.get('databases').get('gg_usearch')
    gg_fasta = conf.get('databases').get('gg_fasta')

    if data_files and data_files.get('16S'):
        input_files = data_files.get('16S').get('input')

        ## Create project directories needed in the following steps.
        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                           project,
                                           '16S')
        (deposition_dir, processing_dir, public_dir) = project_dirs

        base_depo_dir = os.path.abspath(os.path.join(deposition_dir, '..'))

        ## Our md5sums are going to be in the directories with their 
        ## respective files so we'll want to collate them somehow.
        input_dirs = [group[0] for group in itertools.groupby(input_files, 
                                                              os.path.dirname)]
        md5sum_files = chain.from_iterable([find_files(input_dir, '.md5') for 
                                            input_dir in input_dirs])
        merged_md5sum_file = tempfile.mktemp('.md5')
        checksums_file = create_merged_md5sum_file(md5sum_files, 
                                                   merged_md5sum_file)

        ## Validate files to ensure file integrity via md5 checksums
        validated_files = verify_files(workflow, input_files, 
                                       checksums_file)

        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_depo_dir)

        ## Files are generated by the Broad and will live under 
        ## the /seq/picard_aggregation directory so we'll want to symlink
        ## them over to our data deposition directory.
        deposited_files = stage_files(workflow,
                                      validated_files,
                                      deposition_dir,
                                      symlink=True)

        fastq_files = bam_to_fastq(workflow, deposited_files, 
                                   processing_dir)
        
        qc_fasta_outs = quality_control(workflow,
                                        fastq_files,
                                        processing_dir,
                                        args.threads,
                                        1,
                                        150)

        closed_ref_tsv = taxonomic_profile(workflow,
                                           qc_fasta_outs[0],
                                           qc_fasta_outs[1],
                                           qc_fasta_outs[2],
                                           processing_dir,
                                           args.threads,
                                           '0.97',
                                           gg_usearch,
                                           gg_fasta,
                                           gg_tax,
                                           2)
        
        predict_metagenomes_tsv = functional_profile(workflow, 
                                                     closed_ref_tsv, 
                                                     processing_dir)

        ## TODO: We are making public directories and not individual files 
        ## here so I need to be really careful not to make public files that 
        ## we don't want public
        make_files_web_visible(workflow, 
                               [qc_fasta_outs,
                                [closed_ref_tsv],
                                [predict_metagenomes_tsv]])

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
