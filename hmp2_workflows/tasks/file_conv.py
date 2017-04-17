# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.file_conv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains functions that convert between different file types.

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

from biobakery_workflows import utilities as bbutils
from biobakery_workflows.tasks.sixteen_s import convert_to_biom_from_tsv


def bam_to_fastq(workflow, input_files, output_dir, threads = 1):
    """Converts BAM sequence files to a single interleaved FASTQ file using
    the samtools bam2fq utility.

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 Workflow object to append 
            the BAM to FASTQ conversion step to.
        input_files (list): A list containing all BAM files to be converted.
        output_dir (string): The output directory to write converted files too.
        threads (int): The number of threads/cores to use for BAM -> FASTQ 
            conversion.

    Requires:
        samtools v1.12+: Utilities for the Sequence Alignment/Map (SAM) format

    Returns:
        list: A list of the newly-converted FASTQ files.

    Example:
        from anadama2 import Workflow

        from hmp2_workflows.tasks.file_conv import bam_to_fastq

 
        workflow = Workflow()
        fastq_files = bam_to_fastq(workflow,
                                   ['/tmp/fooA.bam', '/tmp/fooB.bam'],
                                   '/seq/ibdmdbd/out_dir'])
    """
    output_files = bbutils.name_files(map(os.path.basename, input_files),
                                      output_dir, 
                                      extension = 'fastq')

    workflow.add_task_group_gridable('samtools bam2fq [depends[0]] > [targets[0]]',
                                     depends = input_files,
                                     targets = output_files,
                                     time = 24*60,
                                     mem = 4096,
                                     cores = threads)

    return output_files


def batch_convert_tsv_to_biom(workflow, tsv_files): 
    """Batch converts tsv files to the biom format. BIOM files will be 
    deposited in the same folder as source TSV files and will carry the 
    same filenames.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        tsv_files (list): A list containing all TSV files to be converted 
            to BIOM format.
    
    Requires:
        Biom v2: A tool for general use formatting of biological data.

    Returns: 
        list: A list containing paths to all converted BIOM files.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import common

        workflow = anadama2.Workflow()

        tsv_files = ['/tmp/foo.tsv', '/tmp/bar.tsv', '/tmp/baz.tsv']
        biom_files = common.batch_convert_tsv_to_biom(workflow, tsv_files)

        print biom_files
        ## ['/tmp/foo.biom', '/tmp/bar.biom', '/tmp/baz.biom']
    """
    biom_files = []

    tsv_fnames = bbutils.sample_names(tsv_files)
    tsv_dir = os.path.dirname(tsv_files[0])

    biom_files = [os.path.join(tsv_dir, biom_fname) for biom_fname in 
                  bbutils.name_files(tsv_fnames, tsv_dir, extension='biom')]

    for (tsv_file, biom_file) in zip(tsv_files, biom_files):
        convert_to_biom_from_tsv(workflow, tsv_file, biom_file)

    return biom_files

