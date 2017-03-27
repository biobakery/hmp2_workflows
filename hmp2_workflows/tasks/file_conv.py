# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.file_conv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains functions that convert between different file types.
"""

from biobakery_workflows import utils as bbutils


def bam_to_fastq(workflow, input_files, output_dir, threads):
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
                                      ext='fastq')

    workflow.add_task_group_gridable('samtools bam2fq [depends[0]] > [target[0]]',
                                     depends = input_files,
                                     targets = output_files,
                                     time = 24*60,
                                     mem = 4096,
                                     cores = threads)

    return output_files
