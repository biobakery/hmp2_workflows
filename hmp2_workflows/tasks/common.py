# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.common
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains common functions that are used across the board in
all HMP2 AnADAMA2 workflows.
"""

from hmp2_workflows import utils


def verify_file_integrity(workflow, input_file_dir, checksums_file):
    """Verifies the integrity of all files found under the supplied directory 
    using md5 checksums. In order for this function to work properly an file 
    contanining md5 checksums must have been generated on the source side of 
    the files and provided when these files were uploaded.


    Args:
        input_file_dir (string): Path to directory containing files to be 
                                 checked.
        checksums_file (string): Path to file containing the source side md5 
                                 checksums.

    Requires:
        None

    Returns:
        list: A list of files that have passed integrity verification

    Example:
        valid_files = verify_input_integrity('/seq/ibdmdb/staging/new_files',
                                             '/seq/ibdmdb/md5_checksums.txt')
    """
    pass
