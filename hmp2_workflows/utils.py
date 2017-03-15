# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils
~~~~~~~~~~~~~~~~~~~~

This module provides utility functions that are used within the AnaDaMa2 
HMP2 workflows.
"""

import os
import subprocess
import tempfile

from itertools import islice


def _chunks(n, iterable):
    """Takes an iterable object and splits it into n "even-sized" chunks.

    Source: http://stackoverflow.com/a/19264525

    Args:
        n (integer): The number of chunks desired.
        iterable (iterable): An iterable object.
    Requires:
        None
    Returns:
        list: A list containing a slice of the iterable object.
    Example:
        chunked_numbers = chunks(2, [1,2,3,4,5,6])
    """
    iterable = iter(iterable)

    while True:
        yield tuple(islice(iterable, n)) or iterable.next()


def _chunk_md5sum_files(checksums_file):
    """ Takes a file containing md5 checksums in the following format:
        
        f36daa1eaade89762474d956833cffd0 *160928_SM-7BCZ1_400.raw
        0e59da3ff5b9af1e35e8013ef938d60f *160928_SM-7BCZJ_397.raw
        84e16a3e70d3b425df494fe38dc54765 *160928_SM-7GYK6_434.raw
    
    and splits it into N chunks each containing a subset of the md5 checksums
    in the file.

    Args:
        checksums_file (string): Path to file containing md5 checksums

    Requires:
        None

    Returns:
        list: A list containing a subset of md5 checksums.

    Example:
        checksums_list = chunk_md5sum_files('/path/to/md5sum.txt')
    """
    ## Stupid waste but we need to know the number of lines to determine an 
    ## efficient chunk size.
    checksums_fh = open(checksums_file)
    num_lines = sum(1 for line in checksums_fh)
    
    checksums_fh.seek(0)
    yield chunks(num_chunks, open(checksums_file))


def stage_data_files(origin_dir, dest_dir, files, symlink=True):
    """Moves data files from the supplied origin directory to the supplied
    destination directory. By default this function will create symlinks
    frmo the deposition directory to the processing directory to save on 
    space used but can be overriden.

    If the symlink parameter is set to False this function will use rsync to
    copy files to the destination directory to ensure file integrity.

    Args:
        origin_dir (string): Origin directory containing files to be moved.
        dest_dir (string): Path to destination directory where files should 
                           be moved.
        symlink (boolean): By default create symlinks from the origin 
                           directory to the destination directory. If set to 
                           False files will be copied using rsync.

    Requires:
        rsync v3.0.6+: A versatile file copying tool.

    Returns:
        list: A list of all files that were successfuly "moved"

    Example:
        staged_files = stage_data_files('/seq/ibdmdb/data_deposition/',
                                        '/seq/ibdmdb/processing/')
    
    """
    for target_dir in [origin_dir, dest_dir]:
        if not os.path.exists(target_dir):
            raise OSError(2, 'No such directory', target_dir)
        if not os.path.isdir(target_dir):
            raise OSError(2, 'Target is not a directory', target_dir)

    if not symlink:
        ## TODO: Add code to rsync copy files over if the symlink flag is 
        ## set to false.
    else:
        def _stage_file(target_file):
            """Handles staging an individual directory from the origin 
            directory to the destination directory. If the file does not exist
            the resulting Exception is squashed to allow downstream processing
            to continue.

            Args:
                target_file (string): The file to be staged.
            Requires:
                None
            Returns:
                None
            Example:
                _stage_file('fileA.bam')
                staged_files = map(_stage_file, files)
            """
            origin_file = os.path.join(origin_dir, target_file)
            dest_file = os.path.join(dest_dir, target_file)

            if not os.path.exists(origin_file):
                raise OSError(2, 'Origin file does not exists', origin_file)
            if not os.path.exists(dest_file):
                raise OSError(2, 'Destination file cannot be staged', 
                              dest_file)

            os.path.symlink(origin_file, dest_file)

        try:
            staged_files = map(_stage_file, files)
        except OSError as oerr
            ## Generally if we end up in here it means that symlinking did 
            ## not work out but we still want to proceed. It would be good to 
            ## log which files failed upstream of this function
            pass
    
    return staged_files
    
