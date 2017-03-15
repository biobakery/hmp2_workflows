# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.common
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains common functions that are used across the board in
all HMP2 AnADAMA2 workflows.
"""

from biobakey_workflow import utils as bbakery_utils
from hmp2_workflows import utils as hmp2_utils


def verify_files(workflow, input_files, checksums_file):
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
        from anadama2 import Workflow

        from hmp2_workflows.common import verify_files

        workflow = Workflow()
        valid_files = verify_files(workflow, 
                                   ['/tmp/fooA.bam', '/tmp/fooB.bam'],
                                   '/tmp/foo_checksums.txt)
    """
    checksums_dict = hmp2_utils.parse_checksums_file(checksums_file)

    for input_file in input_files:
        md5sum = checksums_dict.get(os.path.basename(input_file))
        
        if not md5sum:
            raise KeyError('MD5 checksum not found.', input_file)

        workflow.add_task_gridable('echo "[args[0]]" [depends[0]] | md5sum -c -',
                                   depends =[input_file],
                                   args = [md5sum_str],
                                   time = 24*60,
                                   mem = 1024,
                                   cores = 1)

    ## Kind of wonky but if the workflow doesn't fail than the files we 
    ## passed in should all be valid. Right?
    return input_files


def stage_files(workflow, input_files, target_dir, symlink=False):
    """Moves data files from the supplied origin directory to the supplied
    destination directory. In order to include a file verification check in
    the staging process rsync is used by default to copy files.

    If the symlink parameter is set to True this function will instead create
    symlinks from the origin directory to the target directory.

    An optional parameter may be provided to only stage files with the 
    corresponding extension.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        input_files: A collection of input files to be staged.
        dest_dir (string): Path to destination directory where files should 
                           be moved.
        symlink (boolean): By default create symlinks from the origin 
                           directory to the destination directory. If set to 
                           False files will be copied using rsync.

    Requires:
        rsync v3.0.6+: A versatile file copying tool.

    Returns:
        list: A list of all files that were successfuly staged

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import common

        workflow = anadama2.Workflow()

        staged_files = common.stage_files(workflow, 
                                          ['/tmp/fooA.sam', '/tmp/fooB.sam'],
                                          '/tmp/out_dir')

        workflow.go()
    """
    if not os.path.exists(target_dir):
        raise OSError(2, 'Target directory does not exist', target_dir)

    ## TODO: We need to preserve the file directory structure here because
    ## it tells when the files were received and is used by the website.
    target_files = bbakery_utils.name_files(input_files, target_dir)

    ## TODO: Figure out a better way to handle this rather than creating 
    ## N rsync calls.
    stage_cmd = "rsync -rvz [depends[0]] [targets[0]]"
    if symlink:
        stage_cmd = "ln -s [depends[0]] [targets[0]]"

    workflow.add_task_group(stage_cmd,
                            depends = input_files,
                            targets = target_files)

    return target_files
