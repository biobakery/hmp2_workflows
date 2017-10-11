# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.common
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains common functions that are used across the board in
all HMP2 AnADAMA2 workflows.

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

import itertools
import os
import shutil
import tempfile

from biobakery_workflows import utilities as bb_utils
from hmp2_workflows import utils as hmp_utils


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
    checksums_dict = hmp_utils.misc.parse_checksums_file(checksums_file)

    for input_file in input_files:
        md5sum = checksums_dict.get(os.path.basename(input_file))
        
        if not md5sum:
            raise KeyError('MD5 checksum not found.', input_file)

        workflow.add_task_gridable('echo "[args[0]] *[depends[0]]" | md5sum -c -',
                                   depends = [input_file],
                                   args = [md5sum],
                                   time = 24*60,
                                   mem = 1024,
                                   cores = 1)

    ## Kind of wonky but if the workflow doesn't fail than the files we 
    ## passed in should all be valid. Right?
    return input_files


def stage_files(workflow, input_files, target_dir, delete=False, 
                preserve=False, symlink=False):
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
        preserve (boolean): If set to True preserve the source subdirectory 
            structure on the target side.
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
    target_files = bb_utils.name_files(input_files, target_dir)

    ## TODO: Figure out a better way to handle this rather than creating 
    ## N rsync calls.
    stage_cmd = "remove_if_exists.py [targets[0]] ; rsync -avz [depends[0]] [targets[0]]"
    
    if preserve:
        stage_cmd = stage_cmd.replace('-avz', 
                                      '--rsync-path=\"mkdir -p `dirname '
                                      '[depends[0]]`\" -avz')
    if symlink:
        stage_cmd = "remove_if_exists.py [targets[0]] ; ln -s [depends[0]] [targets[0]]"

    workflow.add_task_group(stage_cmd,
                            depends = input_files,
                            targets = target_files)

    return target_files


def make_files_web_visible(workflow, files):
    """Receives a list of files to be disseminated and ensures that they are 
    visible on the IBDMDB website.

    Args:
        files (list): A list of files that should be made visible on the 
                      IBDMDB website.

    Requires:
        None   

    Returns:
        list: A list of all files that should now be publicly visible on the 
            IBDMDB website.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import common

        workflow = anadama2.Workflow()

        visible_files = common.make_files_web_visible(workflow, 
                                                      ['/tmp/public/fooA.sam'])

        workflow.go()
    """
    ## Making a group of files web visible is as simple as touching a file 
    ## named complete.html in the directories containing the files we want 
    ## to show up on the website
    public_files = itertools.chain.from_iterable(files)
    public_dirs = itertools.groupby(public_files, os.path.dirname)
    complete_files = [os.path.join(item[0], 'complete.html') for item
                      in public_dirs]

    workflow.add_task_group('touch [targets[0]]',
                            depends = map(os.path.dirname, complete_files),
                            targets = complete_files)

    ## Again kinda lazy, but the files passed in will be the ones that are
    ## made public.
    return files


def tar_files(workflow, files, output_tarball, depends, compress=True):
    """Creates a tarball package of the provided files with the given output
    tarball file path.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        files (list): A list of files to package together into a tarball.
        output_tarball (string): The desired output tarball file.
        depends (list): A list of files that can be tied to the files to 
            be tar'd. These dependencies are not tar'd but needed to make 
            sure that this step in the workflow is not run out of order.
        compress (string):             

    Requires:
        None

    Returns:
        string: Path to the newly created tarball file.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import common

        workflow = anadama2.Workflow()

        files_to_tar = ['/tmp/foo.txt', '/tmp/bar.txt']
        tar_file = '/tmp/foo_bar.tar'

        out_tar = common.tar_files(workflow, files_to_tar, tar_file)
    """
    ## Though there may be a cleaner way of doing this we need to create 
    ## a temporary folder to symlink all the files we want to package into 
    ## a tarball to get rid of tar'ing the directory structure as well.
    tmp_dir = tempfile.mkdtemp()
    tar_args = "-hcvf"

    for target_file in files:
        symlink_path = os.path.join(tmp_dir, os.path.basename(target_file))
        os.symlink(target_file, symlink_path)
    
    depend_files = []
    if depends:
        depend_files.extend(depends)
    else:
        depend_files.extend(files)

    if compress:
        tar_args = tar_args + "z"
    
    workflow.add_task('tar -hcvzf [targets[0]] -C [args[0]] .; rm -rf [args[0]]',
                      depends = depend_files,
                      targets = [output_tarball],
                      args = [tmp_dir])

    return output_tarball


def generate_md5_checksums(workflow, files):
    """Generates MD5 checksums for the provided set of files. All checksums 
    are written to a file containing the same name as the input but with the 
    "md5" extension appended.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        files (list): A list of files to package together into a tarball.
        output_tarball (string): The desired output tarball file.

    Requires:
        None

    Returns:
        list: A list of the generated md5 checksum files.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import common

        workflow = anadama2.Workflow()

        files = ['/tmp/foo.txt', '/tmp/bar.txt']

        md5sum_files = common.generate_md5_checksums(workflow, files)
    """
    output_dir = os.path.dirname(files[0])
    checksum_files = bb_utils.name_files(bb_utils.sample_names(files),
                                         output_dir,
                                         extension=".md5")

    workflow.add_task_gridable('md5sum [depends[0]] > [targets[0]]',
                               depends=files,
                               targets=checksum_files)

    return checksum_files