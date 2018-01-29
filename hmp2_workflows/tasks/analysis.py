# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Contains any left-over analysis functions that are not found in the 
biobakery workflows.

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


def generate_ko_files(workflow, genefamilies, output_dir):
    """Derives kegg-orthology files from the provided genefamilies files
    using the humann2_regroup_table utility.

    Args:
        worfklow (anadama2.Workflow): The AnADAMA2 workflow instance.
        genefamilies (list): A list of all genefamilies files from 
            which to derive KO files from.
        output_dir (string): The output directory to write KO files too.            

    Requires:
        None

    Returns:
        string: The path to the merged KOs files and merged normalized KOs files.
    """
    pass