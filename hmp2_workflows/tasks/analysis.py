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

import pandas as pd

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


def generate_relative_abundance_file(workflow, counts_table, output_dir):
    """Computes the relative abundances for the provided tab-delimited counts file 
    (can be taxonomic counts, OTU counts, etc.) and returns the path to the newly 
    created file.

    Args:
        worfklow (anadama2.Workflow): The AnADAMA2 workflow instance.
        counts_table (string): Path to the tab-delimited counts file to generate
            relative abundances from.
        output_dir (string): Output directory to write relative abundance file too.

    Requires:
        None

    Returns:
        string: Path to the newly created relative abundance file.

    Example:
        from hmp2_workflows.tasks.analysis import generate_relative_abundance_file

        tax_profile = "/tmp/tax_counts.tsv"
        output_dir = "/tmp/new"

        rel_abund_file = generate_relative_abundance_file(tax_profile, output_dir)

        print rel_abund_file
        # /tmp/new/tax_counts.rel_abund.tsv
    """    
    input_basename = os.path.splitext(os.path.basename(counts_table))[0]
    abundance_table = os.path.join(output_dir, input_basename.replace('.tsv', '.rel_abund.tsv'))

    def _compute_relative_abundances(task):
        """Take an input counts file and computes a corresponding relative 
        abundances file.
        """
        counts_df = pd.read_table(counts_table, index_col=0)
        abund_df = counts_df / counts_df.sum()
        abund_df.write_csv(abundance_table, sep='\t', index=True)

    workflow.add_task(_compute_relative_abundances,
                      depends=counts_table,
                      target=abundance_table)

    return abundance_table