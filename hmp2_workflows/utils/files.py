# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils.file
~~~~~~~~~~~~~~~~~~~~~~~~~

A collection of functions related around file I/O needed for the HMP2-IBDMDB
AnADAMA2 workflows.

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

import datetime
import os

import pandas as pd

from glob2 import glob

from biobakery_workflows import utilities as bb_utils


def create_project_dirs(directories, project, submit_date, data_type):
    """Creates project directories that are required for raw, intermediate
    and output files produced by HMP2 workflows. The template for directories
    created is in the following format:

            <BASE DIR>/<PROJECT>/<TIME STAMP>/<DATA TYPE>

    Returns a list containing all created directories.

    Args:
        directories (list): A list containing all directories which will 
            have the above template applied to create a project directory.
        project (string): The project to use when creating project 
            directories.
        data_type (string): The data type for these project directories.

    Requires:
        None

    Returns:
        list: A list containing create project directories.

    Example:
        from hmp2_workflows import utils

        base_dirs = ['/tmp/foo', '/tmp/bar', '/tmp/baz']
        project = 'HMP'
        data_type = '16S'

        project_dirs = utils.create_project_directories(base_dirs, 
                                                        project, 
                                                        data_type)

        print project_dirs
        ## ['/tmp/foo/HMP/2017-04-17/16S',
            '/tmp/bar/HMP/2017-04-17/16S',
            '/tmp/baz/HMP/2017-04-17/16S']
    """
    project_dirs = []
    date_stamp = str(submit_date)

    for directory in directories:
        target_dir = os.path.join(directory, project, date_stamp, data_type)
        bb_utils.create_folders(target_dir)
        project_dirs.append(target_dir)

    return project_dirs


def find_files(path, pattern=None, extension=None):
    """Searches a directory in a recurisve fashion for all files. If the 
    extension and/or pattern parameters are provided a search for files 
    that match the provided pattern and/or have the provided extension 
    is done.

    Args:
        path (string): Path to search for files.
        pattern (string): A pattern to search for within the provided path.
            Can include the wildcard character '*'
        extension (string): If provided return files that match this 
            extension.

    Requires:
        None

    Returns:
        list: A list of all files that match the provided criteria.

    Example:
        from hmp2_workflows.utils import file

        path = "/tmp/test_dir"
        ext = ".fastq"
        pattern = "foo*"

        match_files = file.find_files_recursive(path, pattern, ext)
        print match_files

        ## ['/tmp/test_dir/foo.fastq',
        ##  '/tmp/test_dir/foo2.fastq',
        ##  '/tmp/test_dir/fooB.fastq',
        ##  '/tmp/test_dir/fooc.fastq',
        ##  '/tmp/test_dir/football.fastq']
    """
    search_path = path

    if pattern and extension:
        search_path = os.path.join(search_path, pattern + extension)
    elif pattern:
        search_path = os.path.join(search_path, pattern)
    elif extension:
        search_path = os.path.join(search_path, "*" + extension)
    else:
        search_path = os.path.join(search_path, "*")

    files = glob(search_path)
    return files
        

def get_sequence_file_from_gid(gid, broad_storage_path):
    """Given a Broad Project GID for an associated sample do a search
    over the Broad sequence storage to find the associated raw sequence 
    files.

    Args:
        gid (string): The Broad Project GID assocaited with a single sample. 
            These identifiers are unique per sequencing run per sample.
        broad_storage_path (string): The path to the Broad sequencing product
            archive.

    Requires:
        None

    Returns:
        string: The path to the associated sequence file if it exists. None 
            if a sequence file cannot be found.

    Example:
        from hmp2_workflows.utils import misc

        gid = "19047774"
        seq_dir = "/n/broad/seq_archive/"

        seq_file = get_sequence_file_from_gid(gid, seq_dir)
        print seq_file
        
        # /n/broad/seq_archive/19047774/seq.bam
    """
    seq_file = None

    broad_gid_dir = os.path.join(broad_storage_path, gid)
    if os.path.exists(broad_gid_dir):
        seq_file_path = os.path.join(broad_gid_dir, "*", "current")
        seq_files = find_files(seq_file_path, extension=".bam")
        
        if seq_files and len(seq_files) == 1:
            seq_file = seq_files[0]
        else:
            raise OSError('Multiple sequence files exist for GID', gid)
    
    return seq_file

def match_tax_profiles(mtx_fastqs, mtx_col_id,
                       tax_profiles, tax_col_id, 
                       metadata_file,
                       tax_tag='_taxonomic_profile'):
    """Takes two sets of files and attempts to match them together based on 
    the supplied HMP2 metadata file. Test

    Args:
        files_a (list): The first set of files to match against.
        files_a_id (string): Column to search for filename/ID in.
        files_b (list): The second set of files to match against.
        files_b_id (string): Column to search for filename/ID in. 
        data_type (string): Data-type that files "B" are.
        metadata_file (string): Path to a tab-delimited look-up file that 
            could provide a way to match any files up.
        tags (list): Any tags that are attached to files that can be 
            stripped prior to matching.

    Requires:
        None

    Returns:
        list: A list of files in set a that match to files set b; in order
        list: A list of files in set b that match to files in set a; in order

    Example:
        from hmp2_workflows.utils import files

        files_a = ['sampleA.fastq', 'sampleC.fastq', 'sampleD.fastq']
        files_b = ['sampleC_tax.tsv', 'sampleA_tax.tsv', 'sampleD_tax.tsv']

        (matched_files_a, matched_files_b) = files.match_files(files_a, files_b)
        # matched_files_a
        # ['sampleA.fastq', 'sampleC.fastq', 'sampleD.fastq']
        # matched_files_b
        # ['sampleA_tax.tsv', 'sampleC_tax.tsv', 'sampleD_tax.tsv']
    """
    matching_mtx_fastq = []
    matching_tax_profiles = []
 
    # Making the assumption here that the sample ID we will use for lookup in 
    # our metadata file is the filename once we remove the extension
    mtx_sample_names = bb_utils.sample_names(mtx_fastqs)
    mtx_sample_map = dict(zip(mtx_sample_names, mtx_fastqs))    
 
    tax_profiles_fnames = map(os.path.basename, tax_profiles)
    tax_profiles_map = dict(zip(tax_profiles_fnames, tax_profiles))

    metadata_df = pd.read_csv(metadata_file)
    metadata_df_subset = metadata_df[(metadata_df.data_type == "metatranscriptomics") &
                                     (metadata_df[mtx_col_id].isin(mtx_sample_names))]

    for (idx, row) in metadata_df_subset.iterrows():
        mtx_id = row.get(mtx_col_id)
        tax_profile_fname = row.get(tax_col_id) + tax_tag

        if tax_profile_fname in tax_profiles_fnames:
            matching_mtx_fastq.append(mtx_sample_map.get(mtx_id))
            matching_tax_profiles.append(tax_profiles_map.get(tax_profile_fname))

    return (matching_mtx_fastq, matching_tax_profiles)