# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils.misc
~~~~~~~~~~~~~~~~~~~~~~~~~

A collection of miscellaneous functions needed for the HMP2-IBDMDB AnADAMA2
workflows.

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
import re

import yaml

from biobakery_workflows import utilities as bb_utils


def create_merged_md5sum_file(checksum_files, merged_checksum_file):
    """Parses a list of files containing md5checksums for a respective 
    file in the same directory. These files should only contain an md5checksum
    and when merged together will produce a file that resembles the following:

        f36daa1eaade89762474d956833cffd0 *160928_SM-7BCZ1_400.raw
        0e59da3ff5b9af1e35e8013ef938d60f *160928_SM-7BCZJ_397.raw
        84e16a3e70d3b425df494fe38dc54765 *160928_SM-7GYK6_434.raw

    Filenames for md5 sum files are expected to be in the format of 
    <BASENAME>.<BASE EXTENSION>.md5 -- an example is below:

        fooB.bam.md5
        fooC.txt.md5
        fooD.fastq.md5

    Args:
        checkum_files (list): A list containing paths to one or more 
            md5checksum files.
        merged_checksum_file (string): Path to the desired merged checksum
            output file.

    Requires:
        None

    Returns:
        string: Path to the merged checksums file.

    Example:
        from hmp2_workflows import utils

        md5_list = ['/tmp/fooA.md5', '/tmp/fooB.txt.md5']
        merged_md5_file = '/tmp/foo_merged.txt.md5'
    
        utils.create_merged_md5checksum_file(md5_list, merged_md5_file)
    """
    with open(merged_checksum_file, 'w') as merged_checksum_fh:
        for checksum_file in checksum_files:
            ## Double check to make sure our md5 checksum files are in the proper
            ## naming convention.
            file_parts = os.path.basename(checksum_file).split(os.extsep)
            file_exts = file_parts[1:]

            if not file_exts or len(file_exts) <= 1:
                raise ValueError('Checksum file name format is incorrect', 
                                 checksum_file)

            source_file = checksum_file.replace('.md5', '')
            checksum_str = open(checksum_file).readline().strip()

            if len(checksum_str) != 32:
                raise ValueError('MD5 Checksum value is not valid', checksum_str)

            merged_checksum_fh.write("%s *%s\n" % (checksum_str, source_file))
    
    return merged_checksum_file


def parse_cfg_file(config_file, section=None):
    """Parses the provided YAML config file. If a specific section is 
    provided to parse only this section is returned.

    Args:
        config_file (string): Path to YAML config file 
        section (string): Specific section in the supplied config file 
                          to return.
    Requires:
        None
    Returns: 
        dict: A dictionary containing all configuration parameters found
              in the supplied config file.
    """
    stream = file(config_file, 'r')
    config = yaml.load(stream)
    
    if section and not section in config:
        raise KeyError('Section not found in config file', section)

    return config.get(section) if section else config


def parse_checksums_file(checksums_file):
    """Parses a file containing MD5 checksums in the following format:

        f36daa1eaade89762474d956833cffd0 *160928_SM-7BCZ1_400.raw
        0e59da3ff5b9af1e35e8013ef938d60f *160928_SM-7BCZJ_397.raw
        84e16a3e70d3b425df494fe38dc54765 *160928_SM-7GYK6_434.raw
    
    And returns a dictionary keyed on filename with MD5 checksums as values.

    Args:
        checksums_file (string): Path to file containing MD5 checksums

    Requires:
        None

    Returns:
        dict: A dictionary with filenames as keys and MD5 checksums as values.
    
    Example:
        md5_file = "/tmp/file_checksums.txt"
        md5_map = parse_checksums_file(md5_file)
    """        
    md5_map = {}

    with open(checksums_file) as md5_fh:
        for line in md5_fh:
            (checksum, filename) = line.split()
            
            #filename = filename.split(os.extsep, 1)[0].replace('*', '')
            filename = os.path.basename(filename).replace('*', '')

            # Quick check to make sure the checksum is 32 characters long   
            # since MD5 checksums should always be 32 characters long.
            if len(checksum) != 32:
                raise ValueError('MD5 Checksum value is not valid', checksum)

            md5_map[filename] = checksum

    return md5_map            


def get_sample_id_from_fname(filename):
    """Given a sequence file extract the sample ID that is mappable to one of 
    the metadata sources associated with the HMP2 workflows.
    
    Extracting the sample ID from a given filename is going to vary based 
    off of the source of the file. 

    Broad-sourced:

        * ESM5MEBP.bam -> SM-5MEBP
        * CSM5MCVZ.bam -> SM-5MCVZ
        * MSM5LLE3.bam -> SM-5LLE3

    PNNL-sourced:
        
        * 160916_SM-9W3BK_307.raw -> SM-9W3BK
        * 160626-SM-A2J3A-52.raw  -> SM-A2J3A
        * 160923_SM-71T45_394.raw -> SM-71T45

    Args:
        filename (string): The input filename to parse a sample ID out of

    Requires:
        None

    Returns:
        string: A sample ID that is mappable to HMP2 workflow metadata 
            sources.
    """
    ## Just in-case we don't have a basename here; get our basename
    filename = os.path.basename(filename)

    pattern = re.compile('.*[-|_]?(S[M|m]-?[0-9a-zA-Z]{4,5})[-|_]?.*')
    matches = re.match(pattern, filename)
    if matches:
        sample_id = matches.group(1)
        sample_id = (sample_id.replace(sample_id[:2], sample_id[:2] + '-', 1) if '-' not in sample_id
                     else sample_id)
    else:
        sample_id = os.path.splitext(filename)[0]

    return sample_id

def reset_column_headers(data_frame, header_col=0):
    """Resets the column header on a pandas DataFrame using the supplied
    row index as the new column header.

    Args:
        data_frame (pandas.DataFrame): The pandas DataFrame to reset column
            header
        header_col (int): Row to use as new header row.

    Requires:
        None

    Returns:
        pandas.DataFrame: DataFrame with new column headers.                
    """    
    data_frame.rename(columns=data_frame.iloc[header_col], inplace=True)
    data_frame.drop(data_frame.index[header_col], inplace=True)

    return data_frame
