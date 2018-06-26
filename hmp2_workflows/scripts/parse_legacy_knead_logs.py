# -*- coding: utf-8 -*-

"""
parse_legacy_knead_logs.py
~~~~~~~~~~~~~~~~~~~~~~~~~~
Parses legacy KneadData logs to create the KneadData read counts table 
required by the current visualization pipelines.

Copyright (c) 2018 Harvard School of Public Health

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


import argparse
import os

import pandas as pd

from glob2 import glob

LEGACY_TAG_MAP = {
    'Initial number of reads': 'raw',
    'Total reads after trimming': 'trimmed',
    'Total reads after removing those found in reference database': 'decontaminated Homo_sapiens',
    'Total reads after merging results from multiple databases': 'final'
}

def parse_cli_arguments():
    """Parses any command-line arguments passed into the script.

    Args:
        None
    Requires:
        None
    Returns:
        argparse.ArgumentParser: Object containing arguments passed in by user.

    """
    parser = argparse.ArgumentParser('Parses legacy KneadData logs to create '
                                     'read count tables files.')
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Input directory containing KneadData log files.')                                     
    parser.add_argument('-o', '--output-counts-table', required=True,
                        help='Output table file containing all read count statistics.')
    parser.add_argument('-d', '--data-type', default='MGX', choices=['MGX', 'MTX'],
                        help='OPTIONAL. The data type these log files were generated from. '
                        'Will convey a different set of columns between MGX and MTX')
    parser.add_argument('-p', '--pair-identifier', default="r1", 
                        help='OPTIONAL. If working with paired-end sequences the identifier '
                        'to differentiate between mate pair files.')
    
    return parser.parse_args(columns=[''])


def get_all_log_files(input_dir):
    """Scans the provided input directory to retrieve all KneadData log files.

    Args:
        input_dir (string): The input directory to scan for log files.
    Requires:
        None
    Returns:
        list: A list containing the path to all log files to parse.
    """
    return glob(os.path.join(input_dir, "*.log"))


def _add_aux_info_to_tag(tag, fname, pair_identifier):
    """Adds pair identifier information to a KneadData column header tag if we
    are dealing with a paired-end sample and correctly identifies if we are dealing
    with orphan reads.

    Args:
        tag (string): The KneadData category tag 
        fname (string): The file name we are pulling reads counts out of.
        pair_identifier (string): The pair identifier to differentiate between
            mate pair files.
    Requires:
        None
    Returns:
        string: If dealing with a paired-end sample the tag updated with the proper
            pair identifier appended to it.
    """
    pair_identifier_2 = (pair_identifier[0] + str(int(pair_identifier[-1])+1) 
                         if len(pair_identifier) > 1 else str(int(pair_identifier[0])+1))

    if "single" in fname.lower():
        tag = tag + " orphan"

    if pair_identifier in fname.lower():
        tag = tag + " pair 1"
    elif pair_identifier_2 in fname.lower():
        tag = tag + " pair 2"

    return tag


def _parse_log_file(log_file, data_type, pair_identifier):
    """Parses a log file to extract the several read counts at the many steps 
    of KneadData

    Args:
        log_file (string): Path to the log file to parse
        data_type (string): The data type for the sample that produced this 
            log file.
        pair_identifier (string): The pair identifier to differentiate between
            mate pair files.
    Requires:
        None
    Returns:
        dict: A dictionary containing the read count statistics parsed 
            from this log file.
    """
    sample_name = os.path.splitext(os.path.basename(log_file))[0]
    read_stats = dict.setdefault(sample_name, {})

    for line in open(log_file):
        if "number of reads" in line or "Total" in line:
            (_date, legacy_tag_raw, read_count) = line.split(':')

            (legacy_tag, fname) = legacy_tag_raw.split('(')
            fname.replace(')', '').strip()
            legacy_tag = legacy_tag.strip()
            new_tag = LEGACY_TAG_MAP[legacy_tag]

            new_tag = _add_aux_info_to_tag(new_tag, fname, pair_identifier)
            read_stats[sample_name][new_tag] = read_count

    return read_stats


def parse_legacy_knead_logs(log_files, data_type, pair_identifier):
    """Parses the provided list of legacy KneadData logs and recreates 
    the KneadData read counts tables from all logs.

    Args:
        log_files (list): A list of log files to parse.
        data_type (string): The data type these log files were generated from.
        pair_identifier (string): The pair identifier to differentiate between
            mate pair files.
    Requires:
        None
    Returns:
        pandas.DataFrame: A pandas DataFrame containing read count 
            statistics.
    """
    if   data_type == "MGX":
        columns = ['raw pair1', 'raw pair2', 'trimmed pair1', 'trimmed pair2', 
                   'trimmed orphan1', 'trimmed orphan2', 
                   'decontaminated Homo_Sapiens pair1', 
                   'decontaminated Homo_sapiens pair2',
                   'decontaminated Homo_sapiens orphan1',
                   'decontaminated Homo_sapiens orphan2',
                   'final pair1', 'final pair2', 'final orphan1', 'final orphan2']
    elif data_type == "MTX":
        columns = ['raw pair1', 'raw pair2', 'trimmed pair1', 'trimmed pair2', 
                   'trimmed orphan1', 'trimmed orphan2', 
                   'decontaminated Homo_Sapiens pair1', 
                   'decontaminated Homo_sapiens pair2',
                   'decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA pair1',
                   'decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA pair2',
                   'decontaminated human_hg38_refMrna pair1',
                   'decontaminated human_hg38_refMrna pair2',
                   'decontaminated Homo_sapiens orphan1',
                   'decontaminated Homo_sapiens orphan2',
                   'decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA orphan1',
                   'decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA orphan2',
                   'decontaminated human_hg38_refMrna orphan1',
                   'decontaminated human_hg38_refMrna orphan2',
                   'final pair1', 'final pair2', 'final orphan1', 'final orphan2']

    read_counts_df = pd.DataFrame(columns=columns)

    for log_file in log_files:
        log_stats = _parse_log_file(log_file data_type, pair_identifier)

    return read_counts_df


def main(args):
    log_files = get_all_log_files(args.input_dir)
    read_counts_df = parse_legacy_knead_logs(log_files, args.data_type, args.pair_identifier)


if __name__ == "__main__":
    main(parse_cli_arguments()())