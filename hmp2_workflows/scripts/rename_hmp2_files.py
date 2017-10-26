# -*- coding: utf-8 -*-

"""
rename_hmp2_files.py
~~~~~~~~~~~~~~~~~~~~

Renames any HMP2 related files that contains one of our many identifiers 
to another identifier that can be found in the HMP2 metadata file.

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

import argparse

import os
import pandas as pd

import shutil

from glob2 import glob


def parse_cli_arguments():
    """Parses any command-line arguments passed into this script.

    Args:
        None

    Requires:
        None

    Returns:
        argparse.ArgumentParser: argparse object containing the arguments
            passed in by the user.
    """

    parser = argparse.ArgumentParser('Renames HMP2-related files from one '
                                    'identifier to another as found in the '
                                    ' HMP2 metadata file.')
    parser.add_argument('-i', '--input-dir', required=True, 
                        help='Input directory containing files to be '
                        'renamed.')
    parser.add_argument('-e', '--input-extension', required=True,
                        default='.fastq', 
                        help='Input extension to search input directory '
                        'for.')                        
    parser.add_argument('-m', '--metadata-file', required=True, 
                        help='HMP2 metadata file.')                    
    parser.add_argument('-f', '--from-id', required=True,
                        help='The identifier that files are currently using '
                        'for naming.')
    parser.add_argument('-t', '--to-id', required=True, default='External ID',
                        help='Identifier that will be used to rename files.')
    parser.add_argument('-d', '--data-type', required=True, 
                        choices=['metagenomics', 'metatranscriptomics',
                        'viromics', 'proteomics', 'amplicon', 'host_genome',
                        'host_transcriptomics', 'metabolomics', 'methylome',
                        'serology', 'biopsy_16S'],
                        help='The data-type of the files being renamed.')
    parser.add_argument('-o', '--output-dir', 
                        help='OPTIONAL. Output directory to write renamed files '
                        'too.')                    
    parser.add_argument('-s', '--symlink', action='store_true', default=False,
                        help='OPTIONAL. Instead of moving files create a symlink.')
    parser.add_argument('-c', '--copy', action='store_true', default=False, 
                        help='OPTIONAL. Instead of moving files create a copy of '
                        'the file.')
    parser.add_argument('--dry-run', action='store_true', default=False,
                        help='OPTIOANAL. If provided do a dry-run of renaming.')                    

    return parser.parse_args()


def main(args):
    input_files = glob(os.path.join(args.input_dir, "*" + args.input_extension))
    metadata_df = pd.read_csv(args.metadata_file, dtype='str')

    for seq_file in input_files:
        (seq_fname, ext) = os.path.splitext(os.path.basename(seq_file))

        row = metadata_df[(metadata_df[args.from_id] == seq_fname) 
                           & (metadata_df['data_type'] == args.data_type)]

        if row.empty:
            continue 

        rename_id = row.get(args.to_id).values[0]                           
        output_dir = args.output_dir if args.output_dir else os.path.dirname(seq_file)
        rename_file = os.path.join(output_dir, rename_id + ext)

        print "Renaming file %s to %s" % (seq_file, rename_file)

        if not args.dry_run:
            if args.symlink:
                os.symlink(seq_file, rename_file)
            elif args.copy:
                shutil.copy(seq_file, rename_file)
            else:
                shutil.move(seq_file, rename_file)


if __name__ == "__main__":
    main(parse_cli_arguments())