# -*- coding: utf-8 -*-

"""
harmonize_hmp2_pcl.py
~~~~~~~~~~~~~~~~~~~~~

Takes a list of samples that exist in HMP2 TSV file but not in any HMP2 
metadata file and prunes out those columns from the TSV file.
"""

import argparse
import itertools

import pandas as pd


def parse_cli_arguments():
    """Parses any command-line arguments passed into the script.

    Args:
        None
    Requires:
        None
    Returns:
        argparse.ArgumentParser: Object containing arguments passed in by user.

    """
    parser = argparse.ArgumentParser('Removes samples present in HMP2 TSV file '
                                     'but missing from the HMP2 metadata file.')
    parser.add_argument('-i', '--input-tsv-file', required=True,
                        help='HMP2 product TSV file.')
    parser.add_argument('-s', '--samples-to-remove', required=True,
                        help='A list of the samples to be pruned.')
    parser.add_argument('-o', '--output-file', required=True,
                        help='Desired output pruned HMP2 TSV file.')
    
    return parser.parse_args()


def main(args):
    samples_to_remove = [l.strip() for l in open(args.samples_to_remove).readlines()]

    product_df = pd.read_csv(args.metadata_file)
    product_df.drop(samples_to_remove, axis=1, inplace=True)
    product_df.to_csv(args.output_file, sep='\t')


if __name__ == "__main__":
    main(parse_cli_arguments())