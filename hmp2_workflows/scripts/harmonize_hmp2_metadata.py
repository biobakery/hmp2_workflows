# -*- coding: utf-8 -*-

"""
harmonize_hmp2_metadata.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

Takes a list of samples that exist in HMP2 metadata but not in any HMP2 
product files and prunes out those rows from the metadata file.
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
    parser = argparse.ArgumentParser('Removes samples present in HMP2 metadata '
                                     'but missing from any product files.')
    parser.add_argument('-m', '--metadata-file', required=True,
                        help='OPTIONAL: An existing metadata file if '
                        'available.')
    parser.add_argument('-s', '--samples-to-remove', required=True,
                        help='A list of the samples to be pruned.')
    parser.add_argument('-o', '--output-file', required=True,
                        help='Desired output pruned HMP2 metadata file.')
    
    return parser.parse_args()


def main(args):
    samples_to_remove = [l.strip() for l in open(args.samples_to_remove).readlines()]

    metadata_df = pd.read_csv(args.metadata_file)

    for (group, items) in itertools.groupby(samples_to_remove, lambda x: x.split('\t')[-1]):
        samples = [sample.split('\t')[0] for sample in list(items)]
        metadata_df.drop(metadata_df[(metadata_df['External ID'].isin(samples)) & 
                         (metadata_df['data_type'] == group)].index, inplace=True)

    metadata_df.to_csv(args.output_file, index=False)                            


if __name__ == "__main__":
    main(parse_cli_arguments())