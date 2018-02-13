# -*- coding: utf-8 -*-

"""
qc_biopsy_and_blood_samples.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Takes an HMP2 manifest file containing different blood and biopsy samples 
(Host Genome, Host Transcriptomics, 16S biopsy, serology) and attempts to 
map all samples to the StudyTrax clinical metadata.
"""

import argparse
import itertools
import os

import funcy
import pandas as pd

from hmp2_workflows.utils.misc import parse_cfg_file 


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

    parser = argparse.ArgumentParser('Verifies and maps blood and biopsy samples '
                                     'against a StudyTrax metadata sheet.')
    parser.add_argument('-i' ,'--input-manifest', required=True,
                        help='Input MANIFEST file containing samples to map.')
    parser.add_argument('-s', '--studytrax-metadata', required=True,
                        help='StudyTrax clinical metadata sheet.')
    parser.add_argument('-b', '--broad-sample-sheet', required=True,
                        help='Broad sample tracking spreadsheet.')                        
    parser.add_argument('-c', '--config-file', required=True,
                        help='Project configuration file.')                        

    return parser.parse_args()


def check_sample_mapping(samples, clinical_df, mapping_cols):
    """
    Checks the provided columns for the list of sample ID's and returns samples
    that are found. 

    Args:
        samples (list): A list of samples to search for in the clinical metadata
        clinical_df (pandas.DataFrame): DataFrame containing clinical metadata
        mapping_cols (list): The columns in the clinical metadata to search 
            for the provided sample IDs in.

    Requires: 
        None

    Returns:
        set: A set of the sample ID's found in the provided columns in the 
            clinical metadata.                        
    """
    found_sample_ids = []

    for col in mapping_cols:
        found_df = clinical_df[clinical_df[col].isin(samples)]
        found_sample_ids.extend(found_df[col].tolist())

    return set(found_sample_ids)


def main(args):
    manifest = parse_cfg_file(args.input_manifest)
    analysis_cfg = parse_cfg_file(args.config_file)
    mapping_cols = analysis_cfg.get('mapping_columns')

    data_files = manifest.get('submitted_files')
    studytrax_df = pd.read_csv(args.studytrax_metadata, dtype='str')
    broad_sample_df = pd.read_csv(args.broad_sample_sheet, dtype='str')

    # When we have some samples that don't map to their desired columns 
    # we are going to want to run them through every single sample ID column
    # and see if we can find a hit.
    clinical_search_cols = ['st_q13', 'st_q12', 'st_q10', 'st_q4', 'st_q17', 'st_q11']
    broad_search_cols = ['Viromics', 'MbX', 'Proteomics', 'Parent Sample A', 
                         'Parent Sample B', 'DNA/RNA']

    if data_files:
        for (dtype, mapping_cols) in [('MVX', mapping_cols.get('MVX')),
                                      ('HTX', mapping_cols.get('HTX')), 
                                      ('RRBS', mapping_cols.get('RRBS')),
                                      ('SER', mapping_cols.get('SER')),
                                      ('HG', mapping_cols.get('HG'))]:

            if not data_files.get(dtype):
                continue

            samples = set([os.path.splitext(os.path.basename(sample_id))[0] 
                           for sample_id in data_files.get(dtype).get('input')])

            if not samples:
                continue

            mapping_cols = funcy.flatten(mapping_cols)
            found_samples = check_sample_mapping(samples, studytrax_df, mapping_cols)
            print "Correctly mapped %s samples" % len(found_samples)

            missing_samples = samples - found_samples
            print "%s samples were not mapped" % len(missing_samples)
            print

            if len(missing_samples) == 0:
                continue

            # Sometimes we have a combination of characters that were recorded incorrectly 
            # (i.e. 1 becomes I or O becomes 0 so we want to isolate and test for all of these)
            for (orig_char, replace_char) in [('1', 'I'), ('O', '0'), ('I', '1'),
                                              ('0', 'O'), ('SM-', 'SM')]:
                mod_samples = map(lambda s: s.replace(orig_char, replace_char), missing_samples)
                found_samples = check_sample_mapping(mod_samples, studytrax_df, mapping_cols)

                if found_samples:
                    print "Found %s more samples after replacing character %s with %s in sample IDs:" % (len(found_samples), orig_char, replace_char)
                    print "\n".join(found_samples)
                    print

                    found_samples = map(lambda s: s.replace(replace_char, orig_char), found_samples)
                    missing_samples = missing_samples - set(found_samples)

            for col in clinical_search_cols:
                found_samples = check_sample_mapping(missing_samples, studytrax_df, [col])

                if found_samples:
                    print
                    print "Found %s samples for data type %s in incorrect column %s:" % (len(found_samples), dtype, col)
                    print "\n".join(list(found_samples))
                    print

                    missing_samples = missing_samples - found_samples

            # For Viromics data we can reference the Broad tracking sheet.
            if dtype == 'MVX':
                for col in broad_search_cols:
                    found_samples = check_sample_mapping(missing_samples, broad_sample_df, [col])

                    if found_samples:
                        print
                        print "Found %s samples for data type %s in Broad tracking sheet, column %s:" % (len(found_samples), dtype, col)
                        print "\n".join(list(found_samples))
                        print

                    missing_samples = missing_samples - found_samples                        

            print "Final missing samples: %s" % (len(missing_samples))
            print "\n".join(missing_samples) 


if __name__ == "__main__":
    main(parse_cli_arguments())