# -*- coding: utf-8 -*-

"""
make_metadata_human_readable.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Takes a metadata table dump from the IBDMDB website and replaces the coded
column headers with human readable ones.

This is a bit of a stop-gap script as the metadata pipeline will do this
properly.
"""

import argparse

import pandas as pd


def parse_cli_arguments():
    """Parses any command-line arguments passed in by the user."""
    parser = argparse.ArgumentParser('Converts coded metadata column names to '
                                     'human readable ones.')
    parser.add_argument('-m', '--metadata-file', required=True, help='IBDMDB '
                        'metadata table file.')
    parser.add_argument('-d', '--data-dictionary', required=True, help='Data '
                        'dictionary containing coded column names to human '
                        'readable column names look-up')
    parser.add_argument('-o', '--output-file', required=True, help='Path to '
                        'desired output file.')

    return parser.parse_args()


def populate_value_lookup(row, lookup):
    """Grabs all valid values for a given column in the studytrax metadata
    file and populates a lookup map of numeric code to human readable value 
    that will be used in creating our human readable metadata file.
    """
    code_col = row.get('Code')
    raw_values = row.get('Pick Lists  (Value, Missing, Name)')
    values = raw_values.split('\r\n')

    lookup.setdefault(code_col, {})
    for value in values:
        (code_val, _null, human_val) = value.split(',', 2)
        lookup[code_col][str(float(code_val))] = human_val.replace(':', '').strip()


def main(args):
    ## First parse the metadata file 
    metadata_df = pd.read_csv(args.metadata_file, dtype='object')

    ## Then parse the data dictionary
    dictionary_df = pd.read_excel(args.data_dictionary)

    ## Now let's create a dictionary lookup for the coded to human readable
    col_name_lookup = pd.Series(dictionary_df['Variable Name'].values, index=dictionary_df['Code'].values)
    
    value_field_lookup = {}
    dictionary_df[dictionary_df['Pick Lists  (Value, Missing, Name)'].notnull()].apply(populate_value_lookup, axis=1, args=(value_field_lookup,))

    ## Need to account for one pesky case here
    value_field_lookup['diagnosis'] = dict((key.replace('.0', ''), val) for (key, val) in value_field_lookup['diagnosis'].iteritems())

    ## Replace coded values
    metadata_df.replace(value_field_lookup, inplace=True)

    ## Replace some other left-over yes/no fields
    replace_yes_no = {'0': 'No', '1': 'Yes', '0.0': 'No', '1.0': 'Yes'}
    replace_cols = ['bx_q31', 'bx_q33', 'bx_q35', 'i_q3', 'i_q4',
                    'i_q5', 'i_q6', 'i_q7', 'i_q8', 'i_q9', 'i_q10',
                    'i_q11', 'i_q12', 'i_q13', 'i_q14', 'i_q15', 
                    'ic_q1', 'ic_q5', 'ic_q6',
                    'i_q16', 'i_q17', 'i_q18', 'i_q19', 'i_q20', 'i_q21',
                    'i_q22', 'i_q23', 'i_q24', 'i_q25', 'i_q26', 'i_q27',
                    'i_q28', 'i_q29', 'i_q30', 'i_q31', 'i_q32', 'i_q33',
                    'i_q34', 'i_q35', 'i_q36', 'i_q37', 'i_q38', 'i_q39',
                    'i_q40', 'i_q41', 'i_q42', 'i_q43', 'i_q44', 'i_q45',
                    'i_q46', 'i_q47', 'i_q48', 'i_q49', 'i_q50', 'bl_q12',
                    'bl_q14', 'bl_q16', 'st_q1', 'st_q9', 'st_q19', 'st_q21',
                    'st_q23', 'hbi_q2', 'hbi_q9', 'hbi_q10', 'hbi_q11', 'hbi_q12',
                    'hbi_q13', 'hbi_q14', 'hbi_q15', 'hbi_q16', 'sccai_q2', 'sccai_q2',
                    'sccai_q11', 'sccai_q12', 'sccai_q13', 'sccai_q14', 'ses_score2', 
                    'mbs_q15', 'dr_q2', 'dr_q2a', 'dr_q2b', 'dr_q2c', 'dr_q3', 'dr_q4',
                    'dr_q5', 'dr_q6', 'dr_q7']

    map(lambda field: metadata_df[field].replace(replace_yes_no, inplace=True), replace_cols)
                    
    ## Rename and write out new CSV file
    metadata_df.rename(columns=col_name_lookup, inplace=True)
    metadata_df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main(parse_cli_arguments())