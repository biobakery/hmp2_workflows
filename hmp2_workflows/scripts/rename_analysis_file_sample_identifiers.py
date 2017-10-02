"""
rename_analysis_file_sample_identifiers.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Takes any tab-delimited HMP2 product files and renames the header columns
using the type of identifier currently found in the analysis file 
(site_sub_coll, External ID, etc.) to the new identifier specified.

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

import pandas as pd

from hmp2_workflows.utils.misc import parse_cfg_file 


def parse_cli_arguments():
    """Parses any command-line arguments passed in by the user.
    
    Args:
        None

    Requires:
        None

    Returns:
        argparse.ArgumentParser: An ArgumentParser object containing all 
            arguments supplied by the user.
    """
    parser = argparse.ArgumentParser('Replaces column headers in any HMP2 '
                                     'analysis files with a specified identifier '
                                     'that can be found in the HMP2 metadata table.')
    parser.add_argument('-i', '--input-analysis-file', required=True,
                        help='The analysis file to have column identifiers '
                        'replaced.')
    parser.add_argument('-m', '--metadata-file', required=True, 
                        help='HMP2 metadata table that will be used to replace '
                        'column identifiers.')
    parser.add_argument('-c', '--config-file', required=True,
                        help='HMP2 project configuration file.')                        
    parser.add_argument('-oi', '--old-id', required=True,
                        help='Column identifiers field currently in the '
                        'analysis file. Must exist in metadata table.')
    parser.add_argument('-ni', '--new-id', required=True,
                        help='New column identifiers field to replace old '
                        'column identifiers.')
    parser.add_argument('-o', '--output-file', required=True,
                        help='New output file containing new column '
                        'identifiers.')
    parser.add_argument('-d', '--data-type', required=True,
                        help='Data type for these analysis files '
                        '(i.e. metagenomics)')  
    parser.add_argument('--no-tag', default=False, action='store_true',
                        help='OPTIONAL. If set don\'t attempt to parse a '
                        'tag out of the column headers.')                         

    return parser.parse_args()


def get_column_mapping(metadata_df, analysis_cols, old_id, 
                       new_id, dtype, config, no_tag=False):
    """Attempts to map column identifiers 
    """
    col_map = {}    
    ids_not_found = []

    metadata_cols = metadata_df.filter([old_id, new_id])
    replace_strs = config['base']['analysis_col_patterns']

    for sample_id in analysis_cols:
        sample_id_clean = sample_id
        for replace_str in replace_strs:
            sample_id_clean = sample_id_clean.replace(replace_str, '')

        row = metadata_cols[metadata_cols[old_id] == sample_id_clean]

        if not row.empty:
            new_sample_id = row.get(new_id).values[0]
            tags = sample_id.split('_', 1)[-1]

            new_sample_id = new_sample_id if no_tag else new_sample_id + "_" + tags
            col_map[sample_id] = new_sample_id
        else:
            ids_not_found.append((sample_id, sample_id_clean))

    return (col_map, ids_not_found)


def main(args):
    metadata_df = pd.read_csv(args.metadata_file, dtype='str')
    input_analysis_df = pd.read_table(args.input_analysis_file, dtype='str')
    analysis_cols = input_analysis_df.columns.tolist()[1:]
    config = parse_cfg_file(args.config_file)

    metadata_subset_df = metadata_df[metadata_df['data_type'] == args.data_type]

    if (args.old_id not in metadata_df.columns 
        or args.new_id not in metadata_df.columns):
        raise ValueError('Could not find current column identifier or new '
                         'column identifier in HMP2 metadata file.')

    (column_mapping, not_found) = get_column_mapping(metadata_subset_df, 
                                                     analysis_cols, args.old_id,
                                                     args.new_id, args.data_type, 
                                                     config, args.no_tag)

    ## TODO: Deal with the not found IDs here at some point.
    input_analysis_df.rename(columns=column_mapping, inplace=True)

    filter_cols = column_mapping.values()
    filter_cols.insert(0, input_analysis_df.columns[0])
    input_analysis_df = input_analysis_df.filter(filter_cols)
    
    input_analysis_df.to_csv(args.output_file, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main(parse_cli_arguments())
