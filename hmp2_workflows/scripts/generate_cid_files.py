# -*- coding: utf-8 -*-

"""
generate_cid_files.py
~~~~~~~~~~~~~~~~~~~~~

This script parses a CSV file and attempts to generate a cutplace interface 
definition (CID) file to be used in validation.

An example CID file resembles the following:

D,Format,CSV
D,Header,1
,Name,Example,Empty,Length,Type,Rule
F,ProjectID,,,,Integer,
F,sex,,,,Choice,"0,1"

The first two lines define the type of file and which line the header can be 
found on while the lines starting with 'F' are the fields in the CSV file. 
This script attempts to fill in whether or not an empty value is valid in 
the column, the type of column and define a rule when possible.

The produced CID file should be reviewed manually to make sure all defined
fields are satisfactory.
"""

import argparse
import csv

import pandas as pd


def parse_cli_arguments():
    """Parses any command-line arguments passed into the workflow.

    Args:
        None
    Requires:
        None
    Returns:
        argparse.ArgumentParser: Object containing arguments passed in by user.
    """
    parser = argparse.ArgumentParser('Parses a CSV file and attempts to '
                                     'generate a vladiate-ready validation '
                                     'object that can be used by the HMP2 '
                                     'metadata workflow.')
    parser.add_argument('-i', '--input-csv', required=True, 
                        help='Input CSV file to generate a validator from.')
    parser.add_argument('-o', '--output-file', required=True,
                        help='Target output python file')
    parser.add_argument('-l', '--min-uniq-count', default=10,
                        help='The threshold to consider a collection of '
                        'unique values as a cutplace Choice type.')

    return parser.parse_args()


def create_cid_field_definitions(csv_df, min_uniq_count):
    """Iteraters over all columns in the supplied Pandas DataFrame and 
    attempts to guess the best combination of values to generate a 
    well-formed CID file.

    Args:
        csv_df: The dataframe containing the CSV file to generate a CID 
                file from.
        min_uniq_count: The threshold to consider a collection of unique
                        values as cutplace type Choice. If the number of 
                        unique values is greater than this we do not 
                        consider it of type Choice.
    Requires:
        None
    Returns:
        list: A list containing lines for our CID file.
    """
    cid_list = []

    for (col_name, values) in csv_df.iteritems():
        field_type = "Text"
        allow_empty = ""
        rule = ""

        ## Grab all non-nan unique values for further processing
        unique_values = values.dropna().unique()

        ## If we have any 'nan' values in our collection of unique values 
        ## this indicates we allow empty values in this column.
        if values.isnull().any():
            allow_empty = "X"

        ## We need to be careful here since integers are usually interpreted 
        ## as floats with Pandas. So we need to check if the field really is 
        ## a float or an integer
        if values.dtype == "float64" and len(unique_values) > 0:
            field_type = "Decimal"

            ## So in order to check whether or not we actually have a true 
            ## float column we'll want to check if we have anything to the 
            ## right of the decimal of any of our values
            if all(val % 1 == 0.0 for val in unique_values if val != 0.0):
                field_type = "Integer"
                unique_values = map(int, unique_values)
        elif values.dtype == "int64":
            field_type = "Integer"

        ## If our type is Text or Integer and the number of unique values is 
        ## below a threshold we can likely mark this as of type Choice and 
        ## define a rule
        if (field_type in ['Text', 'Integer'] and 
            len(unique_values) > 0 and len(unique_values) <= min_uniq_count):
            field_type = "Choice"
            rule = ",".join('\'%s\'' % str(val) for val in unique_values)
        elif field_type == "Integer" and 'id' not in col_name.lower():
            ## Might be really presumptous here but if we have integer types
            ## and it isn't an ID of any sort we might be able to use a range
            ## rule here.
            min_val = min(unique_values)
            max_val = max(unique_values)
            rule = "%s...%s" % (min_val, max_val)

        cid_list.append(['F', 
                         col_name, 
                         '', 
                         allow_empty,
                         '', 
                         field_type,
                         rule])

    return cid_list


def write_cid_file(cid_list, out_file):
    """Takes the provided list containing CID definition lines and writes out
    a valid CID file to be used with the source CSV provided.

    Args:
        cid_list: List containing lines for a CID file that can be used to 
                  validate the source CSV file provided.
        out_file: Path to the desired output file to write CID too.                  
    Requires:
        None
    Returns:
        None
    """
    with open(out_file, 'w') as csvfile:
        cidwriter = csv.writer(csvfile)
        cidwriter.writerow(['D', 'Format', 'CSV'])
        cidwriter.writerow(['D', 'Header', '1'])
        
        for row in cid_list:
            cidwriter.writerow(row)


def main(args):
    csv_df = pd.read_csv(args.input_csv)
    cid_list = create_cid_field_definitions(csv_df, args.min_uniq_count)
    write_cid_file(cid_list, args.output_file)
    

if __name__ == "__main__":
    main(parse_cli_arguments())
