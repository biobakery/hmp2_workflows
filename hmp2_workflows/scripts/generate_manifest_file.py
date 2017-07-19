# -*- coding: utf-8 -*-

"""
generate_manifest_file.py
~~~~~~~~~~~~~~~~~~~~~~~~~

This script generates a project MANIFEST file for the HMP2 AnADAMA workflows.
MANIFEST files are used to indicate new files that are ready for processing.

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
import yaml


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
    parser = argparse.ArgumentParser('Generates a MANIFEST file used by the '
                                     'HMP2 AnADAMA2 workflows.')
    parser.add_argument('-b', '--broad-data-sheet', required=True,
                        help='Broad data product status spreadsheet. '
                        'Contains entries indicating new files to be '
                        'processed.')
    parser.add_argument('-o', '--output-manifest', required=True,
                        help='Path to desired output manifest file.')
    parser.add_argument('-oi', '--origin-institute', required=True,
                        help='Name of institute submitting new files '
                        'to be processed.')
    parser.add_argument('-oc', '--origin-contact', required=True,
                        help='Contact person for corresponding origin '
                        'institute.')
    parser.add_argument('-oe', '--origin-contact-email', required=True,
                        help='Contact email for contact person.')
    parser.add_argument('-p', '--project-name', dest='project', 
                        required=True,
                        help='Project that sequence files belong too.')

    return parser.parse_args()


def generate_yaml_dict(data_df, origin_institute, origin_contact, 
                       origin_contact_email, project_name):
    """Generates a YAML-dictionary representation that will compose a project 
    MANIFEST file. An example MANIFEST file is shown below:

    ## This is an example manifest file for an HMP2 AnADAMA2 workflow run
    ##
    ## The file contains metadata that identifies the origin of the files,
    ## the date files were generated, submitted, a contact person for the
    ## data and which types of data are present in this batch of files.

    --------

    ################
    #   METADATA   #
    ################

    # First some information about who generated/submitted this data
    origin_institute: Broad Insititute
    origin_contact: Tiffany Poon
    origin_contact_email: tpoon@broadinstitute.org

    project: HMP2
    date_of_creation: 2017-04-17T17:07:00
    date_of_submission: 2017-04-17T17:30:00

    ################
    #     DATA     #
    ################

    # The files present in this submission. These will be processed by a
    # the appropriate AnADAMA2 pipelines.
    submitted_files:
        16S:
            input:
                - /seq/picard_aggregation/G79182/MSM5LLFU/current/MSM5LLFU.bam
    
    --------

    Args: 
        data_df (pandas.DataFrame): DataFrame containing Broad data product 
            tracking information.
        origin_institute (string): Name of institute where these files 
            originated            
        origin_contact (string): Name of contact person for origin institute.
        origin_contact_email (string): Email for origin contact person.

    Requires:
        None

    Returns:
        dict: A dictionary in the YAML MANIFEST format.
    """
    data_dict = {}

    ## TODO: These shouldn't be hard-coded here.
    data_dict['origin_institute'] = origin_institute
    data_dict['origin_contact'] = origin_contact
    data_dict['origin_contact_email'] = origin_contact_email
    
    data_dict['project'] = project_name

    data_dict['submitted_files'] = {}
    for (data_format, data_rows) in data_df.groupby('Format'):
        input_files = [row[0] for row in data_rows.values]
        data_dict['submitted_files'][data_format] = {}
        data_dict['submitted_files'][data_format]['input_files'] = input_files

    return data_dict


def main(args):
    broad_data_df = pd.read_csv(args.broad_data_sheet)
    yaml_file = generate_yaml_dict(broad_data_df, 
                                   args.origin_institute,
                                   args.origin_contact,
                                   args.origin_contact_email,
                                   args.project)

    yaml.dump(yaml_file, 
              open(args.output_manifest, 'w'), 
              default_flow_style=False)


if __name__ == "__main__":
    main(parse_cli_arguments())
