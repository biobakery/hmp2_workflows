# -*- coding: utf-8 -*-

"""
generate_manifest_file_from_file_list.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script generates a project MANIFEST file for the HMP2 AnADAMA workflows.
MANIFEST files are used to indicate new files that are ready for processing 
from a text file containing a list of files, one-per-line.

Each file can be associated with a data-type that is prepended to the start of
each line as shown below in the example:

    $ cat manifest_files.txt
    mgx;/local/data/mbx/sample1.bam
    16s;/local/data/16s/sample2.bam
    /local/data/mgx/sample3.bam

In the case of a file without an associated data-type the data-type supplied
via the parameter "--data-type" will be applied.    

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
import datetime

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
                                     'HMP2 AnADAMA2 workflows from a text '
                                     'file containing a list of data files.')
    parser.add_argument('-i', '--input-file-list', required=True,
                        help='A list of data files to use as the basis for '
                        'a manifest file.')
    parser.add_argument('-m', '--md5sums-file', default=None,
                        help='MD5sums file for all files provided.')
    parser.add_argument('-of', '--output-file-list',
                        help='A list of output files to add to the manifest '
                        'file.')
    parser.add_argument('-oi', '--origin-institute', required=True,
                        help='Name of institute submitting new files '
                        'to be processed.')
    parser.add_argument('-oc', '--origin-contact', required=True,
                        help='Contact person for corresponding origin '
                        'institute.')
    parser.add_argument('-oe', '--origin-contact-email', required=True,
                        help='Contact email for contact person.')
    parser.add_argument('-p', '--project', required=True,
                        help='Project that sequence files belong too.')
    parser.add_argument('-d', '--data-type',
                        choices=['MGX', 'MBX', '16S', 'MTX', 'MPX', 'TX',
                                 'BP', 'HG', 'RRBS', 'SER', 'MVX'],
                        help='A blanket data-type to apply to all files in '
                        'the input file list. Over-riden by any data type '
                        'specified in the input file list.')
    parser.add_argument('-o', '--output-manifest', required=True,
                        help='Path to desired output manifest file.')

    return parser.parse_args()


def parse_file_list(file_list, data_type):
    """Parses an list of data files, one per line, and returns a dictionary
    containing files grouped by data type (one of mgx (Metagenome), 
    mbx (Metabolome), mtx (Metatranscriptome), 16s (16S) or mpx (Proteome)).

    An example of an input file can be seen below:

    ----------

    ## Format is <DATA TYPE>;<PATH TO FILE>
    ##
    ## Can also see one data file per line which will be assigned to data
    ## type supplied in --data-type parameter.

    mgx;/local/data/mbx/sample1.bam
    16s;/local/data/16s/sample2.bam
    /local/data/mgx/sample3.bam

    Args:
        file_list (string): A path to a text file containing a list of 
            data files to generate a manifest file for.
        data_type (string): A data-type to apply to any file in the supplied 
            input file list that does not have an explicit data-type specified.
            Must be one of mbx, mgx, mtx, 16s or prot.

    Require:
        None

    Returns:
        dict: A dictionary containing each set of files grouped with their 
            corresponding data type.                    
    """
    files_dict = {}

    with open(file_list, 'r') as input_fh:
        for line in input_fh:
            file_elts = line.strip().split(';')

            if len(file_elts) > 1:
                files_dict.setdefault(file_elts[0], []).append(file_elts[1])
            else:
                if not data_type:
                    raise ValueError('Must supply default data-type if '
                                     'providing files with no data type:',
                                     line)

                files_dict.setdefault(data_type, []).append(file_elts[0])                                            

    return files_dict


def generate_yaml_dict(data_files, md5sums_file, origin_institute,
                       origin_contact, origin_contact_email, project_name):
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
        data_files (dict): DataFrame containing Broad data product 
            tracking information.
        md5sums_file (string): Path to a file containing md5 checksums for 
            all files to be in our manifest file.
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
    now = datetime.datetime.now()

    data_dict['origin_institute'] = origin_institute
    data_dict['origin_contact'] = origin_contact
    data_dict['origin_contact_email'] = origin_contact_email
    data_dict['project'] = project_name

    data_dict['submitted_files'] = {}
    data_dict['submission_date'] = now.strftime('%Y-%m-%d')

    for (data_type, files) in data_files.iteritems():
        data_dict['submitted_files'].setdefault(data_type, {})
        data_dict['submitted_files'][data_type].setdefault('input', [])
        data_dict['submitted_files'][data_type]['input'] = files.get('input')


        if md5sums_file:
           data_dict['submitted_files'][data_type]['md5sums_file'] = md5sums_file
        if 'output' in files:
            data_dict['submitted_files'][data_type].setdefault('output', [])
            data_dict['submitted_files'][data_type]['output'] = files.get('output')

    return data_dict


def main(args):
    input_files_dict = parse_file_list(args.input_file_list,
                                             args.data_type)

    output_files_dict = None
    if args.output_file_list:
        output_files_dict = parse_file_list(args.output_file_list,
                                            args.data_type)
    files_dict = {}

    for data_type in input_files_dict.keys():
        files_dict.setdefault(data_type, {})
        files_dict[data_type].setdefault('input', input_files_dict[data_type])
      
        if output_files_dict and data_type in output_files_dict:
           files_dict[data_type].setdefault('output', output_files_dict[data_type])

    yaml_file = generate_yaml_dict(files_dict,
                                   args.md5sums_file,
                                   args.origin_institute,
                                   args.origin_contact,
                                   args.origin_contact_email,
                                   args.project)

    yaml.dump(yaml_file, 
              open(args.output_manifest, 'w'), 
              default_flow_style=False)


if __name__ == "__main__":
    main(parse_cli_arguments())
