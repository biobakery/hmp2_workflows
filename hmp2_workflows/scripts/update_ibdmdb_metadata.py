#!/usr/bin/env python

"""
merge_ibdmdb_metadata.py
~~~~~~~~~~~~~~~~~~~~~~~~~

This script takes the multiple sources of metadata used by the IBDMDB 
project and merges them together into one finalized large metadata table.

The primary work done here is normalizing all the ID's used across the board
and joining all the sources together properly.

Currently metadata is gathered from the following sources:

    * StudyTrax metadata - All clinical data; provided by MGH
    * Sample tracking status - Sample status and what sequencing products 
        have been derived from each sample; provided by Broad
    
Additionally miscellaneous metadata is handled on a case-by-case basis 
as these need special handling. Currently the following misc. metadata 
is used in the merging procedure:
    
    * Proteomics metadata - Metadata tied to Proteomics input and output
        files; provided by PNNL

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
import os

import glob2
import numpy as np
import pandas as pd

import biobakery_workflows.utilities as bb_utils

from hmp2_workflows.utils.misc import (parse_cfg_file, 
                                       get_sample_id_from_fname)



def parse_cli_arguments():
    """Parses any command-line arguments passed into the script.

    Args:
        None
    Requires:
        None
    Returns:
        argparse.ArgumentParser: Object containing arguments passed in by user.

    """
    parser = argparse.ArgumentParser('Creates merged IBDMDB metadata file '
                                     'from several sources of metadata.')
    parser.add_argument('-c', '--config', required=True,
                        help='Project configuration containing parameters '
                        'required by this script.')
    parser.add_argument('-m', '--metadata-file',
                        help='OPTIONAL: An existing metadata file if '
                        'available.')
    parser.add_argument('-o', '--output-dir',
                        help='Desired output folder to write updated metadata '
                        'and auxillary file too.')
    parser.add_argument('-f', '--manifest-file',
                        help='Manifest file detailing new or updated data '
                        'files from which metadata should be generated.')
    parser.add_argument('-s', '--studytrax-metadata', required=True,
                        help='StudyTrax clinical metadata.')
    parser.add_argument('-b', '--broad-sample-tracking', required=True,
                        help='Broad Sample Tracking spreadsheet containing '
                        'information about each sample and all products '
                        'generated from each sample.')
    parser.add_argument('-p', '--proteomics-metadata', 
                        help='Any metadata associated with Proteomics data '
                        'supplied by the PNNL.')
    parser.add_argument('-f', '--dump-full-table', default=False,
                        action='store_true')                        
    parser.add_argument('-a', '--refresh-all', action="store_true",
                        help='OPTIONAL: Refresh all metadata from scratch '
                        'before adding files specified in accompanying ' 
                        'manifest.')
    
    return parser.parse_args()


def _add_previous_collection_date(dataframe):
    """For each row in a dataframe takes the previous date of reception
    and adds it as a new column to aid in computation of the week_num and
    interval_days columns in our final metadata table.

    Args:
        dataframe (pandas.DataFrame): A dataframe containing all rows 
            of metadata grouped together for a given subject.
    
    Requires:
        None

    Returns:
        pandas.DataFrame: Modified grouped dataframe containing a new 
            column containing a reference to the previous collection 
            date when available.
    """
    dataframe['prev_coll_date'] = dataframe['Actual Date of Receipt'].shift()
    return dataframe


def get_collection_dates(metadata_df):
    """Retrieves all collection dates for a given subject and returns a 
    dicitonary key'd on subject.

    Args:
        metadata_df (pandas.DataFrame): Broad sample tracking status 

    Requires:
        None

    Returns:
        dict: A dictionary key'd on subject ID containing rows of metadata
            containing collection dates per visit.
    """
    collection_dict = dict((subj, df) for (subj, df) in 
                           metadata_df.ix[:,'Subject':'Actual Date of Receipt']
                           .sort_values(by='Actual Date of Receipt').groupby('Subject'))
    collection_dict = dict((subj, _add_previous_collection_date(df))
                            for (subj, df) in collection_dict.iteritems())

    return collection_dict


def get_all_sequence_files(input_dir, extensions):
    """Scans the given directory (recursively) for all sequence files 
    that will be used to construct a new IBDMDB metadata file.

    Args:
        input_dir (string): Input directory to search for sequence files.
        extensions (list): A list of string extensions to filter files 
            with.

    Requires:
        None

    Returns:
        list: A list containing all sequence files that will be used to 
            construct an IBDMDB metadata table.
    """
    seq_files = [seq_file for seq_file in glob2.glob(os.path.join(input_dir, '**', '*'))
                 if os.path.splitext(os.path.basename(seq_file))[-1] 
                 in extensions]
    
    return seq_files


def get_project_id(row):
    """Populates the 'Project' column in the HMP2 metadata table based off 
    whether or not data already exists or data can be pulled from an 
    auxillary column. 

    Args:
        row (pandas.Series): Row of metadata from an HMP2 metadata table

    Requires:
        None

    Returns:
        pandas.Series: Row of metadata with Project column populated.
    """
    project_id = row.get('Project')

    ## This specific case is applicable to Proteomics data only but the 
    ## function can be expanded to handle other scenarios in the future
    if project_id is None and row.get('Job') is not None:
        project_id = row['Job']

    return project_id


def get_pdo_number(row):
    """Populates the PDO number when it can be obtained from other pieces of 
    metadata.

    Args:
        row (pandas.Series): Row of metadata from an HMP2 metadata table

    Requires:
        None

    Returns:
        string: The corresponding PDO number for the given metadata row
    """
    pdo_num = row.get('PDO Number')
    gid = row.get('Project')

    ## This specific case is applicable to Proteomics data only but the 
    ## function can be expanded to handle other scenarios in the future
    if (not isinstance(gid, float) and '.raw' in gid 
        and row.get('data_type') == 'proteomics'):
        filename = os.path.basename(gid).replace('_', '-')
        batch_num = filename.split('-')[0]
        
        if batch_num.isdigit():
            pdo_num = batch_num

    return pdo_num


def get_data_type(row):
    """Populates the 'data_type' column in the HMP2 metadata table. Attempts
    guess what kind of data type the row of metadata is referencing based off 
    off of the values in other columns.

    Args:
        row (pandas.Series): Row of metadata from an HMP2 metadata table

    Requires:
        None

    Returns:
        pandas.Series: Row of metadata with data_type column populated
    """
    data_type = None

    if not pd.isnull(row.get('Job')):
        data_type = 'proteomics'

    return data_type


def get_metadata_rows(config, studytrax_df, sample_df, proteomics_df, 
                      data_type, sequence_files, pair_identifier):
    """Extracts metadata from the supplied sources of metadata for the
    provided sequence files. 

    Args:
        config (dict): Configuration parameters for metadata
        studytrax_df (pandas.DataFrame): StudyTrax clinical metadata
        sample_df (pandas.DataFrame): Broad sample status metadata
        proteomics_df (pandas.DataFrame): Proteomics metadata.
        data_type (string): Data type of the provided sequence files
        sequence_files (list): A list of sequence files that metadata
            should be pulled for if available.
        pair_identifier (string): If working with paired-end files the 
            identifier to distinguish the first file from its pair.            

    Requires:
        None

    Returns:
        pandas.DataFrame: Slice of metadata for files provided.
    """
    metadata_df = None
    sample_mapping = dict(zip(bb_utils.sample_names(sequence_files, pair_identifier),
                              map(get_sample_id_from_fname, sequence_files)))
    sample_ids = sample_mapping.values()   
    sample_ids = [sid.replace(pair_identifier, '') for sid in sample_ids]

    data_type_mapping = config.get('dtype_mapping')

    ## Grab subset of Broad sample tracking spreadsheet
    sample_subset_df = sample_df[(sample_df['Parent Sample A'].isin(sample_ids)) |
                                 (sample_df['Proteomics'].isin(sample_ids)) |
                                 (sample_df['MbX'].isin(sample_ids)) |
                                 (sample_df['Site/Sub/Coll']).isin(sample_ids)]

    ## TODO: Figure out if we have any samples that did not have aassociated metadata
    metadata_df = sample_subset_df.merge(studytrax_df, 
                                         left_on='Parent Sample A',
                                         right_on='st_q4',
                                         how='left')

    ## We sometimes get a situation where our studytrax metadata is missing 
    ## some of the proteomics sample ID's so we need to make sure we 
    ## replicate them.
    metadata_df.loc[metadata_df['st_q17'].isnull(), 'st_q17'] = metadata_df['Proteomics']
    metadata_df.loc[metadata_df['st_q11'].isnull(), 'st_q11'] = metadata_df['MbX']

    if proteomics_df is not None:   
        ## In order to merge our proteomics data properly we'll need to 
        ## first create a subset of our Broad sample tracking sheet 
        ## that isolates just rows related to Proteomics data (the column 
        ## Proteomics Status should be EXPORTED). 
        sample_filter_df = sample_subset_df[sample_subset_df['Proteomics status'] == 'EXPORTED']
        sample_filter_df = sample_filter_df[['Parent Sample A', 'Proteomics']]

        proteomics_df['sample_ids'] = proteomics_df['Dataset'].replace(sample_mapping)
        proteomics_df['PDO Number'] = proteomics_df['Dataset'].map(lambda did: did.replace('-', '_')
                                                                                  .split('_')[0])
        proteomics_df = sample_filter_df.merge(proteomics_df,
                                               left_on='Proteomics',
                                               right_on='sample_ids',
                                               how='right')
        proteomics_df = proteomics_df.drop('Proteomics', 1)

        metadata_df = metadata_df.merge(proteomics_df,
                                        on='Parent Sample A',
                                        how='left')
    
    metadata_df['data_type'] = data_type_mapping.get(data_type)

    return metadata_df
    

def remove_columns(metadata_df, drop_cols):
    """Removes the specified columns from the provided DataFrame.

    Args:
        metadata_df (pandas.DataFrame): DataFrame to remove columns from
        drop_cols (list): A lost of column names to remove

    Requires:
        None

    Returns:
        pandas.DataFrame: DataFrame with specified columns removed.
    """
    return metadata_df.drop(drop_cols, 1)


def generate_external_id(row):
    """Retrieves or produces the external ID for the given row of metadata. 
    An external ID is an identifier that the contributor of a given sample 
    can identify their samples via.

    Args:
        row (pandas.Series): A row of metadata from our metadata table.

    Requires:
        None

    Returns:
        string: The external ID for the given row of metadata
    """
    stool_id = row['st_q4']

    if pd.isnull(stool_id):
        stool_id = row['Parent Sample A']

    site_sub_coll = row['Site/Sub/Coll']
    return site_sub_coll[0] + stool_id.replace('-', '')


def fix_site_sub_coll_id(row, site_mapping):
    """Adds the SiteName abbreviation to the Site/Sub/Coll ID in the instances
    where it is not present. This abbreviation is necessary to de-dupe rows.

    Args:
        row (pandas.Series): A row of metadata from our metadata table.
        site_mapping (dict): The mapping of Site Name to Site abbreviation

    Requires:
        None

    Returns:
        string: The completed Site/Sub/Coll ID
    """
    site_sub_coll_id = row['Site/Sub/Coll ID']
    site_abbrev = site_sub_coll_id[0]

    if site_abbrev.isdigit():
        site_sub_coll_id = site_mapping.get(row['SiteName']) + site_sub_coll_id

    return site_sub_coll_id


def fill_visit_nums(row):
    """In the case of a missing visit number in a metadata row will attempt
    to parse the visit number from the Site/Sub/Coll ID. This ID should 
    encapsulate the visit number in the following format: XXXXC<VISIT_NUM>
    
    Example: C3010C9

    Args:
        row (pandas.Series): A row of metadata from our metadata table.

    Requires:
        None

    Returns:
        string: The corresponding visit number for the given row.
    """
    site_sub_coll_id = row['Site/Sub/Coll ID']
    visit_num = row['visit_num']

    if np.isnan(visit_num):
        visit_num = site_sub_coll_id.split('C')[-1]
    
    return visit_num


def generate_collection_statistics(metadata_df, collection_dict):
    """Generates the week_num and interval_days columns which contain
    the number of weeks between the past collection date and days between 
    the last collection date respectively.

    Args:
        metadata_df (pandas.DataFrame): DataFrame containing all metadata
        collection_dict (dict): Dictionary containing collection dates for 
            each subject grouped by subject ID.

    Requires:
        None

    Returns:
        pandas.DataFrame: Updated DataFrame with week_num and interval_days
            columns populated for each row.
    """
    def _add_collection_columns(row):
        """Adds the week_num and interval_days columns to the provided row 
        of a dataframe.

        Args:
            row (pandas.Series): A row of the metadata dataframe
    
        Requires:
            None

        Returns:
            pandas.Series: Modified row containing populated week_num and 
                interval_day columns.
        """
        subject_id = row['Subject']
        visit_num = row['Collection #']
        receipt_date = row['Actual Date of Receipt']
        
        collection_dates = collection_dict.get(subject_id)
        initial_visit_date = collection_dates.iloc[0]['Actual Date of Receipt']
        subj_collection_row = collection_dates[collection_dates['Collection #'] 
                                               == visit_num]
        prev_visit_date = subj_collection_row['prev_coll_date'].values[0]

        if pd.isnull(prev_visit_date):
            prev_visit_date = initial_visit_date

        row['week_num'] = (receipt_date - initial_visit_date).days / 7
        row['interval_days'] = (receipt_date - prev_visit_date).days

        return row

    return metadata_df.apply(_add_collection_columns, axis=1)


def reorder_columns(metadata_df, cols_to_move):
    """Re-order's column headers for an HMP2 metadata table to match
    the ordering required for our final output files.

    Args:
        metadata_df (pandas.DataFrame): HMP2 metadata table stored in a
            DataFrame
        cols_to_move (list): A list of columns to move to the front of
            our metadata table.

    Requires:
        None

    Returns:
        pandas.DataFrame: DataFrame with new ordering.
    """
    metadata_cols = metadata_df.columns.tolist()
    metadata_cols = [col for col in metadata_cols 
                     if col not in cols_to_move]
    metadata_cols = cols_to_move + metadata_cols
    metadata_df = metadata_df[metadata_cols]

    return metadata_df


def get_no_sequence_metadata(metadata_df, clinical_df):
    """Generates a dataframe containing any remaiing clinical metadata that 
    is not tied to sequence data. 

    Args:
        metadata_df (pandas.DataFrame): HMP2 metadata table stored in a pandas
            DataFrame.
        clinical_df (pandas.DataFrame): StudyTrax clinical metadata stored in
            pandas DataFrame.

    Requires:
        None

    Returns:
        pandas.DataFrame: DataFrame containing any clinical metadata not 
            associated with a specific sequence file.                    
    """
    sample_ids = [sid[1:3] + '-' + sid[3:] for sid in metadata_df['External ID']]
    clinical_noseq_df = clinical_df[-clinical_df['st_q4'].isin(sample_ids)]

    return clinical_noseq_df


def main(args):
    config = parse_cfg_file(args.config)  

    study_trax_df = pd.read_csv(args.studytrax_metadata)
    broad_sample_df = pd.read_csv(args.broad_sample_tracking,
                                  na_values=['destroyed', 'missed'],
                                  parse_dates=['Actual Date of Receipt'])
    proteomics_df = None
    metadata_df = None
    
    date_today = datetime.date.today()
    metadata_file = os.path.join(args.output_dir, 
                                 'hmp2_metadata_%s.csv' % date_today)
    no_seq_metadata = os.path.join(args.output_dir,
                                   'hmp2_metadata_no_sequence_%s.csv' % date_today)

    ## Before we filter our metadata rows down to just to rows associated
    ## with the files we have present, we'll want a list of all the collection
    ## dates
    collection_dates_dict = get_collection_dates(broad_sample_df)

    if args.proteomics_metadata:
        proteomics_df = pd.read_table(args.proteomics_metadata)

    ## The update procedure either assumes that we have an exisitng metadata
    ## file that we are going to be appending too/updating or that we are 
    ## creating a fresh metadata sheet and will be adding the files in the 
    ## manifest file too it.
    ## TODO: This needs to be re-worked to account for snagging datatypes as as well.
    #if not args.metadata_file or args.refresh_all:
    #    sequence_files.extend(get_all_sequence_files(config.get('deposition_dir'),
    #                                                 config.get('input_extensions')))
    if args.manifest_file:
        manifest = parse_cfg_file(args.manifest_file)
        submitted_files = manifest.get('submitted_files')

        if submitted_files:
            new_metadata = []
            for (dtype, items) in submitted_files.iteritems():
                input_files = items.get('input')
                pair_identifier = items.get('pair_identifier')

                if args.pair_identifier:
                    (input_pair1, input_pair2) = bb_utils.paired_files(input_files,
                                                                       pair_identifier)
                    input_files = input_pair1 if input_pair1 else input_files                                                                       

                new_metadata.append(get_metadata_rows(config,
                                                      study_trax_df, 
                                                      broad_sample_df, 
                                                      proteomics_df,
                                                      dtype,
                                                      input_files,
                                                      args.pair_identifier))
 
            new_metadata_df = pd.concat(new_metadata)
            new_metadata_df['External ID'] = new_metadata_df.apply(generate_external_id, axis=1)

            new_metadata_df['Site/Sub/Coll ID'] = new_metadata_df['Site/Sub/Coll'].map(lambda sid: str(sid))
            new_metadata_df['Site'] = new_metadata_df['SiteName']
            new_metadata_df['Participant ID'] = new_metadata_df['Subject'].map(lambda subj: 'C' + str(subj))
            new_metadata_df['visit_num'] = new_metadata_df['Collection #']
            new_metadata_df['Research Project'] = "ibdmdb"
            new_metadata_df['Project'] = new_metadata_df.apply(get_project_id, axis=1)
            new_metadata_df = generate_collection_statistics(new_metadata_df, 
                                                         collection_dates_dict)

            new_metadata_df = remove_columns(new_metadata_df, config.get('drop_cols'))

    if args.metadata_file:
        metadata_df = pd.read_csv(args.metadata_file)
    
        site_mapping = config.get('site_map')
        metadata_df['Site/Sub/Coll ID'] = metadata_df.apply(fix_site_sub_coll_id, 
                                                            args=(site_mapping,),
                                                            axis=1)
        metadata_df['PDO Number'] = metadata_df.apply(get_pdo_number, axis=1)

        if not new_metadata_df.empty:
            metadata_df = pd.concat([metadata_df, new_metadata_df])
            metadata_df = metadata_df.drop_duplicates(subset=['Site/Sub/Coll ID', 'data_type'], keep='last')
    else:
        metadata_df = new_metadata_df

    metadata_df['visit_num'] = metadata_df.apply(fill_visit_nums, axis=1)
    metadata_df = reorder_columns(metadata_df, config.get('col_order'))
    metadata_df.to_csv(metadata_file, index=False)
    
    no_seq_files_metadata_df = get_no_sequence_metadata(metadata_df, study_trax_df)
    no_seq_files_metadata_df.to_csv(no_seq_metadata, index=False)


if __name__ == "__main__":
    main(parse_cli_arguments())
