# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils.metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A collection of functions related around metadata processing and generation 
related to the HMP2 AnADAMA2 workflows.

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

import numpy as np
import pandas as pd


def get_collection_dates(broad_sample_df):
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
    def _add_previous_collection_date(dataframe):
        """For each row in a dataframe takes the previous date of reception
        and adds it as a new column to aid in computation of the week_num and
        interval_days columns in our final metadata table.
        """
        dataframe['prev_coll_date'] = dataframe['Actual Date of Receipt'].shift()
        return dataframe

    collection_dict = dict((subj, df) for (subj, df) in 
                           broad_sample_df.ix[:,'Subject':'Actual Date of Receipt']
                           .sort_values(by='Actual Date of Receipt').groupby('Subject'))
    collection_dict = dict((subj, _add_previous_collection_date(df))
                            for (subj, df) in collection_dict.iteritems())

    return collection_dict


def add_proteomics_metadata(broad_subset_df, proteomics_metadata, sample_map):
    """Retrieves any corresponding proteomics metadata that should be incldued
    if generating metadata rows for proteomics data.

    Args:
        broad_subset_df (pandas.DataFrame): Broad sample tracking metadata
        proteomics_metadata (string): Proteomics metadata tab-delimited file.
        sample_map (dict): A mapping of filename to a sample ID that can be 
            found in our metadata table.

    Requires:
        None

    Returns:
        pandas.DataFrame: A dataframe containing any corresponding proteomics
            metadata.        
    """
    proteomics_df = pd.read_table(proteomics_metadata)

    sample_filter_df = broad_subset_df[broad_subset_df['Proteomics status'] == 'EXPORTED']
    sample_filter_df = sample_filter_df[['Parent Sample A', 'Proteomics']]

    proteomics_df['sample_ids'] = proteomics_df['Dataset'].replace(sample_map)
    proteomics_df['PDO Number'] = proteomics_df['Dataset'].map(lambda did: did.replace('-', '_')
                                                                                .split('_')[0])
    proteomics_df = sample_filter_df.merge(proteomics_df,
                                            left_on='Proteomics',
                                            right_on='sample_ids',
                                            how='right')
    proteomics_df = proteomics_df.drop('Proteomics', 1)

    return proteomics_df


def add_auxillary_metadata(metadata_df, auxillary_metadata):
    """Supplements the provided DataFrame containing metadata rows
    with auxillary metadata. Any auxillary method must begin with the 
    following site_sub_coll and data_type columns as these act as keys to map
    back to existing metadata rows.

    Args:
        metadata_df (pandas.DataFrame): A DataFrame containing metadata. 
        auxillary_metadata (list): A list of files containing metadata to 
            be appended to the provided metadata table.

    Requires:
        None

    Returns:
        pandas.DataFrame: A dataframe updated with the provided auxillary
            metadata.                    
    """
    for aux_metadata_file in auxillary_metadata:
        aux_metadata_df = pd.read_table(aux_metadata_file)

    return metadata_df


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


def merge_metadata_files(metadata_files, output_file):
    """Merges two or more metadata files together with logic to update 
    any rows and remove duplicate rows.

    Args:
        metadata_files (list): A list containing all files to merge together.
        output_file (string): Path to merged metadata file

    Requires:
        None

    Returns:
        string: Path to merged metadata file.                
    """
    if len(metadata_files) > 1:
        merged_metadata_df = pd.concat(metadata_files)
  
        # If we have duplicate rows we currently are going to drop the last
        # set of rows until update logic is put in place
        merged_metadata_df = merged_metadata_df.drop_duplicates(subset=['Site/Sub/Coll ID', 'data_type'], 
                                                                keep='last')
        merged_metadata_df.to_csv(output_file, index=False)
    else:
        output_file = metadata_files[0]

    return output_file