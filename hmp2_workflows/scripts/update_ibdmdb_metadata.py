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

import sys

import argparse
import collections
import datetime
import math
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
    parser.add_argument('-d', '--biopsy-dates', 
                        help='Spreadsheet containg metadata for any 16S'
                        'biopsy data.')
    parser.add_argument('-as', '--add-all-stool-collections', required=False,
                        action='store_true', default=False,
                        help='OPTIONAL. Add in all stool '
                        'collections to the metadata file regardless of whether  '
                        'or not a data product was generated.')
    parser.add_argument('-a', '--auxillary-metadata', action='append',
                        default=[], help='Additional auxillary metadata '
                        'to use in populating the HMP2 metadata table.')
    
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


def parse_biopsy_dates(biopsy_dates):
    """Parses the supplementary Studytrax biopsy dates to allow for 
    attaching dates to any samples of type 'Screening Colonoscopy' or 
    'Additional Biopsy'

    Args:
        biopsy_dates (string): The path to the file containing biopsy dates.

    Requires:
        None

    Returns:
        dict: A dictionary containing biopsy dates keyed on subject ID and 
            sample type (i.e. Screening Colonoscopy)
    """
    biopsy_dict = {}

    for line in open(biopsy_dates):
        if line.startswith('#'):
            continue

        week_num = ""
        (subj_id, days, interval) = line.strip().split(',')
        if days != "Days Unknown" and not days == "":
            week_num = int(math.floor(int(days)/7))

        biopsy_dict.setdefault(subj_id, {})
        biopsy_dict[subj_id][interval] = week_num

    return biopsy_dict


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
    type_mapping = {'host_transcriptomics': 'HTX',
                    'biopsy_16S': 'BP',
                    'metatranscriptomics': 'MTX',
                    'metagenomics': 'MGX',
                    'viromics': 'MVX',
                    'host_genome': 'HG',
                    'methylome': 'RRBS',
                    'serology': 'SER' }


    ## This specific case is applicable to Proteomics data only but the 
    ## function can be expanded to handle other scenarios in the future
    if pd.isnull(project_id) and not  pd.isnull(row.get('Job')):
        project_id = row['Job']
    else:
        project_id = row.get('Site/Sub/Coll') + '_' + type_mapping.get(row.get('data_type'))

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
    elif row.get('data_type') == 'amplicon' and 'PDO' not in pdo_num:
        pdo_num = "PDO-%s" % pdo_num

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


def get_biopsy_site_sub_coll(row):
    subj_id = row.get('ProjectSpecificID')
    interval_name = row.get('IntervalName')

    if interval_name == "Screening Colonoscopy":
        collection_num = "S"
    elif interval_name == "Additional Biopsy":
        collection_num = "B"
    
    site_sub_coll = subj_id + "C" + collection_num

    if studytrax_df['Site/Sub/Coll'].isin(site_sub_coll):
        coll_num = studytrax_df[studytrax_df['Site/Sub/Coll ID'] == site_sub_coll][-1]  
        new_coll_num = int(coll_num) + 1

    return site_sub_coll


def _get_non_stool_site_sub_coll(row):
    """For any non-stool collection StudyTrax entries generates Site/Sub/Coll
    ID's by using a combination of Project Specific ID, Site Name and 
    IntervalName.
    """
    interval_name = row.get('IntervalName')
    subj_id = row.get('ProjectSpecificID')
    site_name = row.get('SiteName')

    site_mapping = {'Cincinnati': 'H',
                    'Massachusetts General Hospital': 'M',
                    'Emory': 'E',
                    'MGH Pediatrics': 'P',
                    'Cedars-Sinai': 'C'}
    interval_mapping = {'Screening Colonoscopy': 'SC',
                        'Additional Biopsy': 'B',
                        'Baseline (IBD and Healthy)': 'BL'}

    coll_num = "1"
    if "follow-up" in interval_name.lower():
        coll_num = interval_name.lower().replace('follow-up (month ', '').replace(')','')
        interval_name_recode = "FU"
    else:
        interval_name_recode = interval_mapping.get(interval_name)

    return site_mapping.get(site_name) + subj_id + 'C' + interval_name_recode + coll_num


def resolve_dupe_ssc_ids(metadata_df):
    prev_ssc = None
    counter = 1

    ids = metadata_df['Site/Sub/Coll']
    for (idx, row) in metadata_df[ids.isin(ids[ids.duplicated()])].sort_values("Site/Sub/Coll").iterrows():
        ssc_id = row.get('Site/Sub/Coll')
        if ssc_id != prev_ssc:
            prev_ssc = ssc_id
            counter = 1
            continue

        metadata_df.loc[idx, 'Site/Sub/Coll'] = ssc_id[:-1] + str(counter + 1)
        counter = counter + 1
    
    return metadata_df


def add_all_stool_collections(metadata_df, studytrax_df, broad_df):
    """Adds any of the missing stool sample collections for which a product was not
    generated.

    Args:
        metadata_df (pandas.DataFrame): All existing metadata in a pandas DataFrame.
        studytrax_df (pandas.DataFrame): StudyTrax clinical metadata.
        broad_df (pandas.DataFrame): Broad sample spreadsheet metadata.

    Requires:
        None

    Returns:
        pandas.DataFrame: An updated dataframe containing all stool samples that 
        haven't been associated with a data type.
    """
    site_mapping = {'Cincinnati': 'H',
                    'Massachusetts General Hospital': 'M',
                    'Emory': 'E',
                    'MGH Pediatrics': 'P',
                    'Cedars-Sinai': 'C'}

    def _gen_site_sub_coll(row):
        site_abbrev = site_mapping.get(row.get('SiteName'))
        coll_num = row.get('IntervalName').replace('Stool Collection #', '')
        subject_id = row.get('ProjectSpecificID')
        return site_abbrev + subject_id + "C" + coll_num


    ## This is a temporary hack that allows us to loop in all the stool 
    ## samples that were received but did not have a corresponding data point

    ## In order to collect these we need to create a new temporary row that 
    ## is going to be similar to our Site/Sub/Coll ID's
    metadata_stool_df = metadata_df[metadata_df['IntervalName'].str.startswith('Stool')]
    studytrax_stool_df = studytrax_df[studytrax_df['IntervalName'].str.startswith('Stool')]
    broad_subset_df = broad_df.filter(['Site/Sub/Coll', 'Actual Date of Receipt'])

    tmp_stool_ids = ["%s_%s" % (x,y) for (x,y) in zip(metadata_stool_df['ProjectSpecificID'],
                                                      metadata_stool_df['IntervalSequence'])]
    tmp_stool_ids = list(set(tmp_stool_ids))

    studytrax_stool_df['stool_id'] = studytrax_stool_df.apply(lambda row: "%s_%s" % (row['ProjectSpecificID'],
                                                                                     row['IntervalSequence']), axis=1)
    studytrax_stool_df['Site/Sub/Coll ID'] = studytrax_stool_df.apply(_gen_site_sub_coll, axis=1)
    studytrax_noprod_df = studytrax_stool_df[-studytrax_stool_df['stool_id'].isin(tmp_stool_ids)]
    studytrax_noprod_df['data_type'] = "noproduct"
    studytrax_noprod_df['ProjectSpecificID'].astype('int')
    studytrax_noprod_df = studytrax_noprod_df.merge(broad_subset_df, left_on='Site/Sub/Coll ID', right_on='Site/Sub/Coll', how='left')
    studytrax_noprod_df['Site'] = studytrax_noprod_df['SiteName']

    studytrax_noprod_df = studytrax_noprod_df.drop('stool_id', 1)
    studytrax_noprod_df = studytrax_noprod_df.drop('Site/Sub/Coll', 1)
    
    metadata_df = pd.concat([metadata_df, studytrax_noprod_df], ignore_index=True)

    return metadata_df


def _clean_blood_sample_ids(sample_id):
    """Cleans up some of the funny looking sample IDs in the bl_q5 column 
    used to map serology samples to the StudyTrax metadata

    Args:
        sample_id (string): A sample ID in the Stuydtrax bl_q5 column to 
            clean-up

    Requires: None

    Returns: 
        string: Cleaned sample ID.            
    """
    if not pd.isnull(sample_id):
        sample_id = (sample_id.replace(' 1', '')
                                .replace('-1', '')
                                .replace('s1', '')
                                .replace('S1', '')
                                .replace('.1', ''))

    return sample_id


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

    if pair_identifier:
       sample_ids = [sid.replace(pair_identifier, '') for sid in sample_ids]

    sample_ids_techreps = [sid for (k, sid) in sample_mapping.iteritems() if "techrep" in k]
    sample_ids = set(sample_ids) - set(sample_ids_techreps)

    data_type_mapping = config.get('dtype_mapping')

    ## Grab subset of Broad sample tracking spreadsheet
    sample_subset_df = pd.DataFrame()
    if data_type != "HTX":
        sample_subset_df = sample_df[(sample_df['Parent Sample A'].isin(sample_ids)) |
                                    (sample_df['Proteomics'].isin(sample_ids)) |
                                    (sample_df['MbX'].isin(sample_ids)) |
                                    (sample_df['Viromics'].isin(sample_ids)) |
                                    (sample_df['Site/Sub/Coll']).isin(sample_ids)]

    if sample_ids_techreps:
        sample_ids_techreps = [sample_id.replace('_techrep', '') for sample_id 
                               in sample_ids_techreps]
        sample_subset_techreps = sample_df[(sample_df['Parent Sample A'].isin(sample_ids_techreps)) |
                                 (sample_df['Proteomics'].isin(sample_ids_techreps)) |
                                 (sample_df['MbX'].isin(sample_ids_techreps)) |
                                 (sample_df['Viromics'].isin(sample_ids_techreps)) |
                                 (sample_df['Site/Sub/Coll']).isin(sample_ids_techreps)]

        sample_subset_techreps['External ID'] = sample_subset_techreps['Parent Sample A'].map(lambda sid: sid.replace('-', '') + "_TR")
        sample_subset_techreps['External ID'] = sample_subset_techreps.apply(lambda row: row.get('Site/Sub/Coll')[0] + 
                                                                                         row.get('External ID'), axis=1)
        sample_subset_df = pd.concat([sample_subset_df, sample_subset_techreps],
                                     ignore_index=True)

    ## TODO: Figure out if we have any samples that did not have aassociated metadata
    #join_how = 'outer' if full_join else 'left'
    if len(sample_subset_df) == 0:
        other_loc_map = {'0': 'Terminal ileum',
                         '1': 'Neo-ileum',
                         '2': 'Ileocecal Valve',
                         '3': 'Cecum',
                         '4': 'Ascending (right-sided) colon',
                         '5': 'Transverse colon',
                         '6': 'Descending (left-sided) colon',
                         '7': 'Sigmoid Colon',
                         '8': 'Rectum'}

        if data_type == "HTX" or data_type == "RRBS":
        # TODO: Add these to config file
            biopsy_map = {'bx_q5': 'Rectum',
                          'bx_q6': 'Ileum',
                          'bx_q7': 'Other Inflamed',
                          'bx_q9': 'Non-inflamed'}

            new_meta_dfs = []
            for (studytrax_col, location) in biopsy_map.iteritems():
                new_meta_df = studytrax_df[studytrax_df[studytrax_col].isin(sample_ids)]
                new_meta_df['biopsy_location'] = location
                new_meta_df['External ID'] = new_meta_df[studytrax_col].map(lambda sid: sid.replace('-', ''))

                if not new_meta_df['biopsy_location'].empty:
                    if location == "Other Inflamed":
                        new_loc_col = "bx_q8"
                    elif location == "Non-inflamed":
                        new_loc_col = "bx_q10"
                    if not location in ['Rectum', 'Ileum']:
                        new_meta_df['biopsy_location'] = [other_loc_map.get(x) for x in new_meta_df[new_loc_col]]

                new_meta_dfs.append(new_meta_df)

            if data_type == "RRBS":
                blood_df = studytrax_df[studytrax_df['bl_q4'].isin(sample_ids)]
                blood_df['External ID'] = blood_df['bl_q4'].map(lambda sid: sid.replace('-', ''))
                new_meta_dfs.append(blood_df)

            metadata_df = pd.concat([sample_subset_df] + new_meta_dfs, ignore_index=True)
            metadata_df = metadata_df.drop_duplicates(subset=['External ID', 'biopsy_location'], keep='first')
            metadata_df['Site/Sub/Coll'] = metadata_df.apply(_get_non_stool_site_sub_coll, axis=1)
            resolve_dupe_ssc_ids(metadata_df)
        elif data_type == "HG":
            studytrax_col = "bl_q4"
            blood_df = studytrax_df[studytrax_df[studytrax_col].isin(sample_ids)]
            blood_df['External ID'] = blood_df[studytrax_col].map(lambda sid: sid.replace('-', ''))
 
            metadata_df = pd.concat([sample_subset_df, blood_df],
                                    ignore_index=True)
 
            metadata_df['Site/Sub/Coll'] = metadata_df.apply(_get_non_stool_site_sub_coll, axis=1)
            resolve_dupe_ssc_ids(metadata_df)
        elif data_type == "SER":
            new_metadata_dfs = []
            studytrax_cols = {'bl_q5': None}

            # We get some really weird stuff going on in the studytrax mapping column here so let's 
            # clean things up first
            studytrax_df['bl_q5'] = studytrax_df['bl_q5'].map(_clean_blood_sample_ids)

            for (col, label) in studytrax_cols.iteritems():
                new_metadata_df = studytrax_df[studytrax_df[col].isin(sample_ids)]
                new_metadata_df['External ID'] = new_metadata_df[col]
                new_metadata_dfs.append(new_metadata_df)

            metadata_df = pd.concat([sample_subset_df] + new_metadata_dfs, ignore_index=True)
            metadata_df['Site/Sub/Coll'] = metadata_df.apply(_get_non_stool_site_sub_coll, axis=1)
            resolve_dupe_ssc_ids(metadata_df)
        elif data_type == "16SBP":
            biopsy_map = {'bx_q13': 'Rectum',
                          'bx_q14': 'Ileum',
                          'bx_q15': 'Other Inflamed',
                          'bx_q17': 'Non-inflamed'}

            biopsy_dfs = []
            for (studytrax_col, location) in biopsy_map.iteritems():
                biopsy_df = studytrax_df[studytrax_df[studytrax_col].isin(sample_ids)]
                biopsy_df['biopsy_location'] = location
                biopsy_df['External ID'] = biopsy_df[studytrax_col].map(lambda sid: sid.replace('-', ''))

                if not biopsy_df['biopsy_location'].empty:
                    if location == "Other Inflamed":
                        new_loc_col = "bx_q16"
                    elif location == "Non-inflamed":
                        new_loc_col = "bx_q18"
                        
                    if not location in ['Rectum', 'Ileum']:
                        biopsy_df['biopsy_location'] = [other_loc_map.get(x) for x in biopsy_df[new_loc_col]]

                biopsy_dfs.append(biopsy_df)

            metadata_df = pd.concat([sample_subset_df] + biopsy_dfs, ignore_index=True)
            metadata_df.drop_duplicates(subset=['External ID', 'biopsy_location'], keep='first')
            metadata_df['Site/Sub/Coll'] = metadata_df.apply(_get_non_stool_site_sub_coll, axis=1)
            resolve_dupe_ssc_ids(metadata_df)

    else:
        metadata_df = sample_subset_df.merge(studytrax_df,
                                             left_on='Parent Sample A',
                                             right_on='st_q4',
                                             how='left')

        ## We sometimes get a situation where our studytrax metadata is missing 
        ## some of the proteomics sample ID's so we need to make sure we 
        ## replicate them.
        metadata_df.loc[metadata_df['st_q17'].isnull(), 'st_q17'] = metadata_df['Proteomics']
        metadata_df.loc[metadata_df['st_q11'].isnull(), 'st_q11'] = metadata_df['MbX']
        metadata_df.loc[metadata_df['st_q12'].isnull(), 'st_q12'] = metadata_df['Viromics']

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
            metadata_df['External ID'] = None                                            
    
    ## Now if we have techreps in our samples we need to add them in.
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
    metadata_cols = metadata_df.columns.tolist()
    to_drop = set(metadata_cols).intersection(set(drop_cols))

    return metadata_df.drop(to_drop, 1)


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
    blood_id = row['bl_q4']
    external_id  = row.get('External ID')
    base_id = None

    if row.get('data_type') != "noproduct":
        if not pd.isnull(external_id):
            base_id = external_id
        elif not pd.isnull(stool_id):
            base_id = stool_id
        elif not pd.isnull(blood_id):
            base_id = blood_id
        else:
            raise Exception("Could not generate External ID:", row)

        site_sub_coll = row['Site/Sub/Coll ID']
        row['External ID'] = site_sub_coll[0] + base_id.replace('-', '')
        #eturn site_sub_coll[0] + base_id.replace('-', '')
        
    return row


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
    data_type  = row['data_type']

    if pd.isnull(visit_num) and data_type not in ['host_transcriptomics', 'host_genome',
                                                  'biopsy_16S', 'methylome']:
        visit_num = site_sub_coll_id.split('C')[-1]
    
    return visit_num


def generate_collection_statistics(metadata_df, collection_dict, biopsy_dates=None):
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
        data_type = row['data_type']

        ## Sticking this in here for the time being...
        participant_id = row['Site/Sub/Coll ID'][:5]
        row['Participant ID'] = participant_id
        subject_id = int(row['Participant ID'][1:])
        
        if data_type in ['host_transcriptomics', 'biopsy_16S', 'methylome']:
            interval_name = row['IntervalName']
            
            if "Baseline" in interval_name:
                row['week_num'] = "0"
            elif not  "Follow" in interval_name:
                row['week_num'] = biopsy_dates[str(subject_id)][interval_name]
        elif data_type != 'host_genome':
            visit_num = row['visit_num']
            receipt_date = row['Actual Date of Receipt']

            ## Sticking this in here for the time being...
            participant_id = row['Site/Sub/Coll ID'].rsplit('C', 1)[0]
            row['Participant ID'] = participant_id

            if pd.isnull(row['Site']):
                row['Site'] = row['SiteName']
            elif pd.isnull(row['SiteName']):
                row['SiteName'] = row['Site']

            if pd.isnull(row['week_num']) and not pd.isnull(receipt_date):
                #print row
                collection_dates = collection_dict.get(subject_id)
                #print "DEBUG:", collection_dates, visit_num
                initial_visit_date = collection_dates.iloc[0]['Actual Date of Receipt']
                subj_collection_row = collection_dates[collection_dates['Collection #']
                                                       == int(visit_num)]
                prev_visit_date = subj_collection_row['prev_coll_date'].values[0]

                if pd.isnull(prev_visit_date):
                    prev_visit_date = initial_visit_date

                row['week_num'] = (receipt_date - initial_visit_date).days / 7
                row['interval_days'] = (receipt_date - prev_visit_date).days

        return row

    metadata_df = metadata_df.apply(_add_collection_columns, axis=1)
    metadata_df = add_biopsy_visit_num(metadata_df)

    return metadata_df


def add_biopsy_visit_num(metadata_df):
    """
    """
    for (idx, row) in metadata_df[metadata_df['data_type'].isin(['host_transcriptomics', 'biopsy_16S', 'methylome'])].iterrows():
        week_num = row.get('week_num')
        visit_num = row.get('visit_num')
        participant_id = row.get('Participant ID')

        if not pd.isnull(visit_num):
            continue

        if pd.isnull(week_num) or week_num == '':
            visit_num = 1
        else:
            ## Do we have an exact match somewhere?
            match = metadata_df.loc[(metadata_df['Participant ID'] == participant_id) &
                                    (metadata_df['week_num'] == week_num) &
                                    (-metadata_df['data_type'].isin(['host_transcriptomics', 'biopsy_16S', 'methylome']))]

            if not match.empty:
                visit_num = match.get('visit_num').values[0]
            else:
                ## If not then we are going to try to get as close as possible here...
                subject_rows = metadata_df.loc[(metadata_df['Participant ID'] == participant_id) &
                                               (-metadata_df['data_type'].isin(['host_transcriptomics', 'biopsy_16S', 'methylome']))]
                subject_rows['week_num'] = pd.to_numeric(subject_rows['week_num'])
                matches = subject_rows.iloc[(subject_rows['week_num']-float(week_num)).abs().argsort()[:2]]

                if matches.empty:
                    visit_num = 1
                else:
                    visit_num = matches['visit_num'].values[0]

        row['visit_num'] = visit_num
        metadata_df.loc[idx] = row
    
    return metadata_df


def fix_site_name(metadata_nosite_df):
    """
    """
    site_mapping = {'H': 'Cincinnati',
                    'M': 'Massachusetts General Hospital',
                    'E': 'Emory',
                    'P': 'MGH Pediatrics',
                    'C': 'Cedars-Sinai'}

    def _fix_site_name(row):
        site_sub_coll = row.get('Site/Sub/Coll ID')
        site_abbrev = site_sub_coll[0]
        row['SiteName'] = site_mapping.get(site_abbrev)
        row['Site'] = row['SiteName']
        return row

    metadata_nosite_df = metadata_nosite_df.apply(_fix_site_name, axis=1)

    return metadata_nosite_df


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


def add_baseline_metadata_values(metadata_df, studytrax_df, baseline_columns):
    """Retrieves metadata from the Studytrax clinical metadata that are 
    applicable and useful to specific stool data.

    Args:
        metadata_df (pandas.DataFrame): HMP2 metadata table stored in a 
            pandas DataFrame.
        studytrax_df (pandas.DataFrame): StudyTrax clinical metadata.
        baseline_columns (list): A list of columns to pull target baseline
            metadata columns from.

    Requires:
        None

    Returns:
        pandas.DataFrame: DataFrame containing metadata with baseline columns
            filled in.                    
    """
    metadata_df.set_index('ProjectSpecificID', inplace=True)

    for col in baseline_columns:
        baseline_df = studytrax_df[(studytrax_df['IntervalName'] == 'Baseline (IBD and Healthy)') &
                                   (studytrax_df[col].notnull())].filter(['ProjectSpecificID', col])
        baseline_df[col].apply(lambda x: x.split("(")[0])
        baseline_df['ProjectSpecificID'] = pd.to_numeric(baseline_df['ProjectSpecificID'])
        baseline_df.set_index('ProjectSpecificID', inplace=True)
        metadata_df.update(baseline_df)

    metadata_df.reset_index(inplace=True)
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

    study_trax_df = pd.read_csv(args.studytrax_metadata, dtype='str')
    broad_sample_df = pd.read_csv(args.broad_sample_tracking,
                                  na_values=['destroyed', 'missed'],
                                  parse_dates=['Actual Date of Receipt'])
    proteomics_df = None
    metadata_df = None
    new_metadata_df = None

    date_today = datetime.date.today()
    metadata_file = os.path.join(args.output_dir,
                                 'hmp2_metadata_%s.csv' % date_today)

    ## Before we filter our metadata rows down to just to rows associated
    ## with the files we have present, we'll want a list of all the collection
    ## dates
    collection_dates_dict = get_collection_dates(broad_sample_df)

    biopsy_date_map = None
    if args.proteomics_metadata:
        proteomics_df = pd.read_table(args.proteomics_metadata)
    if args.biopsy_dates:
        biopsy_date_map = parse_biopsy_dates(args.biopsy_dates)

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

                if pair_identifier:
                    (input_pair1, input_pair2) = bb_utils.paired_files(input_files,
                                                                       pair_identifier)
                    input_files = input_pair1 if input_pair1 else input_files

                else:
                    new_metadata.append(get_metadata_rows(config,
                                                          study_trax_df, 
                                                          broad_sample_df, 
                                                          proteomics_df,
                                                          dtype,
                                                          input_files,
                                                          pair_identifier))
 
            new_metadata_df = pd.concat(new_metadata, ignore_index=True)

            #new_metadata_df[new_metadata_df['External ID'].isnull()] = None
            new_metadata_df['Site/Sub/Coll ID'] = new_metadata_df['Site/Sub/Coll'].map(lambda sid: str(sid))
            #new_metadata_df['Participant ID'] = new_metadata_df['Subject'].map(lambda subj: 'C' + str(subj))
            if 'Collection #' in new_metadata_df.columns:
                new_metadata_df['visit_num'] = new_metadata_df['Collection #']
            new_metadata_df['Project'] = new_metadata_df.apply(get_project_id, axis=1)
            new_metadata_df['ProjectSpecificID'] = pd.to_numeric(new_metadata_df['ProjectSpecificID'])
            new_metadata_df['Site'] = new_metadata_df['SiteName']
            new_metadata_df = new_metadata_df.apply(generate_external_id, axis=1)

            new_metadata_df = remove_columns(new_metadata_df, config.get('drop_cols'))

    if args.metadata_file:
        metadata_df = pd.read_csv(args.metadata_file, parse_dates=['Actual Date of Receipt'])
    
        site_mapping = config.get('site_map')
        metadata_df['Site/Sub/Coll ID'] = metadata_df.apply(fix_site_sub_coll_id,
                                                            args=(site_mapping,),
                                                            axis=1)
        metadata_df['PDO Number'] = metadata_df.apply(get_pdo_number, axis=1)

        if new_metadata_df and not new_metadata_df.empty:
            metadata_df = pd.concat([metadata_df, new_metadata_df], ignore_index=True)
            metadata_df = metadata_df.drop_duplicates(subset=['External ID', 'Site/Sub/Coll ID', 'data_type'], keep='last')
    else:
        metadata_df = new_metadata_df

    metadata_df[metadata_df['External ID'].isnull()] = metadata_df[metadata_df['External ID'].isnull()].apply(generate_external_id, axis=1)

    if args.auxillary_metadata:
        for aux_file in args.auxillary_metadata:
            supp_df = pd.read_table(aux_file)
            supp_columns = supp_df.columns.tolist()

            idx_offset = 1
            if 'data_type' in supp_columns:
                join_id = supp_columns[:1] + ['data_type']
                idx_offset = 2
            else:
               join_id = supp_columns[0]

            ## We need to do this in two stages. If the columns already exist
            ## here we want to update them. If they do not exist we append
            ## them.
            metadata_cols = metadata_df.columns.tolist()
            new_cols = set(supp_columns[idx_offset:]) - set(metadata_cols)
            existing_cols = set(supp_columns[idx_offset:]).intersection(metadata_cols)

            if new_cols:
                supp_new_df = supp_df.filter(items=supp_columns[:idx_offset] + list(new_cols))
                metadata_df = metadata_df.merge(supp_new_df, how='left', on=join_id)

            if existing_cols:
                supp_existing_df = supp_df.filter(items=supp_columns[:idx_offset] + list(existing_cols))
                metadata_df.set_index(join_id, inplace=True)
                supp_existing_df.set_index(join_id, inplace=True)

                metadata_df.update(supp_existing_df)
                metadata_df.reset_index(inplace=True)

    if args.add_all_stool_collections:
        metadata_df = add_all_stool_collections(metadata_df, study_trax_df, broad_sample_df)

    metadata_df['Actual Date of Receipt'] = pd.to_datetime(metadata_df['Actual Date of Receipt'])
    metadata_df['visit_num'] = metadata_df.apply(fill_visit_nums, axis=1)

    metadata_df['hbi_score'] = pd.to_numeric(metadata_df['hbi_score'])
    if 'Site' in metadata_df.columns.tolist():
        metadata_df['SiteName'] = metadata_df['Site']
    else:
        metadata_df['Site'] = metadata_df['SiteName']

    ## Couple small remaining changes
    metadata_df.ix[metadata_df.hbi_score > 900, 'hbi_score'] = None 
    metadata_df.ix[metadata_df.consent_age > 150, 'consent_age'] = None 
    metadata_df['total_reads'].loc[metadata_df['total_reads'].astype('str').str.startswith('PDO')] = None 
    metadata_df['Research Project'] = "ibdmdb"

    metadata_df = generate_collection_statistics(metadata_df, collection_dates_dict, biopsy_date_map)
    metadata_df = add_baseline_metadata_values(metadata_df, study_trax_df, config.get('baseline_cols'))

    metadata_df[metadata_df['SiteName'].isnull()] = fix_site_name(metadata_df[metadata_df['SiteName'].isnull()])
    metadata_df = reorder_columns(metadata_df, config.get('col_order'))
    metadata_df.drop(['Site'], 1, inplace=True)

    metadata_df = metadata_df.sort_values(['data_type', 'Participant ID', 'visit_num'])
    metadata_df.to_csv(metadata_file, index=False)


if __name__ == "__main__":
    main(parse_cli_arguments())
