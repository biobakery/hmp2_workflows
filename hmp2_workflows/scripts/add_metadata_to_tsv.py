# -*- coding: utf-8 -*-

"""
add_metadata_to_tsv.py
~~~~~~~~~~~~~~~~~~~~~~

A stand-alone script to add metadata to an existing tab-delimited or 
comma-separated analysis file.

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

# Running on Harvard resouces we run into some issues so giving this hacky
# solution a shot.

import sys

import argparse
import funcy
import os

import pandas as pd

import biobakery_workflows.utilities as bb_utils

from numpy import nan as Nan

from hmp2_workflows.utils.misc import (get_sample_id_from_fname,
                                       parse_cfg_file,
                                       reset_column_headers)


def parse_cli_arguments():
    """
    Parses any command-line arguments passed into this script.

    Args:
        None

    Requires:
        None

    Returns:
        argparse.ArgumentParser: Object containing all command-line args.
    """
    parser = argparse.ArgumentParser('Adds metadata to the provided analysis '
                                     'file.')
    parser.add_argument('-i', '--input-file', required=True, action='append',
                        help='Input analysis file (tab or comm delimited).')
    parser.add_argument('-m', '--metadata-file', required=True,
                        help='HMP2 metadata file.')
    parser.add_argument('-c', '--config-file', required=True,
                        help='Config file containing project config parameters.')
    parser.add_argument('-d', '--data-type', choices=['MGX', 'MTX', 'HTX',
                                                      'MVX', 'MPX', 'MBX',
                                                      '16S', 'RRBS', 'SER',
                                                      '16SBP'],
                        required=True,
                        help='Data type of data used to generate analysis '
                        'file.')
    parser.add_argument('-t', '--id-col', required=True,
                        help='Target ID column to match up analysis samples '
                        'to samples found in the metadata.')
    parser.add_argument('--drop-missing-cols', action='store_true', default=False,
                        help='OPTIOANL. If set to True drop any columns that '
                        'do not map to our metadata.')                        
    parser.add_argument('-s', '--supplement', default=[], action='append',
                        help='Any supplementary metadata to be added. '
                        'should be keyed on the same ID provided in the '
                        '--id-col parameter.')


    return parser.parse_args()


def add_metadata_to_tsv(analysis_files, metadata_file,
                        data_type,
                        id_col, col_replace, 
                        drop_cols,
                        target_cols=[],
                        supplement=[]):
    """Adds metadata to the top of a tab-delimited file. This function is
    meant to be called on analysis files to append relevant metadata to the 
    analysis output found in the file. An example can be seen below:

        
        sample  Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
        Age 87  78  3   2   32  10  39  96
        Cohort  Healthy Healthy Healthy Healthy IBD IBD IBD IBD
        Favorite_color  Yellow  Blue    Green   Yellow  Green   Blue    Green 
        Height  60  72  63  67  71  65  61  64
        Sex 0   1   0   1   1   0   1   0
        Smoking 0   0   1   0   1   1   1   0
        Star_Trek_Fan   1   1   0   0   1   0   0   1
        Weight  151 258 195 172 202 210 139 140
        Bacteria    1   1   1   1   1   1   1   1
        Bacteria|Actinobacteria|Actinobacteria  0.0507585   0.252153    0.161725   

    Args:
        analysis_files (list): Target TSV's to add metadata too
        metadata_file (string): The path to the metadata file to pull from.
        id_col (string): The column name in the supplied metadata file to 
            attempt to subset on using ID's from the analysis file.
        col_replace (list): A list of string fragments that should be searched 
            for and replaced in either of the column headers of the analysis 
            or metadata files.
        target_cols (list): A list of columns to filter the metadata file on.
        drop_cols (boolean): If True drop any columns that don't map to our 
            metadata file.
        supplement (list): Any additional metadata files to integrate into 
            analysis files. 

    Requires:
        None

    Returns: 
        list: A list containing the path to all modified files.
    """
    metadata_df = pd.read_csv(metadata_file, dtype='str')

    col_offset = -1    
    metadata_rows = None
    na_rep = ""

    def _add_metadata_to_tsv(analysis_file, pcl_out):
        if analysis_file.endswith('.csv'):
            analysis_df = pd.read_csv(analysis_file, dtype='str', header=None)
        else:
            analysis_df = pd.read_table(analysis_file, dtype='str', header=None)

        pcl_metadata_df = None
        header = True
            
        # Going to make the assumption that the next row following our PCL 
        # metadata rows is the row containing the ID's that we will use to merge
        # the analysis file with our metadata file and we can use these same 
        # ID's to merge the PCL metadata rows into the larger metadata file.
        if metadata_rows:
            pcl_metadata_df = analysis_df[:metadata_rows+1]
            header = None

            offset_cols = range(0, col_offset+1)
            pcl_metadata_df.drop(pcl_metadata_df.columns[offset_cols[:-1]],
                                    axis=1,
                                    inplace=True)

            pcl_metadata_df = pcl_metadata_df.T.reset_index(drop=True).T
            pcl_metadata_df.xs(metadata_rows)[0] = id_col
            
            pcl_metadata_df = pcl_metadata_df.T
            pcl_metadata_df = reset_column_headers(pcl_metadata_df)

            analysis_df.drop(analysis_df.index[range(0,metadata_rows)], inplace=True)
            analysis_df.rename(columns=analysis_df.iloc[0], inplace=True)
        else:
            analysis_df = reset_column_headers(analysis_df)

        sample_ids = analysis_df.columns.tolist()[col_offset+1:]
            
        if len(sample_ids) == 1:
            raise ValueError('Could not parse sample ID\'s:',
                             sample_ids)

        if col_replace:
            new_ids = sample_ids
            for replace_str in col_replace:
                new_ids = [sid.replace(replace_str, '') if not pd.isnull(sid)
                           else sid for sid in new_ids]

            if new_ids != sample_ids:
                sample_ids_map = dict(zip(sample_ids, new_ids))
                sample_ids = new_ids
    
                analysis_df.rename(columns=sample_ids_map, inplace=True)

        subset_metadata_df = metadata_df[(metadata_df.data_type == data_type) &
                                         (metadata_df[id_col].isin(sample_ids))]
        mapping_sample_ids = subset_metadata_df[id_col].tolist()                                         

        if supplement:
            for aux_file in supplement:
                aux_metadata_df = pd.read_table(aux_file, dtype='str')
                aux_metadata_cols = aux_metadata_df.columns.tolist()
                join_id = aux_metadata_cols[0]

                aux_metadata_df = aux_metadata_df[aux_metadata_df[join_id].isin(sample_ids)]

                ## We need to do this in two stages. If the columns already exist
                ## here we want to update them. If they do not exist we append
                ## them.
                subset_metadata_cols = subset_metadata_df.columns.tolist()
                new_cols = set(aux_metadata_cols[1:]) - set(subset_metadata_cols)
                existing_cols = set(aux_metadata_cols[1:]).intersection(subset_metadata_cols)

                if new_cols:
                    aux_metadata_new_df = aux_metadata_df.filter(items=aux_metadata_cols[:1] +
                                                                 list(new_cols))
                    subset_metadata_df = pd.merge(subset_metadata_df, aux_metadata_new_df,
                                                  how='left', on=join_id)

                if existing_cols:
                    aux_metadata_existing_df = aux_metadata_df.filter(items=aux_metadata_cols[:1] +
                                                                      list(existing_cols))
                    subset_metadata_df.set_index(join_id, inplace=True)
                    aux_metadata_existing_df.set_index(join_id, inplace=True)

                    subset_metadata_df.update(aux_metadata_existing_df)
                    subset_metadata_df.reset_index(inplace=True)

        if pcl_metadata_df and not pcl_metadata_df.empty:
            subset_metadata_df = pd.merge(subset_metadata_df, pcl_metadata_df,
                                          how='left', on=id_col)

        if target_cols:
            target_cols.insert(0, id_col)
            subset_metadata_df = subset_metadata_df.filter(target_cols)

        subset_metadata_df = subset_metadata_df.drop_duplicates(subset=[id_col], keep='first')

        subset_metadata_df = subset_metadata_df.T
        subset_metadata_df = reset_column_headers(subset_metadata_df)
        subset_metadata_df = subset_metadata_df.reset_index()
        subset_metadata_df.fillna('NA', inplace=True)

        _col_offset = col_offset-1 if col_offset != -1 else col_offset
        col_name = analysis_df.columns[_col_offset+1]

        col_name = '' if col_name == "index" else col_name
        subset_metadata_df.rename(columns={'index': col_name}, inplace=True)

        analysis_df.index = analysis_df.index + len(subset_metadata_df.index)
        
        ## add an empty row to make things easier on the R loader scripts
        empty_df = pd.DataFrame([[""]*len(subset_metadata_df.columns)], columns=subset_metadata_df.columns)
        #subset_metadata_df = subset_metadata_df.append(empty_row, ignore_index=True)

        analysis_metadata_df = pd.concat([subset_metadata_df,
                                          analysis_df], axis=0)

        if drop_cols:
            analysis_metadata_df.filter(mapping_sample_ids)

        analysis_metadata_df.fillna('NA', inplace=True)
        analysis_metadata_df = pd.concat([analysis_metadata_df[:len(subset_metadata_df)],
                                          empty_df,
                                          analysis_metadata_df[len(subset_metadata_df):]], axis=0)
        analysis_metadata_df = analysis_metadata_df[analysis_df.columns]

        analysis_metadata_df.to_csv(pcl_out,
                                    index=False,
                                    header=header,
                                    sep='\t',
                                    na_rep=na_rep)

    output_folder = os.path.dirname(analysis_files[0])
    pcl_files = bb_utils.name_files(analysis_files,
                                    output_folder,
                                    extension="pcl.tsv")

    for (a_file, p_file) in zip(analysis_files, pcl_files):
        _add_metadata_to_tsv(a_file, p_file)

    return pcl_files


def main(args):
    config = parse_cfg_file(args.config_file, section=args.data_type)
    analysis_col_patterns = config.get('analysis_col_patterns')
    target_metadata_cols = funcy.flatten(config.get('target_metadata_cols'))
    data_type_label = config.get('data_type_mapping').get(args.data_type)

    output_pcl = add_metadata_to_tsv(args.input_file,
                                     args.metadata_file,
                                     data_type_label,
                                     args.id_col,
                                     analysis_col_patterns,
                                     args.drop_missing_cols,
                                     target_metadata_cols,
                                     supplement=args.supplement)


if __name__ == "__main__":
    main(parse_cli_arguments())
