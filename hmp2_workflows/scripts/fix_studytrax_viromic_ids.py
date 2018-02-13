#!/usr/bin/env python

import argparse

import pandas as pd


parser = argparse.ArgumentParser('Fixes mis-mapped viromics IDs in the '
                                 'Studytrax clinical data')
parser.add_argument('-i', '--input-samples')
parser.add_argument('-s', '--studytrax-metadata')
parser.add_argument('-b', '--broad-tracking-sheet')
parser.add_argument('-o', '--output-file')

args = parser.parse_args()

studytrax_df = pd.read_csv(args.studytrax_metadata)
broad_df = pd.read_csv(args.broad_tracking_sheet)

with open(args.input_samples) as samples_fh:
    for line in samples_fh:
        sample_id = line.strip()

        # First check if we have a straight matching into our viromics columns
        if any(studytrax_df['st_q12'].isin([sample_id])):
            continue
        # If that fails check to make sure that the viromics ID isn't in our sample storage tube
        elif any(studytrax_df['st_q13'].isin([sample_id])):
            studytrax_row = studytrax_df[studytrax_df['st_q13'] == sample_id]
            storage_id = studytrax_row['st_q12'].values[0]
            parent_sample_id = studytrax_row['st_q4'].values[0]

            broad_row = broad_df[broad_df['Parent Sample A'] == parent_sample_id]

            if broad_row['Storage'].values[0] == storage_id:
                # Swap em!
                studytrax_row['st_q13'] = storage_id
                studytrax_row['st_q12'] = sample_id
            elif broad_row['DNA/RNA'].values[0] == storage_id:
                # This is now more complicated... 
                #       
                # The storage ID contains the DNA/RNA ID which should be in the st_q10 col
                # and likely the ID in the st_q10 col is the storage ID so we need to do 
                # a bunch of swaps here.
                dna_rna_id = storage_id
                storage_id = studytrax_row['st_q10'].values[0]

                if (broad_row['Storage'].values[0] == storage_id):
                    studytrax_row['st_q10'] = dna_rna_id
                    studytrax_row['st_q12'] = sample_id
                    studytrax_row['st_q13'] = storage_id

            else:
                raise ValueError("Could not identify sample ID in st_q17: %s" % sample_id)

            studytrax_df[studytrax_df['st_q13'] == sample_id] = studytrax_row
        else:
            # Now we are assuming that the st_q12 column is empty for this row which means we need 
            # to populate it from the Broad tracking spreadsheet.
            broad_df_row = broad_df[broad_df['Viromics'] == sample_id]
            parent_sample_id = broad_df_row['Parent Sample A'].values[0]
            subject_id = broad_df_row['Subject'].values[0]

            studytrax_row = studytrax_df[studytrax_df['st_q4'] == parent_sample_id]
            if pd.isnull(studytrax_row['st_q12']).values[0]:
                studytrax_patient = studytrax_row['ProjectSpecificID'].values[0]

                if studytrax_patient == subject_id:
                    studytrax_row['st_q12'] = sample_id
            else:
                raise ValueError('st_q12 column is not empty for studytrax row: %s' % studytrax_row['st_q12'])

            studytrax_df[studytrax_df['st_q4'] == parent_sample_id] = studytrax_row

studytrax_df.to_csv(args.output_file, index=False)
