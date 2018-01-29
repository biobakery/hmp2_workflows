# -*- coding: utf-8 -*-

"""
generate_biopsy_sample_counts.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parses an HMP2 metadata file and a StudyTrax clinical file to retrieve the 
unique biopsy sample count. 

Biopsy samples are unique in the sense that they don't have a traditional 
collection number like stool samples do but we can take advantage of the 
IntervalSequence identifier in the clinical metadata to collapse each 
data product down to a specific patient + visit.
"""

import argparse

import pandas as pd


def _generate_biopsy_id(row):
    subject_id = row.get('ProjectSpecificID')
    interval_seq = row.get('IntervalSequence')
    
    return "%s_%s" % (subject_id, interval_seq)


parser = argparse.ArgumentParser()

parser.add_argument('-m', '--hmp2-metadata-file')
parser.add_argument('-s', '--studytrax-metadata-file')

args = parser.parse_args()

metadata_df = pd.read_csv(args.hmp2_metadata_file)
studytrax_df = pd.read_csv(args.studytrax_metadata_file)

all_biopsy_ids = []

for dtype in ['host_transcriptomics', 'biopsy_16S', 'methylome']:
    filtered_df = metadata_df[(metadata_df['data_type'] == dtype) & 
                               ~(metadata_df['IntervalName'].str.startswith('Stool')) & 
                               ~(metadata_df['IntervalName'].str.startswith('Follow-up')) &
                               ~(metadata_df['IntervalName'].str.startswith('Follow-Up')) &
                               ~(metadata_df['IntervalName'].str.startswith('Baseline'))]
    biopsy_ids = list(filtered_df.apply(_generate_biopsy_id, axis=1))
    print "Data type: %s, Number of samples: %s" % (dtype, len(biopsy_ids))

    all_biopsy_ids.extend(biopsy_ids)

print
print "Total biopsy samples: %s" % len(all_biopsy_ids)


#unique_biopsy_ids = set(all_biopsy_ids)
#print
#print "Unique biopsy samples %s" % len(unique_biopsy_ids)