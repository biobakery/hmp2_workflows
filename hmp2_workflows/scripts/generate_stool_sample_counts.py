# -*- coding: utf-8 -*-

"""
generate_stool_sample_counts.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parses an HMP2 metadata file and a StudyTrax clinical file to retrieve the 
unique sample count. 

This sample count is all unique site/sub/coll ID's for all stool sample based
data products plus the number of samples that have a fecalcal measurement that 
are not in the set of stool based data products.
"""

import argparse

import pandas as pd

SITE_MAPPING = {'Cincinnati': 'H',
                'Massachusetts General Hospital': 'M',
                'Emory': 'E',
                'MGH Pediatrics': 'P',
                'Cedars-Sinai': 'C'}


def _generate_site_sub_coll_id(row):
    patient_id = row['ProjectSpecificID']
    site_abbrev = SITE_MAPPING[row['SiteName']]
    stool_coll_num = row['IntervalName'].replace('Stool Collection #', '')

    return site_abbrev + str(patient_id) + "C" + stool_coll_num
                                

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--hmp2-metadata-file')
parser.add_argument('-s', '--studytrax-metadata-file')

args = parser.parse_args()

metadata_df = pd.read_csv(args.hmp2_metadata_file)
studytrax_df = pd.read_csv(args.studytrax_metadata_file)

site_sub_coll_ids = []

for stool_dtype in ['metagenomics', 'metabolomics', 'stool_16S', 'viromics', 'proteomics', 'metatranscriptomics']:
    samples = metadata_df[(metadata_df['data_type'] == stool_dtype) & ~(metadata_df['Project'].str.endswith('BP'))]['site_sub_coll'].tolist()
    print "Data type: %s, Number of samples: %s" % (stool_dtype, len(samples))

    site_sub_coll_ids.extend(samples)

print
print "Total non-unique stool samples: %s" % len(site_sub_coll_ids)

site_sub_coll_ids = set(site_sub_coll_ids)
print "Unique stool samples (not including fecalcal): %s" % len(site_sub_coll_ids)

studytrax_stool_df = (studytrax_df[(studytrax_df['IntervalName'].str.startswith('Stool')) &
                                   (studytrax_df['st_q16'].notnull())])

print
print "Total number of fecalcal samples: %s" % len(studytrax_stool_df)
studytrax_site_sub_coll = set(studytrax_stool_df.apply(_generate_site_sub_coll_id, axis=1))

print "Unique fecalcal samples: %s" % len(studytrax_site_sub_coll)

print "Number fecalcal samples with no stool equivalent: %s" % (len(studytrax_site_sub_coll - site_sub_coll_ids))

site_sub_coll_ids = site_sub_coll_ids.union(studytrax_site_sub_coll)

print
print "Number of unique stool samples: %s" % len(site_sub_coll_ids)