import sys
import csv
import json

from toolz import pluck

metadata_json = ("/seq/ibdmdb/data_deposition/HMP2"
                 "/Metadata/json/joined_subcoll.json")

metadata_names = ("/seq/ibdmdb/data_deposition/HMP2"
                  "/Metadata/json/hospital_mapping.csv")

def load_csv(fname):
    with open(fname) as f:
        reader = csv.reader(f)
        lines = list(reader)
    return lines


hosp_map = dict(pluck([1,2], load_csv(metadata_names)))

with open(metadata_json) as f:
    data = json.load(f)

# headers = [ hosp_map.get(h, h) for h in data['headers'] ]
headers = data['headers']

csv.writer(sys.stdout).writerows([headers]+data['data'])

