import csv
from collections import OrderedDict
from os.path import basename

from . import fname_search

default_kvs = lambda : OrderedDict([
    ("sample_name",""),
    ("collection_timestamp",""),
    ("physical_location","Broad Insitute 415 Main Street Cambridge, MA 02142"),
    ("taxon_id","9606"),
    ("description",""),
    ("scientific_name","Homo sapiens"),
    ("sample_type","Stool"),
    ("physical_specimen_remaining","TRUE"),
    ("dna_extracted","TRUE"),
    ("latitude",""),
    ("longitude",""),
    ("host_subject_id",""),
    ("preprocessed_name", ""),
    ("has_extracted_data", "TRUE"),
    ("required_sample_info_status", "completed"),
    ("has_physical_specimen", "TRUE"),
])

latlon = {"C": (34.07533, -118.39253), #c.sinai
          "H": (39.13787, -84.50358),  #cin. children's
          "M": (42.36279, -71.06807),  #mass. general
          "E": (33.79322, -84.32333)}  #emory


def output(out_dicts, outtsv_fname):
    with open(outtsv_fname, 'w') as f:
        print >> f, "\t".join(default_kvs().iterkeys())
        for d in out_dicts:
            print >> f, "\t".join(d.itervalues())


def create(fna_fnames, incsv_fname, outtsv_fname, skip_header=True,
           sample_name_idx=18, subject_idx=17, fname_prefix_idx=16,
           collection_time_idx=13, description_idx=11):
    def _out_dicts():
        with open(incsv_fname, 'r') as in_f:
            reader = csv.reader(in_f)
            if skip_header:
                next(reader)
            for line in reader:
                d = default_kvs()
                d['sample_name'] = line[sample_name_idx]
                d['latitude'], d['longitude'] = map(
                    str, latlon.get(line[subject_idx][0], ("", "")) )
                d['preprocessed_name'] = basename(
                    fname_search(fna_fnames, line[fname_prefix_idx]))
                d['collection_timestamp'] = line[collection_time_idx]+" 00:00"
                d['description'] = line[description_idx]
                d['host_subject_id'] = line[subject_idx]
                yield d

    return output(_out_dicts(), outtsv_fname)
            

