import csv
from collections import namedtuple
from itertools import count

import cutlass
from toolz import groupby, first

from .. import matcher


def load(session, sample_ids):
    if not hasattr(sample_ids, "__iter__"): # exclude strings
        sample_ids = (sample_ids,)
    sample_ids = set(list(sample_ids))
    for i in count(1):
        res = session.get_osdf().oql_query(
            "ihmp", '"subject"[node_type]', page=i
        )['results']

        for data in res:
            if not sample_ids:
                break
            if data['meta']['rand_subject_id'] in sample_ids:
                yield cutlass.Subject.load_subject(data)
                sample_ids.remove(data['meta']['rand_subject_id'])
        if not sample_ids:
            break
            

def fields(fname):
    with open(fname, 'rb') as f:
        reader = csv.reader(f)
        Subject = namedtuple("Subject", next(reader), rename=True)
        for row in reader:
            yield Subject._make(row)

rename_cache = dict()
def _rename(s, valid_choices):
    global rename_cache
    s = s.lower()
    if s in rename_cache:
        return rename_cache[s]
    else:
        results = matcher.closest(s, valid_choices)
        if len(results) > 1:
            results = [r for r in results if r[1] == s]
            if not results:
                raise ValueError("Not sure how to rename `%s'. Choices: %s"%(
                    s, [r[1] for r in results]))
        distance, ret = results[0]
        if distance / (len(s)-1.) > 0.75:
            raise ValueError("Variable rename probably failed: %s -> %s" %(
                s, ret))
        rename_cache[s] = ret
        return ret
        
def merge(new, into):
    into.gender = _rename(new.gender, cutlass.Subject.valid_genders)
    into.race = _rename(new.race, cutlass.Subject.valid_races)
    return into

_racemap = {
    'Black or African American': "african_american",
    'American Indian or Alaska Native': "american_indian_or_alaska_native",
    'Other': "ethnic_other",
    'More than one race': "ethnic_other",
    'White': "caucasian",
    'Asian': "asian",
    '': "ethnic_other"
}

def parse_race(record):
    if not record.race:
        return "ethnic_other"
    if bool(record.hispa.strip()) and "not" not in record.hispa.lower():
        return "hispanic_or_latino"
    k = _rename(record.race, _racemap.keys())
    return _racemap.get(k)

def parse_record(record, study, subject_cache):
    if record[0] in subject_cache:
        s = subject_cache[record[0]] 
    else:
        s = cutlass.Subject()
    try:
        s.gender = _rename(record.sex, cutlass.Subject.valid_genders)
    except ValueError:
        s.gender = "unknown"
    s.race = parse_race(record)
    s.rand_subject_id = record[0]
    s.links['participates_in'] = [study.id]

    if record[8].strip().lower() == "yes":
        s.tags.append("subject_withdrew")
    if record[9].strip().lower() == "yes":
        s.tags.append("terminated_by_investigator")
    return s
    
def from_file(fname, study):
    records = map(first, groupby(0, fields(fname)).itervalues())
    subject_cache = dict([ (s.rand_subject_id,s) for s in study.subjects() ])
    for record in records:
        yield parse_record(record, study, subject_cache)


    
def sync(study, input_fname, delete_missing=False):
    db_subjs = dict([ (s.rand_subject_id, s) for s in study.subjects() ])
    for in_subj in fields(input_fname):
        if in_subj[0] not in db_subjs:
            yield parse_record(in_subj)
        else:
            yield merge(in_subj, db_subjs.pop(in_subj[0]))
    if delete_missing:
        for s in db_subjs.itervalues():
            s.delete()

