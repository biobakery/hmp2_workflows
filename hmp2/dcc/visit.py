
import cutlass
import dateutil.parser
from itertools import chain
from toolz import groupby

from .subject import fields
from .project import default_mixs_dict

def can_parse_record(record, date_idx=12):
    return bool(record[date_idx].strip())

def parse_record(subject, record, i, date_idx=12, etoh_samplid_idx=13,
                 h2o_sampleid_idx=14):
    # add visit_cache to argument
    v = cutlass.Visit()
    v.visit_id = "{}_{}".format(subject.rand_subject_id, i)
    v.visit_number = i
    date = dateutil.parser.parse(
        record[date_idx]
    ).strftime(cutlass.Visit.date_format)
    v.date = date
    v.interval = 0 if i == 1 else 1
    v.links['by'] = [subject.id]
    ret = [v]
    if bool(record[etoh_samplid_idx].strip()):
        s = cutlass.Sample()
        s.body_site = "stool"
        s.fma_body_site = "FMA:64183" # Feces
        s.tags.append(record[etoh_samplid_idx].strip())
        d = default_mixs_dict()
        d['biome'] = "ENVO:00009003"
        d['collection_date'] = date
        d['material'] = "ethanol"
        d['env_package'] = "Human-associated"
        d['feature'] = "ENVO:00002003" # Feces
        d['geo_loc_name'] = "USA"
        d['lat_lon'] = "+42.363664 -71.069230"
        d['body_product'] = "stool"
        s.mixs = d
        ret.append(s)
    if bool(record[h2o_sampleid_idx].strip()):
        s = cutlass.Sample()
        s.body_site = "stool"
        s.fma_body_site = "FMA:64183" # Feces
        s.tags.append(record[h2o_sampleid_idx].strip())
        d = default_mixs_dict()
        d['biome'] = "ENVO:00009003"
        d['collection_date'] = date
        d['material'] = "ENVO:00005791" # sterile water
        d['env_package'] = "Human-associated"
        d['feature'] = "ENVO:00002003" # Feces
        d['geo_loc_name'] = "USA"
        d['lat_lon'] = "+42.363664 -71.069230"
        d['body_product'] = "stool"
        s.mixs = d
        ret.append(s)
        
    return ret


def from_file(fname, subjects, nest=True):
    records = groupby(0, fields(fname))
    # visit_cache = dict([(s.rand_subject_id,s) for s in subjects() ])
    for subject in subjects:
        visit_records = records.get(subject.rand_subject_id, None)
        if not visit_records:
            pass
        groups = iter( parse_record(subject, rec, i)
                       for i, rec in enumerate(visit_records, 1)
                       if can_parse_record(rec) )
        if nest:
            yield groups
        else:
            for group in groups:
                for obj in group:
                    yield obj

class SaveError(ValueError):
    def __init__(self, validation_errors, *args, **kwargs):
        self.validation_errors = validation_errors
        return super(SaveError, self).__init__(*args, **kwargs)

def save_all(visit_groups):
    visit_groups = map(list, visit_groups)
    visit_groups = filter(None, visit_groups) # remove empty lists
    visit_groups = list(chain.from_iterable(visit_groups)) # one visit per item
    valids = list()
    valid_errors = list()
    sample_groups = list()
    visits = list()
    for g in visit_groups:
        v, samples = g[0], g[1:]
        valid = v.is_valid()
        valids.append(valid)
        if not valid:
            valid_errors.append(v.validate())
        else:
            valid_errors.append(None)
        visits.append(v)
        sample_groups.append(samples)
    if not all(valids):
        raise SaveError(valid_errors)
    for v in visits:
        v.save()
    valids = list()
    valid_errors = list()
    for v, g in zip(visits, sample_groups):
        group_valids = list()
        group_errors = list()
        for sample in g:
            sample.links['collected_during'] = [v.id]
            valid = sample.is_valid()
            group_valids.append(valid)
            if valid:
                group_errors.append(None)
            else:
                group_errors.append(sample.validate())
        valids.append(group_valids)
        valid_errors.append(group_errors)
    if not all(map(all, valids)):
        raise SaveError(valid_errors)
    for g in sample_groups:
        for sample in g:
            sample.save()

    
    


