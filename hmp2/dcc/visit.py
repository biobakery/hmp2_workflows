
import cutlass
import dateutil.parser
from itertools import chain
from toolz import groupby

from .subject import fields
from .project import default_mixs_dict

class settings:
    water = "ENVO:00005791" # sterile water
    ethanol = "ethanol"
    biome = "ENVO:00009003"
    feature = "ENVO:00002003" # Feces
    gen_loc_name = "USA"
    lat_lon = "+42.363664 -71.069230"
    body_product = "stool"
    env_package = "Human-associated"
    fma_body_site = "FMA:64183" # Feces
    body_site = "stool"


def can_parse_record(record, date_idx=12):
    return bool(record[date_idx].strip())


def parse_record(subject, record, i, visit_cache, date_idx=12,
                 etoh_samplid_idx=13, h2o_sampleid_idx=14):
    # add visit_cache to argument
    date = dateutil.parser.parse(
        record[date_idx]
    ).strftime(cutlass.Visit.date_format)
    if date in visit_cache:
        v = visit_cache[date]
    else:
        v = cutlass.Visit()
    v.visit_id = "{}_{}".format(subject.rand_subject_id, i)
    v.visit_number = i
    v.date = date
    v.interval = 0 if i == 1 else 1
    v.links['by'] = [subject.id]
    ret = [v]
    samples = groupby(lambda s: s.mixs['material'] , v.samples())
    if bool(record[etoh_samplid_idx].strip()):
        if settings.ethanol in samples:
            s = samples[settings.ethanol][0]
        else:
            s = cutlass.Sample()
        s.body_site = settings.body_site
        s.fma_body_site = settings.fma_body_site
        s.tags.append(record[etoh_samplid_idx].strip())
        d = default_mixs_dict()
        d['biome'] = settings.biome
        d['collection_date'] = date
        d['material'] = settings.ethanol
        d['env_package'] = settings.env_package
        d['feature'] = settings.feature
        d['geo_loc_name'] = settings.gen_loc_name
        d['lat_lon'] = settings.lat_lon
        d['body_product'] = settings.body_product
        s.mixs = d
        ret.append(s)
    if bool(record[h2o_sampleid_idx].strip()):
        if settings.water in samples:
            s = samples[settings.water][0]
        else:
            s = cutlass.Sample()
        s.body_site = settings.body_site
        s.fma_body_site = settings.fma_body_site
        s.tags.append(record[h2o_sampleid_idx].strip())
        d = default_mixs_dict()
        d['biome'] = settings.biome
        d['collection_date'] = date
        d['material'] = settings.water
        d['env_package'] = settings.env_package
        d['feature'] = settings.feature
        d['geo_loc_name'] = settings.gen_loc_name
        d['lat_lon'] = settings.lat_lon
        d['body_product'] = settings.body_product
        s.mixs = d
        ret.append(s)
        
    return ret


def from_file(fname, subjects, nest=True):
    records = groupby(0, fields(fname))
    for subject in subjects:
        visit_records = records.get(subject.rand_subject_id, None)
        if not visit_records:
            continue
        visit_cache = dict([ (v.date, v) for v in subject.visits() ])
        groups = iter( parse_record(subject, rec, i, visit_cache)
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

    
    


