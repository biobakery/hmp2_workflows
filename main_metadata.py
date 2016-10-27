import os
import re
import sys
import csv
import json

from toolz import pluck, groupby
from dateutil.parser import parse as dateparse

#settings

valid_data_types = ("amplicon", "proteomics", 
                    "metagenomics", "metatranscriptomics")

CUSTOM_HEADERS = ['Project',                'External ID',      'Participant ID', 
                  'Site/Sub/Coll ID',       'data_type',        'week_num', 
                  'interval_days',          'visit_num',        'Research Project',
                  'PDO Number',             'GSSR IDs',         'Product', 
                  'LCSET',                  'Aggregated Lanes', 'WR ID', 
                  '# Lanes in Aggregation', 'Total Reads']

studytrax_names = ("/seq/ibdmdb/data_deposition/HMP2"
                   "/Metadata/json/hospital_mapping.csv")

studytrax_rows = ("/seq/ibdmdb/data_deposition/HMP2"
                   "/Metadata/json/hospital.csv")

receipt_map_json = ("/seq/ibdmdb/data_deposition/HMP2"
                    "/Metadata/json/receipt_date_map.json")

fecalcal_csv = ("/seq/ibdmdb/data_deposition/HMP2"
                "/Metadata/json/fecalcal.csv")

age_map_csv = ("/seq/ibdmdb/data_deposition/HMP2"
               "/Metadata/json/age_map.csv")

melanie_ordering = ("/seq/ibdmdb/data_deposition/HMP2/Metadata/json/"
                    "melanie_order.csv")

prx_upload_dir = "/seq/ibdmdb/data_deposition/HMP2/Proteomics/1633/rschwager/"

mtx_pilot1 = ("/seq/ibdmdb/data_deposition/HMP2"
               "/Metadata/json/Aggregated Picard Report for IBDMDB TagSeq_wCollabSampleID2.csv")

wgs_pilot1 = ("/seq/ibdmdb/data_deposition/HMP2"
               "/Metadata/json/IBDMDB_PDO-5457_5458_WGS_SeqComplete_2.19.2015_wCollabSampleID2.csv")

sxs_pilot1 = ("/seq/ibdmdb/data_deposition/HMP2"
               "/Metadata/json/IBDMDB_16S_PDO-4852_SeqComplete.csv")

wgs_pilot2 = ("/seq/ibdmdb/data_deposition/HMP2"
               "/Metadata/json/IBDMDB_PilotII_WGS_SeqComplete_4.27.2016.csv")

prx_pilot1 = ("/seq/ibdmdb/data_deposition/HMP2/Metadata/json/"
              "proteomics.csv")

# end settings


def err(msg, *args, **kwargs):
    print >> sys.stderr, msg.format(*args, **kwargs)

    
def load_csv(fname):
    with open(fname) as f:
        reader = csv.reader(f)
        lines = list(reader)
    return lines


# setup
with open(receipt_map_json) as f:
    receipt_map = json.load(f)
receipt_map.pop("Site/Sub/Coll")
receipt_map = groupby(lambda k: re.search(r'^[A-Z](\d{4})C\d+', k[0]).group(1),
               receipt_map.items())
for sub in receipt_map:
    l = sorted(receipt_map[sub], key=lambda v: int(v[0].split("C")[1]))
    l = [ (a, dateparse(b)) for a, b in l
          if "/" in b]
    if not l:
        continue
    top = l[0][1]
    prev = top
    ret = {}
    for i, (subcoll, d) in enumerate(l):
        visit_no = subcoll.split("C")[1]
        ret[visit_no] = [(d-top).days/7,
                         (d-prev).days,
                         int(visit_no)]
        prev = d
    receipt_map[sub] = ret
    

fecalcal_map = {}
for k, v in pluck([-1, 5], load_csv(fecalcal_csv)[1:]):
    match = re.search(r'^[A-Z](\d+)C(\d+)', k)
    if not match:
        continue
    fecalcal_map[match.groups()] = v

age_map = dict( load_csv(age_map_csv) )

mel_ordering = load_csv(melanie_ordering)

st_rows = load_csv(studytrax_rows)
headers = st_rows[0]
studytrax_grp = groupby(headers.index("ProjectSpecificID"), st_rows[1:])
idx = headers.index("IntervalName")
for k in studytrax_grp:
    if not k:
        continue
    g = groupby(lambda i: re.sub(r'.*#(\d+).*', r'\1', i[idx]), 
                studytrax_grp[k])
    studytrax_grp[k] = g


prx_filemap = groupby(
    lambda f: re.sub(r'.*(.M.[A-Z0-9]{5}).*', r'\1', f),
    os.listdir(prx_upload_dir)
    )

# end setup

def lookup_date(sub, coll):
    if sub not in receipt_map:
        err("Unable to locate receipt date for subject {}.",
            sub)
        return None, None, None
    if coll not in receipt_map[sub]:
        err("No sample numbered `{}' yet received for subject {}.",
            coll, sub)
        return None, None, None
    return receipt_map[sub][coll]


def lookup_fecalcal(sub, coll):
    if (sub, coll) in fecalcal_map:
        return fecalcal_map[sub,coll]
    err("Unable to find fecalcal measurement for subject {} collection {}",
        sub, coll)


def lookup_age(sub):
    if sub in age_map:
        return age_map[sub]
    err("Unable to find age for subject {}.", sub)


def lookup_studytrax(sub, coll):
    if sub not in studytrax_grp:
        err("Unable to locate subject {} in studytrax dump.",
            sub)
        return None
    if coll not in studytrax_grp[sub]:
        err("Unable to locate collection `{}' for subject {}.",
            coll, sub)
        return None
    return studytrax_grp[sub][coll][0]
    

def lookup_prx(lims):
    if lims in prx_filemap:
        hits = prx_filemap[lims]
        if len(hits) > 1:
            msg = "Multiple proteomics files for lims label {}.".format(lims)
            err(msg)
        return hits
    err("Unable to locate proteomics file for lims label {}.",
        lims)

    
def create_row(gnum_or_path, lims_id, subject_id, collection_num, 
               data_type, studytrax_row,
               pdo_num=None, gssr_ids=None, product=None,
               lcset=None, aggregated_lanes=None, wr_id=None, 
               num_lanes_in_aggregation=None, total_reads=None):
    if data_type not in valid_data_types:
        raise ValueError("Data type "+data_type+" not recognized.")

    week_num, interval_days, visit_num = lookup_date(subject_id, collection_num)
    fecalcal = lookup_fecalcal(subject_id, collection_num)
    age_at_dx = lookup_age(subject_id)
    return ([ gnum_or_path, lims_id.replace("-", ""), "C"+subject_id, 
              subject_id+"C"+collection_num, data_type, week_num, interval_days, 
              visit_num, "ibdmdb", pdo_num, gssr_ids, product, lcset, 
              aggregated_lanes, wr_id, num_lanes_in_aggregation, total_reads ]
            +studytrax_row[:5]+[age_at_dx]+studytrax_row[5:10]+[fecalcal]
            +studytrax_row[10:] )


    
if __name__ == '__main__':
    rows = []
    for r in load_csv(mtx_pilot1)[1:]:
        subcoll = re.sub(r".*?(\d+C\d+)", r'\1', r[3])
        sub, coll = subcoll.split("C")
        studytrax_row = lookup_studytrax(sub, coll)
        if not studytrax_row:
            continue
        rows.append( create_row(r[0], r[1], sub, coll, "metatranscriptomics",
                                lookup_studytrax(sub, coll), pdo_num=r[6], gssr_ids=r[7],
                                lcset=r[10], product=r[9], aggregated_lanes=r[11], 
                                wr_id=r[12], num_lanes_in_aggregation=r[13], total_reads=r[14])
            )

    for i, r in enumerate(load_csv(wgs_pilot1)[1:]):
        if len(r) < 14:
            err("row {} too short for file: {}",
                i, wgs_pilot1)
            continue
        subcoll = re.sub(r".*?(\d+C\d+)", r'\1', r[3])
        sub, coll = subcoll.split("C")
        studytrax_row = lookup_studytrax(sub, coll)
        if not studytrax_row:
            continue
        rows.append( create_row(r[0], r[1], sub, coll, "metagenomics",
                                lookup_studytrax(sub, coll), pdo_num=r[5], gssr_ids=r[6],
                                product=r[7], lcset=r[8], aggregated_lanes=r[9], 
                                wr_id=r[10], num_lanes_in_aggregation=r[11], total_reads=r[12])
            )


    for r in load_csv(sxs_pilot1)[1:]:
        subcoll = re.sub(r".*?(\d+C\d+)", r'\1', r[3])
        sub, coll = subcoll.split("C")
        studytrax_row = lookup_studytrax(sub, coll)
        if not studytrax_row:
            continue
        rows.append( create_row(r[0], r[1], sub, coll, "amplicon",
                                lookup_studytrax(sub, coll), pdo_num="4852", gssr_ids=None,
                                product=None, lcset=None, aggregated_lanes=None, 
                                wr_id=None, num_lanes_in_aggregation=None, total_reads=r[8])
            )


    for r in load_csv(wgs_pilot2)[1:]:
        subcoll = re.sub(r".*?(\d+C\d+)", r'\1', r[3])
        sub, coll = subcoll.split("C")
        studytrax_row = lookup_studytrax(sub, coll)
        if not studytrax_row:
            continue
        rows.append( create_row(r[0], r[1], sub, coll, "metagenomics",
                                lookup_studytrax(sub, coll), pdo_num=r[7], gssr_ids=None,
                                product=r[8], lcset=r[9], aggregated_lanes=r[10], 
                                wr_id=None, num_lanes_in_aggregation=r[11], total_reads=r[12])
            )
        
    
    for r in load_csv(prx_pilot1)[1:]:
        subcoll = re.sub(r".*?(\d+C\d+)", r'\1', r[9])
        sub, coll = subcoll.rstrip("B").split("C")
        prx_fn = lookup_prx(r[0])
        if not prx_fn:
            continue
        for fn in prx_fn:
            prx_path = os.path.join(prx_upload_dir, fn)
            studytrax_row = lookup_studytrax(sub, coll)
            if not studytrax_row:
                continue
            rows.append( create_row(prx_path, r[1], sub, coll, "proteomics",
                                    lookup_studytrax(sub, coll)) )


    grp = groupby(0, rows)
    ordered = []
    for gid in mel_ordering:
        if gid[0] in grp:
            ordered.append(grp.pop(gid[0])[0])
    ordered.extend([r[0] for r in grp.values()])

    json.dump(
        {"headers": (CUSTOM_HEADERS+
                     st_rows[0][:5]+
                     ["age_at_consent"]+
                     st_rows[0][5:10]+
                     ["fecalcal_ng_ml"]+
                     st_rows[0][10:]),
         "data": ordered},
        sys.stdout
        )
