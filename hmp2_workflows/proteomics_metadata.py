import json
from toolz import groupby
import re
import csv

def load_csv(fname):
    with open(fname) as f:
        reader = csv.reader(f)
        lines = list(reader)
    return lines


def write_csv(fname, data):
    with open(fname, 'w') as out:
        csv.writer(out).writerows(data)


with open("joined_subcoll.json") as f:
    data = json.load(f)

headers = data['headers']
data = data['data']

hdata = load_csv("hospital.csv")
hheaders = hdata[0]
hdata = hdata[1:]

grp = groupby(hheaders.index("ProjectSpecificID"), hdata)
idx = hheaders.index("IntervalName")
for k in grp:
    g = groupby(lambda i: re.sub(r'.*#(\d+).*', r'\1', i[idx]), grp[k])
    grp[k] = g

pdata = load_csv("proteomics.csv")
pheaders = pdata[0]
pdata = pdata[1:]

pmap = dict([ (r[0], (r[1], re.sub(r'B$', '', r[-1]))) for r in pdata])

for k in pmap:
    l = list(pmap[k])
    l[1] = re.sub(r'^[A-Z]', '', l[1])
    pmap[k] = tuple(l)

results = {}
for k, fn in filemap.iteritems():
    sids = pmap[k]
    row = { "pid": k,
            "fname": fn,
            "sid": sids[0] }
    results[sids[1]] = row


with open("receipt_date_map.json") as f:
    receipt_map = json.load(f)

receipt_map.pop("Site/Sub/Coll")
from dateutil.parser import parse as dateparse
rgrp = groupby(lambda v: re.search(r'^[A-Z](\d{4})C\d+', v[0]).group(1),
               receipt_map.items())
for sub in rgrp:
    l = sorted(rgrp[sub], key=lambda v: int(v[0].split("C")[1]))
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
    rgrp[sub] = ret
        
age_map = dict( load_csv("../1623/eaandrews/Liz_Melanie_age_at_consent_6.6.2016_16.06.07_14.22.20.csv") )

fecalcal_map = {}
for k, v in pluck([-1, 5], load_csv("fecalcal.csv")[1:]):
    match = re.search(r'^[A-Z](\d+)C(\d+)', k)
    if not match:
        continue
    fecalcal_map[match.groups()] = v


path = "/seq/ibdmdb/data_deposition/HMP2/Proteomics/1633/rschwager/"
prs = []
for k, v in results.iteritems():
    subcoll = re.sub(r'^[A-Z]', '', k)
    sub, coll = subcoll.split("C")
    hrow = grp[sub].get(coll, [None])[0]
    if not hrow:
        print "can't find subcoll: ", subcoll
        continue
    packed = rgrp.get(sub, {}).get(coll, (None, None, None))
    week_no, int_days, visit_num = packed
    prs.append( [ path+v['fname'],
                  v['pid'].replace("-", ''),
                  "C"+sub,
                  subcoll,
                  "proteomics",
                  week_no, int_days, visit_num,
                  "ibdmdb",
                  "", None, None, None,
                  None, None, None, None]
                +hrow[:5]+[ age_map.get(sub, None) ]
                +hrow[5:10]+[ fecalcal_map.get((sub,coll), None) ]
                +hrow[10:]
    )

        
with open("joined_subcoll.json", 'w') as f:
    data.extend(prs)
    json.dump({"headers": headers, "data": data}, f)
