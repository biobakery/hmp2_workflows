import os
import re
import sys
import json
import operator
from hashlib import md5
from glob import glob

from toolz import groupby

from all_pipeline import output_dirs as dirs

dietmap = [
    ('"Tea or coffee no sugar and no sugar replacement"',                                  'dr_q9' ),
    ('"Soft drinks, tea or coffee with sugar (corn syrup, maple syrup, cane sugar, etc)"', 'dr_q10'),
    ('"Diet soft drinks, tea or coffee with sugar (Stevia, Equal, Splenda etc)"',          'dr_q11'),
    ('"Fruit juice (orange, apple, cranberry, prune etc.)"',                               'dr_q12'),
    ('"Water"',                                                                            'dr_q13'),
    ('"Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)"',               'dr_q14'),
    ('"Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)"',   'dr_q15'),
    ('"Dairy (milk, cream, ice cream, cheese, cream cheese)"',                             'dr_q16'), 
    ('"Probiotic"',                                                                        'dr_q17'),
    ('"Fruits (no juice) (Apples, raisins, bananas, oranges, strawberries, blueberries"',  'dr_q18'),
    ('"Vegetables (salad, tomatoes, onions, greens, carrots, peppers, green beans, etc)"', 'dr_q19'),
    ('"Beans (tofu, soy, soy burgers, lentils, Mexican beans, lima beans etc)"',           'dr_q20'),
    ('"Whole grains (wheat, oats, brown rice, rye, quinoa, wheat bread, wheat pasta)"',    'dr_q21'),
    ('"Starch (white rice, bread, pizza, potatoes, yams, cereals, pancakes, etc.)"',       'dr_q22'),
    ('"Eggs"',                                                                             'dr_q23'),
    ('"Processed meat (other red or white meat such as lunch meat, ham, salami, bologna"', 'dr_q24'),
    ('"Red meat (beef, hamburger, pork, lamb)"',                                           'dr_q25'),
    ('"White meat (chicken, turkey, etc.)"',                                               'dr_q26'),
    ('"Shellfish (shrimp, lobster, scallops, etc.)"',                                      'dr_q27'),
    ('"Fish (fish nuggets, breaded fish, fish cakes, salmon, tuna, etc.)"',                'dr_q28'),
    ('"Sweets (pies, jam, chocolate, cake, cookies, etc.)"',                               'dr_q29'),
]

first = operator.itemgetter(0)
last = operator.itemgetter(-1)

def load(fn):
    with open(fn) as f:
        a,b = [row.split('\t') for row in map(str.strip, f)]
    return dict(zip(a,b))

try:
    mapfile, dietfile = sys.argv[1:]
except Exception as e:
    print >> sys.stderr, str(e)
    print >> sys.stderr, "Usage: {} map.txt diet.txt".format(sys.argv[0])
    sys.exit(1)

files = []
for dir in (dirs.metagenomics.taxprof, dirs.amplicon.taxprof):
    files.extend(glob(os.path.join(dir, "*.tsv")))

rows = []
for f in files:
    biom = os.path.abspath(f).replace(".tsv", '.biom')
    if not os.path.exists(biom):
        print >> sys.stderr, "no biom for ", f
        continue
    try:
        with open(biom) as biomf:
            b = json.load(biomf)
        if -0.001 < sum(map(last, b['data'])) < 0.001:
            continue
    except:
        continue
    row = load(f)
    row['biom'] = biom
    rows.append(row)
    
grp = groupby(lambda o: re.sub(r'\D+', '', o['Participant ID']), rows)
with open(mapfile, 'w') as mf, open(dietfile, 'w') as df:
    print >> mf, "\t".join(("#SampleID", "DataFile", "Hash"))
    print >> df, "\t".join(["Project Specific Id"]+map(first, dietmap))
    for sid, samples in grp.iteritems():
        md5sum = md5(sid).hexdigest()[:8]
        for i, sample in enumerate(sorted(samples, key=lambda r: int(r['visit_num']))):
            m = [sid+"."+str(i), sample['biom'], md5sum+"."+str(i)]
            d = [sid]+[str(sample.get(k, '')) for  _, k in dietmap]
            print >> mf, "\t".join(m)
            print >> df, "\t".join(d)


