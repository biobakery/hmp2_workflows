import sys
import csv
import operator
import optparse
from functools import partial
from itertools import imap

from toolz import join, concat
from toolz import unique
from toolz.curried import get

def output(records, to_file=sys.stdout):
    close_end = False
    if type(to_file) is str:
        to_file = open(to_file, 'wb')
        close_end = True
    w = csv.writer(to_file, quoting=csv.QUOTE_MINIMAL)
    for r in records:
        w.writerow(list(r))
    if close_end:
        to_file.close()


def fields(fname):
    with open(fname, 'rb') as f:
        for row in csv.reader(f):
            yield row

def main(fname_a, fname_b, join_idx_a_func, join_idx_b_func,
         output_file, do_output=True):
    a_rows = fields(fname_a)
    b_rows = fields(fname_b)
    headers = concat([ next(a_rows), next(b_rows) ])
    joined = imap(concat, join(join_idx_a_func, a_rows,
                               join_idx_b_func, b_rows))
    if do_output:
        output(concat([[headers], joined]), output_file)
    else:
        return headers, joined
        

def join_16s(six_fname, other_fname, outfile, six_idx=0, other_idx=14):
    def _a(row):
        try:
            return row[six_idx][3:]
        except IndexError:
            pass
    def _b(row):
        try:
            return row[other_idx][3:]
        except IndexError:
            pass
    headers, joined = main(six_fname, other_fname, _a, _b, outfile,
                           do_output=False)
    joined = map(list, joined)
    output(concat([[headers], unique(joined, key=get(-2))]), outfile)
    

join_wgs = partial(join_16s, six_idx=1, other_idx=14)

    
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-1", "--join_idx_a", type=int, default=0)
    parser.add_option("-2", "--join_idx_b", type=int, default=1)
    opts, args = parser.parse_args()
    kwargs = dict( join_idx_a_func=operator.itemgetter(opts.join_idx_a),
                   join_idx_b_func=operator.itemgetter(opts.join_idx_b),
                   output_file=sys.stdout)
    sys.exit( main(*args, **kwargs) )
