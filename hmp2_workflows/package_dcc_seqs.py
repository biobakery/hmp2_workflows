import os
from itertools import izip
from gzip import GzipFile
from glob import glob

import anadama.reporters
from anadama import RunContext
from anadama.sge import SGEPowerup
from anadama.deps import HugeFileDependency as HFD
from anadama.util import fname
from anadama.helpers import system
from anadama.helpers import rm_r

from bunch import Bunch

from all_pipeline import output_dirs as prod_dirs


output_dirs = Bunch(
    metagenomics = "/seq/ibdmdb/tmp/package_dcc_seqs/wgs",
    metatranscriptomics = "/seq/ibdmdb/tmp/package_dcc_seqs/mtx",
    amplicon = "/seq/ibdmdb/tmp/package_dcc_seqs/16s",
)

def fastq_pair_split(fastq, outr1, outr2):
    def _actually_split(task):
        with open(fastq, 'r') as f, \
             GzipFile(outr1, 'w') as o1, \
             GzipFile(outr2, 'w') as o2:
            reads = izip(*[iter(f)]*4)
            for read in reads:
                if read[0].endswith("/1\n"):
                    o1.write("".join(read))
                elif read[0].endswith("/2\n"):
                    o2.write("".join(read))
    return _actually_split


def split_compact(fastq, gztar):
    tmpdir = fname.rmext(gztar, all=True)
    head, tail = os.path.split(tmpdir)
    outr1 = fname.mangle(gztar, dir=tmpdir, ext="fastq.gz", tag="1")
    outr2 = fname.mangle(gztar, dir=tmpdir, ext="fastq.gz", tag="2")
    return [ fastq_pair_split(fastq, outr1, outr2),
             system(["tar", "c", tail], stdout_clobber=gztar, working_dir=head),
             rm_r([tmpdir]) ]

def main():
    ctx = RunContext(grid_powerup = None)# SGEPowerup(queue="short", tmpdir="/seq/ibdmdb/tmp/anasge"))
    reporter = anadama.reporters.LoggerReporter(ctx, "DEBUG")
    for datatype in output_dirs:
        for f in glob(prod_dirs[datatype].cleanseq+"/*.fastq"):
            d = HFD(f)
            ctx.already_exists(d)
            gztar = fname.mangle(f, dir=output_dirs[datatype], ext="tar")
            ctx.grid_add_task(split_compact(f, gztar),
                              targets=HFD(gztar),
                              depends=d,
                              name="Split and compact "+f,
                              mem=100, cores=1, time=30)
            ctx.grid_add_task("md5sum {depends[0]} | awk '{{print $1;}}' > {targets[0]}",
                              depends=HFD(gztar),
                              targets=fname.mangle(gztar, ext="md5"),
                              mem=100, cores=1, time=30)
    ctx.go(reporter=reporter, n_parallel=2, n_grid_parallel=2)

if __name__ == '__main__':
    main()
