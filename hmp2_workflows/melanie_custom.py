from glob import glob
from all_pipeline import (
    knead,
    metaphlan2,
    humann2,
    h2ntdb,
    h2pdb,
    kneaddb_deps,
    kneaddb,
    conf,
    bin
    )
from anadama import RunContext
from anadama.backends import LevelDBBackend
from anadama.reporters import LoggerReporter
from anadama.sge import SGEPowerup
from anadama.deps import HugeFileDependency as HFD
from anadama.deps import GlobDependency
from anadama.util import fname

be = LevelDBBackend("melanie_db", False)

dogdir = "/broad/geverslab_microbiome/melanie/HMP2_March2016/data/TechRep_DogStool/bam"
inbams = glob(dogdir+"/*.bam")

ctx = RunContext(storage_backend=be,
                 grid_powerup=SGEPowerup(queue="long", tmpdir="/seq/ibdmdb/tmp/anasge"))
reporter = LoggerReporter(ctx, loglevel_str='debug', 
                          logfile="/seq/ibdmdb/melanie_run.log")
outdir = "/seq/ibdmdb/melanie"
dogdb = "/seq/ibdmdb/centos6/var/lib/CanFam3.1"
dogdb_dep = GlobDependency(dogdb+"/*")
ctx.already_exists(dogdb_dep)
ctx.already_exists(*bin.values())
ctx.already_exists(h2ntdb, h2pdb)
ctx.already_exists(*kneaddb_deps.values())
ctx.already_exists(*map(HFD, conf.values()))

for bam in inbams:
    bamdep = HFD(bam)
    ctx.already_exists(bamdep)
    cleanfq = fname.mangle(bam, dir=outdir, ext="fastq")
    knead(ctx, bamdep, cleanfq,
          fname.mangle(bam, dir=outdir, tag="clean", ext="log"),
          dbs=[dogdb, kneaddb.six], db_deps=[dogdb_dep, kneaddb_deps.six])
    wgs_taxtxt = fname.mangle(bam, dir=outdir, tag="tax", ext="txt")
    metaphlan2(ctx, cleanfq, 
               fname.mangle(bam, dir=outdir, ext="biom"),
               wgs_taxtxt, 
               fname.mangle(bam, dir=outdir, tag="strain", ext="bam"),
               mem=3500, time=120)
    humann2(ctx, cleanfq, wgs_taxtxt, 
            fname.mangle(bam, dir=outdir, tag="humann2", ext="tar.bz2"))

ctx.go(reporter=reporter, n_grid_parallel=35, n_parallel=2)
