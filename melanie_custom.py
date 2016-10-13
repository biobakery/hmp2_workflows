from glob import glob
from all_pipeline import (
    knead,
    metaphlan2,
    humann2,
    h2ntdb,
    h2pdb,
    kneaddb_deps,
    conf,
    bin
    )
from anadama import RunContext
from anadama.reporters import LoggerReporter
from anadama.sge import SGEPowerup
from anadama.deps import HugeFileDependency as HFD
from anadama.util import fname

dogdir = "/broad/geverslab_microbiome/melanie/HMP2_March2016/data/TechRep_DogStool/bam"
inbams = glob(dogdir+"/*.bam")

ctx = RunContext(grid_powerup=SGEPowerup(queue="long", tmpdir="/seq/ibdmdb/tmp/anasge"))
reporter = LoggerReporter(ctx, loglevel_str='debug', 
                          logfile="/seq/ibdmdb/melanie_run.log")
outdir = "/seq/ibdmdb/melanie"
ctx.already_exists(*bin.values())
ctx.already_exists(h2ntdb, h2pdb)
ctx.already_exists(*kneaddb_deps.values())
ctx.already_exists(*map(HFD, conf.values()))

for bam in inbams:
    bamdep = HFD(bam)
    ctx.already_exists(bamdep)
    cleanfq = fname.mangle(bam, dir=outdir, ext="fastq")
    knead(ctx, bamdep, cleanfq,
          fname.mangle(bam, dir=outdir, tag="clean", ext="log"))
    wgs_taxtxt = fname.mangle(bam, dir=outdir, tag="tax", ext="txt")
    metaphlan2(ctx, cleanfq, 
               fname.mangle(bam, dir=outdir, ext="biom"),
               wgs_taxtxt, 
               fname.mangle(bam, dir=outdir, tag="strain", ext="bam"))
    humann2(ctx, cleanfq, wgs_taxtxt, 
            fname.mangle(bam, dir=outdir, tag="humann2", ext="tar.bz2"))

ctx.go(reporter=reporter, n_grid_parallel=35, n_parallel=2)
