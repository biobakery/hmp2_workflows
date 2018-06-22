#!/usr/bin/env python

"""
mgx_assembly_gene_calling.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Runs assembly utilizing MEGAHIT and gene calling using Prodigal on the provided
metagenomic dataset.

Copyright (c) 2018 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
"""


import itertools
import os
import shutil
import tempfile

from itertools import chain

from anadama2 import Workflow

from glob2 import glob



def parse_cli_arguments():
    """Parses any command-line arguments passed into the workflow.

    Args:
        None
    Requires:
        None
    Returns:
        anadama2.Workflow: The workflow object for this pipeline
        anadama2.cli.Configuration: Arguments passed into this workflow.
    """
    workflow = Workflow(version='0.1', description='A workflow to assemble '
                        'metagenomic data and run a gene caller on the '
                        'resulting contigs.')
    workflow.add_argument('contaminant-db', desc='KneadData DNA contaminants database.')
    workflow.add_argument('file-extension', desc='Extension of input files to '
                          'assemble and gene call on.', default='.fastq.gz')
    workflow.add_argument('threads', desc='number of threads/cores for each '
                          'task to use', default=1)
    workflow.add_argument('memory', desc='The amount of memory to use for each '
                          'assembly job. Provided in GB', default='10240')


    return workflow


def main(workflow):
    args = workflow.parse_args()

    sequence_files = glob(os.path.join(args.input, "*%s" % args.file_extension))
    samples = [os.path.basename(s).split(os.extsep)[0] for s in sequence_files]
    ## MEGAHIT will work a bit better if we are working with paired-end
    ## data here so let's split things apart
    split_dir = os.path.join(args.output, 'deinterleave')
    workflow.add_task('mkdir -p [targets[0]]',
                      targets=split_dir)


    sorted_dir = os.path.join(args.output, 'sorted')
    workflow.add_task('mkdir -p [targets[0]]',
                      depends=split_dir,
                      targets=sorted_dir)

    ## It seems like most of the data we are dealing with are not actual
    ## proper deinterleaved sequences so we'll need to sort them beforehand
    sorted_seqs = [os.path.join(sorted_dir, "%s.sorted.fastq.gz" % s) for s in samples]

    for (in_seq, out_seq) in zip(sequence_files, sorted_seqs):
        sample_name = os.path.basename(in_seq).split(os.extsep)[0]
        temp_dir = os.path.join(sorted_dir, "%s.tmp" % sample_name)
    
        workflow.add_task('mkdir -p [targets[0]]',
                          depends=[sorted_dir],
                          targets=[temp_dir])

        workflow.add_task_gridable('zcat [depends[0]] | paste - - - - | '
                                   'sort -T [args[1]] -k1,1 -S [args[0]] | tr \'\t\' \'\n\' | '
                                   'pigz --best -p 4 > [targets[0]]',
                                   depends=[in_seq, sorted_dir, temp_dir],
                                   targets=out_seq,
                                   args=['10G', temp_dir],
                                   cores='4',
                                   mem=args.memory,
                                   time='120')

        # Delete our temp directories from the sort step
        workflow.add_task('rm -rf [depends[0]]',
                          depends=[temp_dir, out_seq])


    ## Drop out any reads in our interleaved file that do not have a matching paired
    ## sequence
    matched_seqs = [s.replace('.sorted.fastq', '.paired.fastq') for s in sorted_seqs]
    workflow.add_task_group_gridable('seqtk dropse [depends[0]] | pigz --best -p 4 > [targets[0]]',
                                     depends=sorted_seqs,
                                     targets=matched_seqs,
                                     cores='4',
                                     mem='4096',
                                     time='60')

    for (sorted_seq, matched_seq) in zip(sorted_seqs, matched_seqs):
        workflow.add_task('rm [depends[0]]',
                          depends=[sorted_seq, matched_seq])
    
    # And finally peel off the orphans into a separate file
    for (raw_seq, matched_seq) in zip(sorted_seqs, matched_seqs):
       orphan_seq = raw_seq.replace('.paired.fastq', '.orphans.fastq')
       workflow.add_task_gridable('extract_orphans.sh [depends[0]] [depends[1]] [args[0]]',
                                  depends=[raw_seq, matched_seq],
                                  targets=[orphan_seq],
                                  args=[split_dir])

    split_files = []
    for in_seq in matched_seqs:
        seq_base = os.path.basename(in_seq).split(os.extsep)[0]
        f_seq = os.path.join(split_dir, "%s_R1.fastq.gz" % seq_base)
        r_seq = os.path.join(split_dir, "%s_R2.fastq.gz" % seq_base)
        split_files.append((f_seq, r_seq))

        workflow.add_task_gridable('seqtk seq -1 [depends[0]] | pigz --best -p 4 > [targets[0]] && '
                                   'seqtk seq -2 [depends[0]] | pigz --best -p 4 > [targets[1]]',
                                   depends=[in_seq, split_dir],
                                   targets=[f_seq, r_seq],
                                   cores='4',
                                   mem='4096',
                                   time='120')



    ## We need to run KneadData on our sequences first.
    qc_out_dir = os.path.join(args.output, 'qc')
    workflow.add_task('mkdir -p [targets[0]]',
                      depends=split_dir,
                      targets=qc_out_dir,
                      name="mkdir_qc")

    cleaned_seqs = []
    unmatched_seqs = []
    for (f_seq, r_seq) in split_files:
        seq_base = os.path.basename(f_seq).split(os.extsep)[0]
        f_seq_cleaned = os.path.join(qc_out_dir, '%s_kneaddata_paired_1.fastq' % seq_base)
        r_seq_cleaned = os.path.join(qc_out_dir, '%s_kneaddata_paired_2.fastq' % seq_base)
        f_seq_unmatched = os.path.join(qc_out_dir, '%s_kneaddata_unmatched_1.fastq' % seq_base)
        r_seq_unmatched = os.path.join(qc_out_dir, '%s_kneaddata_unmatched_2.fastq' % seq_base)

        workflow.add_task_gridable('kneaddata --input [depends[0]] --input [depends[1]] '
                                   '--reference-db [args[0]] --output [args[1]] --serial',
                                   depends=[f_seq, r_seq],
                                   targets=[f_seq_cleaned, r_seq_cleaned, f_seq_unmatched, r_seq_unmatched],
                                   args=[args.contaminant_db, qc_out_dir],
                                   cores=args.threads,
                                   mem='8192',
                                   time='180')

        cleaned_seqs.append((f_seq_cleaned, r_seq_cleaned))
        unmatched_seqs.append((f_seq_unmatched, r_seq_unmatched))

    ## Now run MEGAHIT
    assembly_dir = os.path.join(args.output, 'assembly')
    workflow.add_task('mkdir -p [targets[0]]',
                      depends=qc_out_dir,
                      targets=assembly_dir)

    megahit_contigs = []
    for (cleaned_seqs, unmatched_seqs) in zip(cleaned_seqs, unmatched_seqs):
        seq_base = os.path.basename(cleaned_seqs[0]).split(os.extsep)[0].replace('_R1_kneaddata_paired_1', '')
        megahit_contig_dir = os.path.join(assembly_dir, seq_base)
        megahit_contig = os.path.join(megahit_contig_dir, '%s.contigs.fa' % seq_base)

        ## MEGAHIT needs memory in a byte format so let's take care of thata
        float_mem = float(args.memory) * 1000000

        workflow.add_task('mkdir -p [targets[0]]',
                          depends=assembly_dir,
                          targets=megahit_contig_dir)

        workflow.add_task_gridable('megahit -1 [depends[0]] -2 [depends[1]] '
                                   '-r [depends[2]],[depends[3]] -t [args[0]] '
                                   '-m [args[1]] -o [targets[0]] --out-prefix [args[2]]',
                                   depends=cleaned_seqs + unmatched_seqs,
                                   targets=[megahit_contig_dir, megahit_contig],
                                   args=[args.threads, float_mem, seq_base],
                                   cores=args.threads,
                                   mem=args.memory,
                                   time='240')

        megahit_contigs.append(megahit_contig)

    ## As per Damians workflow let's filter out some of the smaller contigs
    filtered_contigs = []
    for contig in megahit_contigs:
        contig_base = os.path.basename(contig).split(os.extsep)[0]
        contig_dir = os.path.dirname(contig)
        filtered_contig_file = os.path.join(contig_dir, '%s.min500.contigs.fa' % contig_base)

        workflow.add_task_gridable('cat [depends[0]] | awk -v var="[args[0]]" '
                                   '\'{ if(substr($0,0,1) == ">") {header=substr($0, 2,length($0))} '
                                   'else {seq=$0; if(length($0) >= 500) {print ">"var"_"header"\\n"seq}} }\''
                                   ' > [targets[0]]',
                                   depends=[contig],
                                   targets=[filtered_contig_file],
                                   time='60',
                                   mem='4096',
                                   cores=args.threads,
                                   args=[contig_base])

        filtered_contigs.append(filtered_contig_file)


    annotation_dir = os.path.join(args.output, 'annotation')
    workflow.add_task('mkdir -p [targets[0]]',
                      depends=assembly_dir,
                      targets=annotation_dir)

    ## And finally Prodigal
    for contig in filtered_contigs:
        contig_base = os.path.basename(contig).split(os.extsep)[0]
        gff_file = os.path.join(annotation_dir, '%s.gff' % contig_base)
        cds_file = os.path.join(annotation_dir, '%s.fna' % contig_base)
        cds_aa = os.path.join(annotation_dir, '%s.faa' % contig_base)
        stderr_log = os.path.join(annotation_dir, '%s.stderr.log' % contig_base)
        stdout_log = os.path.join(annotation_dir, '%s.stdout.log' % contig_base)

        workflow.add_task_gridable('prodigal -p meta -i [depends[0]] '
                                   '-f gff -o [targets[0]] -d [targets[1]] '
                                   '-a [targets[2]] '
                                   '2> [args[0]] > [args[1]]',
                                   depends=[contig],
                                   targets=[gff_file, cds_file, cds_aa],
                                   args=[stderr_log, stdout_log],
                                   cores=args.threads,
                                   mem='8192',
                                   time='120')

    workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
