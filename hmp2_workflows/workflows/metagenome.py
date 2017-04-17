# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.metagenome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles HMP2 WGS metagenome data and metadata files.

Copyright (c) 2017 Harvard School of Public Health

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
import tempfile

from anadama2 import Workflow

from biobakery_workflows.tasks.shotgun import (quality_control, 
                                               taxonomic_profile,
                                               functional_profile)

from biobakery_workflows.utilities import find_files

from hmp2_workflows.tasks.common import (verify_files, stage_files,
                                         create_project_dirs,
                                         tar_files,
                                         make_files_web_visible)
from hmp2_workflows.tasks.file_conv import (bam_to_fastq,
                                            batch_convert_tsv_to_biom)
from hmp2_workflows.utils import (parse_cfg_file, 
                                  create_merged_md5sum_file)


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
    workflow = Workflow(version='0.1', description='A workflow to handle HMP2 '
                        'WGS data.', remove_options=['input', 'output'])
    workflow.add_argument('--manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('--config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')

    return (workflow, workflow.parse_args())


def main(workflow, args):
    conf = parse_cfg_file(args.config_file, section='wgs')

    contaminate_db = conf.get('database').get('knead_dna')
    project = conf.get('project')

    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')


    if data_files and data_files.get('WGS'):
        (input_files, output_files) = data_files.get('WGS').values()

        ## Create all the project directories we will need in the following 
        ## steps. These will be in the layout of:
        ##
        ##      <BASE DIRECTORY>/<PROJECT>/<TIME STAMP>/<DATE TYPE>
        ##
        ##      e.x. /tmp/data_deposition/HMP2/2017-04-14/WGS/
        project_dirs = create_project_dirs([conf.get('deposition_dir'),
                                            conf.get('processing_dir'),
                                            conf.get('public_dir')],
                                           project,
                                           'WGS')
        (deposition_dir, processing_dir, public_dir) = project_dirs

        ## We'll want the base deposition directory as well to drop our 
        ## MANIFEST file into for book-keeping purposes.
        base_depo_dir = os.path.abspath(os.path.join(deposition_dir, '..'))

        ## Our md5sums are going to be in the directories with their 
        ## respective files so we'll want to collate them somehow.
        ## TODO: Gather up all md5sums for files and create an md5sums file
        input_dirs = [group[0] for group in itertools.grouby(input_files, 
                                                             os.path.dirname)]
        md5sum_files = find_files(input_dirs, ext='.md5')
        merged_md5sum_file = tempfile.mktemp('.md5')
        checksums_file = create_merged_md5sum_file(md5sum_files, 
                                                   merged_md5sum_file)

        ## Validate files to ensure file integrity via md5 checksums
        validated_files = verify_files(workflow, input_files, 
                                       checksums_file, args.threads)

        manifest_file = stage_files(workflow,
                                    [args.manifest_file],
                                    base_depo_dir)

        ## WGS files are generated by the Broad and will live under 
        ## the /seq/picard_aggregation directory so we'll want to symlink
        ## them over to our data deposition directory.
        deposited_files = stage_files(workflow,
                                      validated_files,
                                      deposition_dir,
                                      symlink=True)

        fastq_files = bam_to_fastq(workflow, deposited_files, 
                                   processing_dir, args.threads)

        (cleaned_fastqs, read_counts) = quality_control(workflow, 
                                                        fastq_files, 
                                                        processing_dir,
                                                        args.threads,
                                                        contaminate_db)
        
        ## Generate taxonomic profile output. Output are stored in a list 
        ## and are the following:   
        ##
        ##      * Merged taxonomic profile 
        ##      * Individual taxonomic files
        ##      * metaphlan2 SAM files 
        tax_profile_outputs = taxonomic_profile(workflow,
                                                cleaned_fastqs,
                                                public_dir,
                                                args.threads)

        ## Generate functional profile output using humann2. Outputs are the 
        ## the following:
        ## 
        ##      * Merged gene families
        ##      * Merged relative abundances
        ##      * Merged pathway abundances
        func_profile_outputs = functional_profile(workflow,
                                                  cleaned_fastqs,
                                                  public_dir,
                                                  args.threads,
                                                  tax_profile_outputs[1])

        ## BIOM files are also distributed on the IBDMDB website so we need 
        ## to make sure we round those all up and move them to the public 
        ## directory
        biom_files = batch_convert_tsv_to_biom(workflow, tax_profile_outputs[1])
        tax_biom_files = stage_files(workflow,
                                     biom_files,
                                     public_dir)
        
        ## The functional files are a little trickier. We need to get 
        ## each of the individual files and package them up into a tarball
        norm_genefamilies = find_files(os.path.join(public_dir, 'genes'),
                                       'relab.tsv')
        norm_ecsfiles = find_files(os.path.join(public_dir, 'ecs'),
                                  '.relab.tsv')
        norm_path_files = find_files(os.path.join(public_dir, 'pathways'), 
                                     '.relab.tsv')
        
        func_tar_files = []
        for (genefam_file, ecs_file, path_file) in zip(norm_genefamilies,
                                                       norm_ecsfiles,
                                                       norm_path_files):
            sample_name = os.path.splitext(os.path.basename(genefam_file))
            tar_path = os.path.join(public_dir, "%s_humann2.tar" % sample_name)
            func_tar_file = tar_files(workflow, 
                                      [genefam_file, ecs_file, path_file],
                                      tar_path)
            func_tar_files.append(func_tar_file)

        ## TODO: Add static webpage visualization generation in here
        make_files_web_visible(workflow, 
                               cleaned_fastqs,
                               tax_profile_outputs[0:1],
                               tax_biom_files,
                               func_tar_files)

        workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
