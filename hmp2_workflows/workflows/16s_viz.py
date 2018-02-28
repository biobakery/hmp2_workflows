# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.16s_vis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles visualzation of output from the HMP2 
16S metagenome pipeline. 

This workflow generates a static HTML summary page that summarizes the
results for easier consumption

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


import os

from glob2 import glob

from anadama2 import Workflow

from hmp2_workflows import document_templates
from biobakery_workflows import utilities, files


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
    workflow = Workflow(version='1.0', description='A workflow to handle visualization '
                        'of HMP2 16S data.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('metadata-file', desc='Accompanying metadata file '
                           'for the provided data files.', default=None)
    workflow.add_argument('source', desc='The source of the output files generated. '
                          '[biobakery, CMMR]', default='biobakery', 
                          choices=['biobakery', 'CMMR'])


    return workflow


def main(workflow):
    args = workflow.parse_args()

    ## Because we accept files analyzed by Baylor here or our own files analyzed via the biobakery worfklow 
    ## derivation we need to be able to handle either set of data. This is specified by the source parameter 
    ## provided to the viz script.
    templates = []
    vars = {}

    # This file should exist in either scenario
    eestats_table = glob(os.path.join(args.input + "/**/",
        files.SixteenS.file_info['eestats2'].keywords.get('names')))[0]

    if args.source == 'biobakery':
        otu_table = glob(os.path.join(args.input + '/**/', 
            files.SixteenS.file_info['otu_table_closed_reference'].keywords.get('names')))[0]
        otu_table_open = glob(os.path.join(args.input + '/**/',
            files.SixteenS.file_info['otu_table_open_reference'].keywords.get('names')))[0]
        read_counts_table = glob(os.path.join(args.input + "/**/",
            files.SixteenS.file_info['read_count_table'].keywords.get('names')))[0]
        centroid_fasta = glob(os.path.join(args.input + "/**/",
            files.SixteenS.file_info['msa_nonchimera'].keywords.get('names')))[0]
        centroid_closed_fasta = glob(os.path.join(args.input + "/**/",
            files.SixteenS.file_info['msa_closed_reference'].keywords.get('names')))[0]
        centroid_closed_fasta = glob(os.path.join(args.input + "/**/",
            files.SixteenS.file_info['msa_closed_reference'].keywords.get('names')))[0]
        log_file = files.Workflow.path('log', args.input)

        dependencies = [otu_table, otu_table_open, read_counts_table, eestats_table,
                        centroid_fasta, centroid_closed_fasta, log_file]
        vars.update({
            'log': log_file,
            'read_count_table': read_counts_table,
        })

        templates.append(document_templates.get_template('header_16S'))
        templates.append(document_templates.get_template('quality_control_16S_CMMR'))
    elif args.source == 'CMMR':
        otu_table = glob(os.path.join(args.input + '/**/' + 'OTU_Table.tsv'))
        centroid_fasta = glob(os.path.join(args.input + '/**/', 'CentroidInformation.fa'))
        dependencies = [otu_table, eestats_table, centroid_fasta]

        templates.append(document_templates.get_template('header_16S_CMMR'))
        templates.append(document_templates.get_template('quality_control_16S_CMMR'))

    templates.append(document_templates.get_template('taxonomy_16S'))
    templates.append(document_templates.get_template('footer'))

    vars.update({
        'summary_title': "HMP2: 16S Data Summary Report",
        'otu_table': otu_table,
        'eestats_table': eestats_table
    })

    doc_task = workflow.add_document(
        templates = templates,
        depends = dependencies,
        targets = workflow.name_output_files("summary.html"),
        vars = vars,
        table_of_contents = False
    )

    workflow.add_task("sed -i -e 's/ xmlns=\"http:\/\/www.w3.org\/1999\/xhtml\" lang=\"\" xml:lang=\"\"//' [depends[0]]",
                      depends = workflow.name_output_files('summary.html'),
                      targets = workflow.name_output_files('summary.html'))

    workflow.add_task("sed -i -e '/<header>/,/<\/header>/d' [depends[0]]",
                      depends = workflow.name_output_files('summary.html'),
                      targets = workflow.name_output_files('summary.html'))

    workflow.add_task("sed -i -e '/img.* alt/! s/img/img alt=\"\"/' [depends[0]];",
                      depends = workflow.name_output_files('summary.html'),
                      targets = workflow.name_output_files('summary.html'))

    workflow.add_task("rm [depends[0]]",
                      depends = workflow.name_output_files('summary.html-e'))

    workflow.add_task('tidy -q -i --wrap 0 -m [depends[0]]',
                      depends = workflow.name_output_files('summary.html'),
                      targets = workflow.name_output_files('summary.html'))

    workflow.go()


if __name__ == "__main__":
    main(parse_cli_arguments())
