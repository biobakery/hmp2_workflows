# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.wmgx_vis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles visualzation of output from the HMP2 
WGS metagenome pipeline. 

This workflow generates a static HTML summary page that summarizes the
results for easier consumption

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
                        'of HMP2 WGS metagenome data.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')
    workflow.add_argument('metadata-file', desc='Accompanying metadata file '
                          'for the provided data files.', default=None)


    return workflow


def main(workflow):
    args = workflow.parse_args()

    taxonomic_profile = glob(os.path.join(args.input + '/**/', 
        files.ShotGun.file_info['taxonomic_profile'].keywords.get('names')))[0]
    dna_read_counts = glob(os.path.join(args.input + '/**/',
        files.ShotGun.file_info['kneaddata_read_counts'].keywords.get('names')))[0]
    pathabundance = glob(os.path.join(args.input + "/**/",
        files.ShotGun.file_info['pathabundance_relab'].keywords.get('names')))[0]
    read_counts = glob(os.path.join(args.input + "/**/",
        files.ShotGun.file_info['humann2_read_counts'].keywords.get('names')))[0]
    feature_counts = glob(os.path.join(args.input + "/**/",
        files.ShotGun.file_info['feature_counts'].keywords.get('names')))[0]

    templates = []
    templates.append(document_templates.get_template('header'))
    templates.append(document_templates.get_template('wmgx'))
    templates.append(document_templates.get_template('quality_control_dna'))
    templates.append(document_templates.get_template('taxonomy'))
    templates.append(document_templates.get_template('functional_dna'))
    templates.append(document_templates.get_template('footer'))

    doc_task = workflow.add_document(
        templates = templates,
        depends = [taxonomic_profile],
        targets = workflow.name_output_files("summary.html"),
        vars = {
            'summary_title': "HMP2: Metagenomes Data Summary Report",
            'metadata_file': args.metadata_file,
            'taxonomic_profile': taxonomic_profile,
            'dna_read_counts': dna_read_counts,
            'dna_pathabundance': pathabundance,
            'read_counts': read_counts,
            'feature_counts': feature_counts
        },
        table_of_contents = False
    )

    workflow.add_task("sed -i -e 's/ xmlns=\"http:\/\/www.w3.org\/1999\/xhtml\" lang=\"\" xml:lang=\"\"//' [depends[0]]",
                      depends = workflow.name_output_files('summary.html'),
                      targets = workflow.name_output_files('summary.html'))

    workflow.add_task("sed -i -e '/<header>/,/<\/header>/d' [depends[0]]",
                      depends = workflow.name_output_files('summary.html'),
                      targets = [workflow.name_output_files('summary.html'), 
                                 workflow.name_output_files('summary.html-e')])

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
