# -*- coding: utf-8 -*-

"""
hmp2_workflows.workflows.upload_dcc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An AnADAMA2 workflow that handles uploading metadata and data files to 
Data Coordination Center (DCC). 

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
import tempfile

import cutlass
import pandas as pd

from anadama2 import Workflow

from hmp2_workflows.utils.misc import (parse_cfg_file, parse_checksums_file, 
                                       create_merged_md5sum_file) 
from biobakery_workflows.utilities import find_files

from hmp2_workflows.utils import dcc
from hmp2_workflows.tasks.dcc import upload_data_files


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
    workflow = Workflow(version='0.1', description='A workflow to handle '
                        'uploading metadata and data files to the Data '
                        'Coordination Center (DCC)', 
                        remove_options=['input', 'output'])
    workflow.add_argument('manifest-file', desc='Manifest file containing '
                          'files to process in this workflow run.')
    workflow.add_argument('metadata-file', desc='Accompanying metadata file '
                           'for the provided data files.', default=None)
    workflow.add_argument('broad-data-sheet', desc='Master sample data sheet '
                          'tracking all samples and their reception dates.')
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')

    return workflow


def main(workflow):
    args = workflow.parse_args()
    conf = parse_cfg_file(args.config_file, section='dcc')
    
    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')

    metadata_df = pd.read_csv(args.metadata_file)
    coll_df = pd.read_csv(args.broad_data_sheet)
    md5sums_map = {}

    ## For every set of data files we're going to want to iterate over each  
    ## section and process the data files one section at a time, handling 
    ## the unique pieces of metadata for each data type as needed.
    if data_files:
        username = conf.get('username')
        password = conf.get('password')
        session = cutlass.iHMPSession(username, password)

        dcc_objs = []
        dcc_project = dcc.get_project(conf, session)
        dcc_study = dcc.get_or_update_study(conf, 
                                            session,
                                            dcc_project.id)
        dcc_subjects = dcc.group_osdf_objects(dcc_study.subjects(),
                                              'rand_subject_id')

        for data_type in data_files:
            dtype_metadata = conf.get(data_type.lower())

            md5sums_file = data_files.get(data_type).get('md5sums_file')
            if md5sums_file:
                md5sums_map.update(parse_checksums_file(md5sums_file))
            else:
                ## If the files we are working with are provided by the Broad
                ## we need to gather up all the md5sums into one file.
                input_files = data_files['data_type']['input']
                input_dirs = [group[0] for group in itertools.groupby(input_files, 
                                                                      os.path.dirname)]
                md5sum_files = chain.from_iterable([find_files(input_dir, '.md5') for 
                                                    input_dir in input_dirs])
                merged_md5sum_file = tempfile.mktemp('.md5')
                checksums_file = create_merged_md5sum_file(md5sum_files, 
                                                           merged_md5sum_file)
                md5sums_map.update(parse_checksums_file(checksums_file))

            ## Getting data files from different souces means that the 
            ## identifier we use to map a file to a piece of metadata may 
            ## be different. We need to account for this and create a map 
            ## of the 'universal ID' back to the specific file it references
            id_col = conf['metadata_id_mappings'][data_type]
            sample_map = dcc.create_sample_map(data_type,
                                               data_files[data_type]['input'])
            sample_ids = sample_map.keys()
            sample_metadata_df = metadata_df[(metadata_df[id_col].isin(sample_ids)) &
                                             (metadata_df['data_type'] == data_type)]

            ## Once we've mapped to our samples to our metadata we need to 
            ## grab another set of ID's so we can map to the Broad sample 
            ## tracking sheet.
            tracking_map_col = conf.get('metadata_to_tracking_mapping')
            tracking_id_col = conf.get('tracking_id_col')

            tube_ids = sample_metadata_df[tracking_map_col]
            samples_coll_df = coll_df[coll_df[tracking_id_col].isin(tube_ids)]
            samples_coll_df = samples_coll_df.set_index(tracking_id_col)
            samples_coll_df.index.names = [None]

            for (subject_id, metadata) in sample_metadata_df.groupby(['Participant ID']):
                dcc_subject = dcc.create_or_update_subject(dcc_subjects,
                                                           subject_id.replace('C', ''),
                                                           dcc_study.id,
                                                           metadata,
                                                           conf)
                dcc_visits = dcc.group_osdf_objects(dcc_subject.visits(),
                                                    'visit_number')
                
                for (idx, row) in metadata.iterrows():
                    sample_coll_row = samples_coll_df.xs(row.get(tracking_map_col))

                    dcc_visit = dcc.create_or_update_visit(dcc_visits, 
                                                           row.get('visit_num'),
                                                           dcc_subject.id, 
                                                           row)
                    dcc_sample = dcc.create_or_update_samples(dcc_visit,
                                                              sample_coll_row,
                                                              row)

                    data_file = sample_map.get(dcc_sample.name)
                    if data_file:
                        if   data_type == "proteomics": 
                            dcc_prep = dcc.create_or_update_microbiome_prep(dcc_sample,
                                                                            conf,
                                                                            row)
                        elif data_type == "MTX":
                            dcc_prep = dcc.create_or_update_wgs_dna_prep(dcc_sample,
                                                                         conf,
                                                                         row)
                        elif data_type == "WGS":
                            dcc_prep = dcc.create_or_update_wgs_dna_prep(dcc_sample,
                                                                         conf,
                                                                         row)
                        elif data_type == "16S":
                            dcc_prep = dcc.create_or_update_16s_dna_prep(dcc_sample,
                                                                         conf,
                                                                         row)

                        dcc_objs.extend([dcc_project, dcc_study, dcc_subject,
                                         dcc_visit, dcc_sample, dcc_prep])
     
            dcc_files = upload_data_files(workflow, dtype_metadata, dcc_objs)
        

if __name__ == "__main__":
    main(parse_cli_arguments())