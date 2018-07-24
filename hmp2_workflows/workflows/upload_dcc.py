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

import itertools
import logging
import os
import sys
import tempfile

import cutlass
import numpy as np
import pandas as pd

from anadama2 import Workflow

from hmp2_workflows.utils.misc import (parse_cfg_file, parse_checksums_file, 
                                       create_merged_md5sum_file) 
from biobakery_workflows.utilities import find_files

from hmp2_workflows.utils import dcc
from hmp2_workflows.tasks.dcc import upload_data_files


def set_logging():
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)

set_logging()


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
    workflow.add_argument('baseline-metadata-file', desc='Metadata file '
                          'containing baseline visit metadata per subject.')                           
    workflow.add_argument('config-file', desc='Configuration file '
                          'containing parameters required by the workflow.')

    return workflow


def main(workflow):
    args = workflow.parse_args()
    conf = parse_cfg_file(args.config_file)
    data_type_mapping = conf.get('datatype_mapping')

    manifest = parse_cfg_file(args.manifest_file)
    data_files = manifest.get('submitted_files')

    metadata_df = pd.read_csv(args.metadata_file)
    baseline_metadata_df = pd.read_csv(args.baseline_metadata_file)
    md5sums_map = {}

    ## For every set of data files we're going to want to iterate over each  
    ## section and process the data files one section at a time, handling 
    ## the unique pieces of metadata for each data type as needed.
    if data_files:
        username = conf.get('username')
        password = conf.get('password')
        session = cutlass.iHMPSession(username, password, ssl=False)

        dcc_objs = []
        dcc_project = dcc.get_project(conf, session)
        dcc_study = dcc.crud_study(conf, 
                                   session,
                                   dcc_project.id)
        dcc_subjects = dcc.group_osdf_objects(dcc_study.subjects(),
                                              'rand_subject_id')
        dcc_subjects = dcc.crud_subjects(dcc_subjects, dcc_study, baseline_metadata_df, conf)

        for data_type in data_files:
            dtype_metadata = conf.get(data_type)

            input_files = data_files.get(data_type, {}).get('input')
            output_files = data_files.get(data_type, {}).get('output')
            md5sums_file = data_files.get(data_type).get('md5sums_file')
            file_tags = data_files.get(data_type).get('tags', [])

            if md5sums_file:
                md5sums_map.update(parse_checksums_file(md5sums_file))
            else:
                raise ValueError("MD5 checksums file is required.")

            ## Getting data files from different souces means that the 
            ## identifier we use to map a file to a piece of metadata may 
            ## be different. We need to account for this and create a map 
            ## of the 'universal ID' back to the specific file it references
            id_cols = conf['metadata_id_mappings'][data_type]
            seq_fname_map = dcc.create_seq_fname_map(data_type, input_files, tags=file_tags)
 
            sample_ids = seq_fname_map.keys()
            dtype_name = data_type_mapping.get(data_type)
            id_col = conf['metadata_id_mappings'][data_type]

            sample_metadata_df = metadata_df[(metadata_df[id_col].isin(sample_ids)) &
                                             (metadata_df['data_type'] == dtype_name)]

            ## Just in case there are more samples
            missing_samples = set(sample_ids) - set(sample_metadata_df[id_col].tolist())

            ## In our proteomics dataset we occasionally see two datasets tied to the same
            ## sample so we need to do a little extra work to figure out which dataset
            ## we are working with.
            is_proteomics = True if data_type == "proteomics" else False

            ## Add an extra column to our sample metadata with the corresponding sequence file 
            ## for easy access later on.
            #sample_metadata_df['seq_file'] = None
            sample_metadata_df = sample_metadata_df.apply(dcc.map_sample_id_to_file,
                                                          args=(id_col, seq_fname_map, 
                                                                is_proteomics),
                                                          axis=1)
            #sample_metadata_df = sample_metadata_df.dropna(axis=0, subset=['seq_file'])
            
            output_files_map = None
            if output_files:
                ## Do a bunch of stuff here since we have output files
                output_files_map = dcc.create_output_file_map(data_type, output_files, tags=file_tags)
    
            for (subject_id, metadata) in sample_metadata_df.groupby(['Participant ID']):
                dcc_subject = dcc_subjects.get(subject_id[1:])
                if dcc_subject:
                    dcc_subject = dcc_subject[0]
                else:
                    raise ValueError('Could not find Subject object for subject ID %s' % subject_id)                        

                dcc_visits = dcc.group_osdf_objects(dcc_subject.visits(),
                                                    'visit_id')
                
                for (idx, row) in metadata.iterrows():
                    dcc_visit = dcc.crud_visit(dcc_visits, 
                                               row['visit_num'],
                                               dcc_subject.id,
                                               data_type,
                                               row,
                                               conf)
                    dcc_visits.setdefault(dcc_visit.visit_id, []).append(dcc_visit)

                    dcc_samples = dcc.group_osdf_objects(dcc_visit.samples(),
                                                         'name')
                    dcc_sample = dcc.crud_sample(dcc_samples,
                                                 row.get('site_sub_coll'),
                                                 dcc_visit.id, 
                                                 conf,
                                                 row)

                    input_dcc_objs = []
                    if data_type == "MBX": 
                        dcc_prep = dcc.crud_host_assay_prep(dcc_sample, 
                                                            conf.get('data_study'),
                                                            data_type,
                                                            conf.get(data_type),
                                                            row)

                        ## MBX is a curious case in that we have multiple raw files (and multiple outputs) so 
                        ## we'll need to handle all these files (in another loop.) I need to refactor this 
                        ## to be way more elegant in the future.
                        ##
                        ## This is going to be a bit trickier than other data types because of that 
                        input_metabolomes = row.filter(like='metabolome')

                        for (metabolome_type, metabolome_file) in input_metabolomes.iteritems():
                            metabolome_fname = os.path.basename(metabolome_file[0])
                            analysis_type = metabolome_type.split('_', 1)[-1]

                            dcc_seq_obj = dcc.crud_metabolome(dcc_prep,
                                                              metabolome_file[0],
                                                              md5sums_map.get(metabolome_fname),
                                                              dcc_sample.name,
                                                              conf.get('data_study'),
                                                              dtype_metadata,
                                                              row)
                            input_dcc_objs.append(dcc_seq_obj)
                    elif data_type == "MPX":
                        #url_param = '_raw_url'
                        dcc_prep = dcc.crud_microb_assay_prep(dcc_sample,
                                                                conf.get('data_study'),
                                                                data_type,
                                                                dtype_metadata,
                                                                row)
                        dcc_seq_obj = dcc.crud_proteome(dcc_prep,
                                                        file_md5sum,
                                                        dcc_sample.name,
                                                        dtype_metadata,
                                                        row) 
                    elif data_type == "HTX":
                        dcc_prep = dcc.crud_host_seq_prep(dcc_sample,
                                                            conf.get('data_study'),
                                                            data_type,
                                                            dtype_metadata,
                                                            row)
                        dcc_seq_obj = dcc.crud_host_tx_raw_seq_set(dcc_prep,
                                                                    file_md5sum,
                                                                    dcc_sample.name,
                                                                    conf.get(data_type),
                                                                    row)
                    elif data_type == "HG":
                        dcc_prep =  dcc.crud_host_seq_prep(dcc_sample,
                                                            conf.get('data_study'),
                                                            data_type,
                                                            dtype_metadata,
                                                            row)
                        dcc_seq_obj = dcc.crud_host_wgs_raw_seq_set(dcc_prep,
                                                                    file_md5sum,
                                                                    dcc_sample.name,
                                                                    dtype_metadata,
                                                                    row)
                    elif data_type == "MTX":
                        dcc_prep = dcc.crud_wgs_dna_prep(dcc_sample,
                                                            conf.get('data_study'),
                                                            data_type,
                                                            dtype_metadata,
                                                            row)
                        dcc_seq_obj = dcc.crud_microb_transcriptomics_raw_seq_set(dcc_prep,
                                                                                    file_md5sum,
                                                                                    dcc_sample.name,
                                                                                    dtype_metadata,
                                                                                    row)
                    elif data_type == "MGX":
                        dcc_prep = dcc.crud_wgs_dna_prep(dcc_sample,
                                                            conf.get('data_study'),
                                                            data_type,
                                                            dtype_metadata,
                                                            row)
                        dcc_seq_obj = dcc.crud_wgs_raw_seq_set(dcc_prep,
                                                                file_md5sum,
                                                                dcc_sample.name,
                                                                dtype_metadata,
                                                                row)
                    elif data_type == "MVX":
                        dcc_prep = dcc.crud_wgs_dna_prep(dcc_sample,
                                                            conf.get('data_study'),
                                                            data_type,
                                                            dtype_metadata,
                                                            row) 
                        dcc_seq_obj = dcc.crud_wgs_raw_seq_set(dcc_prep,
                                                                file_md5sum,
                                                                dcc_sample.name,
                                                                dtype_metadata,
                                                                row,
                                                                private=True)
                    elif data_type == '16SBP':
                        raw_seq_set_fname = os.path.basename(row.get('16S_raw_seq_set'))
                        trimmed_seq_set_fname = os.path.basename(row.get('16S_trimmed_seq_set'))

                        dcc_prep = dcc.crud_sixs_dna_prep(dcc_sample,
                                                          conf.get('data_study'),
                                                          data_type,
                                                          dtype_metadata,
                                                          row)
                        dcc_raw_seq_set = dcc.crud_sixs_raw_seq_set(dcc_prep,
                                                                    md5sums_map.get(raw_seq_set_fname),
                                                                    dtype_metadata,
                                                                    row)
                        dcc_trimmed_seq_set = dcc.crud_sixs_trimmed_seq_set(dcc_raw_seq_set,
                                                                            md5sums_map.get(trimmed_seq_set_fname),
                                                                            dtype_metadata,
                                                                            row)
                        input_dcc_objs.extend([dcc_raw_seq_set, dcc_trimmed_seq_set])
                    elif data_type == 'RRBS':
                        raw_epigenetics_seq_set = row.get('host_epigenetics_raw_seq_set')
                        raw_epigenetics_fname = os.path.basename(raw_epigenetics_seq_set[0])

                        dcc_prep = dcc.crud_host_seq_prep(dcc_sample,
                                                            conf.get('data_study'),
                                                            data_type,
                                                            dtype_metadata,
                                                            row)
                        dcc_seq_obj = dcc.crud_host_epigenetics_raw_seq_set(session,
                                                                            dcc_prep,
                                                                            raw_epigenetics_seq_set,
                                                                            md5sums_map.get(raw_epigenetics_fname),
                                                                            dcc_sample.name,
                                                                            conf.get('data_study'),
                                                                            conf.get(data_type),
                                                                            row)
                    elif data_type == 'SER':
                        dcc_prep = dcc.crud_host_assay_prep(dcc_sample, 
                                                            conf.get('data_study'),
                                                            data_type,
                                                            conf.get(data_type),
                                                            row)
                        dcc_seq_obj = dcc.crud_serology(session,
                                                        dcc_prep, 
                                                        file_md5sum,
                                                        dcc_sample.name,
                                                        conf.get('data_study'),
                                                        dtype_metadata,
                                                        row)

                    if len(input_dcc_objs) == 0:
                        input_dcc_objs.append(dcc_seq_obj)

                    ##uploaded_files = upload_data_files(workflow, input_dcc_objs)

                    ## The only output type currently supported are AbundanceMatrices 
                    ## so those are the only we will work with. Short-sided and 
                    ## ugly but can re-work this later.
                    if output_files_map and row.get('External ID') in output_files_map:
                        seq_out_files = output_files_map.get(row.get('External ID'))
                        dcc_output_objs = []

                        for (output_ftype, output_file) in seq_out_files.iteritems():
                            output_filename = os.path.basename(output_file[0])
                            output_md5sum = md5sums_map.get(output_filename)


                            if not output_md5sum:
                                raise ValueError("Could not find md5sum for file", output_filename)

                            if data_type == "MBX" and file_tags:
                                output_base = os.path.splitext(output_filename)[0]

                                ## MBX data is a bit tricky since we can have multiple inputs and outputs
                                ## that need to be threaded together.
                                analysis_type = output_base.split('_', 1)[-1]
                                dcc_parent_obj = next((p for p in input_dcc_objs if analysis_type in p.urls[0]), None)
                            else:
                               dcc_parent_obj = input_dcc_objs[-1]

                            ## We need a special case here when dealing with Host Genomes... 
                            ## TODO: Clean this up to make this a lot better...
                            if data_type == "HG":
                                dcc_output_obj  = dcc.crud_host_variant_call(session,
                                                                             dcc_parent_obj,
                                                                             output_file[0],
                                                                             output_md5sum,
                                                                             conf.get('data_study'),
                                                                             dtype_metadata,
                                                                             row)
                            else:
                                dcc_output_obj = dcc.crud_abundance_matrix(session,
                                                                           dcc_parent_obj,
                                                                           output_file[0],
                                                                           output_md5sum,
                                                                           dcc_sample.name,
                                                                           conf.get('data_study'),
                                                                           dtype_metadata,
                                                                           row)

                            dcc_output_objs.append(dcc_output_obj)

                        uploaded_file = upload_data_files(workflow, dcc_output_objs)


if __name__ == "__main__":
    main(parse_cli_arguments())
