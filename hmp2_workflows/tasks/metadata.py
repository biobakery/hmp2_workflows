# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains functions required to update, validate, refresh and 
disseminate HMP2 metadata  via the IBDMDB website and the DCC.

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

import pandas as pd

from biobakery_workflows import utilities as bb_utils


def validate_metadata_file(workflow, input_file, validation_file):
    """Validates an HMP2 metadata file using the cutplace utility. 
    A valid cutplace interface definition file must exist for the provided 
    input metadata file.

    Args: 
        workflow (anadama2.Workflow): The workflow object.
        input_file (string): Path to input metadata file in CSV format.
        validation_file (string): Path to cutplace intrface definition file 
            used to validate the provided input file.
 
    Requires:
        cutplace: Python module/binary to facilitate validation.

    Returns:
        boolean: True or False if the file is valid.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import metadata

        workflow = anadama2.Workflow()

        is_valid = metadata.validate_file(workflow, '/tmp/hmp2.csv', 
                                          '/tmp/hmp2.cid')

        workflow.go()
    """
    if not os.path.exists(input_file):
        raise OSError(2, 'Input CSV file does not exist', input_file)
    if not os.path.exists(validation_file):
        raise OSError(2, 'Input CID file does not exists', validation_file)

    workflow.add_task_gridable('cutplace [depends[0] [depends[1]]',
                               depends=[input_file, validation_file])


def generate_sample_metadata(workflow, data_type, samples, metadata_file, 
                             output_dir, id_column = 'External ID'):
    """Generates a series of individual metadata files in CSV format 
    from the provided merged metadata file. Each of the provided samples
    has a metadata file generated to accompany any product files generated 
    by the analysis pipelines.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        data_type (string): The data type of the provided samples. One of 
            either 'metageonimcs', 'proteomices', 'amplicon'.
        samples (list): A list of all samples that will have individual 
            metadata files written.
        metadata_file (string): Path to the merged metadata file.
        id_column (string): The ID column to attempt to map sample names to 
            in the merged metadata file. Default set to "External ID" but 
            can change depending on the data type.
        output_dir (string): Path to output directory to write each 
            sample metadata file too.

    Requires:
        None

    Returns:
        list: A list containing the path to all sample metadata files created.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import metadata

        workflow = anadama2.Workflow()
        
        samples = ['sampleA', 'sampleB']
        metadadta_file = '/tmp/merged_metadata.csv'
        output_dir = '/tmp/metadata'

        metadata_files = metadata.generate_sample_metadata(workflow, 
                                                           'metagenomics',
                                                           samples, 
                                                           metadata_file, 
                                                           output_dir)
        print metadata_files
        ## ['/tmp/metadata/sampleA.csv', '/tmp/metadata/sampleB.csv']
    """
    metadata_df = pd.read_csv(metadata_file)
    
    output_metadata_files = bb_utils.name_files(samples, 
                                                output_dir, 
                                                extension = 'csv',
                                                subfolder = 'metadata',
                                                create_folder = True)
    sample_metadata_dict = dict(zip(samples, output_metadata_files))

    def _workflow_gen_metadata(task):
        metadata_subset = metadata_df.loc[(metadata_df[id_column].isin(samples)) &
                                          (metadata_df['data_type'] == data_type)]
    
        if metadata_subset.empty:
            raise ValueError('Could not find metadata associated with samples.',
                             ",".join(samples))

        for (sample_id, row) in metadata_subset.iterrows():
            sample_metadata_file = sample_metadata_dict.get(row[id_column])
            metadata_subset.xs(sample_id).to_csv(sample_metadata_file, index=False)
    
    workflow.add_task(_workflow_gen_metadata,
                      targets = output_metadata_files,  
                      depends = [metadata_file],
                      name = 'Generate sample metadata')

    return sample_metadata_dict.values()

