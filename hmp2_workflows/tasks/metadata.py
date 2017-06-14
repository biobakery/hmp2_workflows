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
import shutil
import tempfile

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


def add_metadata_to_tsv(workflow, analysis_files, metadata_file,
                        id_col, col_replace, target_cols=[]):
    """Adds metadata to the top of a tab-delimited file. This function is
    meant to be called on analysis files to append relevant metadata to the 
    analysis output found in the file. An example can be seen below:

        
        sample  Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
        Age 87  78  3   2   32  10  39  96
        Cohort  Healthy Healthy Healthy Healthy IBD IBD IBD IBD
        Favorite_color  Yellow  Blue    Green   Yellow  Green   Blue    Green 
        Height  60  72  63  67  71  65  61  64
        Sex 0   1   0   1   1   0   1   0
        Smoking 0   0   1   0   1   1   1   0
        Star_Trek_Fan   1   1   0   0   1   0   0   1
        Weight  151 258 195 172 202 210 139 140
        Bacteria    1   1   1   1   1   1   1   1
        Bacteria|Actinobacteria|Actinobacteria  0.0507585   0.252153    0.161725   

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 workflow object.
        analysis_files (list): Target TSV's to add metadata too
        metadata_file (string): The path to the metadata file to pull from.
        id_col (string): The column name in the supplied metadata file to 
            attempt to subset on using ID's from the analysis file.
        col_replace (list): A list of string fragments that should be searched 
            for and replaced in either of the column headers of the analysis 
            or metadata files.
        target_cols (list): A list of columns to filter the metadata file on.
        

    Requires:
        None

    Returns: 
        list: A list containing the path to all modified files.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import metadata

        workflow = anadama2.Workflow()
    
        target_cols = ['age', 'sex', 'smoking']
        col_replace = ['_taxonomic_profile', '_functional_profile']
        out_files = metadata.add_metadata_to_tsv(workflow, ['/tmp/metaphlan2.out'], 
                                                 'External ID',
                                                 col_replace,
                                                 '/tmp/metadata.tsv',
                                                 target_cols)

        print out_files
        ## ['/tmp/metaphlan2.out']
    """
    metadata_df = pd.read_csv(metadata_file, dtype='str')
    
    def _workflow_add_metadata_to_tsv(task):
        analysis_file = task.depends[0].name
        pcl_out = task.targets[0].name

        analysis_df = pd.read_table(analysis_file, dtype='str', index_col=0)
     
        # Big assumption that our sample IDs are in the first line of our 
        # tab delimited analysis file. Could end poorly. Need a better
        # way to handle this.
        sample_ids = analysis_df.columns.tolist()
        if sample_ids == 1:
            # Something went wrong here.
            raise ValueError('Could not parse sample ID\'s:', 
                             sample_ids)
        if col_replace:
            new_ids = sample_ids
            for replace_str in col_replace:
                new_ids = [sid.replace(replace_str, '') for sid in new_ids]

            if new_ids != sample_ids:
                sample_ids_map = dict(zip(sample_ids, new_ids))
                sample_ids = new_ids

                analysis_df.rename(columns=sample_ids_map, inplace=True)

        subset_metadata_df = metadata_df[metadata_df[id_col].isin(sample_ids)]

        if target_cols:
            target_cols.insert(0, id_col)
            subset_metadata_df = subset_metadata_df.filter(target_cols)

        subset_metadata_df = subset_metadata_df.T
        subset_metadata_df.rename(columns=subset_metadata_df.iloc[0], inplace=True)
        subset_metadata_df.drop(subset_metadata_df.index[0], inplace=True)

        row_order = subset_metadata_df.index.tolist() + analysis_df.index.tolist()
        analysis_metadata_df = pd.concat([analysis_df, subset_metadata_df], 
                                         axis=0)
        analysis_metadata_df = analysis_metadata_df.reindex(row_order)
    
        analysis_metadata_df.to_csv(pcl_out, sep='\t', na_rep="NA")

    output_folder = os.path.dirname(analysis_files[0])
    pcl_files = bb_utils.name_files(analysis_files, 
                                    output_folder, 
                                    extension="pcl")

    workflow.add_task_group(_workflow_add_metadata_to_tsv,
                            depends = analysis_files,
                            targets = pcl_files,
                            time = 1*60,
                            mem = 1024,
                            cores = 1,
                            name =  "Generate analysis PCL output file")


    return pcl_files
