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
import tempfile

import funcy
import pandas as pd

from anadama2.tracked import Container

import hmp2_workflows.utils.metadata as m_utils

from biobakery_workflows import utilities as bb_utils

from hmp2_workflows.utils.misc import (get_sample_id_from_fname, 
                                       reset_column_headers)

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

    workflow.add_task_gridable('cutplace [depends[0]] [[depends[1]]',
                               depends=[input_file, validation_file])


def generate_metadata_file(workflow, config, data_files, studytrax_metadata, 
                           broad_tracking_sheet, proteomics_metadata=None,
                           auxillary_metadata=[]):
    """Generates a metadata file for the set of sequence files provided 
    using the StudyTrax metadata sheet, Broad Sample Tracking data 
    sheet and any auxillary metadata files provided.

    Any auxillary metadata files provided should have the column to merge
    into the StuyTrax and Broad tracking sheets in the first column.

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 workflow object.
        config (dict): A dictionary containing metadata configuration 
            parameters used in metadata generation.
        data_files (list): The set of files to retrieve metadata for. 
            These files should be in format [<DATATYPE>: [FILES]].
        studytrax_metadata (string): Path to StudyTrax clinical metadata 
            file.
        broad_tracking_sheet (string): Path to Broad sample tracking 
            spreadsheet file.
        proteomics_metadata (string): Path to any proteomics metadata to 
            accompany PNNL proteomics data.
        auxillary_metadat (list): Any auxillary metadata files to pull
            metadata and integrate into an HMP2 metadata file.                                  
    """
    metadata_files = []
    temp_dir = tempfile.mkdtemp(prefix='generate_metadata_file')

    def _generate_metadata_file(task):
        input_files = [seq_file.name for seq_file in task.depends[:-3]]
        studytrax_metadata = task.depends[-4].name
        broad_sample_sheet = task.depends[-3].name
        auxillary_metadata = task.depends[-1].name if task.depends[-1].name != '/dev/null' else None
        metadata_out_file = task.targets[0].name

        data_type_map = config.get('dtype_mapping')

        studytrax_df = pd.read_csv(studytrax_metadata)
        broad_sample_df = pd.read_csv(broad_sample_sheet, 
                                      na_values=['destroyed', 'missed'],
                                      parse_dates=['Actual Date of Receipt'])

        collection_dates_dict = m_utils.get_collection_dates(broad_sample_df)

        if pair_identifier:
            (input_pair1, input_pair2) = bb_utils.paired_files(input_files, pair_identifier)
            input_files = input_pair1 if input_pair1 else input_files

        sample_mapping = dict(zip(bb_utils.sample_names(input_files, pair_identifier),
                               map(get_sample_id_from_fname, input_files)))        
        sample_ids = [sid.replace(pair_identifier, '') for sid in sample_mapping.values()]

        sample_subset_df = broad_sample_df[(broad_sample_df['Parent Sample A'].isin(sample_ids)) |
                                           (broad_sample_df['Proteomics'].isin(sample_ids)) |
                                           (broad_sample_df['MbX'].isin(sample_ids)) |
                                           (broad_sample_df['Site/Sub/Coll']).isin(sample_ids)]        

        metadata_df = sample_subset_df.merge(studytrax_df,
                                            left_on='Parent Sample A',
                                            right_on='st_q4',
                                            how='left')

        ## We sometimes get a situation where our studytrax metadata is missing
        ## some of the proteomics sample ID's so we need to make sure we replicate
        ## them
        metadata_df.loc[metadata_df['st_q17'].isnull(), 'st_q17'] = metadata_df['Proteomics']
        metadata_df.loc[metadata_df['st_q11'].isnull(), 'st_q11'] = metadata_df['MbX']
        metadata_df['data_type'] = data_type_map.get(data_type)

        if proteomics_metadata:
            proteomics_df = m_utils.add_proteomics_metadata(sample_subset_df, 
                                                            proteomics_metadata,
                                                            sample_mapping)
            metadata_df = metadata_df.merge(proteomics_df,
                                            on='Parent Sample A',
                                            how='left')


        metadata_df['External ID'] = new_metadata_df.apply(generate_external_id, axis=1)

        metadata_df['Site/Sub/Coll ID'] = metadata_df['Site/Sub/Coll'].map(lambda sid: str(sid))
        metadata_df['Site'] = metadata_df['SiteName']
        metadata_df['Participant ID'] = metadata_df['Subject'].map(lambda subj: 'C' + str(subj))
        metadata_df['visit_num'] = metadata_df['Collection #']
        metadata_df['Research Project'] = config.get('research_project')
        metadata_df['Project'] = metadata_df.apply(m_utils.get_project_id, axis=1)
        metadata_df = generate_collection_statistics(metadata_df,
                                                     collection_dates_dict)
        metadata_df = metadata_df.drop(config.get('drop_cols'), axis=1, inplace=True)
 
        if auxillary_metadata:
            ## Auxillary metadata are columns that will be added into our
            ## existing metadata rows. 
            metadata_df = m_utils.add_auxiliary_metadata(metadata_df,auxillary_metadata)

        metadata_df.to_csv(metadata_out_file, index=False)


    for (data_type, items) in data_files.iteritems():
        metadata_file = tempfile.mkstemp(prefix='%s_metadata_' % data_type,
                                         suffix='.csv',
                                         dir=temp_dir)[-1]
        pair_identifier = items.get('pair_identifier', '')
        sequence_files = items.get('input')

        # Setup the group of dependencies the task needs
        metadata_dependencies = [studytrax_metadata, broad_tracking_sheet]

        if auxillary_metadata:
           metadata_dependencies.extend(auxillary_metadata)
        else:
            metadata_dependencies.append('/dev/null')

        workflow.add_task(_generate_metadata_file,
                                   targets=[metadata_file],
                                   depends=sequence_files + 
                                           metadata_dependencies,
                                   time="1*60",
                                   memory="2*1024",
                                   cores=1,
                                   name="Generate metadata file")

        metadata_files.append(metadata_file)                                   

    return metadata_files


def generate_sample_metadata(workflow, data_type, in_files, metadata_file, 
                             output_dir, id_column = 'External ID'):
    """Generates a series of individual metadata files in CSV format 
    from the provided merged metadata file. Each of the provided samples
    has a metadata file generated to accompany any product files generated 
    by the analysis pipelines.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        data_type (string): The data type of the provided samples. One of 
            either 'metageonimcs', 'proteomices', 'amplicon'.
        in_files (list): A list of files that should have corresponding 
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
    samples = bb_utils.sample_names(in_files)

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
                      targets=output_metadata_files,  
                      depends=in_files + [metadata_file],
                      name='Generate sample metadata')

    return sample_metadata_dict.values()


def add_metadata_to_tsv(workflow, analysis_files, metadata_file, dtype,
                        id_col, col_replace=None, col_offset=-1, 
                        metadata_rows=None, target_cols=None, 
                        aux_files=None, na_rep=""):
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
        dtype (string): Data type of files for which metadata is being refreshed
            to include.
        id_col (string): The column name in the supplied metadata file to 
            attempt to subset on using ID's from the analysis file.
        col_replace (list): A list of string fragments that should be searched 
            for and replaced in either of the column headers of the analysis 
            or metadata files.
        col_offset (int): In certain situations a series of metadata columns
            will be present prior to columns containing analysis results.
            In these cases an offset needs to be provided for proper creation 
            of PCL files.
        metadata_rows (int): If our analysis file already contains some 
            metadata files at the top of the file (in effect already a PCL
            file) this parameter indicates how many rows of metadata exist.
        target_cols (list): A list of columns to filter the metadata file on.
        aux_files (list): Any additional metadata files to integrate into 
            analysis files. 
        na_rep (string): String representation for any empty cell in our 
            PCL file. Defaults to an empty string.

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
    metadata_df = pd.read_csv(metadata_file, dtype='str',
                              parse_dates=['date_of_receipt'])
    
    def _workflow_add_metadata_to_tsv(task):
        analysis_file = task.depends[0].name
        pcl_out = task.targets[0].name

        analysis_df = pd.read_csv(analysis_file, dtype='str', header=None)
        pcl_metadata_df = None
        header = True
            
        # Going to make the assumption that the next row following our PCL 
        # metadata rows is the row containing the ID's that we will use to merge
        # the analysis file with our metadata file and we can use these same 
        # ID's to merge the PCL metadata rows into the larger metadata file.
        if metadata_rows:
            pcl_metadata_df = analysis_df[:metadata_rows+1]
            header = None

            offset_cols = range(0, col_offset+1)
            pcl_metadata_df.drop(pcl_metadata_df.columns[offset_cols[:-1]], 
                                    axis=1,
                                    inplace=True)

            pcl_metadata_df = pcl_metadata_df.T.reset_index(drop=True).T
            pcl_metadata_df.xs(metadata_rows)[0] = id_col
            
            pcl_metadata_df = pcl_metadata_df.T
            pcl_metadata_df = reset_column_headers(pcl_metadata_df)

            analysis_df.drop(analysis_df.index[range(0,metadata_rows)], inplace=True)
            analysis_df.rename(columns=analysis_df.iloc[0], inplace=True)
        else:
            analysis_df = hmp2_utils.misc.reset_column_headers(analysis_df)

        sample_ids = analysis_df.columns.tolist()[col_offset+1:]
            
        if len(sample_ids) == 1:
            raise ValueError('Could not parse sample ID\'s:', 
                             sample_ids)

        if col_replace:
            new_ids = sample_ids
            for replace_str in col_replace:
                new_ids = [sid.replace(replace_str, '') if not pd.isnull(sid) 
                           else sid for sid in new_ids]

            if new_ids != sample_ids:
                sample_ids_map = dict(zip(sample_ids, new_ids))
                sample_ids = new_ids
    
                analysis_df.rename(columns=sample_ids_map, inplace=True)

        subset_metadata_df = metadata_df[(metadata_df.data_type == dtype) &
                                         (metadata_df[id_col].isin(sample_ids))]

        if aux_files:
            for aux_file in aux_files:
                aux_metadata_df = pd.read_table(aux_file, dtype='str')
                aux_metadata_cols = aux_metadata_df.columns.tolist()
                join_id = aux_metadata_cols[0]

                ## We need to do this in two stages. If the columns already exist
                ## here we want to update them. If they do not exist we append
                ## them.
                subset_metadata_cols = subset_metadata_df.columns.tolist()
                new_cols = set(aux_metadata_cols[1:]) - set(subset_metadata_cols)
                existing_cols = set(aux_metadata_cols[1:]).intersection(subset_metadata_cols)

                if new_cols:
                    aux_metadata_new_df = aux_metadata_df.filter(items=aux_metadata_cols[:1] + 
                                                                 list(new_cols))
                    subset_metadata_df = pd.merge(subset_metadata_df, aux_metadata_new_df, 
                                                  how='left', on=join_id)

                if existing_cols:
                    aux_metadata_existing_df = aux_metadata_df.filter(items=aux_metadata_cols[:1] + 
                                                                      list(existing_cols))
                    subset_metadata_df.set_index(join_id, inplace=True)
                    aux_metadata_existing_df.set_index(join_id, inplace=True)

                    subset_metadata_df.update(aux_metadata_existing_df)
                    subset_metadata_df.reset_index(inplace=True)

        if not pcl_metadata_df.empty:
            subset_metadata_df = pd.merge(subset_metadata_df, pcl_metadata_df,
                                          how='left', on=id_col)

        if target_cols:
            target_cols.insert(0, id_col)
            subset_metadata_df = subset_metadata_df.filter(target_cols)

        subset_metadata_df = subset_metadata_df.T
        subset_metadata_df = reset_column_headers(subset_metadata_df)
        subset_metadata_df = subset_metadata_df.reset_index()
        subset_metadata_df.fillna('NA', inplace=True)

        _col_offset = col_offset-1 if col_offset != -1 else col_offset
        col_name = analysis_df.columns[_col_offset+1]

        col_name = '' if col_name == "index" else col_name
        subset_metadata_df.rename(columns={'index': col_name}, inplace=True)

        analysis_df.index = analysis_df.index + len(subset_metadata_df.index)

        analysis_metadata_df = pd.concat([subset_metadata_df,
                                          analysis_df], axis=0)
        analysis_metadata_df = analysis_metadata_df[analysis_df.columns]
        analysis_metadata_df.to_csv(pcl_out, 
                                    index=False, 
                                    header=header, 
                                    sep='\t',
                                    na_rep=na_rep)

    output_folder = os.path.dirname(analysis_files[0])
    pcl_files = bb_utils.name_files(analysis_files, 
                                    output_folder, 
                                    extension="pcl.tsv")

    # Because of how YAML inherits lists we'll need to see if we can't 
    # flatten this list out. 
    target_cols = funcy.flatten(target_cols)

    workflow.add_task_group(_workflow_add_metadata_to_tsv,
                            depends=analysis_files,
                            targets=pcl_files,
                            time="1*60 if ( file_size('depends[0]]') < 1 else 2*60",
                            mem="4*1024 if ( file_size('depends[0]]') < 1 else 3*12*1024" ,
                            cores=1,
                            name="Generate analysis PCL output file")

    return pcl_files