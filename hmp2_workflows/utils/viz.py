# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils.viz
~~~~~~~~~~~~~~~~~~~~~~~~~

A collection of utility functions needed for HMP2 AnADAMA2 visualization 
workflows.

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

import json

import pandas as pd

from operator import itemgetter


def convert_table_to_datatables_json(table_file, output_dir):
    """Converts a tab-delimited text file to a JSON file that can be read 
    by the jquery DataTables library.

    Args:
        table_file (string): Path to the table file to be converted.
        output_dir (string): Path to the output directory to write the 
            converted JSON file.

    Requires:
        None

    Returns:
        string: Path to the converted JSON file.

    Example:
        from hmp2_workflows.utils import viz

        species_counts_tbl = "/tmp/species.tsv"

        viz.table_to_datetables_json(species_counts_tbl, "/tmp")
    """
    table_basename = os.path.splitext(os.path.basename(table_file))[0]
    output_json_file = os.path.join(output_dir, '%s.json' % table_basename)

    table_df = pd.read_table(table_file)
    table_df.columns = table_df.columns.str.replace('# ', '')
    table_json = json.loads(table_df.to_json(orient='records'))

    with open(output_json_file, 'w') as json_out:
        json.dump(table_json, json_out)
    
    return output_json_file


def _generate_group_barplot_json(table_file, sort_on, transpose=False, rescale_vals=True):
    """Generates a JSON string that can be consumed by the Plotly JS library
    to produce a grouped barplot.

    Args:
        table_file (string): Path to the tab-delimited table file to be 
            parsed and converted to a JSON object.
        sort_on (string): Key to sort plot data on.
        transpose (boolean): Whether or not to transpose the dataframe
            once loaded. Depending on the table passed in this may need 
            to be done.       
        rescale_vals (boolean): Rescale values from 0-1 to 0-100 if True.

    Requires: 
        None

    Returns:
        string: A JSON representation of the table data.        
    """
    traces = []

    table_df = pd.read_table(table_file, index_col=0)
    
    if rescale_vals:
        table_df = table_df * 100

    table_df = table_df.sort_values(by=sort_on) if sort_on else table_df
    table_df = table_df.T if transpose else table_df
    
    table_split = json.loads(table_df.to_json(orient='split'))

    for (idx, col) in enumerate(table_split.get('columns')):
        traces.append({'name': col, 
                       'type': 'bar',
                        'x': table_split.get('index'), 
                        'y': map(itemgetter(idx), table_split.get('data'))})

    return traces


def convert_table_to_plotly_barplot_json(table_file, output_dir, sort_on=None, 
                                         hide_x_axis=False, plot_type='group', 
                                         relab_100_scale=True,
                                         xaxis_font_size='10', legend_order='reverse'):
    """Converts a tab-delimited text file to a JSON file that can be read 
    by the Plotly js library to generate dynamic charts.

    Args:
        table_file (string): Path to the table file to be converted.
        output_dir (string): Path to the output directory to write the 
            converted JSON file.
        sort_on (string): Key to sort plot data on [Default: None]
        hide_x_axis (boolean): If True hide x-axis and display the number of 
            x elements [Default: False]
        plot_type (string): The type of barplot to produce. [group, stacked]
        relab_100_scale (boolean): If True rescale relab'd values from 0-1 to
            0-100 [Default: True]
        xaxis_font_size (string): Size of x-axis label font [Default: 10]
        legend_order (string): Ordering of plot legend [Default: reverse]

    Requires:
        None

    Returns:
        string: Path to the converted JSON file.

    Example:
        from hmp2_workflows.utils import viz

        species_counts_tbl = "/tmp/species.tsv"

        viz.table_to_plotly_json(species_counts_tbl, "/tmp", type="group")
    """
    plot_json = {}

    table_basename = os.path.splitext(os.path.basename(table_file))[0]
    output_json_file = os.path.join(output_dir, '%s_plotly.json' % table_basename)
    output_json_file = output_json_file.replace('table', 'plot')

    transpose = False if plot_type == "group" else True
    plot_traces = _generate_group_barplot_json(table_file, sort_on, transpose, rescale_vals=True)

    plot_json['data'] = plot_traces
    
    plot_json['layout'] = {'barmode': plot_type}
    plot_json['layout']['xaxis'] = {}
    plot_json['layout']['yaxis'] = {'title': 'Number of Reads'}

    plot_json['layout']['legend'] = {'traceorder': legend_order}

    if hide_x_axis:
        plot_json['layout']['xaxis']['showticklabels'] = False

        ## If we are doing this we're going to want to add a label indicating 
        ## how many samples are in this dataset.
        plot_json['layout']['xaxis']['title'] = "Number of Samples: %s" % len(plot_traces[0]['x'])
    else:
        plot_json['layout']['xaxis']['title'] = "Samples"
        plot_json['layout']['xaxis']['tickfont'] = {'size': xaxis_font_size}

    with open(output_json_file, 'w') as json_out:
        json.dump(plot_json, json_out)

    return output_json_file


def generate_viral_read_proportions_file(viral_read_counts, out_dir):
    """Takes viral read counts as generated by the Baylor CMMR viromics analysis 
    pipeline and generate the proportion of reads that are viral out of the raw 
    and filtered readsets.

    Args:
        viral_read_counts (string): Path to the file containing read count data 
            as produced by viromics analysis pipeline.
        out_dir (string): Path to write viral read proportions file.

    Requires:
        None

    Returns:
        string: A path to a file containing viral read proportions.                    

    Example:
        from hmp2_workflows.utils import viz

        viral_read_counts = "/tmp/viral_read_counts.tsv"

        viz.table_to_plotly_json(viral_read_counts, "/tmp/output")
    """
    out_file = os.path.join(out_dir, "viral_read_proportions.tsv")
    viral_read_counts_df = pd.read_table(viral_read_counts)

    viral_read_counts_df['raw_proportion'] = viral_read_counts_df['Viral Reads'] / viral_read_counts_df['Raw Reads']
    viral_read_counts_df['trim_proportion'] = viral_read_counts_df['Viral Reads'] / viral_read_counts_df['Trimmed Reads']

    viral_read_counts_df.drop('Viral Reads', axis=1, inplace=True)
    viral_read_counts_df.drop('Trimmed Reads', axis=1, inplace=True)
    viral_read_counts_df.drop('Raw Reads', axis=1, inplace=True)

    viral_read_counts_df.to_csv(out_file, sep='\t', index=False)

    return out_file


def filter_taxonomic_profiles(tax_profile, metadata_file, output_folder, min_reads_threshold=1000000):
    """Filters taxonomic profiles based on a number of reads threshhold as well 
    as screening out any samples whos abundance can be defined by just one bug.

    Args:
        tax_profile (string): Path to taxonomic profile file.
        metadata_file (string): Path to project metadata.
        output_folder (string): Path to write filtered taxonomic profile file.
        min_reads_threshold (int): The minimum number of reads needed for a sample to be 
            included in downstream analysis (DEFAULT: 1000000)

    Requires:
        None

    Returns:
        pandas.DataFrame: Project metadata filtered to only samples passing criteria.
        string: Path to filtered taxonomic profile file.

    Example:
        import pandas as pd

        from hmp2_workflows.utils import viz

        taxonomic_profile = "/tmp/hmp2_tax_profiles.tsv"
        output_folder = "/home/carze/"
        metadata_file = "/home/carze/hmp2_metadata.csv"

        viz.filter_taxonomic_profiles(taxonomic_profile, metadata, output_folder)
    """
    filtered_tax_file = os.path.join(output_folder, "filtered_taxonomic_profiles.tsv")
    taxonomy_df = pd.read_table(tax_profile)
    metadata_df = pd.read_csv(metadata_file)

    ## If a metadata file is provided we should have access to number of reads here so we can 
    ## execute the same filtering steps from the manuscript
    valid_samples = metadata_df[metadata_df['filtered_reads'] >= 1000000]['External ID'].values.tolist()

    taxonomy_species_df = taxonomy_df[(taxonomy_df['#SampleID'].str.contains('s__')) &
                                      -(taxonomy_df['#SampleID'].str.contains('t__'))]
    weird_samples_series = (taxonomy_species_df == 100).any()
    weird_samples = weird_samples_series[weird_samples_series == True].index.tolist()
    weird_samples = [s.replace('_taxonomic_profile', '') for s in weird_samples]
    valid_samples = set(valid_samples) - set(weird_samples)

    taxonomy_df = taxonomy_df.filter(['#SampleID'] + 
                                     ['%s_taxonomic_profile' % s for s in valid_samples])
    taxonomy_df.to_csv(filtered_tax_file, sep='\t', index=False)

    metadata_df = metadata_df[metadata_df['External ID'].isin(valid_samples)]

    return (metadata_df, filtered_tax_file)