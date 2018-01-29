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


def _generate_group_barplot_json(table_file, transpose=False):
    """Generates a JSON string that can be consumed by the Plotly JS library
    to produce a grouped barplot.

    Args:
        table_file (string): Path to the tab-delimited table file to be 
            parsed and converted to a JSON object.
        transpose (boolean): Whether or not to transpose the dataframe
            once loaded. Depending on the table passed in this may need 
            to be done.            

    Requires: 
        None

    Returns:
        string: A JSON representation of the table data.        
    """
    traces = []

    table_df = pd.read_table(table_file, index_col=0)
    table_df = table_df.T if transpose else table_df
    
    table_split = json.loads(table_df.to_json(orient='split'))

    for (idx, col) in enumerate(table_split.get('columns')):
        traces.append({'name': col, 
                       'type': 'bar',
                        'x': table_split.get('index'), 
                        'y': map(itemgetter(idx), table_split.get('data'))})

    return traces


def convert_table_to_plotly_barplot_json(table_file, output_dir, plot_type='group'):
    """Converts a tab-delimited text file to a JSON file that can be read 
    by the Plotly js library to generate dynamic charts.

    Args:
        table_file (string): Path to the table file to be converted.
        output_dir (string): Path to the output directory to write the 
            converted JSON file.
        plot_type (string): The type of barplot to produce. [group, stacked]            

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
    output_json_file = os.path.join(output_dir, '%s.json' % table_basename)
    output_json_file = output_json_file.replace('table', 'plot')

    transpose = False if plot_type == "group" else True
    plot_traces = _generate_group_barplot_json(table_file, transpose)

    plot_json['data'] = plot_traces
    plot_json['layout'] = {'barmode': plot_type}

    with open(output_json_file, 'w') as json_out:
        json.dump(plot_json, json_out)

    return output_json_file