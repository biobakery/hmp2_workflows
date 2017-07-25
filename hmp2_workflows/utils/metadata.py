# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils.metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A collection of functions related around metadata processing and generation 
related to the HMP2 AnADAMA2 workflows.

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


def get_collection_dates(broad_sample_df):
    """Retrieves all collection dates for a given subject and returns a 
    dicitonary key'd on subject.

    Args:
        metadata_df (pandas.DataFrame): Broad sample tracking status 

    Requires:
        None

    Returns:
        dict: A dictionary key'd on subject ID containing rows of metadata
            containing collection dates per visit.
    """
    def _add_previous_collection_date(dataframe):
        """For each row in a dataframe takes the previous date of reception
        and adds it as a new column to aid in computation of the week_num and
        interval_days columns in our final metadata table.
        """
        dataframe['prev_coll_date'] = dataframe['Actual Date of Receipt'].shift()
        return dataframe

    collection_dict = dict((subj, df) for (subj, df) in 
                           broad_sample_df.ix[:,'Subject':'Actual Date of Receipt']
                           .sort_values(by='Actual Date of Receipt').groupby('Subject'))
    collection_dict = dict((subj, _add_previous_collection_date(df))
                            for (subj, df) in collection_dict.iteritems())

    return collection_dict