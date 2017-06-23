# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.dcc
~~~~~~~~~~~~~~~~~~~~~~~~

This module contains tasks related to uploading metadata and data files 
to the HMP2 Data Coordination Center (DCC) using the cutlass API.

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


def upload_data_files(workflow, metadata, dcc_coll):
    """Transfers the provided iHMP OSDF object to the DCC using the cutlass
    API and aspera. This task is meant to upload in parallel to the DCC to
    account for the large amount of files that are present in the IBDMDB 
    datasets.

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 workflow object.
        metadata (pandas.DataFrame): Dataframe containing accompanying metadata
            for all files being submitted to the DCC
        dcc_coll (list): A list of lists containing all the DCC objects 
            that have been submitted and the sequence file (and accompanying 
            cutlass object) that will be submitted.

    Requires:
        None

    Returns:
        list: A list of the final sequence files submitted to the DCC.
    """
    uploaded_files = []

    #def _dcc_upload(task):
    def _dcc_upload(dcc_objects):
        """Invokes upload of sequencing product(s) to the DCC
        making using of Cutlass' aspera transfer functionality.

        Args: 
            task (anadama2.Task): AnADAMA2 Task object.

        Requires: 
            None

        Returns:
            None
        """
        seq = dcc_objects[-1]
        
        if seq:
            success = seq.save()
            if not success:
                raise ValueError('Saving sequence to DCC failed: %s' % seq.sample_name)
            else:
                uploaded_files.append(seq.sample_name)
    
    for dcc_objects in dcc_coll:
        if not dcc_objects[-1]:
            continue
     
        _dcc_upload(dcc_objects)
        # workflow.add_task(_dcc_upload,
        #                  depends = seq_files,
        #                  name = "DCC Upload %s" % dcc_objects[-1].sample_name,
        #                  time = 2*60,
        #                  mem = 1024,
        #                  cores = 1)

    return uploaded_files
