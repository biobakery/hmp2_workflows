# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains functions required to update, validate, refresh and 
disseminate HMP2 metadata  via the IBDMDB website and the DCC.
"""

import os

import cutplace
import cutplace.errors

from biobackery_workflow import utils as b_utils
from hmp2_workflows import utils as h_utils


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

    workflow.add_task( 
   
