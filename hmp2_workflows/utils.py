# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils
~~~~~~~~~~~~~~~~~~~~

This module provides utility functions that are used within the AnaDaMa2 
HMP2 workflows.
"""

import os
import subprocess
import tempfile

from itertools import islice

import yaml


def parse_cfg_file(config_file, section=None):
    """Parses the provided YAML config file. If a specific section is 
    provided to parse only this section is returned.

    Args:
        config_file (string): Path to YAML config file 
        section (string): Specific section in the supplied config file 
                          to return.
    Requires:
        None
    Returns: 
        dict: A dictionary containing all configuration parameters found
              in the supplied config file.
    """
    stream = file(config_file, 'r')
    config = yaml.load(stream)
    
    if not section in config:
        raise KeyError('Section not found in config file', section)

    return config.get(section) if section else config


def parse_checksums_file(checksums_file):
    """Parses a file containing MD5 checksums in the following format:

        f36daa1eaade89762474d956833cffd0 *160928_SM-7BCZ1_400.raw
        0e59da3ff5b9af1e35e8013ef938d60f *160928_SM-7BCZJ_397.raw
        84e16a3e70d3b425df494fe38dc54765 *160928_SM-7GYK6_434.raw
    
    And returns a dictionary keyed on filename with MD5 checksums as values.

    Args:
        checksums_file (string): Path to file containing MD5 checksums

    Requires:
        None

    Returns:
        dict: A dictionary with filenames as keys and MD5 checksums as values.
    
    Example:
        md5_file = "/tmp/file_checksums.txt"
        md5_map = parse_checksums_file(md5_file)
    """        
    md5_map = {}

    with open(checksums_file) as md5_fh:
        for line in md5_fh:
            (checksum, filename) = md5_fh.split('\t')        
            
            filename = os.path.basename(filename)
            checksum_str = "%s  %s" 

            if "*" in filename:
                filename = filename.replace('*', '')
                checksum_str = "%s *%s"

            # Quick check to make sure the checksum is 32 characters long   
            # since MD5 checksums should always be 32 characters long.
            if len(checksum) != 32:
                raise ValueError('MD5 Checksum value is not valid', checksum)

            md5_map[filename] = checksum_str % (checksum, filename)

    return md5_map            
