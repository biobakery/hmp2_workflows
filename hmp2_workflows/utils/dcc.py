# -*- coding: utf-8 -*-

"""
hmp2_workflows.utils.dcc
~~~~~~~~~~~~~~~~~~~~~~~~

A collection of functions that are utilized when uploading data files and 
metadata to the Data Coordination Center (DCC).

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

import hashlib
import itertools
import operator
import os
import tempfile

import cutlass


def group_osdf_objects(osdf_collection, group_by_field):
    """Retrieves all subjects from the provided OSDF Study object. Subjects 
    are returned keyed on the rand_subject_id field which is mapable back to 
    a metadata table.

    Args:
        osdf_collection (cutlass.*): A collection of OSDF objects (Study,
            Subject, Sample etc.)
        group_by_field (string): The field (present in the OSDF objects)
            that the collection will be grouped by.

    Requires:
        None

    Returns:
        dictionary: A dictionary containing all OSDF Subject keyed on 
            the provided group by field.
    """
    return itertools.groupby(operator.attrgetter(group_by_field), 
                             osdf_collection)


def get_fields_to_update(new_metadata, osdf_object):
    """Compares existing values for metadata in the provided iHMP OSDF 
    object (Study, Subject, Sample, etc) to the provided metadata dictionary
    passed in. Any fields where values differ are returned to be updated.

    Args:
        new_metadata (dict): A dictionary containing proposed new metadata
        osdf_object (cutlass.*): An OSDF object that is to be updated using 
            the provided metadata.

    Requires:
        None

    Returns:
        list: A list of fields which will be updated.
    """
    required_fields = new_metadata.keys()

    ## Nested paramters in the tags and mixs fields need to be handled
    ## separately from our other fields. In this case the decision to replace
    ## the contents on a whole is made instead of checking each value and 
    ## updating piece meal
    if 'tags' in required_fields:
        required_fields.remove('tags')

        tags_new = sorted(new_metadata.get('tags'))
        tags_curr = sorted(osdf_object.tags)
        updated_fields.append('tags') if tags_new != tags_curr

    if 'mixs' in required_fields:
        required_fields.remove('mixs')

        mixs_new = sorted(new_metadata.get('mixs').values())
        mixs_curr = sorted(osdf_object.mixs.values())
        updated_fields.append('mixs') if mixs_new != mixs_curr

    updated_fields.extend([key for key in new_metadata.keys() 
                           if new_metadata.get(key) != osdf_object.get(key)])
    
    return updated_fields


def is_metadata_updates(cutlass_obj, metadata_new):
    """Given an existing cutlass object, attempt to figure out if  new
    metadata being passed contains any new information that constitutes 
    an update to the object.

    Args:
        dcc_object (cutlass.*): A cutlass object of type Study, Subject, 
            Sample, Visit, or file type (Proteome, WgsDnaPrep etc.).
        metadata_new (pyyaml.Loader): Metadata to be validated (against 
            objects metadata).

    Requires:
        None

    Returns:
        boolean: True or False depending on whether or not the metadata 
            provided constitutes an update to the cutlass object.
    """
    fields_to_check = metadata_new.keys()

    ## Create 
    req_metadata_new = map(lambda key: metadata_new.get(key, None),
                           fields_to_check)
    req_metadata_obj = map(lambda key: getattr(key, cutlass_obj),
                           fields_to_check)

    ## Because of mixs objects we need to also append metadata there 
    req_metadata_obj.extend(map(lambda key: getattr(key, cutlass_obj.mixs),
                                fields_to_check))
    
    req_metadata_new_md5 = hashlib.md5(repr(req_metadata_new)).hexdigest()
    req_metadata_obj_md5 = hashlib.md5(repr(req_metadata_obj)).hexdigest()

    return req_metadata_new_md5 == req_metadata_obj_md5


def get_project(conf, session):
    """Retrieves iHMP OSDF project provided a project ID to search for.
    
    Args:
        conf (pyyaml.Loader): A python object representation of our YAML 
            configuration object containing parameters and metadata required
            for the DCC upload process.
        session (cutlass.iHMPSession): Session object that represents a 
            a connection to the iHMP OSDF instance.

    Requires:
        None

    Returns: 
        cutlass.Project: An iHMP OSDF Project object containing 
            metadata for the specific project.

    Example:
        import cutlass
        import yaml

        from hmp2_workflows.utils import dcc

        session = cutlass.iHMPSession('user', 'pass')
        config = yaml.load('/tmp/foo/bar.conf')

        project = dcc.create_or_update_project(config, session)
    """
    namespace = conf.get('namespace')
    project_id = conf.get('project_id')

    osdf = session.get_osdf()
    query = '{ "query": { "match": { "meta.name": "%s" } } }' % project_id
    query_res = osdf.query_all_pages('ihmp', query)        
 
    if query_res.get('search_result_total') != 1:
        raise ValueError('Could not find existing project: %s' % project_id)

    dcc_project_id = query_res.get('id')

    return cutlass.Project.load(dcc_project_id)


def get_or_update_study(conf, session, project_id):
    """Retrieves an iHMP OSDF study using the provided study metadata. If any
    study metadata has been updated since the prior submission the updated 
    metadata will submitted.

    Args:
        conf (pyyaml.Loader): A python object representation of a YAML 
            configuration file. Contains parameters and metadata required
            for the DCC upload process.
        session (cutlass.iHMPSession): Session object that represents a 
            connection to the iHMP OSDF instance.
        project_id (string): iHMP OSDF Project ID acting as a parent to 
            the existing or to be created Study.

    Requires:
        None

    Returns:
        cutlass.Study: An iHMP OSDF Study object in a python representation.

    Example:
        import cutlass
        import yaml

        from hmp_workflows.utils import dcc
        from hmp_workflows.utils.misc import parse_cfg_file

        session - cutlass.iHMPSession('user', 'pass')
        config = parse_cfg_file('/tmp/foo/bar.conf')

        study = dcc.create_or_update_study(config, session)
    """
    namespace = conf.get('namespace')
    study_metadata = conf.get('study')
    study_id = study_metadata.get('name')

    osdf = session.get_osdf()
    query = '{ "query": { "match": { "meta.name": "%s" } } }' % study_id
    query_res = osdf.query_all_pages(namespace, query)

    if query_res.get('search_result_total') != 1:
        raise ValueError('Could not find existing study: %s' % study_id)

    study_id = query_res.get('id')
    study = cutlass.Study.load(dcc_study_id)

    fields_to_update = get_fields_to_update(study_metadata, study)
    map(lambda key: setattr(study, key, study_metadata.get(key)),
        fields_to_update)
        
    study.links['part_of'] = project_id if fields_to_update
        
    if study.is_valid():
        success = study.save()
        raise ValueError('Saving study %s failed.' % study.name) if not success
    else:
        raise ValueError('Study validation failed: %s' % study.validate())

    return study


def create_or_update_subject(subjects, metadata_subject_id, study_id, 
                             metadata):
    """Creates an iHMP OSDF Subject object if it does not exist or updates 
    an existing Subject object with new metadata when present.

    Args:
        subjects (list): A list of OSDF Subject objects
        metadata_subject_id (string): The subject ID (pulled from the metadata
            table) to check if an OSDF Subject exists for.
        study_id (string): OSDF Study ID that supplied subject should be 
            linked too.
        metadata (panda.Series): All metadata for one subject/sample combo

    Requires:
        None

    Returns:
        cutlass.Subject: The created or updated OSDF Subject object.
    """
    subject = subjects.get(metadata_subject_id)

    if not subject: 
        subject = Cutlass.Subject()
    
    req_metadata = {}
    req_metadata['rand_subject_id'] = metadata_subject_id
    req_metadata['gender'] = metadata['Sex'].lower()
    req_metadata['race'] = metadata['race']
    req_metadata['tags']['diagnosis'] = metadata['Diagnosis'].lower().replace(' ', '_')
    req_metadata['tags']['age_at_dx'] = metadata['age_at_dx']
    req_metadata['tags']['highest_education'] = metadata['Education Level']

    fields_to_upadte = get_fields_to_update(req_metadata, subject)
    map(lambda key: setattr(subject, key, req_metadata.get(key)),
        fields_to_update)

    subject.links['participates_in'] = [study_id]

    if subject.is_valid():
        success = subject.save()
        raise ValueError('Saving subject %s failed.' % subject.rand_subject_id)
    else:
        raise ValueError('Subject validation failed: %s' % subject.validate())
    
    return subject


def create_or_update_visit(visits, subject, metadata):
    """Creates an iHMP OSDF Visit object if it does not exist or updates an
    already existing Visit object with the provided metadata.

    Args:
        visits (list): A list of exisitng OSDF Visit objects.
        subject (cutlass.Subject): The OSDF Subject that the 
    """
    pass
