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


import importlib
import itertools
import operator
import os
import tempfile

import cutlass
import numpy as np
import pandas as pd

from cutlass.mixs import MIXS
from cutlass.mims import MIMS


def _convert(value, type_):
    """Casts the provided value to the specified type.

    Args:
        value (string): The value to be casted.
        type (string): The type to cast the supplied value too.

    Requires:
        None

    Returns:
        <TYPE>: Value casted to the supplied type.
    """
    module = importlib.import_module('__builtin__')
    cls_ = getattr(module, type_)

    return cls_


def required_mixs_dict():
    """Generates the default sample mixs dictionary needed when submitting a
    Cutlass Sample object.

    Args:
        None

    Requires:
        None

    Returns:
        dictionary: A dictionary containing the fields required for a 
            Cutlass.Sample mixs entry.
    """
    return dict([ (k, v()) for k, v 
                   in MIXS._fields.iteritems() 
                   if k in MIXS.required_fields()])


def required_mims_dict():
    """Generates the default sample mims dictionary needed when submitting a
    Cutlass Sample object.

    Args:
        None

    Requires:
        None

    Returns:
        dictionary: A dictionary containing the fields required for a 
            Cutlass.Sample mixs entry.
    """
    return dict([ (k, v()) for k, v 
                   in MIMS._fields.iteritems() 
                   if k in MIMS.required_fields()])    


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
    grouped_objs = itertools.groupby(osdf_collection, 
                                     operator.attrgetter(group_by_field))

    grouped_dict = {}
    for (key, group) in grouped_objs:
        if hasattr(key, "__iter__"):
            key = key[0]

        grouped_dict[key] = list(group)            

    return grouped_dict


def create_seq_fname_map(data_type, data_files, tags=[]):
    """Creates a mapping of the sequences files to sample identifiers 
    derived from the data file name.

    Args:
        data_type (string): The file type for the files provided in the 
            data_files parameter.
        data_files (list): Filenames to map back to sample identifiers. 
            Depending on the source of the file (Broad, PNNL etc.) naming 
            schemes for files are different and not to be handled on a case
            by case basis.
        tag (list): A list of tags that should be removed from sample 
 
    Requires:
        None

    Returns:
        dict: A dictionary key'd on sample identifier (derived from 
            file name) and values consisting of paths to the corresponding 
            data file.
    """
    sample_id_map = {}

    ## We're going to make a bit of a dangerous assumption here. Currently we
    ## are seeing samples being named two ways:
    ##
    ##  1.) Broad: MSM5LLHX.bam
    ##  2.) PNNL: 160513-SM-AHYMJ-14.raw
    ## 
    ## In both cases here we the sample name is prefixed by 'SM' and may or 
    ## may not have hyphens or underscores separating it. The Broad samples 
    ## are also prefixed by the originating center represented by one 
    ## one character in front of our 'SM' prefix.
    for data_file in data_files:
        ## With proteomics datasets we sometimes get 'pool' files that 
        ## can be ignored for the time being.
        if 'pool' in data_file:
            continue

        (file_name, ext) = os.path.basename(data_file).split(os.extsep, 1)

        if data_type == 'proteomics':
            sample_id = "%s-%s" % ('SM', 
                                   file_name.replace('_', '-').split('-')[2])
        else:
            sample_id = file_name
            
            ## TODO: This shouldn't be hardcoded
            sample_id = sample_id.replace('_taxonomic_profile', '')
    
        for tag in tags:
            sample_id = sample_id.replace(tag, '')

        sample_id_map[sample_id] = data_file

    return sample_id_map
    

def map_sample_id_to_file(row, id_col, fname_map, is_proteomics, out_col='seq_file'):
    """Given a a row from a pandas DataFrame, map the sample identifier 
    in the metadata to a file map and return the corresponding file for 
    the given row.

    Args:
        row (pandas.Series): A row from a pandas Dataframe containing HMP2
            project metadata.
        id_col (string): Column to pull the metadata sample identifier 
            from.
        seq_fname_map (dict): The dictionary map used to map metadata sample 
            identifiers to their corresponding sequence files.
        is_proteomics (boolean): Whether or not we are dealing with proteomics
            data and if we need to account for multiple sequence files 
            attached to the same sample

    Requires:
        None

    Returns:
        string: The path to the corresponding sequence file assocaited
            with a given metadata row.
    """
    sample_id = row.get(id_col)
    sample_file = fname_map.get(sample_id)
    pdo_number = row.get('PDO Number')

    if is_proteomics and str(pdo_number) in sample_file:
        row[out_col] = fname_map.get(sample_id)
    else:
        row[out_col] = sample_file

    return row


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
    updated_fields = []
    required_fields = new_metadata.keys()

    ## Nested paramters in the tags and mixs fields need to be handled
    ## separately from our other fields. In this case the decision to replace
    ## the contents on a whole is made instead of checking each value and 
    ## updating piece meal
    if 'tags' in required_fields:
        required_fields.remove('tags')

        tags_new = sorted(new_metadata.get('tags'))
        tags_curr = sorted(osdf_object.tags)
        if tags_new != tags_curr:
            updated_fields.append('tags')

    if 'mixs' in required_fields:
        if not osdf_object.mixs:
            updated_fields.append('mixs')
        else:
            required_fields.remove('mixs')
            mixs_new = sorted(new_metadata.get('mixs').values())
            mixs_curr = sorted(osdf_object.mixs.values())

            if mixs_new != mixs_curr:
                updated_fields.append('mixs')

    if 'checksums' in required_fields:
        local_files = [(k,v) for (k,v) in new_metadata.iteritems()
                       if 'local_' in k and 'tmp' not in v]

        ## This is super ugly but going to operate under the assumption that 
        ## we only have on raw file per DCC object.
        if len(local_files) > 1:
            raise ValueError('More than one file associated with this object: %s' 
                             % ",".join([lfile for (k, lfile) in local_files]))

        (local_file_field, local_file) = local_files[0]
        required_fields.remove('checksums')
        required_fields.remove(local_file_field)

        ## Most of the objects we deal with will be single files so we only
        ## need to check a single checksum
        osdf_checksum = osdf_object.checksums.get('md5') if osdf_object.checksums else None
        new_checksum = new_metadata['checksums']['md5']

        if new_checksum != osdf_checksum:
            updated_fields.extend(['checksums', local_file_field])

    updated_fields.extend([key for key in required_fields
                           if new_metadata.get(key) != getattr(osdf_object, key)])

    return np.unique(updated_fields).tolist()


def _get_host_assay_prep_abund_matrices(session, prep_id):
    """Returns an iterator of all AbundanceMatrix nodes connnected to the
    provided HostAssayPrep node.

    Args:
        session (cutlass.Session): The current OSDF session object.
        prep_id (string): The prep ID to retrieve any associated AbundanceMatrix
            nodes from.

    Requires:
        None

    Returns:
        Iterator: An iterator containing all children connected to the 
            supplied OSDF object.                    
    """
    query = session.get_odsf().oql_query
    linkage_query = ('"abundance_matrix"[node_type] && '
                     '"{}"[linkage.computed_from]'.format(prep_id))

    from cutlass.AbundanceMatrix import AbundanceMatrix
    from cutlass.HostAssayPrep import HostAssayPrep

    for page_no in itertools.count(1):
        res = query(WgsRawSeqSet.namespace, linkage_query,
                    page=page_no)
        res_count = res['result_count']

        for doc in res['results']:
            yield AbundanceMatrix.load_abundance_matrix(doc)

        res_count -= len(res['results'])

        if res_count < 1:
            break


def _get_wgs_raw_seq_set_abund_matrices(session, seq_set_id):
    """
    Returns an iterator of all AbundanceMatrix nodes connected to 
    the WgsRawSeqSet ID provided.

    Args:
        osdf (cutlass.Session): The current OSDF session
        seq_set_id (string): The sequence set ID from which to grab 
            any associated abundance matrices.

    Requires:
        None

    Returns:
        iterator: An iterator containing all found abundance matrices.                            

    """
    if seq_set_id:
        query = session.get_osdf().oql_query
        linkage_query = ('"abundance_matrix"[node_type] && ' 
                         '"{}"[linkage.computed_from]'.format(seq_set_id))

        from cutlass.AbundanceMatrix import AbundanceMatrix
        from cutlass.WgsRawSeqSet import WgsRawSeqSet            

        for page_no in itertools.count(1):
            res = query(WgsRawSeqSet.namespace, linkage_query,
                        page=page_no)
            res_count = res['result_count']

            for doc in res['results']:
                yield AbundanceMatrix.load_abundance_matrix(doc)

            res_count -= len(res['results'])

            if res_count < 1:
                break


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
    query_resp = osdf.query_all_pages(namespace, query)        
 
    if query_resp.get('search_result_total') != 1:
        raise ValueError('Could not find existing project: %s' % project_id)
    
    query_res = query_resp.get('results')[0]
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
    query = ('{ "query": { "match": { "meta.name": { "query": "%s", "operator"'
            ': "and" } } } }' % study_id)
    query_resp = osdf.query_all_pages(namespace, query)

    if query_resp.get('search_result_total') != 1:
        raise ValueError('Could not find existing study: %s' % study_id)

    query_res = query_resp.get('results')[0]
    study_id = query_res.get('id')
    study = cutlass.Study.load(study_id)

    fields_to_update = get_fields_to_update(study_metadata, study)
    map(lambda key: setattr(study, key, study_metadata.get(key)),
        fields_to_update)
        
    if fields_to_update: 
        study.links['part_of'] = [project_id]
        
        if study.is_valid():
            success = study.save()
            if not success: 
                raise ValueError('Saving study %s failed.' % study.name) 
        else:
            raise ValueError('Study validation failed: %s' % study.validate())

    return study


def _get_medications_list(metadata, med_cols):
    """Parses the provided metadata collection to retrieve a list of 
    medications a subject is currently taking.

    Args:
        metadata (pandas.DataFrame): A collection of technical and clinical
            metadata for a given subject.
        med_cols (list): A list of column names that indicate a type of 
            medication a subject may be currently taking.

    Requires:
        None

    Returns:
        string: A string containing a list of medications a subject is 
            currently taking.        
    """
    medications = None

    all_meds_df = metadata.filter(med_cols)
    curr_meds = all_meds_df.apply(lambda val: val == "Current")
    medications = ",".join([med_name for (med_name, val) in curr_meds.iteritems() if val])

    return medications


def _crud_subject_attribute(subject, metadata, study_id, conf):
    """Creates an iHMP OSDF SubjectAttribute object if it does not exist or
    updates an already existing object with the provided metadata.

    Args:
        subject (cutlass.Subject): The Subject object that this SubjectAttribute 
            object should be linked too.
        metadata (pandas.DataFrame): A collection of metadata that will be used 
            to instantiate or update the SubjectAttribute object
        study_id (string): The study ID this SubjectAttribte object is 
            associated with.
        conf (dict): Python representation of YAML configuration file that 
            contains 

    Requires:
        None

    Returns:
        cutlass.SubjectAttribute: The created or updated OSDF 
            SubjectAttribute object.
    """
    subject_attrs = subject.attributes()
    sa_col_map = conf.get('col_map')

    subject_attr = next(subject_attrs, None)
    if not subject_attr:
        subject_attr = cutlass.SubjectAttribute()
 
    req_metadata = {}
    req_metadata = dict((k, metadata.get(sa_col_map.get(k))) for k in 
                        sa_col_map.keys())
    req_metadata['rx'] = _get_medications_list(metadata, conf.get('meds_cols'))
    req_metadata['study'] = study_id
    req_metadata = dict((k,v) for k, v in req_metadata.iteritems() if v and not pd.isnull(v))


    fields_to_update = get_fields_to_update(req_metadata, subject_attr)
    map(lambda key: setattr(subject_attr, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        subject_attr.links['associated_with'] = [subject.id]

        if subject_attr.is_valid():
            success = subject_attr.save()
            if not success: 
                raise ValueError('Saving subject attribute for subject %s failed.', 
                                  subject.rand_subject_id)
        else:
            raise ValueError('Subject attribute validationg failed: %s' % 
                             subject_attr.validate())                   

    return subject_attr


def crud_subjects(subjects, study, baseline_metadata, conf):
    """Creates iHMP OSDF Subjects objects if they do not exist or 
    updates an existing Subject object with new metadata when present.

    Args:
        subjects (list): A list of OSDF Subject objects
        study_id (cutlass.Studdy): OSDF Study that supplied subject should be 
            linked too.
        baseline_metadata (panda.DataFrame): Collection of metadata containing 
            baseline metadata for all subjects involved in the project.
        conf (dict): A python dictionary representation of the YAML 
            configuration file containing metadata parameters needed for 
            the DCC upload.

    Requires:
        None

    Returns:
        list: A list of OSDF Subject objects key'd on subject ID.
    """
    for (idx, row) in baseline_metadata.iterrows():
        subject_id = str(row.get('ProjectSpecificID'))

        if subject_id in subjects:
            subject = subjects.get(subject_id)[0]
        else:
            subject = cutlass.Subject()

        race_map = conf.get('race_map')

        req_metadata = {}
        req_metadata['tags'] = []
        req_metadata['rand_subject_id'] = subject_id
        
        sex = row.get('sex').lower() if not pd.isnull(row.get('sex')) else None
        req_metadata['gender'] = sex

        race = row.get('race') if not pd.isnull(row.get('race')) else None
        req_metadata['race'] = race_map.get(race)

        fields_to_update = get_fields_to_update(req_metadata, subject)
        map(lambda key: setattr(subject, key, req_metadata.get(key)),
            fields_to_update)

        if fields_to_update:
            subject.links['participates_in'] = [study.id]

            if subject.is_valid():
                success = subject.save()
                if not success:
                    raise ValueError('Saving subject %s failed.' % 
                                    metadata_subject_id)
            
                subjects[subject_id] = subject
            else:
                raise ValueError('Subject validation failed: %s' % 
                                subject.validate())
    
        subject_attr = _crud_subject_attribute(subject, row, 
                                               conf.get('data_study'),
                                               conf.get('subject_attribute'))

    return subjects


def crud_visit(visits, visit_num, subject_id, metadata, conf):
    """Creates an iHMP OSDF Visit object if it does not exist or updates an
    already existing Visit object with the provided metadata.

    Args:
        visits (list): A list of existing OSDF Visit objects.
        visit_num (string): The visit number to search for or create a new 
            Visit object for.
        subject_id (string): Subject ID that this visit is/should be 
            associated with.
        metadata (pandas.Series): Metadata that should be associated with 
            this Visit.
        conf (dict): A python dictionary representation of the YAML 
            configuration file containing metadata parameters needed for 
            the DCC upload.

    Requires:
        None

    Returns:
        cutlass.Visit: The created or updated OSDF Visit object.
    """
    visit = visits.get(visit_num)
    
    ## TODO: Figure out what's going on with Visit Attributes
    if not visit:
        visit = cutlass.Visit()
    else:
        visit = visit[0]
        
    req_metadata = {}
    req_metadata['visit_number'] = visit_num
    req_metadata['visit_id'] = "%s_%s" % (metadata.get('ProjectSpecificID'), 
                                          visit_num)
    req_metadata['interval'] = (int(metadata['interval_days']) 
                                if not np.isnan(metadata['interval_days']) else 0)
    ## This is hard-coded to meet HIPAA compliance.
    req_metadata['date'] = "2000-01-01"     

    fields_to_update = get_fields_to_update(req_metadata, visit)
    map(lambda key: setattr(visit, key, req_metadata.get(key)), 
        fields_to_update)

    if fields_to_update:
        visit.links['by'] = [subject_id]

        if visit.is_valid():
            success = visit.save()
            if not success:
                raise ValueError('Saving visit %s failed.' % visit_num)

        else:
            raise ValueError('Visit validation failed: %s' % visit.validate())

    visit_attr = _crud_visit_attribute(visit, metadata, conf)

    return visit


def _get_visit_attr_metadata(metadata, visit_attr_conf, req_metadata):
    """Retrieves all visit attribute metadata from the provided 
    metadata pandas DataFrame using the mappings found in the 
    supplied dictionary.

    Args:
        metadata (pandas.DataFrame): Collection of metadata that will be used
            to populate VisitAttribute object.
        visit_attr_conf (dict): Dictionary containing pieces of metadata
            required by VisitAttribute object.    
        req_metadat (dict): Dictionary housing all required metadata.            

    Requires:
        None

    Returns:
        dict: Dictionary containing metadata required by VisitAttribute
            object.                            
    """
    visit_attr_cols = visit_attr_conf.get('col_map')

    col_dict = dict(zip(visit_attr_cols.keys(), [metadata.get(k) for k in visit_attr_cols.values()]))
 
     ## Remove any keys that have None as value
    req_metadata.update(dict((k,v) for k,v in col_dict.iteritems() if not pd.isnull(v)))

    ## Out of our metadata dict we will have a bunch of string values 
    ## that we will need to convert to the proper types for OSDF
    va_attrs = cutlass.VisitAttribute.__dict__.get('_VisitAttribute__dict')

    for (key, val) in req_metadata.iteritems():
        if key.startswith('disease'):
            continue

        if key not in va_attrs:
            raise ValueException('Invalid key for VisitAttribute objects:', key)
        
        cls = va_attrs.get(key)[0]
        req_metadata[key] = cls(val)

    return req_metadata


def _crud_visit_attribute(visit, metadata, conf):
    """Creates an iHMP OSDF VisitAttribute object if it does not exist or
    updates an already existing object with the provided metadata.

    Args:
        sample (cutlass.Visit): The Visit object that this VisitAttribute 
            object should be linked too.
        metadata (pandas.Series): A collection of metadata that will be used 
            to instantiate or update the Sample Attribute Object
        conf (dict): Python representation of YAML configuration file that 
            contains 

    Requires:
        None

    Returns:
        cutlass.SampleAttribute: The created or updated OSDF Sample Attribute
            object.
    """
    visit_attrs = visit.visit_attributes()

    ## We should only ever have one VisitAttr object associated with a Visit
    ## object so updating shoudl be a little aeasier.
    visit_attr = next(visit_attrs, None)
    if not visit_attr:
        visit_attr = cutlass.VisitAttribute()

    visit_attr_conf = conf.get('visit_attribute')
    disease_map = conf.get('disease_map')
    disease_desc_map = visit_attr_conf.get('desc_map')

    req_metadata = {}
    #req_metadata['study'] = conf.get('data_study')
    req_metadata['comment'] = ""

    diagnosis = metadata.get('diagnosis')
    if diagnosis != "nonIBD":
        req_metadata['disease_name'] = disease_map.get(diagnosis)
        req_metadata['disease_description'] = disease_desc_map.get(diagnosis)

        if diagnosis in ['UC', 'CD']:
            req_metadata['disease_study_status'] = "affected"
        elif diagnosis == "nonIBD":
            req_metadata['disease_study_status'] = "not_affected"
        else:
            req_metadata['disease_study_status'] = "unknown"

        xref_map = visit_attr_conf.get('xref_map').get(diagnosis)
        if xref_map:
            req_metadata['disease_ontology_id'] = xref_map.get('DO')
            req_metadata['disease_mesh_id'] = xref_map.get('MESH')
            req_metadata['disease_nci_id'] = xref_map.get('NCI')
            req_metadata['disease_umls_concept_id'] = xref_map.get('UML')

    ## None of our subjects have cancer
    req_metadata['cancer'] = "No"

    if metadata.get('IntervalName') == "Screening Colonoscopy":
        req_metadata['colonoscopy'] = True

    req_metadata = _get_visit_attr_metadata(metadata,
                                            visit_attr_conf, 
                                            req_metadata)

    req_metadata['tags'] = []

    if req_metadata.get('hbi_total'):
        req_metadata['hbi'] = True
    if req_metadata.get('sccai_total'):
        req_metadata['sccai'] = True

    fields_to_update = get_fields_to_update(req_metadata, visit_attr)
    map(lambda key: setattr(visit_attr, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        visit_attr.subtype = conf.get('data_study')
        visit_attr.links['associated_with'] = [visit.id]

        if visit_attr.is_valid():
            success = visit_attr.save()
            if not success: 
                raise ValueError('Saving visit attribute for visit %s failed.', 
                                  visit.id)
        else:
            raise ValueError('Visit attribute validation failed: %s' % 
                             visit_attr.validate())                   

    return visit_attr


def _get_biopsy_location(metadata, body_site_map):
    """Extracts biopsy location for a given sample by examining the 
    four biopsy columns available in the HMP2 metadata table. 

    Args:
        metadata (pandas.DataFrame): DataFrame containing sample metadata

    Requires:
        None

    Returns:
        string: The location of the biopsy                
    """
    biopsy_location = metadata.get('biopsy_location')
    
    if biopsy_location == "Other Inflamed":
        biopsy_location = metadata.get('Location of inflamed RNA sample:')
    elif biopsy_location == "Non-inflamed":
        biopsy_location = metadata.get('Location of non-inflamed DNA/RNA sample:')                

    return body_site_map.get(biopsy_location)


def _crud_sample_attribute(sample, metadata, conf):
    """Creates an iHMP OSDF Sample Attribute object if it does not exist or
    updates an already existing object with the provided metadata.

    Args:
        sample (cutlass.Sample): The Sample object that this Sample Attribute 
            object should be linked too.
        metadata (pandas.Series): A collection of metadata that will be used 
            to instantiate or update the Sample Attribute Object
        conf (dict): Python representation of YAML configuration file that 
            contains 

    Requires:
        None

    Returns:
        cutlass.SampleAttribute: The created or updated OSDF Sample Attribute
            object.
    """
    sample_attrs = sample.sampleAttributes()

    ## Sample Attributes are a bit tricky to handle in that they are a 
    ## colleciton associated with a sample and not key'd up on an specific
    ## field so its more difficult to keep track of when we are dealing with
    ## updating an existing SampleAttribute object or we need to create a new 
    ## one.
    ##
    ## For the purposes of the IBD project each sample will only have one 
    ## SampleAttribute object associated with it so this simplifies the process
    ## of checking whether or not we are updating the object or creating a new 
    ## object.
    sample_attr = next(sample_attrs, None)
    if not sample_attr:
        sample_attr = cutlass.SampleAttribute()
    
    req_metadata = {}
    req_metadata['study'] = conf.get('data_study')
    req_metadata['fecalcal'] = str(metadata['fecalcal'])
    
    fields_to_update = get_fields_to_update(req_metadata, sample_attr)
    map(lambda key: setattr(sample_attr, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        sample_attr.links['associated_with'] = [sample.id]

        if sample_attr.is_valid():
            success = sample_attr.save()
            if not success: 
                raise ValueError('Saving sample attribute for sample %s failed.', 
                                  sample.id)
        else:
            raise ValueError('Sample attribute validationg failed: %s' % 
                             sample_attr.validate())                   

    return sample_attr


def crud_sample(samples, sample_id, visit_id, conf, metadata):
    """Creates an iHMP OSDF Sample object if it doesn't exist or updates 
    an already exisiting Sample object with the provided metadata.

    Args:
        samples (list): A list of exisiting OSDF Sample objects.
        sample_id (string): The sample ID to search for or create a new 
            Sample object for.
        visit_id (string): Visit ID that this sample is/should be associated
            with.
        conf (dict): Python representation of YAML configuration file that 
            contains 
        metadata (pandas.Series): Metadata that is associated with this 
            Sample.

    Requires:
        None

    Returns:
        cutlass.Subject: The created or updated OSDF Subject object.
    """
    sample = samples.get(sample_id)
    sample_conf = conf.get('sample')

    body_site_map = sample_conf.get('body_site_map')
    fma_body_site_map = sample_conf.get('fma_body_map')

    if not sample:
        sample = cutlass.Sample()
    else:
        sample = sample[0]
    
    req_metadata = {}
    req_metadata['mixs'] = required_mixs_dict()
    req_metadata['name'] = sample_id

    if metadata['data_type'] == 'host_transcriptomics':
        req_metadata['body_site'] = _get_biopsy_location(metadata, body_site_map)
    else:        
        req_metadata['body_site'] = "stool"

    req_metadata['fma_body_site'] = fma_body_site_map.get(req_metadata['body_site'])
    req_metadata['mixs'].update(sample_conf.get('mixs'))

    fields_to_update = get_fields_to_update(req_metadata, sample)
    map(lambda key: setattr(sample, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        sample.links['collected_during'] = [visit_id]

        if sample.is_valid():
            success = sample.save()
            if not success:
                raise ValueError('Saving sample % failed.' % sample_id)
        else:
            raise ValueError('Sample validation failed; %s' % sample.validate())

    sample_attr = _crud_sample_attribute(sample, metadata, conf)
    
    return sample


def crud_wgs_dna_prep(sample, study_id, dtype_abbrev, conf, metadata):
    """Creates an iHMP OSDF WGS DNA Prep object if it doesn't exist or 
    updates an already existing Prep object with the provided metadata. 

    Args:
        sample (cutlass.Sample): The Sample object that the Prep should be 
            associated with.
        study_id (string): The study ID this microbiome assay prep is
            assocaited with.
        dtype_abbrev (string): Data type abbreviation to be used in constructing 
            unique prep ID.            
        conf (dict): Python dictionary representation of the project YAML
            configuration containing project metadata.
        metadata (pandas.Series): The metadata that is assocaited with this 
            Sample/Prep.

    Requires:
        None

    Returns:
        cutlass.WgsDNAPrep: The created or updated OSDF 
            MicrobiomeAssayPrep object.
    """
    prep_id = "%s_%s" % (metadata.get('External ID'), dtype_abbrev)

    wgs_dna_preps = group_osdf_objects(sample.wgsDnaPreps(),
                                       'prep_id')
    wgs_dna_prep = wgs_dna_preps.get(prep_id)

    if not wgs_dna_prep:
        wgs_dna_prep = cutlass.WgsDnaPrep()
    else:
        wgs_dna_prep = wgs_dna_prep[0]
    
    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}
    req_metadata.update(conf.get('assay'))
    req_metadata['mims'] = required_mims_dict()

    ## Fill in the remaining pieces of metadata needed from other sources
    req_metadata['prep_id'] = prep_id
    req_metadata['mims'].update(conf.get('mims'))
    req_metadata['mims']['collection_date'] = metadata['date_of_receipt']

    fields_to_update = get_fields_to_update(req_metadata, wgs_dna_prep)
    map(lambda key: setattr(wgs_dna_prep, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        wgs_dna_prep.links['prepared_from'] = [sample.id]

        if wgs_dna_prep.is_valid():
            success = wgs_dna_prep.save()
            if not success:
                raise ValueError('Saving WGS DNA prep %s failed.' % 
                                 req_metadata.get('prep_id'))
        else:
            raise ValueError('WGS DNA prep validation failed: %s' % 
                             wgs_dna_prep.validate())
    
    return wgs_dna_prep


def crud_host_assay_prep(sample, study_id, dtype_abbrev, conf, metadata):
    """Creates an iHMP OSDF Host Assay Prep object if it doesn't exist or 
    updates an already existing Prep object with the provided metadata. 

    Args:
        sample (cutlass.Sample): The Sample object that the Prep should be 
            associated with.
        study_id (string): The study ID this microbiome assay prep is
            assocaited with.
        dtype_abbrev (string): The abbreviation for the data type that 
            this prep represents (i.e. MGX, MTX, 16S)            
        conf (dict): Python dictionary representation of the project YAML
            configuration containing project metadata.
        metadata (pandas.Series): The metadata that is assocaited with this 
            Sample/Prep.

    Requires:
        None

    Returns:
        cutlass.HostAssayPrep: The created or updated OSDF HostAssayPrep 
            object.
    """
    prep_id = "%s_%s" % (metadata.get('External ID'), dtype_abbrev)

    host_assay_preps = group_osdf_objects(sample.hostAssayPreps(),
                                          'prep_id')
    host_assay_prep = host_assay_preps.get(prep_id)

    if not host_assay_prep:
        host_assay_prep = cutlass.HostAssayPrep()
    else:
        host_assay_prep = host_assay_prep[0]
    
    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}
    req_metadata.update(conf.get('assay'))

    ## Fill in the remaining pieces of metadata needed from other sources
    req_metadata['prep_id'] = prep_id
    req_metadata['comment'] = "IBDMDB"

    fields_to_update = get_fields_to_update(req_metadata, host_seq_prep)
    map(lambda key: setattr(host_seq_prep, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        host_assay_prep.links['prepared_from'] = [sample.id]

        if host_assay_prep.is_valid():
            success = host_assay_prep.save()
            if not success:
                raise ValueError('Saving host seq prep %s failed.' % 
                                 req_metadata.get('prep_id'))
        else:
            raise ValueError('Host seq prep validation failed: %s' % 
                             host_assay_prep.validate())
    
    return host_assay_prep


def crud_host_seq_prep(sample, study_id, conf, metadata):
    """Creates an iHMP OSDF Host Seq Prep object if it doesn't exist or 
    updates an already existing Prep object with the provided metadata. 

    Args:
        sample (cutlass.Sample): The Sample object that the Prep should be 
            associated with.
        study_id (string): The study ID this microbiome assay prep is
            assocaited with.
        conf (dict): Python dictionary representation of the project YAML
            configuration containing project metadata.
        metadata (pandas.Series): The metadata that is assocaited with this 
            Sample/Prep.

    Requires:
        None

    Returns:
        cutlass.HostSeqPrep: The created or updated OSDF 
            HostSeqPrep object.
    """
    host_seq_preps = group_osdf_objects(sample.hostSeqPreps(),
                                        'prep_id')
    host_seq_prep = host_seq_preps.get(metadata['External ID'])

    if not host_seq_prep:
        host_seq_prep = cutlass.HostSeqPrep()
    else:
        host_seq_prep = host_seq_prep[0]
    
    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}
    req_metadata.update(conf.get('assay'))

    ## Fill in the remaining pieces of metadata needed from other sources
    req_metadata['prep_id'] = metadata['External ID']
    req_metadata['comment'] = "IBDMDB"

    fields_to_update = get_fields_to_update(req_metadata, host_seq_prep)
    map(lambda key: setattr(host_seq_prep, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        host_seq_prep.links['prepared_from'] = [sample.id]

        if host_seq_prep.is_valid():
            success = host_seq_prep.save()
            if not success:
                raise ValueError('Saving host seq prep %s failed.' % 
                                 req_metadata.get('prep_id'))
        else:
            raise ValueError('Host seq prep validation failed: %s' % 
                             host_seq_prep.validate())
    
    return host_seq_prep


def crud_microbiome_prep(sample, study_id, conf, metadata):
    """Creates an iHMP OSDF Microbiome Assay Prep object if it doesn't exist or 
    updates an already existing Prep object with the provided metadata. 

    Args:
        sample (cutlass.Sample): The Sample object that the Prep should be 
            associated with.
        study_id (string): The study ID this microbiome assay prep is
            assocaited with.
        conf (dict): Python dictionary representation of the project YAML
            configuration containing project metadata.
        metadata (pandas.Series): The metadata that is assocaited with this 
            Sample/Prep.

    Requires:
        None

    Returns:
        cutlass.MicrobiomeAssayPrep: The created or updated OSDF 
            MicrobiomeAssayPrep object.
    """
    microbiome_preps = group_osdf_objects(sample.microbAssayPreps(),
                                          'prep_id')
    microbiome_prep = microbiome_preps.get(metadata['Project'])

    if not microbiome_prep:
        microbiome_prep = cutlass.MicrobiomeAssayPrep()
    else:
        microbiome_prep = microbiome_prep[0]
    
    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}
    req_metadata.update(conf.get('assay'))

    ## Fill in the remaining pieces of metadata needed from other sources
    req_metadata['prep_id'] = metadata['Project']
    req_metadata['pride_id'] = 'PRIDE ID'
    req_metadata['sample_name'] = sample.name
    req_metadata['study'] = study_id
    req_metadata['comment'] = "IBDMDB"

    fields_to_update = get_fields_to_update(req_metadata, microbiome_prep)
    map(lambda key: setattr(microbiome_prep, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        microbiome_prep.links['prepared_from'] = [sample.id]

        if microbiome_prep.is_valid():
            success = microbiome_prep.save()
            if not success:
                raise ValueError('Saving microbiome prep %s failed.' % 
                                 req_metadata.get('prep_id'))
        else:
            raise ValueError('Microbiome assay prep validation failed: %s' % 
                             microbiome_prep.validate())
    
    return microbiome_prep


def crud_proteome(prep, md5sum, sample_id, conf, metadata):
    """Creates an iHMP OSDF Proteome object if it doesn't exist or updates
    an already existing Proteome object with the provided metadta.

    This function is different from the other create_or_update_* functions 
    in that the object is not saved but will instead be passed off to 
    AnADAMA2 in an attempt to parallelize the upload to the DCC.

    Args:
        prep (cutlass.MicrobAssayPrep): The Microbiome Assay Prep object that
            this Proteome object will be assocaited with.
        md5sum (string): md5 checksum for the associated sequence file.
        sample_id (string): Sample ID assocaited with this Proteome            
        conf (dict): Config dictionary containing some "hard-coded" pieces of
            metadata assocaited with all Proteomes.
        metadata (pandas.Series): Metadata associated with this Proteome.

    Requires:
        None

    Returns:
        cutlass.Proteome: The Proteome object that needs to be 
            saved.
    """
    raw_file_name = os.path.splitext(os.path.basename(metadata.get('seq_file')))[0]

    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}

    proteomes = group_osdf_objects(prep.proteomes(),
                                   'raw_url')
    proteomes = dict((os.path.splitext(os.path.basename(k))[0], v) for (k,v) 
                      in proteomes.items())
    
    proteome = proteomes.get(raw_file_name)

    if not proteome:
        proteome = cutlass.Proteome()

        ## TODO: Talk to Rick about these dummy files and how we might replace them
        ## with real files.
        req_metadata['local_peak_file'] = tempfile.NamedTemporaryFile(delete=False).name
        req_metadata['local_result_file'] = tempfile.NamedTemporaryFile(delete=False).name
        req_metadata['local_other_file'] = tempfile.NamedTemporaryFile(delete=False).name
    else:
        proteome = proteome[0]

    req_metadata.update(conf.get('proteome'))

    req_metadata['pride_id'] = "NA"
    req_metadata['sample_name'] = sample_id
    req_metadata['checksums'] = { "md5": md5sum }

    req_metadata['local_raw_file'] = metadata.get('seq_file')

    fields_to_update = get_fields_to_update(req_metadata, proteome)

    map(lambda key: setattr(proteome, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        proteome.links['derived_from'] = [prep.id]

        if not proteome.is_valid():
            raise ValueError('Proteome validation failed: %s' % 
                             proteome.validate())
    else:
        ## Really if there are no changes we don't want to save this object so 
        ## we return None to be handled downstream.
        proteome = None            

    return proteome        


def crud_host_tx_raw_seq_set(prep, md5sum, sample_id, conf, metadata):
    """Creates an iHMP OSDF HostTranscriptomicsRawSeqSet  object if it 
    doesn't exist or updates an already existing object with the provided metadta.

    This function is different from the other create_or_update_* functions 
    in that the object is not saved but will instead be passed off to 
    AnADAMA2 in an attempt to parallelize the upload to the DCC.

    Args:
        prep (cutlass.HostSeqPrep): The HostSeqPrep object that this 
            HostTranscriptomicsRawSeqSet object will be associated with.
        md5sum (string): md5 checksum for the associated sequence file.
        sample_id (string): Sample ID assocaited with this transcriptome           
        conf (dict): Config dictionary containing some "hard-coded" pieces of
            metadata assocaited with all transcriptomes
        metadata (pandas.Series): Metadata associated with this transcriptome

    Requires:
        None

    Returns:
        cutlass.HostTranscriptomicsRawSeqSet: The host transcriptome object
            to be saved.
    """
    raw_file_name = os.path.splitext(os.path.basename(metadata.get('seq_file')))[0]

    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}

    transcriptome = group_osdf_objects(prep.derivations(),
                                       'urls')
    transcriptome = dict((os.path.splitext(os.path.basename(k))[0], v) for (k,v) 
                          in transcriptome.items())
    
    transcriptome = transcriptome.get(raw_file_name)

    if not transcriptome:
        transcriptome = cutlass.HostTranscriptomicsRawSeqSet()
    else:
        transcriptome = transcriptome[0]

    req_metadata.update(conf.get('host_transcriptome'))
    req_metadata['checksums'] = { "md5": md5sum }

    req_metadata['local_raw_file'] = metadata.get('seq_file')
    req_metadata['size'] = os.path.getsize(raw_file_name)

    req_metadata['tags'] = []
    biopsy_location = metadata.get('biopsy_location')
    req_metadata['tags'].append('biopsy_location: %s' % biopsy_location)

    fields_to_update = get_fields_to_update(req_metadata, transcriptome)
    map(lambda key: setattr(transcriptome, key, req_metadata.get(key)),
        fields_to_update)

    if fields_to_update:
        transcriptome.links['sequenced_from'] = [prep.id]

        if not transcriptome.is_valid():
            raise ValueError('Host transcriptome validation failed: %s' % 
                             transcriptome.validate())
    else:
        ## Really if there are no changes we don't want to save this object so 
        ## we return None to be handled downstream.
        proteome = None            

    return transcriptome        


def crud_wgs_raw_seq_set(prep, md5sum, sample_id, conf, metadata):
    """Creates an iHMP OSDF WGSRawSeqSet object if it doesn't exist or 
    updates an already existing object with the provided metadta.

    This function is different from the other create_or_update_* functions 
    in that the object is not saved but will instead be passed off to 
    AnADAMA2 in an attempt to parallelize the upload to the DCC.

    Args:
        prep (cutlass.WgsDnaPrep): The WgsDnaPrep object that this 
            WGSRawSeqSet object will be associated with.
        md5sum (string): md5 checksum for the associated sequence file.
        sample_id (string): Sample ID assocaited with this transcriptome           
        conf (dict): Config dictionary containing some "hard-coded" pieces of
            metadata assocaited with all transcriptomes
        metadata (pandas.Series): Metadata associated with this transcriptome

    Requires:
        None

    Returns:
        cutlass.WgsRawSeqSet: The WGS sequence data to be saved.
    """
    raw_file_name = os.path.splitext(os.path.basename(metadata.get('seq_file')))[0]

    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}

    metagenomes= group_osdf_objects(prep.child_seq_sets(),
                                    'urls')
    metagenomes = dict((os.path.splitext(os.path.basename(k))[0], v) for (k,v) 
                          in metagenomes.items())
    
    metagenome = metagenomes.get(raw_file_name)

    if not metagenome:
        metagenome = cutlass.WgsRawSeqSet()
    else:
        metagenome = metagenome[0]

    req_metadata.update(conf.get('metagenome'))
    req_metadata['local_file'] = metadata.get('seq_file')
    req_metadata['size'] = os.path.getsize(metadata.get('seq_file'))
    req_metadata['checksums'] = { "md5": md5sum }

    fields_to_update = get_fields_to_update(req_metadata, metagenome)
    map(lambda key: setattr(metagenome, key, req_metadata.get(key)),
        fields_to_update)

    metagenome.updated = False
    if fields_to_update:
        metagenome.updated = True
        metagenome.links['sequenced_from'] = [prep.id]

        if not metagenome.is_valid():
            raise ValueError('WGS raw seq set validation failed: %s' % 
                             metagenome.validate())

    return metagenome


def crud_microb_transcriptomics_raw_seq_set(prep, md5sum, sample_id, conf, metadata):
    """Creates an iHMP OSDF MicrobTranscriptomicsRawSeqSet object if it 
    doesn't exist or updates an already existing object with the provided
    metadta.

    This function is different from the other create_or_update_* functions 
    in that the object is not saved but will instead be passed off to 
    AnADAMA2 in an attempt to parallelize the upload to the DCC.

    Args:
        prep (cutlass.WgsDnaPrep): The WgsDnaPrep object that this 
            WGSRawSeqSet object will be associated with.
        md5sum (string): md5 checksum for the associated sequence file.
        sample_id (string): Sample ID assocaited with this transcriptome           
        conf (dict): Config dictionary containing some "hard-coded" pieces of
            metadata assocaited with all transcriptomes
        metadata (pandas.Series): Metadata associated with this transcriptome

    Requires:
        None

    Returns:
        cutlass.MicrobTranscriptomicsRawSeqSet: The Microbe Transcriptomics object.
    """
    raw_file_name = os.path.splitext(os.path.basename(metadata.get('seq_file')))[0]

    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}

    metatranscriptomes = group_osdf_objects(prep.child_seq_sets(),
                                            'urls')
    metatranscriptomes = dict((os.path.splitext(os.path.basename(k))[0], v) for (k,v) 
                             in metatranscriptomes.items())
    
    metatranscriptome = metatranscriptomes.get(raw_file_name)

    if not metatranscriptome:
        metatranscriptome = cutlass.MicrobTranscriptomicsRawSeqSet()
    else:
        metatranscriptome = metatranscriptome[0]

    req_metadata.update(conf.get('metatranscriptome'))
    req_metadata['local_file'] = metadata.get('seq_file')
    req_metadata['size'] =  os.path.getsize(raw_file_name)
    req_metadata['checksums'] = { "md5": md5sum }

    fields_to_update = get_fields_to_update(req_metadata, metatranscriptome)
    map(lambda key: setattr(metatranscriptome, key, req_metadata.get(key)),
        fields_to_update)

    metatranscriptome.updated = False

    if fields_to_update:
        metatranscriptome.updated = True
        metatranscriptome.links['sequenced_from'] = [prep.id]

        if not metatranscriptome.is_valid():
            raise ValueError('Microbe Transcritpome validation failed: %s' % 
                             metatranscriptome.validate())

    return metatranscriptome


def crud_abundance_matrix(session, dcc_parent, md5sum, sample_id, 
                          study_name, conf, metadata):
    """Creates or updates an iHMP OSDF AbundanceMatrix object.

    Handles abundance matrices in both BIOM and tab-delimited format.

    Args:
        session (cutlass.Session): The OSDF session instance.
        dcc_parent (cutlass.<SEQ OR ASSAY PREP OBJECTS>): Any OSDF sequence set object 
            or an assay prep from which an abundance matrice may be derived.
             (i.e. WgsRawSeqSet or HostAssayPrep)
        md5sum (string): md5 checksum for the associated sequence file.
        sample_id (string): Sample ID assocaited with this transcriptome           
        conf (dict): Config dictionary containing some "hard-coded" pieces of
            metadata assocaited with all transcriptomes
        metadata (pandas.Series): Metadata associated with this transcriptome

    Requires:
        None

    Returns:
        cutlass.AbundanceMatrix: The abundance matrice to be saved.
    """
    abund_file = metadata.get('out_file')
    abund_fname = os.path.splitext(os.path.basename(abund_file))[0]

    data_type = metadata.get('data_type')

    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}

    if "HostAssayPrep" in dcc_parent.__module__:    
        abund_matrixs = group_osdf_objects(_get_host_assay_prep_abund_matrices(session, dcc_parent.id),
                                           '_urls')
    elif "WgsRawSeqSet" in dcc_parent.__module__:
        abund_matrices = group_osdf_objects(_get_wgs_raw_seq_set_abund_matrices(session, dcc_parent.id),
                                            '_urls')

    abund_matrices = dict((os.path.splitext(os.path.basename(k))[0], v) for (k,v) 
                           in abund_matrices.items())
    abund_matrix = abund_matrices.get(abund_fname)

    if not abund_matrix:
        abund_matrix = cutlass.AbundanceMatrix()
    else:
        abund_matrix = abund_matrix[0]

    req_metadata.update(conf.get('abundance_matrix'))
    req_metadata['local_file'] = metadata.get('out_file')
    req_metadata['size'] = os.path.getsize(abund_file)
    req_metadata['checksums'] = { "md5": md5sum }
    req_metadata['study'] = study_name

    if abund_file.endswith('biom'):
        req_metadata['format'] = "biom"
        req_metadata['format_doc'] = "http://biom-format.org/"
    elif abund_file.endswith('.tsv'):
        req_metadata['format'] = "tbl"
    else:
        raise ValueError("Unknown abundance matrix type:", raw_file_name)

    ## This is really ugly but just about the cleanest way I know to do this..
    if "taxonomic_profile" in abund_file:
        if data_type == "metagenomics":
            req_metadata['matrix_type'] = "wgs_community"
        elif data_type == "amplicon":
            req_metadata['matrix_type'] = "16s_community"
    elif "path" in abund_file or "gene" in abund_file:
        if data_type == "metatranscriptomics":
            req_metadata['matrix_type'] = "microb_metatranscriptome"
        else:            
            req_metadata['matrix_type'] = "wgs_functional"
    elif data_type == "host_transcriptome":
        req_metadata['matrix_type'] = "host_transcriptome"
    elif data_type == "proteomics":
        req_metadata['matrix_type'] = "microb_proteomic"
    elif data_type == "metabolomics":
        req_metadata['matrix_type'] = "microb_metabolome"

    fields_to_update = get_fields_to_update(req_metadata, abund_matrix)
    map(lambda key: setattr(abund_matrix, key, req_metadata.get(key)),
        fields_to_update)

    abund_matrix.updated = False

    if fields_to_update:
        abund_matrix.updated = True
        abund_matrix.links['computed_from'] = [dcc_parent.id]

        if not abund_matrix.is_valid():
            raise ValueError('Abundance matrix validation failed: %s' %
                             abund_matrix.validate())
    else:
        abund_matrix.urls = abund_matrix._urls

    return abund_matrix

def crud_viral_seq_set(prep, md5sum, sample_id, conf, metadata):
    """Creates an iHMP OSDF ViralSeqSet object if it doesn't exist or updates 
    an already existing object with the provided metadta.

    This function is different from the other create_or_update_* functions 
    in that the object is not saved but will instead be passed off to 
    AnADAMA2 in an attempt to parallelize the upload to the DCC.

    Args:
        prep (cutlass.WgsDnaPrep): The WgsDnaPrep object that this 
            ViromicsSetSeq object will be associated with.
        md5sum (string): md5 checksum for the associated sequence file.
        sample_id (string): Sample ID assocaited with this transcriptome           
        conf (dict): Config dictionary containing some "hard-coded" pieces of
            metadata assocaited with all transcriptomes
        metadata (pandas.Series): Metadata associated with this transcriptome

    Requires:
        None

    Returns:
        cutlass.ViralSeqSet: The ViralSeqSet object to be saved.
    """
    raw_file = metadata.get('seq_file')
    raw_file_name = os.path.splitext(os.path.basename(raw_file))

    ## Setup our 'static' metadata pulled from our YAML config
    req_metadata = {}

    virome = group_osdf_objects(prep.child_seq_sets(),
                                'urls')
    virome = dict((os.path.splitext(os.path.basename(k))[0], v) for (k,v) 
                  in virome.items())
    
    virome = virome.get(raw_file_name)

    if not virome:  
         virome = cutlass.ViralSeqSet()
    else:
         virome = virome[0]

    req_metadata.update(conf.get('virome'))
    req_metadata['checksums'] = { "md5": md5sum }
    req_metadata['local_raw_file'] =  raw_file

    req_metadata['tags'] = []

    fields_to_update = get_fields_to_update(req_metadata, virome)
    map(lambda key: setattr(virome, key, req_metadata.get(key)),
        fields_to_update)

    virome.updated = False
    if fields_to_update:
        virome.updated = True
        virome.links['sequenced_from'] = [prep.id]

        if not virome.is_valid():
            raise ValueError('Virome validation failed: %s' % 
                             virome.validate())

    return virome       

