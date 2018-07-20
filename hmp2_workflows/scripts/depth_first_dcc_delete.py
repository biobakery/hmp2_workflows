# -*- coding: utf-8 -*-

"""
depth_first_dcc_delete.py
~~~~~~~~~~~~~~~~~~~~~~~~~

Does a depth-first traversal/delete of the document tree out of the OSDF 
based off an OQL query to provide the starting point for the delete.

Example OQL queries are:

 - '"abundance_matrix"[node_type] && "wgs_community"[meta.matrix_type]'
 - '"visit"[node_type]'

Extreme caution should be taken when using the script as if a specific enough 
OQL query is not provided a large amount of documents could be deleted.

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


import argparse
import importlib
import itertools
import json
import types

import anytree
import cutlass


def parse_cli_arguments():
    """Parses any command-line arguments passed into the workflow.

    Args:
        None
    Requires:
        None
    Returns:
        argparse.ArgumentParser: Object containing arguments passed in by user.
    """
    parser = argparse.ArgumentParser('Deletes a cache of OSDF documents in a '
                                     'depth-first manner using a specific OQL '
                                     'query as the basis.')
    parser.add_argument('-u', '--username', required=True,
                        help='DCC username.')
    parser.add_argument('-p', '--password', required=True,
                        help='DCC password.')
    parser.add_argument('-s', '--study-id', required=True, 
                        default='52d8c92f2d3660b9add954d544a0216e',
                        help='The study ID from which to cascade delete down.')
    parser.add_argument('-t', '--node-type-filter', 
                        help='OPTIONAL. Filter nodes to delete by a node type.')
    parser.add_argument('-q', '--oql-query',
                        help='OQL query to establish the basis from which '
                        'to do a deletion.')
    parser.add_argument('-d', '--dry-run', action='store_true', default=False, 
                        help='Perform a dry-run deletion and list which '
                        'nodes will be delete.')
    parser.add_argument('--delete-root', action='store_true', default=False,
                        help='Delete the root node when this flag is specified.')

    return parser.parse_args()


def build_osdf_tree(study_id):
    """
    Builds a tree structure to contain all our targeted objects from the OSDF.

    Args:
        study_id (string): The study ID to act as the root node from which 
            all children nodes are retrieved.

    Requires:
        None

    Returns:
        anytree.Node: The root tree node of our constructed tree structure.
    """
    ## We want a hash lookup map so that we can map OSDF ID -> OSDF object so
    ## can quickly build our tree without some messy looping.
    osdf_lookup_map = {}

    def _update_osdf_tree(osdf_obj):
        """
        Updates the tree containing OSDF nodes and makes use of lookup table
        set the proper parent node.

        Args:
            osdf_obj (cutlass.*): One of the cutlass objects associted with
                the provided study.
            
        Requires: 
            None

        Returns: 
            None
        """
        ## We should only ever have one parent for each of our objects so 
        ## this is naieve but ok.
        parent_id = osdf_obj.links.values()[0][0]
        parent_node = osdf_lookup_map.get(parent_id)

        if not parent_node:
            print "WARNING: Could not find parent node for following object:", osdf_obj.to_json()
        else:

            osdf_node = anytree.Node(osdf_obj.id, osdf=osdf_obj, type=json.loads(osdf_obj.to_json()).get('node_type'), parent=parent_node)
            return osdf_node

    study_obj = cutlass.Study.load(study_id)
    study_node = anytree.Node("root", osdf=study_obj, type='study')
    
    subjects = list(study_obj.subjects())
    subject_nodes = [anytree.Node(s.id, osdf=s, parent=study_node, type='subject') for s in subjects]
    osdf_lookup_map.update({s.name: s for s in subject_nodes})

    subject_attrs = [list(s.attributes()) for s in subjects]
    subject_attrs = list(itertools.chain.from_iterable(subject_attrs))
    subject_attr_nodes = map(_update_osdf_tree, subject_attrs)
    osdf_lookup_map.update({sa.name: sa for sa in subject_attr_nodes})

    visits = [list(s.visits()) for s in subjects]
    visits = list(itertools.chain.from_iterable(visits))
    visit_nodes = map(_update_osdf_tree, visits)
    osdf_lookup_map.update({v.name: v for v in visit_nodes})

    visit_attrs = [list(v.visit_attributes()) for v in visits]
    visit_attrs = list(itertools.chain.from_iterable(visit_attrs))
    visit_attr_nodes = map(_update_osdf_tree, visit_attrs)
    osdf_lookup_map.update({va.name: va for va in visit_attr_nodes})

    samples = [list(v.samples()) for v in visits]
    samples = list(itertools.chain.from_iterable(samples))
    sample_nodes = map(_update_osdf_tree, samples)
    osdf_lookup_map.update({sp.name: sp for sp in sample_nodes})

    sample_attrs = [list(s.sampleAttributes()) for s in samples]
    sample_attrs = list(itertools.chain.from_iterable(sample_attrs))
    sample_attr_nodes = map(_update_osdf_tree, sample_attrs)
    osdf_lookup_map.update({sa.name: sa for sa in sample_attr_nodes})

    preps = [list(s.preps()) for s in samples]
    preps = list(itertools.chain.from_iterable(preps))
    prep_nodes = map(_update_osdf_tree, preps)
    osdf_lookup_map.update({p.name: p for p in prep_nodes})

    seq_sets = [list(c) for p in preps for c in p.children() if p.children()]
    seq_sets = list(itertools.chain.from_iterable(seq_sets))
    seq_sets = [ss for ss in seq_sets if not isinstance(ss, types.GeneratorType)]
    seq_set_nodes = map(_update_osdf_tree, seq_sets)
    osdf_lookup_map.update({ss.name: ss for ss in seq_set_nodes})

    products = [list(c) for ss in seq_sets for c in ss.children() if ss.children()]
    products = list(itertools.chain.from_iterable(products))
    products = [p for p in products if not isinstance(p, types.GeneratorType)]
    product_nodes = map(_update_osdf_tree, products)
    osdf_lookup_map.update({po.name: po for po in product_nodes})

    ## Sometimes we have another round of products we need to account for here...
    products2 = [list(c) for p in products for c in p.children() if p.children()]
    products2 = list(itertools.chain.from_iterable(products2))
    products2 = [p for p in products2 if not isinstance(p, types.GeneratorType)]
    product_nodes2 = map(_update_osdf_tree, products2)
    osdf_lookup_map.update({po.name: po for po in product_nodes2})

    return study_node


def filter_osdf_tree(root_node, node_type):
    """Filters an existing OSDF tree by a specific node_type.

    Args:
        root_node (anytree.Node): The root node of the tree to filter upon
        node_type (string): THe type of node to filter the tree down by.

    """
    filtered_nodes = anytree.search.findall(study_node, 
                                            filter_=lambda node: node.type in 
                                                node_type)

    filtered_root_node = anytree.Node('root', osdf=root_ndoe.osdf, type='study')

    for filtered_node in filtered_nodes:
        path = filtered_node.path[1:]
        
        parent_node = filtered_root_node
        for node in path:
            new_node = anytree.Node(node.osdf.id, osdf=node.osdf, type=node.type, parent=parent_node)
            parent_node = new_node 

    return filtered_root_node


def delete_nodes(root_node, dry_run, delete_root, stop_node="root"):
    """
    Cascade deletes OSDF nodes in a depth-first manner.

    Args:
        root_node (anytree.Node): The root node of the tree to delete from.
        dry_run (boolean): True/False to delete or print out the nodes to 
            be deleted.
        stop_node (string): Name of the node to stop deletion on. Defaults to 
            the root node but can be any OSDF ID to stop on.
        delete_root (boolean): If the stop_node parameter is set to 'root'
            this parameter can be passed to indicate we want to delete the 
            root node as well.

    Requires:
        None

    Returns:
        None
    """
    deleted = []
    failed_delete = []

    for node in anytree.PostOrderIter(root_node):
        if dry_run:
            if not node.name == "root":
                print "DELETING NODE:", node
        else:
            osdf_obj = node.osdf

            if not node.name == "root":
                print "DELETING NODE:", node
                res = osdf_obj.delete()
               
                if not res:
                    print "FAILED TO DELETE NODE:", node
                    failed_delete.append(osdf_obj)

    if failed_delete:
        print "WARNING: The following OSDF nodes were not deleted:" + "\n".join(failed_delete)

    if delete_root and stop_node == "root":
        print "DELETING ROOT NODE:", root_node

        if not dry_run:
            root_node.osdf.delete()


def main(args):
    session = cutlass.iHMPSession(args.username, args.password, ssl=False)
    osdf = session.get_osdf()

    root_node = build_osdf_tree(args.study_id)

    if args.node_type_filter:
        root_node = filter_osdf_tree(root_node, args.node_type_filter)
            
    delete_nodes(root_node, args.dry_run, args.delete_root)


if __name__ == "__main__":
    main(parse_cli_arguments())
