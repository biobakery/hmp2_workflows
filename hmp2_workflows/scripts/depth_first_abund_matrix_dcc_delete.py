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
    parser.add_argument('-q', '--oql-query', required=True,
                        help='OQL query to establish the basis from which '
                        'to do a deletion.')
    parser.add_argument('-s', '--stop-obj-class', default='study',
                        help='The cutlass object to stop deletion on. '
                        'defaults to the Study class.')
    parser.add_argument('-d', '--dry-run', action='store_true', default=False, 
                        help='Perform a dry-run deletion and list which '
                        'nodes will be delete.')

    return parser.parse_args()


def main(args):
    session = cutlass.iHMPSession('cesar.arze', 'fishnet-socket-dolphin-can', ssl=False)
    osdf = session.get_osdf()

    raw_abund_matrices = osdf.oql_query('ihmp', '"abundance_matrix"[node_type] && "host_transcriptome"[meta.matrix_type]')
    
    abundance_matrices = map(cutlass.AbundanceMatrix.load, [m.get('id') for m in raw_abund_matrices.get('results')])
    host_raw_seq_sets = map(cutlass.HostTranscriptomicsRawSeqSet.load,
                            [m.links.get('computed_from')[0] for m in abundance_matrices])
    host_seq_preps = map(cutlass.HostSeqPrep.load,
                         [seq_set.links.get('sequenced_from')[0] for seq_set in host_raw_seq_sets])
    samples = map(cutlass.Sample.load,
                  [prep.links.get('prepared_from')[0] for prep in host_seq_preps])
    sample_attrs = map(cutlass.Sample.sampleAttributes, samples)
    visits = map(cutlass.Visit.load,
                 [sample.links.get('collected_during')[0] for sample in samples])
    visit_attrs = map(cutlass.Visit.visit_attributes, visits)

    host_tx_tree = zip(visits, visit_attrs, samples, sample_attrs, host_seq_preps, host_raw_seq_sets, abundance_matrices)

    for (visit, visit_attr_g, sample, sample_attr_g, seq_prep, seq_set, abundance_matrix) in host_tx_tree:
        print "DELETING the following objects:"
        print "\t - %s" % visit.visit_id
        print "\t - %s" % sample.name
        print "\t - %s" % seq_prep.prep_id
        print "\t - %s" % seq_set.comment
        print "\t - %s" % abundance_matrix._urls
        print
        
        if not args.dry_run:
            sample_attr = next(sample_attr_g, None)
            visit_attr = next(visit_attr_g, None)

            abundance_matrix.delete()
            seq_set.delete()
            seq_prep.delete()
            if sample_attr:
                sample_attr.delete()
            else:
                print "DEBUG: No sample_attr"
            sample.delete()
            if visit_attr:
                visit_attr.delete()
            else:
                print "DEBUG: No visit_attr"
            visit.delete()


if __name__ == "__main__":
    main(parse_cli_arguments())
