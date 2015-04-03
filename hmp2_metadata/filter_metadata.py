from glob import glob
from doit import get_var
import csv
import os
import sys

def get_spreadsheet():
    spreadsheet = get_var('s', None)
    if not spreadsheet:
        print >> sys.stderr, "no spreadsheet given "
        return None, None, None
    tsvdir = os.path.dirname(spreadsheet)
    tsvfile = os.path.splitext(os.path.basename(spreadsheet))[0]
    tsvfile = os.path.join(tsvdir, tsvfile)
    tsvfile = tsvfile + ".tsv"
    filteredfile = tsvfile + "_filtered" + ".tsv"
    #print >> sys.stderr, "filteredfile: " + filteredfile
    return spreadsheet, tsvfile, filteredfile

def task_make_tsv():

    spreadsheet, tsvfile, filteredfile = get_spreadsheet()
    if not spreadsheet:
        return  {
            "actions": [("date")],
            "verbosity": 2,
            "file_dep": [],
#            "targets": [tsvfile]
        }

    if os.path.splitext(spreadsheet)[1] == ".xlsx":
        return  {
            "actions": [("xlsx2csv -d tab \"" + str(spreadsheet) + "\" \"" + str(tsvfile) + "\"")],
            "verbosity": 2,
            "file_dep": [str(spreadsheet)],
            "targets": [tsvfile]
            }
    else:
        return  {
            "actions": [("date")],
            "verbosity": 2,
            "file_dep": [],
            "targets": [tsvfile]
        }


def task_filter_dates():

    def filter(tsvfile):
        """ create symlink to bamfile """
        tsvdir = os.path.dirname(spreadsheet)
        with open(tsvfile, "r") as input:
            columns = input.readline().split("\t")
            indexes = [i for i,x in enumerate(columns) if "Date" in x]

        with open(tsvfile, "r") as input:
            with open(filteredfile, "w") as output:
                for line in input:
                    line = line.rstrip()
                    for i, col in enumerate(line.split("\t")):
                        if i in indexes:
                            continue
                        output.write("%s\t" % str(col))
                    output.write("\n")

    spreadsheet, tsvfile, filteredfile = get_spreadsheet()

    return {
        "actions": [(filter, [ str(tsvfile) ])],
        "verbosity": 2,
        "file_dep": [ str(tsvfile) ],
        "targets": [filteredfile]
    }

def task_publish_filtered_metadata():

    PROJECT_METADATA_FILE="/seq/ibdmdb/public/HMP2/Metadata/hmp2_project_metadata.tsv"
    spreadsheet, tsvfile, filteredfile = get_spreadsheet()

    return {
        "actions": ['ln -s %s ' % filteredfile + " %s " % PROJECT_METADATA_FILE],
        "verbosity": 2,
        "file_dep": [ str(filteredfile) ],
        "targets": [PROJECT_METADATA_FILE]
    }

DOIT_CONFIG = {
   'default_tasks': ["make_tsv", "filter_dates", "publish_filtered_metadata"],
   'pipeline_name': "Broad Internal Data Build",
   'continue': True
}
