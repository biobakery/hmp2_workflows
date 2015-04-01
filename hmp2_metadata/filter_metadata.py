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
    csvdir = os.path.dirname(spreadsheet)
    csvfile = os.path.splitext(os.path.basename(spreadsheet))[0]
    csvfile = os.path.join(csvdir, csvfile)
    csvfile = csvfile + ".csv"
    filteredfile = csvfile + "_filtered" + ".csv"
    print >> sys.stderr, "filteredfile: " + filteredfile
    return spreadsheet, csvfile, filteredfile

def task_make_csv():

    spreadsheet, csvfile, filteredfile = get_spreadsheet()
    if not spreadsheet:
        return  {
            "actions": [("date")],
            "verbosity": 2,
            "file_dep": [],
#            "targets": [csvfile]
        }

    if os.path.splitext(spreadsheet)[1] == "xlsx":
        return  {
            "actions": [("xlsx2csv -d tab " + str(spreadsheet) + " " + str(csvfile))],
            "verbosity": 2,
            "file_dep": [str(spreadsheet)],
            "targets": [csvfile]
            }
    else:
        return  {
            "actions": [("date")],
            "verbosity": 2,
            "file_dep": [],
            "targets": [csvfile]
        }


def task_filter_dates():

    def filter(csvfile):
        """ create symlink to bamfile """
        csvdir = os.path.dirname(spreadsheet)
        with open(csvfile, "r") as input:
            with open(filteredfile, "w") as output:
                for line in input:
                    #print >> sys.stderr, "line: " + str(line)
                    cols = line.split(',')
                    #import pdb;pdb.set_trace()
                    fline = cols[:4] + cols[5:7] + cols[8:13] + cols[14:16] + cols[22:] 
                    #print >> sys.stderr, "fline: " + str(fline)
                    data="\r\n%s" % str(fline)
                    output.write(data)

    spreadsheet, csvfile, filteredfile = get_spreadsheet()

    return {
        "actions": [(filter, [ str(csvfile) ])],
        "verbosity": 2,
        "file_dep": [ str(csvfile) ],
        "targets": [filteredfile]
    }


DOIT_CONFIG = {
   'default_tasks': ["make_csv", "filter_dates"],
   'pipeline_name': "Broad Internal Data Build",
   'continue': True
}
