from glob import glob
from doit import get_var
import csv
import os
import sys

#print >> sys.stderr, "spreadsheet: " + str(spreadsheet)
#csvdir = os.path.dirname(spreadsheet)
#print >> sys.stderr, "csvdir: " + csvdir
#csvfile = os.path.splitext(os.path.basename(spreadsheet))[0]
#csvfile = os.path.join(csvdir, csvfile)
#csvfile = csvfile + ".csv"
#print >> sys.stderr, "csvfile: " + csvfile

def get_spreadsheet():
    spreadsheet = get_var('s', None)
    if not spreadsheet:
        print >> sys.stderr, "no spreadsheet given "
        return None, None
    csvdir = os.path.dirname(spreadsheet)
    #print >> sys.stderr, "csvdir: " + csvdir
    csvfile = os.path.splitext(os.path.basename(spreadsheet))[0]
    csvfile = os.path.join(csvdir, csvfile)
    csvfile = csvfile + ".csv"
    return spreadsheet, csvfile

def task_make_csv():

    spreadsheet, csvfile = get_spreadsheet()
    if os.path.splitext(spreadsheet)[1] == "xlsx":
        return  {
            "actions": [("xlsx2csv " + str(spreadsheet) + " " + str(csvfile))],
            "verbosity": 2,
            "file_dep": [str(spreadsheet)],
            "targets": [csvfile]
            }
    else:
        return  {
            "actions": [("date")],
            "verbosity": 2,
            "file_dep": [str(spreadsheet)],
            "targets": [csvfile]
        }


def task_make_links():

    def get_bamdirs(csvfile):
        """ yield path to individual bamfiles """
        with open(csvfile, 'r') as f:
            firstline = next(f)
            projectIdx = firstline.index('Project')
            extIDIdx = firstline.index('External ID')
            if not projectIdx or not extIDIdx:
                print >> sys.stderr, "Can't process spreadsheet - missing Project or ExternalID"

            for line in csv.reader(f):
                #print >> sys.stderr, "line: " + str(line)
                path = "/seq/picard_aggregation/" 
                try:
                    path += line[projectIdx]
                except IndexError:
                    path = "/tmp/na"
                path += "/"
                try:
                    path += line[extIDIdx] 
                except IndexError: 
                    path += "/nothere.txt"
                path += "/current/" 
                yield path
    
    def link_bamfiles(csvfile):
        """ create symlink to bamfile """
        csvdir = os.path.dirname(spreadsheet)
        for bamdir in get_bamdirs(csvfile):
            print >> sys.stderr, "bamdir: " + bamdir
            for bamfile in glob(bamdir + "/*.bam"):
                #print >> sys.stderr, "found bamfile: " + bamfile
                bamfilename = os.path.basename(bamfile)
                bamfilelink = os.path.join(csvdir, bamfilename)
           
                #print >> sys.stderr, "found bamfile: " + bamfile
                if not os.path.exists(bamfilelink):
                    print >> sys.stderr, bamfilelink + " -> " + bamfile
                    os.symlink(bamfile, bamfilelink)

    spreadsheet, csvfile = get_spreadsheet()

    return {
        "actions": [(link_bamfiles, [ str(csvfile) ])],
        "verbosity": 2,
        "file_dep": [ str(csvfile) ],
        "targets": []
    }


DOIT_CONFIG = {
   'default_tasks': ["make_csv", "make_links"],
   'pipeline_name': "Broad Internal Data Build",
   'continue': True
}
