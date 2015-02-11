from glob import glob
from doit import get_var
import os
import sys

#here = "/seq/ibdmdb/centos6/tmp/broad_internal_data"
#myspreadsheet = ("/seq/ibdmdb/data_deposition/HMP2/Metadata/1450/kbayer/" +
#   "RISK_WGS_PDO-2276_SeqCompleteLaneLevel_07.25.2014_14.12.12_20.38.30.xlsx")

spreadsheet = get_var('s', None)
if not spreadsheet:
    print >> sys.stderr, "no spreadsheet given "
    exit
    
print >> sys.stderr, "spreadsheet: " + str(spreadsheet)

csvdir = os.path.dirname(spreadsheet)
print >> sys.stderr, "csvdir: " + csvdir
csvfile = os.path.splitext(os.path.basename(spreadsheet))[0]
csvfile = os.path.join(csvdir, csvfile)
csvfile = csvfile + ".csv"

print >> sys.stderr, "csvfile: " + csvfile

def task_make_csv():
    return  {
       "actions": [("xlsx2csv " +spreadsheet+ " " +csvfile)],
       "verbosity": 2,
       "file_dep": [spreadsheet],
       "targets": [csvfile]
    }

def task_make_links():

    def get_bamdirs(csvfile):
        """ yield path to individual bamfiles """
        with open(csvfile, 'r') as f:
            next(f)
            for line in f:
                print >> sys.stderr, "line: " + str(line)
                path = "/seq/picard_aggregation/" 
                try:
                    path += line.split(',')[14]
                except IndexError:
                    path = "/tmp/na"
                path += "/"
                try:
                    path += line.split(',')[3] 
                except IndexError: 
                    path += "/nothere.txt"
                path += "/current/" 
                yield path
    
    def link_bamfiles(csvfile):
        """ create symlink to bamfile """
        for bamdir in get_bamdirs(csvfile):
            print >> sys.stderr, "bamdir: " + bamdir
            for bamfile in glob(bamdir + "/*.bam"):
                print >> sys.stderr, "found bamfile: " + bamfile
                bamfilename = os.path.basename(bamfile)
                bamfilelink = os.path.join(csvdir, bamfilename)
           
                print >> sys.stderr, "found bamfile: " + bamfile
                print >> sys.stderr, "bamfilelink: " + bamfilelink
                if not os.path.exists(bamfilelink):
                    os.symlink(bamfile, bamfilelink)


    return {
        "actions": [(link_bamfiles, [csvfile])],
        "verbosity": 2,
        "file_dep": [csvfile],
        "targets": []
    }


DOIT_CONFIG = {
   'default_tasks': ["make_csv", "make_links"],
   'pipeline_name': "Broad Internal Data Build",
   'continue': True
}
