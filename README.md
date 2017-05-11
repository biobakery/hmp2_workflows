# HMP2 Workflows

HMP2 workflows is a collection of workflows that handle processing several 
different data types assocaited with the IBDMDB project. 

These workflows are built using AnADAMA2 which allows for parallel task 
execution locally and in a grid compute environment.     

##

[TOC]

----

## Requirements

The following Python modules are required for the project:

* [anadama2](https://bitbucket.org/biobakery/anadama2) *HEAD*
* [biobakery\_workflows](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home) *HEAD*
* [pandas](http://pandas.pydata.org/) *0.19.2*
* [pyyaml](http://pyyaml.org/) *3.12*
* [biom-format](http://biom-format.org/) *2.1.5*
* [cutlass](https://github.com/ihmpdcc/cutlass)

Additionally the following software should be installed to launch a production
copy of the website:

* [kneaddata](https://bitbucket.org/biobakery/kneaddata)
* [HUMAnN2](https://bitbucket.org/biobakery/humann2)
* [MetaPHlAn2](https://bitbucket.org/biobakery/metaphlan2)
* [samtools](http://www.htslib.org/)
* [USEARCH](http://www.drive5.com/usearch/)
* [PICRUSt](http://picrust.github.io/picrust/) *1.1*
* [Clustal Omega](http://www.ebi.ac.uk/Tools/msa/clustalo/)

Installing all python modules in a python virtualenv is recommended to create 
an isolated environment free of any python or python module version mismatching.

* [virtualenv](https://virtualenv.pypa.io/en/stable/) >= *1.10.1*

## Installation

Currently installation of the HMP2 workflows is only available through 
download/clone from the bitbucket repository.

1. Download packaged repository [hmp2_workflows.zip](https://bitbucket.org/biobakery/hmp2_workflows/get/cc70eb41860b.zip) OR 
clone repository: `git clone git@bitbucket.org:biobakery/hmp2_workflows.git`
    a. If archive file is downloaded extract: `unzip hmp2_workflows.zip`
    b. cd hmp2_workflows/
2. Install required python modules using pip via requirements.txt: `pip install -r requirements.txt`
    a. OPTIONAL: If using virtualenv, initialize a python virtualenv for project: `virtualenv /home/carze/.venvs/hmp2_biobakery_workflows`
        1. Activate vritualenv before installing any modules: `source /home/carze/.venvs/hmp2_biobaker_workflows/bin/activate`
3. Install software required by biobakery_workflows needed for analysis
    a. Instructions for each piece of software should be available with their distributions.

### Automating Workflow Execution    

## How to Run

## Configuration


## Manifest files

Each of the workflows primarily relies on a manifest file that describes the files that are to be processed. An example
manifest file can be found below:

```
## This is an example manifest file for an HMP2 AnADAMA2 workflow run
##
## The file contains metadata that identifies the origin of the files,
## the date files were generated, submitted, a contact person for the
## data and which types of data are present in this batch of files.

################
#   METADATA   #
################

# First some information about who generated/submitted this data
origin_institute: PNNL
origin_contact: Richard White
origin_contact_email: richard.white@pnnl.gov

project: HMP2
date_of_creation: 2017-04-17T17:07:00
date_of_submission: 2017-04-17T17:30:00

################
#     DATA     #
################

# The files present in this submission. These will be processed by a
# the appropriate AnADAMA2 pipelines.
submitted_files:
    proteomics:
        md5sums_file: /data/ibdmdb/upload/HMP2/Proteomics/test.md5sums.txt
        input:
            - /data/ibdmdb/upload/HMP2/Proteomics/160626-SM-A62DJ-47.raw
            - /data/ibdmdb/upload/HMP2/Proteomics/160618-SM-AIG7A-51.raw
            - /data/ibdmdb/upload/HMP2/Proteomics/160624-SM-73BO3-149.raw
```

The **METADATA** section provides several pieces of metadata that inform each workflow the source 
of data files, which project the files belong too and when the files were created and submitted.

The **DATA** section provides paths to the actual files broken up into sections for each data type as well as 
whether or not an accompanying file containing md5 checksums for all submitted files is available.

## Workflows

HMP2 Workflows contains a collection of workflows that handle processing, analysis and dissemination of several data types 
associated with the IBDMDB project. The analysis workflows leverage the existing biobakery\_workflows heavily and wrap 
the steps in these workflows with integrity checks, formatting and pushing data to the proper locations.

### Metagenomics

Executes several steps from the biobakery\_workflows [Whole Metagenome Shotgun workflow](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!whole-metagenome-shotgun-wmgx).

### Metatranscriptomics

Executes several steps from the biobakery\_workflows [Metatranscriptomics workflow](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!whole-metagenome-and-metatranscriptome-shotgun-wmgx_wmtx).

### 16S

Executes several steps from the biobakery\_workflows [16S workflow](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!16s-rrna-16s)

### Proteomics

This workflow handles processing any proteomics data submitted to the IBDMDB. Proteomics data suubmitted 
to the IBDMDB is received process and with accompanying output files so this workflow works to verify
the integrity of the uploaded files and disseminates them to the proper locations to make them 
available on the IBDMDB website.

#### Inputs
1. Manifest file
    a. Containing proteomics data files to upload
2. Metadata file containing metadata associated with data files in the manifest file.
3. Proteomics peptide roll-up output files assocaited with each spectra file.

#### Outputs
*None*

#### Example Run
```
python proteomics.py --manifest-file /data/ibdmdb/upload/MANIFEST.hmp2-proteomics.yaml 
                     --config-file config/analysis.yaml
```

### DCC Upload

This workflow handles uploading data files and their accompanying metadata to the iHMP Data Coordinating Center (DCC) OSDF instance. 
The upload procedure is handled using the iHMP's cutlass tool which leverages the Aspera file transfer protocol to transfer the 
larger data files.

#### Inputs
1. Manifest file
    a. Containing data files to upload OR 
    b. Data files associated with updated metadata
2. Metadata file containing metadata assocaited with data files in manifest file.

#### Outputs
*None*

#### Example Run
```
python upload_dcc.py --metadata-file hmp2_project_metadata_2016-10-15.hr.csv 
                     --broad-data-sheet IBDMDB_AliquotsAsOf_8.24.2016_updatestatus_1.17.2017_newsamplestatus.csv 
                     --config-file analysis.yaml 
                     --manifest-file /data/ibdmdb/upload/MANIFEST.hmp2-proteomics.yaml
```

### SRA Upload
*TODO*

### Metadata Refresh

This workflow refreshes the master metadata file using clinical metadata from the StudyTrax database, a data file tracking spreadsheet from the Broad Institute, 
a Broad LabOps data generation status file and auxillary metadata files for specific data types (i.e. Proteomics). These files coupled with a set of data files 
that have corresponding new or updated metadata can be used to trigger a metadata refresh

#### Inputs
1. Manifest file
2. StudyTrax metadata table
3. Broad sample tracking spreadsheet
4. Broad LabOps data generation status file
5. OPTIONAL: Auxillary metadata tables

#### Outputs
1. Refreshed project metadata table

#### Example Run
```
python metadata.py --studytrax-metadata /data/ibdmdb/metadata/studytrax.csv
                   --broad-data-sheet /data/ibdmdb/metadata/broad_sample_tracking.csv
                   --tableau-data-status /data/ibdmdb/metadata/labsop_status.csv
                   --config-file /data/ibdmdbd/config/analysis.yaml
```
