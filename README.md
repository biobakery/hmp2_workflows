# HMP2 Workflows

HMP2 workflows is a collection of workflows that handle processing several 
different data types assocaited with the IBDMDB project. 

These workflows are built using AnADAMA2 which allows for parallel task 
execution locally and in a grid compute environment.     
    

[TOC]
* * * 
### Requirements

The following Python modules are required for the project:
* anadama2 HEAD
* biobakery\_workflows HEAD
* pandas 0.19.2
* pyyaml 3.12
* biom-format 2.1.5

Additionally the following software should be installed to launch a production
copy of the website:

* kneaddata
* HUMAnN2
* MetaPHlAn2
* samtools
* USEARCH
* PICRUSt v1.1.
* Clustal Omega

Installing all python modules in a python virtualenv is recommended to create 
an isolated environment free of any python or python module version mismatching.

* virtualenv >= 1.10.1

### Installation

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
