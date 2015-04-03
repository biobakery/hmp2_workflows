#!/bin/sh
##
## HMP2 Broad specific script deployed runs both anpan and mibc_tm with sane defaults
## 
##

broad_data=/seq/ibdmdb/centos6/tmp/hmp2_workflows/broad_internal_data/link_data.py

function mk_map {

  if [ -f "${InputFile}/map.txt" ]; then
    echo "map.txt file present."
  else
    files=`cd "${InputFile}" && ls | grep -v -e \.log -e metadata -e map.txt -e mibc_products | awk -F_ '{print $1}'`
    cat > "${InputFile}/map.txt" << EOF
#SampleID
${files}
EOF
  fi
}

function mk_metadata {

  if [ -f "${InputFile}/metadata.txt" ]; then
    echo "metadata.txt file present."
  else
    if [ $1 == "WGS" ]; then
      sixteenS="false"
    else
      sixteenS="true"
    fi
    files=`cd "${InputFile}" && ls | grep -v -e \.log -e map.txt -e mibc_products | xargs | tr ' ' '	' | sed 's/ $//'`
    cat > "${InputFile}/metadata.txt" << EOF
pi_first_name	Curtis
pi_last_name	Huttenhower
pi_contact_email	chuttenh@hsph.harvard.edu
lab_name	Huttenhower
researcher_first_name	data
researcher_last_name	deposition
researcher_contact_email	kbayer@hsph.harvard.edu
study_title	HMP2
study_description	16S description
sample_type	microbial
filename	${files}
collection_start_date	2013-01-01
collection_end_date	2013-01-01
submit_to_insdc	false
16s_data	$sixteenS
reverse_primer	CCGTCAATTCMTTTRAGT
platform	454
visualize	true
EOF
  fi
}

function mk_broaddata {
 
  spreadsheet=`ls ${InputFile}/*.xlsx 2>/dev/null | sed -n '1p'` 
  if [ -f "${spreadsheet}" ]; then
    echo "processing ${spreadsheet}"
    anadama run -f "${broad_data}" s="${spreadsheet}"
    return 1
  else
    return 0
  fi
}

function call_mibc 
{
  cd /seq/ibdmdb/centos6/tmp/tm
  echo "anpan dag --project='${InputFile}' -a --tmpfiledir=/seq/ibdmdb/centos6/tmp/tm | mibc_tm -t ${Type} -g ${Governor} -p 8888"
  anpan dag --project="${InputFile}" -a --tmpfiledir=/seq/ibdmdb/centos6/tmp/tm | mibc_tm -t ${Type} -g ${Governor} -p 8888
}

if [ $# -ne 3 ]; then
  echo "run_mibc.sh directory [local|lsf] governor"
  exit -1
else
  InputFile=$1
  Type=$2
  Governor=$3
fi

if [ -z "${InputFile}" ]; then
  echo "Error: InputFile env variable empty."
  exit -1
fi


data=`echo "${InputFile}" | awk -F'/' '{print $6}'`

if [ -z $data ]; then
  echo "Error: data variable empty."
  echo "InputDir: ${InputDir}"
  exit -1
fi

echo "data: $data"

case $data in

  WGS)
    mk_broaddata 
    if [ $? == 0 ]; then
      mk_map
      mk_metadata $data
      call_mibc
    fi
    ;;
  16S)
    mk_broaddata 
    if [ $? == 0 ]; then
      mk_metadata $data
      call_mibc
    fi
    ;;
  *)
    echo "Unknown data: $data in InputDir: ${InputDir}"
    echo "To support this new datatype, update script new_data.sh"
    exit -1
    ;;
esac
