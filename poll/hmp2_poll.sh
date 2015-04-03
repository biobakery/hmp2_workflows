#!/bin/sh
##
## Script searches given filesystem(s) for any new or updated file
## in the time range of [since the last time the script ran - to -
## 5 minutes ago].  If file(s) are found, it uses curl to notify
## jenkins giving it the containing directory
## as a parameter.
##
## Note: The time range attempts to prevent files that are still
## being uploaded from being processed before they are fully
## loaded onto the filesystem.
##

# Settable Parameters
TIME_DELAY_SEC=900 # 15 min
TIME_DELAY_MIN=15 # 15 min
DATA_DIR=/seq/ibdmdb/data_deposition
ANADAMA_DIR=/seq/ibdmdb/centos6
TMP_DIR=/seq/ibdmdb/centos6/tmp

##

function call_anadama {

  #echo "mibc_build dag --project='$1' -a --tmpfiledir=${ANADAMA_DIR}/tmp/tm | mibc_tm -t lsf -g 10"
  echo "mibc_run.sh $1 lsf 10"
  mibc_run.sh $1 lsf 10
  mail -s "mibc_run.sh $1 called... " kbayer@broadinstitute.org <<EOF
no body.
EOF

  mibc_run.sh $1 lsf 10 &
  #mibc_build dag --project="$1" -a --tmpfiledir=${ANADAMA_DIR}/tmp/tm | mibc_tm -t lsf -g 10

}

function init_timestamps {

  # only create initial timestamp if it doesn't exist
  if [ ! -f ${DATA_DIR}/.timestamp.empty ]; then
    echo "creating .timestamp.empty"
    touch ${DATA_DIR}/.timestamp.empty
  fi
  # uptime timestamp prior to processing so we don't loose any new events
  # during processing itself - we need them all.
  touch -d `date --date="-${TIME_DELAY_SEC} seconds" +%H:%M:%S` ${DATA_DIR}/.next_timestamp.empty

}

function finalize_timestamps {
  mv ${DATA_DIR}/.next_timestamp.empty ${DATA_DIR}/.timestamp.empty
}

# Source project specific configuration parameters
if [ -f ${ANADAMA_DIR}/bin/activate ]; then
  . ${ANADAMA_DIR}/bin/activate
else
  echo "${ANADAMA_DIR}/bin/activate file not found"
  exit -1
fi

# Are we still running (from the last time?) if so die now
if [ -f ${DATA_DIR}/.next_timestamp.empty ]; then
  echo "Still running last poll; exiting now"
  exit 0
fi

init_timestamps

# We want all files newer than the last run of this script AND older than TIME_DELAY minutes ago
# (to prevent processing the file while it is still in-flight).
#
# After we get that set of files, strip out any that match our DATA_FIND_GREP terms
# Then run each file through dirname to get its directory
# then sort and find only the unique directory entries

files=`find ${DATA_DIR} \
       \( -name mibc_products -o -name anpan_products \) -prune -o \
       -newer ${DATA_DIR}/.timestamp.empty \
       -and ! -newer ${DATA_DIR}/.next_timestamp.empty \
       -type f \
       -print 2>/dev/null | grep -v -e '\.log$' -e '\.empty'`

for file in ${files}
do
  if [[ ${file} == "" ]]; then
    continue
  fi
  echo "file: ${file}"
  dir=`dirname "${file}"`
  directories=`echo "${directories}" && echo "${dir}"`
done

directories=`echo "${directories}" | sort | uniq`

# loop through all directories with new content that have come through quiet period...
for dir in ${directories}
do
  echo "testing ${dir}"
  if test "`find ${dir} -maxdepth 1 -type d -mmin +${TIME_DELAY_MIN}`"; then
    call_anadama "${dir}"
  else
    echo "dir: ${dir} is newer than ${TIME_DELAY} seconds - skipping"
  fi
done

finalize_timestamps

