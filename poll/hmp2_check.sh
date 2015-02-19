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

# Source project specific configuration parameters
if [ -f ${ANADAMA_DIR}/bin/activate ]; then
  . ${ANADAMA_DIR}/bin/activate
else
  echo "${ANADAMA_DIR}/bin/activate file not found"
  exit -1
fi

function call_anadama {

  echo "mibc_run.sh $1 lsf 10"
  #mibc_run.sh $1 lsf 10 &
  #mibc_build dag --project="$1" -a --tmpfiledir=${ANADAMA_DIR}/tmp/tm | mibc_tm -t lsf -g 10

}

# We want all files newer than the last run of this script AND older than TIME_DELAY minutes ago
# (to prevent processing the file while it is still in-flight).
#
# After we get that set of files, strip out any that match our DATA_FIND_GREP terms
# Then run each file through dirname to get its directory
# then sort and find only the unique directory entries

touch -d `date --date="-${TIME_DELAY_SEC} seconds" +%H:%M:%S` /tmp/.next_timestamp.empty

files=`find ${DATA_DIR} \
       -name mibc_products -prune -o \
       -newer ${DATA_DIR}/.timestamp.empty \
       -and ! -newer /tmp/.next_timestamp.empty \
       -type f \
       -print`
      #grep -v -e '/tmp/' -e '\.log$' -e '\.txt$' -e '\.empty'`

for file in $files
do
  echo "file: $file"
done

echo "number of files in window:"
echo ${files} | wc -l

for file in ${files}
do
  if [[ ${file} == "" ]]; then
    continue
  fi
  dir=`dirname "${file}"`
  echo "dir: ${dir}"
  directories=`echo "${directories}" && echo "${dir}"`
done

for dir in "${directories}"
do
  echo "directories: $dir"
done

echo "number of unique directories in window:"
directories=`echo "${directories}" | sort | uniq`
echo ${directories} | wc -l

for dir in ${directories}
do
  echo "testing ${dir}"
  if test "`find ${dir} -maxdepth 1 -type d -mmin +${TIME_DELAY_MIN}`"; then
    call_anadama "${dir}"
  else
    echo "dir: ${dir} is newer than ${TIME_DELAY_SEC} seconds - skipping"
  fi
done
