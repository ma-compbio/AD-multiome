#!/bin/bash


usage() {
  echo "Usage: $(basename $0) [-t threads] [-s dedup_Cpp_script] -p pair_file -o output_file_prefix" 1>&2
  echo "" 1>&2
  echo "  -s should be the path to dedup_Hi-C_cpp binary program" 1>&2
  echo "" 1>&2
  exit 1
}

################################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":t:s:p:o:" opt; do
  case $opt in
  t)
    NTHREADS=$OPTARG
    ;;
  s)
    DEDUP_CPP=$OPTARG
    ;;
  p)
    PAIR_FILE=$OPTARG
    ;;
  o)
    OUTPUT_FILE_PREFIX=$OPTARG
    ;;
  \?)
   echo "Invalid option: -$OPTARG" >&2
   usage
   ;;
  [?])
   usage
   ;;
  :)
   echo "Option -$OPTARG requires an argument." >&2
   echo "" >&2
   usage
   ;;
  esac
done

if [[ -z $NTHREADS ]]; then
  echo "Thread number not specified, use 1" >&2
  NTHREADS=1
fi

if [[ -z $DEDUP_CPP ]]; then
  echo "dedup_Hi-C_cpp not specified" >&2
  usage
fi

if [[ -z $PAIR_FILE ]]; then
  echo "Pair file not specified" >&2
  usage
fi

if [[ -z $OUTPUT_FILE_PREFIX ]]; then
  echo "Output file not specified" >&2
  usage
fi  

OUTPUT_nodup="${OUTPUT_FILE_PREFIX}_nodup.pairs.gz"
OUTPUT_nodup_wo1k="${OUTPUT_FILE_PREFIX}_nodup_wo1k.pairs.gz"

pigz -cdp"$NTHREADS" "$PAIR_FILE" | grep -v 'chrUn\|random\|EBV' | grep -P 'UU|UR|RU' | cut -f2-7,11- | "$DEDUP_CPP" 500 500 9216 25 | \
  grep -vP '\tDD\t' | ( tee >(pigz -6p2 > "${OUTPUT_nodup}") ) | \
  awk '$2!=$4 || $5-$3>=1000' | pigz -6p2 > "${OUTPUT_nodup_wo1k}"

