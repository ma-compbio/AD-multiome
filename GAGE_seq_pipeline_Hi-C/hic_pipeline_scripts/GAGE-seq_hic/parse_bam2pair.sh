#!/bin/bash

usage() {
  echo "Usage: $(basename $0) [-t threads] [-a assembly] -w walk_policy -b bam_file -r rescue_script -c chrom_size -p pair_file" 1>&2
  echo "" 1>&2
  echo "  walk_policy should be one of wowalk, linear, complete, walkonly" 1>&2
  echo "" 1>&2
  exit 1
}

################################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":t:a:w:b:r:c:p:" opt; do
  case $opt in
  t)
    NTHREADS=$OPTARG
    ;;
  a)
    ASSEMBLY=$OPTARG
    ;;
  w)
    WALK_POLICY=$OPTARG
    ;;
  b)
    BAM_FILE=$OPTARG
    ;;
  r)
    RESCUE_SCRIPT=$OPTARG
    ;;
  c)
    CHROM_SIZE=$OPTARG
    ;;  
  p)
    PAIR_FILE=$OPTARG
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

if [[ -z $ASSEMBLY ]]; then
  echo "Assembly not specified, use hg38" >&2
  ASSEMBLY="hg38"
fi

if [[ -z $WALK_POLICY ]]; then
  echo "Walk policy not specified" >&2
  usage
fi

if [[ -z $BAM_FILE ]]; then
  echo "Bam file not specified" >&2
  usage
fi

if [[ -z $RESCUE_SCRIPT ]]; then
  echo "Path to rescue sscript not specified" >&2
  usage
fi

if [[ -z $CHROM_SIZE ]]; then
  echo "Chrom size file not specified" >&2
  usage
fi

if [[ -z $PAIR_FILE ]]; then
  echo "Pair file not specified" >&2
  usage
fi

# ----  pairtools  ----
# run pairtools
# simplify its output
# flip and mark reads
function simplify {
  awk -v OFS='\t' '{
    NF = 11
    p = match($9,"BC:Z:")+5; $11 = substr($9,p,match(substr($9,p),"\x19")-1)
    $9="."; $10="."
    print
  }'
}
function flip {
  awk -v OFS='\t' -F'[\t]' '
    BEGIN { chra2i["!"] = 0 } ARGIND==1 { chra2i[$1] = FNR ; next }
    { if ( chra2i[$2] > chra2i[$4] || ( $2 == $4 && $3 > $5 ) ) {
      print $1 OFS $4 OFS $5 OFS $2 OFS $3 OFS $7 OFS $6 OFS substr($8,2,1) substr($8,1,1) OFS $10 OFS $9 OFS $11 OFS "F"
    } else { print $0 OFS "O" } }
  ' $CHROM_SIZE -
}
function keep_walk {
  awk -v OFS='\t' -F'[\t]' '
    {
      if (qname == $1) { cnt += 1; print last_line }
      else { if ( cnt > 1 ) { print last_line } cnt = 1 }
      qname = $1; last_line = $0
    }
    END { if ( cnt > 1 ) print last_line }
  '
}

if [[ $WALK_POLICY == "wowalk" ]]; then
  function f { pairtools parse --walks-policy mask \
    --no-flip --min-mapq=10 --assembly=$ASSEMBLY --chroms-path=$CHROM_SIZE | \
    grep -v '^#' | simplify | flip | grep '' --line-buffered ; }
elif [[ $WALK_POLICY == "linear" ]]; then
  function f { pairtools parse --walks-policy all \
    --no-flip --min-mapq=10 --assembly=$ASSEMBLY --chroms-path=$CHROM_SIZE | \
    grep -v '^#' | simplify | flip | grep '' --line-buffered ; }
elif [[ $WALK_POLICY == "complete" ]]; then
  function f { pairtools parse --walks-policy all \
    --no-flip --min-mapq=10 --assembly=$ASSEMBLY --chroms-path=$CHROM_SIZE | \
    grep -v '^#' | simplify | python $RESCUE_SCRIPT | flip | grep '' --line-buffered ; }
elif [[ $WALK_POLICY == "walkonly" ]]; then
  function f { pairtools parse --walks-policy all \
    --no-flip --min-mapq=10 --assembly=$ASSEMBLY --chroms-path=$CHROM_SIZE | \
    grep -v '^#' | simplify | keep_walk | python $RESCUE_SCRIPT | flip | grep '' --line-buffered ; }
else echo "Not Implemented"; fi


# 
if [[ $BAM_FILE = "/dev/stdin" ]]; then
  f | gzip -6 > "$PAIR_FILE"
else
  samtools view -h@"$NTHREADS" "$BAM_FILE" | f | gzip -6 > "$PAIR_FILE"
  #echo "samtools view -h@$NTHREADS $BAM_FILE | f | gzip -6 > $PAIR_FILE"
fi
