#!/bin/bash


usage() {
  echo "Usage: $(basename $0) -b bam_file -p pair_file -d pair_nodup -w pair_nodup_wo1k -o output_file_prefix" 1>&2
  echo "" 1>&2
  echo "  -b bam file" 1>&2
  echo "  -p pair file" 1>&2
  echo "  -d pair file without duplicates" 1>&2
  echo "  -w pair file without duplicates and with distance > 1k" 1>&2
  echo "" 1>&2
  exit 1
}

###############################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":b:p:d:w:o:" opt; do
  case $opt in
  b)
    BAM_FILE=$OPTARG
    ;;
  p)
    PAIR_FILE=$OPTARG
    ;;
  d)
    PAIR_FILE_NODUP=$OPTARG
    ;;
  w)
    PAIR_FILE_WO1K=$OPTARG
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

if [[ -z $BAM_FILE ]]; then
  echo "Bam file not specified" >&2
  usage
fi

if [[ -z $PAIR_FILE ]]; then
  echo "Pair file not specified" >&2
  usage
fi

if [[ -z $PAIR_FILE_NODUP ]]; then
  echo "Pair file without duplicates not specified" >&2
  usage
fi

if [[ -z $PAIR_FILE_WO1K ]]; then
  echo "Pair file without duplicates and with distance > 1k not specified" >&2
  usage
fi

if [[ -z $OUTPUT_FILE_PREFIX ]]; then
  echo "Output file not specified" >&2
  usage
fi

# ----  pair type  ----
gzip -cd $PAIR_FILE | sed 's/chr/hg38_chr/g' | cut -f8 | awk '{n[$0]+=1} END{for(k in n) print n[k]"\t"k}' | sort -k2,2n | (tee \
    >(grep 'N'               | cut -f1 | paste -sd'+' | bc | (printf "1\t"; cat) | grep '' --line-buffered) \
    >(grep 'M' | grep -v 'N' | cut -f1 | paste -sd'+' | bc | (printf "2\t"; cat) | grep '' --line-buffered) \
    >(grep 'WW'              | cut -f1 | paste -sd'+' | bc | (printf "3\t"; cat) | grep '' --line-buffered) \
    >(grep 'UU\|UR\|RU'      | cut -f1 | paste -sd'+' | bc | (printf "4\t"; cat) | grep '' --line-buffered) \
  > /dev/null) | sort -k1,1n | cut -f2 | tr '\n' '\t' | \
  (printf "pair type - one or both Null / one or both Multiple mapping / Walk / Unique/Rescue: \t"; cat)
# ----  read pair distances  ----
function f {
  awk -v OFS='\t' '{if ($1==$4) print $1 OFS $2 OFS $3 OFS $5 OFS $6 ; else print -2}' | \
  awk -v OFS='\t' '{if (NF>1) if ($2==$4) print $3 OFS $5 OFS $1 ; else print -1 OFS $1 ; else print $0}' | \
  awk -v OFS='\t' '{if (NF>2) print $2-$1 OFS $3 ; else print $0}' | \
  awk '
  BEGIN { n_human_0=0; n_human_1=0; n_human_2=0; n_human_3=0; n_human_4=0; n_human_5=0; n_human_inter=0;
    n_inter_spec=0; n_mouse_0=0; n_mouse_1=0; n_mouse_2=0; n_mouse_3=0; n_mouse_4=0; n_mouse_5=0; n_mouse_inter=0 }
  $1 < 0 {
    if ($1==-2) n_inter_spec += 1;
    else if ( $2 == "hg38" ) n_human_inter += 1; else n_mouse_inter += 1 ;
    next }
  $1 < 1000    { if ( $2 == "hg38" ) n_human_0 += 1; else n_mouse_0 += 1 ; next }
  $1 < 20000   { if ( $2 == "hg38" ) n_human_1 += 1; else n_mouse_1 += 1 ; next }
  $1 < 50000   { if ( $2 == "hg38" ) n_human_2 += 1; else n_mouse_2 += 1 ; next }
  $1 < 100000  { if ( $2 == "hg38" ) n_human_3 += 1; else n_mouse_3 += 1 ; next }
  $1 < 1000000 { if ( $2 == "hg38" ) n_human_4 += 1; else n_mouse_4 += 1 ; next }
               { if ( $2 == "hg38" ) n_human_5 += 1; else n_mouse_5 += 1 ; next }
  END { print n_human_0 "\t" n_human_1 "\t" n_human_2 "\t" n_human_3 "\t" n_human_4 "\t" n_human_5 "\t" n_human_inter \
    "\t" n_inter_spec "\t" n_mouse_0 "\t" n_mouse_1 "\t" n_mouse_2 "\t" n_mouse_3 "\t" n_mouse_4 "\t" n_mouse_5 "\t" n_mouse_inter }
'
}
# cols 2-5,8: mm10_chr7	12941065	mm10_chr7	24611449  UU
# after `tr`: mm10	chr7	12941065	mm10	chr7	24611449
cat \
  <(gzip -cd $PAIR_FILE       | sed 's/chr/hg38_chr/g' | cut -f2-5,8 | grep -P 'UU|UR|RU' | cut -f1-4 | tr _ '\t' | f) \
  <(gzip -cd $PAIR_FILE_NODUP | sed 's/chr/hg38_chr/g' |cut -f2-5                                    | tr _ '\t' | f) | \
  paste -sd'\t' | (printf "Interaction type - all / nodup: \t"; cat)
# ----  duplicate enrichment  ----
gzip -cd $PAIR_FILE_NODUP      | cut -f10 | awk '{n[$0] += 1} END {for (i in n) print i "\t" n[i]}' | \
  sort -k1,1n > "stat_dupHisto_${OUTPUT_FILE_PREFIX}.out"
gzip -cd $PAIR_FILE_WO1K | cut -f10 | awk '{n[$0] += 1} END {for (i in n) print i "\t" n[i]}' | \
  sort -k1,1n > "stat_dupHisto_wo1k_${OUTPUT_FILE_PREFIX}.out"
