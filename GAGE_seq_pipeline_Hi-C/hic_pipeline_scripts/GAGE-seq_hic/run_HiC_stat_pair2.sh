#!/bin/bash

library=$1

cd $library

bam="${library}.bam"
pairs="${library}.pairs.gz"
pairs_nodup="${library}_nodup.pairs.gz"
pairs_nodup_wo1k="${library}_nodup_wo1k.pairs.gz"
# ----  pair type  ----
gzip -cd $pairs | cut -f8 | awk '{n[$0]+=1} END{for(k in n) print n[k]"\t"k}' | sort -k2,2n | (tee \
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
  <(gzip -cd $pairs       | cut -f2-5,8 | grep -P 'UU|UR|RU' | cut -f1-4 | tr _ '\t' | f) \
  <(gzip -cd $pairs_nodup | cut -f2-5                                    | tr _ '\t' | f) | \
  paste -sd'\t' | (printf "Interaction type - all / nodup: \t"; cat)
# ----  duplicate enrichment  ----
gzip -cd $pairs_nodup      | cut -f10 | awk '{n[$0] += 1} END {for (i in n) print i "\t" n[i]}' | \
  sort -k1,1n > "stat_dupHisto_${library}_${species}_${genome}.out"
gzip -cd $pairs_nodup_wo1k | cut -f10 | awk '{n[$0] += 1} END {for (i in n) print i "\t" n[i]}' | \
  sort -k1,1n > "stat_dupHisto_wo1k_${library}_${species}_${genome}.out"
