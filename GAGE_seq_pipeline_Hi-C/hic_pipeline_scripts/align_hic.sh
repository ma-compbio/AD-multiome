#!/bin/bash

# ---- path ----
path2rawdata="/work/magroup/xinyuelu/AD-proj/data"
path2results="/work/magroup/xinyuelu/AD-proj/results"
path2scripts_alignment="/work/magroup/xinyuelu/AD-proj/scripts_alignment"
path2scripts_results="/work/magroup/xinyuelu/AD-proj/scripts_results"
path2genome="/work/magroup/xinyuelu/genome/hg38_mm10"

# # -------- prepare file structure --------

cd $path2results

mkdir -p alignment_hic alignment_hic_complete alignment_hic_linear alignment_hic_wowalk alignment_hic_walkonly
mkdir -p alignment_rna
mkdir -p results_hic results_hic_complete results_hic_linear results_hic_wowalk results_hic_walkonly
mkdir -p results_rna

for lib in $(cut -f1 "../library_list.txt")
do  
  for dir in alignment_hic* results_hic*
  do
      mkdir -p "$dir/$lib"
  done
done

for lib in $(cut -f2 "../library_list.txt")
do  
  for dir in alignment_rna* results_rna*
  do
      mkdir -p "$dir/$lib"
  done
done

# -------- align HiC --------
cd alignment_hic
species="all"
genome="hg38_mm10"
for lib in $(cut -f1 "../raw/library_list.txt" | sed -n '1,5p')
do
  pwd
  echo "$lib"
  cd $lib
  date
  OPENBLAS_NUM_THREADS=1 $path2scripts_alignment/run_HiC_align.sh "$lib" "$species" "$genome" 1> "HiC_align_${lib}.out" 2> "HiC_align_${lib}.err"
  date
  cd ..
  sleep 5
done
cd ..

# -------- HiC merge pairs --------
# for walk_policy in complete linear wowalk walkonly
for walk_policy in complete
do
  echo "$walk_policy"
  p=16

  cd "alignment_hic_${walk_policy}"
  pwd
  ls
  date
  for lib in $(cut -f1 "/work/magroup/xinyuelu/AD-proj/library_list.txt")
  do
  {
    pigz -cdp8      "202406/${lib}/${lib}_all_hg38_mm10.pairs.gz"
    pigz -cdp8      "202408/${lib}/${lib}_all_hg38_mm10.pairs.gz"
  } | pigz -6p$p >  "combined/${lib}/${lib}_all_hg38_mm10.pairs.gz"
  done
  
  cd ..
done

# -------- HiC dedup --------
# for walk_policy in complete linear wowalk walkonly
for walk_policy in complete
do
  cd "alignment_hic_${walk_policy}/combined"
  date
  for lib in $(cut -f1 "/work/magroup/xinyuelu/AD-proj/library_list.txt")
  do
    cd ${lib}
    $path2scripts_alignment/run_HiC_dedup3.sh ${lib} all hg38_mm10
    cd ..
  done
  date
  cd ../..
done


