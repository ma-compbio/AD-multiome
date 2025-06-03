#!/bin/bash

usage() {
  echo "Usage: $(basename $0) -b bam_file " 1>&2
  echo "" 1>&2
  echo "  -b bam file" 1>&2
  echo "" 1>&2
  exit 1
}

###############################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":b:o:" opt; do
  case $opt in
  b)
    BAM_FILE=$OPTARG
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

# ----  # of un-/unique/multi/chimeric reads  ----
awk_script_compressSam='
  { if (and($2,0x100)) printf "S" ; else printf "P" }
  { if ($5>=30)        printf "H" ; else printf "L" }
  { if ($0~/\tXA:Z:/)  printf "M" ; else printf "U" }
  { if ($0~/\tSA:Z:/)  printf "C" ; else printf "L" }
  { if (and($2,0x40 )) printf "1" ; else printf "2" }
  { printf substr($3,1,match($3,"_")-1) }
  { print "" }
'
awk_script_stat='
  { species_r = substr($0,5) }
  substr($0,1,1) == "P"                                                                            { n_primary       [species_r] += 1 }
  substr($0,1,1) == "P" && substr($0,2,1) == "H"                                                   { n_30            [species_r] += 1 }
  substr($0,1,1) == "P" &&                          substr($0,3,1) == "M"                          { n_multi         [species_r] += 1 }
  substr($0,1,1) == "P" && substr($0,2,1) == "H" && substr($0,3,1) == "M"                          { n_multi30       [species_r] += 1 }
  substr($0,1,1) == "P" &&                          substr($0,4,1) == "C"                          { n_chimeric      [species_r] += 1 }
  substr($0,1,1) == "P" && substr($0,2,1) == "H" && substr($0,4,1) == "C"                          { n_chimeric30    [species_r] += 1 }
                                                    substr($0,4,1) == "C"                          { n_chimeric_all  [species_r] += 1 }
                           substr($0,2,1) == "H" && substr($0,4,1) == "C"                          { n_chimeric30_all[species_r] += 1 }
  substr($0,1,1) == "P" &&                          substr($0,4,1) == "L" && substr($0,3,1) == "U" { n_unique        [species_r] += 1 }
  substr($0,1,1) == "P" && substr($0,2,1) == "H" && substr($0,4,1) == "L" && substr($0,3,1) == "U" { n_unique30      [species_r] += 1 }
  END {
    print "mapped mm10 R1 = \t"n_primary["1mm10"]"\t"n_30["1mm10"]"\t"n_multi["1mm10"]"\t"n_multi30["1mm10"]"\t"n_chimeric["1mm10"]"\t"n_chimeric30["1mm10"]"\t"n_chimeric_all["1mm10"]"\t"n_chimeric30_all["1mm10"]"\t"n_unique["1mm10"]"\t"n_unique30["1mm10"]
    print "mapped mm10 R2 = \t"n_primary["2mm10"]"\t"n_30["2mm10"]"\t"n_multi["2mm10"]"\t"n_multi30["2mm10"]"\t"n_chimeric["2mm10"]"\t"n_chimeric30["2mm10"]"\t"n_chimeric_all["2mm10"]"\t"n_chimeric30_all["2mm10"]"\t"n_unique["2mm10"]"\t"n_unique30["2mm10"]
    print "mapped hg38 R1 = \t"n_primary["1hg38"]"\t"n_30["1hg38"]"\t"n_multi["1hg38"]"\t"n_multi30["1hg38"]"\t"n_chimeric["1hg38"]"\t"n_chimeric30["1hg38"]"\t"n_chimeric_all["1hg38"]"\t"n_chimeric30_all["1hg38"]"\t"n_unique["1hg38"]"\t"n_unique30["1hg38"]
    print "mapped hg38 R2 = \t"n_primary["2hg38"]"\t"n_30["2hg38"]"\t"n_multi["2hg38"]"\t"n_multi30["2hg38"]"\t"n_chimeric["2hg38"]"\t"n_chimeric30["2hg38"]"\t"n_chimeric_all["2hg38"]"\t"n_chimeric30_all["2hg38"]"\t"n_unique["2hg38"]"\t"n_unique30["2hg38"]
  }
'
samtools view -@8 -F0x4 "$BAM_FILE" | cut -f1-9,12- | sed 's/chr/hg38_chr/g' | \
  awk "$awk_script_compressSam" | tqdm --mininterval 1 | awk "$awk_script_stat" | \
  cut -f2- | tr '\n' '\t' | sed 's/\t$/\n/' | \
  (printf "mapping type - mm10R1 / mm10R2 / hg38R1 / hg38R2: \t" ; cat)
