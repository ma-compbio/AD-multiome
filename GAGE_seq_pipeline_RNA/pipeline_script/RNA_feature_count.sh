#!/bin/bash


usage() {
  echo "Usage: $(basename $0) [-t threads] -b nodup_bam -g gene_gtf -l library_name -o output_folder" 1>&2
  echo "" 1>&2
  echo "  gene_gtf is the gene annotation file " 1>&2
  echo "  library_name is the name of the library" 1>&2
  echo "  output_folder is the folder where the output files will be stored" 1>&2   
  echo "" 1>&2
  exit 1
}

################################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":t:b:g:l:o:" opt; do
  case $opt in
  t)
   NTHREADS=$OPTARG
   ;;
  b)
    NODUP_BAM=$OPTARG 
    ;;
  g)
    GENE_GTF=$OPTARG
   ;;
  l)
    LIBRARY_ID=$OPTARG
    ;;
  o)
    OUTDIR=$OPTARG
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

################################
### Check COMMAND LINE ARGS ####
################################

shift $((OPTIND-1))

if [ -z $NTHREADS ]; then
  if [ -z $NUMTHREADS ]; then
    NTHREADS=16
  else
    NTHREADS=$NUMTHREADS
  fi
fi

if [ -z $NODUP_BAM ]; then
  echo "must define a path to no duplicate bam file -b"
fi

if [[ ! -f $NODUP_BAM ]]; then
  echo "nodup bam file $NODUP_BAM does not exist"
fi


if [[ -z $GENE_GTF ]]; then
  echo "must define a path to the gene gtf annotation file -g"
fi

if [[ ! -f $GENE_GTF ]]; then
  echo "gene gtf file $GENE_GTF does not exist"
fi

if [[ -z $LIBRARY_ID ]]; then
  echo "must define a library name -l"
fi

if [[ -z $OUTDIR ]]; then
  echo "must define an output folder -o"
fi

################################
### Hard code some shared var ##
################################
# samtools only work with 1.9 or 1.8
samtools="/jet/home/yzhang38/Software/samtools-1.21/bin/bin/samtools"
featurecounts="/jet/home/yzhang38/Software/subread-2.0.8-Linux-x86_64/bin/featureCounts"

################################
### gene quantification ########
################################
out_feature_count="${OUTDIR}/${LIBRARY_ID}_human_hg38_fc.out"
${featurecounts} -T 8 -F GTF -a $GENE_GTF -t gene -g gene_id -s 1 --minOverlap 1 --fracOverlap .5 \
  --ignoreDup -o "${out_feature_count}" -R BAM "$NODUP_BAM"
mv $NODUP_BAM.featureCounts.bam "${OUTDIR}/${LIBRARY_ID}_human_hg38.featureCounts.bam"
bamfc="${OUTDIR}/${LIBRARY_ID}_human_hg38.featureCounts.bam"

################################
### rna stat ###################
################################
# ---- # overall mapping stat
echo "flagstat:"
$samtools flagstat "$NODUP_BAM" #| tr '\n' '\t'
echo ""

# ---- # of reads per species ----
function f { awk '
  { species = "hg38"; isDup = and($1,0x400)!=0;
    for(i=3;i<=NF;i++) if(substr($i,1,2)=="NH") { isUnique = $i=="NH:i:1" ; break }
    n_mapped[species isUnique] += 1;
    if (!isDup) n_nodup[species isUnique] += 1; }
  END { print n_mapped["hg381"] "\t" n_mapped["hg380"] "\t" n_nodup["hg381"]  "\t" n_nodup["hg380"] }
' | grep '' --line-buffered ; }
$samtools view --no-PG -@4 -F0x104 "$NODUP_BAM" | cut -f2,3,12- | f | \
 ( printf "Mapping per species - mapped unique | mapped multi | nodup unique | nodup multi: \n"; cat )

# ---- per cell stat ----
echo "per cell stat"
out_cell_stat="${OUTDIR}/statPerCell_mapping_${LIBRARY_ID}_human_hg38.out"
$samtools view --no-PG -@8 -F0x504 -d "NH:1" "$NODUP_BAM" | cut -f3,16 | awk '
  BEGIN { chr[0] = "A"; chr[1] = "B"; chr[2] = "C"; chr[3] = "D"; chr[4] = "E"; chr[5] = "F"; chr[6] = "G"; chr[7] = "H";
   for (i=0;i<8;i++) for (j=1;j<=12;j++) for (k=0;k<8;k++) for (l=1;l<=12;l++) { w = chr[i] j "," chr[k] l ;
    n_human[w] = 0} }
  { wellID = substr($2,6); n_human[wellID] += 1  }
  END { for ( w in n_human ) print w "\t" n_human[w]}
' > $out_cell_stat 
echo ""

# --- report gene count 
out_RNAcount="${OUTDIR}/${LIBRARY_ID}_human_hg38_counts.txt"
echo "gene count: see $out_RNAcount"
$samtools view --no-PG -@4 -F0x504 -d "NH:1" -u "$bamfc" | $samtools view --no-PG -@4 -d "XS:Assigned" | cut -f3,16,21 | tr ':' '\t' | cut -f1,4,7 | sed 's/\.[[:digit:]]\+$//' > $out_RNAcount
echo ""

# --- mapped to genes ----
echo "uniquely mapped to genes | mapped to multiple genes | unassigned | unassigned due to low overlapping"
count=$($samtools view --no-PG -@4 -F0x504 -d "NH:1" -u "$bamfc" | $samtools view --no-PG -@4 -d "XS:Assigned" | wc -l)
echo $count
count=$($samtools view --no-PG -@4 -F0x504 -d "NH:1" -u "$bamfc" | $samtools view --no-PG -@4 -d "XS:Unassigned_Ambiguity" | wc -l)
echo $count
count=$($samtools view --no-PG -@4 -F0x504 -d "NH:1" -u "$bamfc" | $samtools view --no-PG -@4 -d "XS:Unassigned_NoFeatures" | wc -l)
echo $count
count=$($samtools view --no-PG -@4 -F0x504 -d "NH:1" -u "$bamfc" | $samtools view --no-PG -@4 -d "XS:Unassigned_Overlapping_Length" | wc -l)
echo $count
