#!/bin/bash

usage() {
  echo "Usage: $(basename $0) [-t threads] -i genome_index_folder -f first_fastq -s second_fastq -l library_name -o output_folder" 1>&2
  echo "" 1>&2
  echo "  bwaIndex can be a path to a bwa index prefix or a tarball of an bwa index" 1>&2
  echo "  genome_index_folder is the folder where the STAR index files are stored" 1>&2
  echo "  first_fastq is the path to the first fastq file" 1>&2
  echo "  second_fastq is the path to the second fastq file" 1>&2
  echo "  library_name is the name of the library" 1>&2
  echo "  output_folder is the folder where the output files will be stored" 1>&2   
  echo "" 1>&2
  exit 1
}

################################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":t:i:f:s:l:o:" opt; do
  case $opt in
  t)
   NTHREADS=$OPTARG
   ;;
  i)
   INDEXFOLDER=$OPTARG
   ;;
  f)
    FASTQ_R1=$OPTARG
    ;;
  s)
    FASTQ_R2=$OPTARG
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

if [[ -z $INDEXFOLDER ]]; then
  echo "must define a path to a STAR genome annotation folder -i"
fi

if [[ -z $FASTQ_R1 ]]; then
  echo "must define a path to the first fastq file -f"
fi

if [[ ! -f $FASTQ_R1 ]]; then
  echo "fastq file $FASTQ_R1 does not exist"
fi

if [[ -z $FASTQ_R2 ]]; then
  echo "must define a path to the second fastq file -s"
fi

if [[ ! -f $FASTQ_R2 ]]; then
  echo "fastq file $FASTQ_R2 does not exist"
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

adaptor_path="/ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/annotation/rna_adaptor"
adaptor_file_bc1="$adaptor_path/adaptor_cDNA_bc1"
adaptor_file_bc2="$adaptor_path/adaptor_cDNA_bc2"
adaptor_file_l0="$adaptor_path/adaptor_cDNA_l0"
adaptor_file_l1="$adaptor_path/adaptor_cDNA_l1"
adaptor_file_l2="$adaptor_path/adaptor_cDNA_l2"
# simplify read names
function trimSeqName { awk -v OFS='\t' 'NR==1 { idx=1; for(c=0;c<4;idx++) if(substr($0,idx,1)==":") c++; } NR%4==1 { $0 = "@"substr($0,idx) } 1' ; }
barcode2keep=""
# ----  for STAR alignment  ----
genomeDir=${INDEXFOLDER}
genomeLoad="LoadAndKeep"
STARtool="/jet/home/yzhang38/Software/STAR_2.7.11a/STAR"
# samtools only work with 1.9 or 1.8
samtools="/jet/home/yzhang38/Software/samtools-1.9/bin/samtools"

################################
### Demultiplexed reads ########
################################
cd ${OUTDIR}

script_demultiplex="/ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/pipeline_script/demultiplex.py"
OUT_FASTQ_R1="${OUTDIR}/${LIBRARY_ID}_demulti_R1.fastq"
OUT_FASTQ_R2="${OUTDIR}/${LIBRARY_ID}_demulti_R2.fastq"
python ${script_demultiplex} \
  --barcode=$adaptor_file_bc2 --pos1="slice(10, 18)" --pos2="[]" --mode=1 \
  --barcode=$adaptor_file_l2  --pos1="slice(18, 48)" --pos2="[]" --mode=6 \
  --barcode=$adaptor_file_bc1 --pos1="slice(48, 56)" --pos2="[]" --mode=1 \
  --barcode=$adaptor_file_l1  --pos1="slice(56, 71)" --pos2="[]" --mode=6 \
  --infile1  <(gzip -cd "${FASTQ_R1}" | trimSeqName) \
  --infile2  <(gzip -cd "${FASTQ_R2}" | trimSeqName) \
  --outfile1 >(grep -a . > "${OUT_FASTQ_R1}") \
  --outfile2 >(grep -a . > "${OUT_FASTQ_R2}") \
  --barcode2keep="$barcode2keep" 

echo "Demultiplexing finished"
echo "Demultiplexed reads are stored in ${OUT_FASTQ_R1} and ${OUT_FASTQ_R2}"
echo ""

################################
#####     STAR alignment #######
################################
echo "Start alignment"

out_bam="${OUTDIR}/${LIBRARY_ID}_human_hg38.bam"
$STARtool --runThreadN=8 --genomeLoad $genomeLoad --genomeDir="$genomeDir" --outSAMtype SAM --outStd SAM \
  --outSAMunmapped Within --outFileNamePrefix="${library}_${species}_${genome}_" --readFilesIn=<(
    paste -d'\0' "${OUT_FASTQ_R2}" <(awk 'NR%4==2{printf "_" substr($0,1,10) "\n\n\n\n"}' "${OUT_FASTQ_R1}")
) | awk -v OFS='\t' '{split($1,bcumi,"_"); $1=bcumi[1]; print $0 "\tBC:Z:" bcumi[2] "," bcumi[4] "\tRX:Z:" bcumi[6]}' | \
  ${samtools} view -hb@8 > "${out_bam}"

echo "Alignment finished"