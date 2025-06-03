#!/bin/bash


#SBATCH -p RM-shared
#SBATCH -t 2-00:00:00
#SBATCH --ntasks-per-node=8
##SBATCH --mem=128Gb
#SBATCH --job-name bwa_index
#SBATCH --output log_bwa_index-%J.txt

# activate conda
source activate ad_project

# goto location
cd /jet/home/yzhang38/bighive/shared/AD_project/genome

# script starting from here
#bwa index hg38.fa bwa_hg38_index

samtools faidx GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

