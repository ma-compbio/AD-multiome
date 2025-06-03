#!/bin/bash

#SBATCH -p ma1
#SBATCH -w oven-0-25
#SBATCH -t 2-00:00:00
#SBATCH --ntasks-per-node=64
#SBATCH --mem=80Gb
#SBATCH --job-name STAR_build
#SBATCH --output log-STAR_build_%J.txt

# activate conda env
source ~/mambaforge/etc/profile.d/conda.sh
#module load cuda-12.2
conda activate ad_rna

# go to working folder
cd /work/magroup/shared/4DN_project/ad/pipeline_RNA/annotation

# run scripts
genome_fasta=/work/magroup/shared/4DN_project/ad/pipeline_RNA/genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
out_STAR_folder=/work/magroup/shared/4DN_project/ad/pipeline_RNA/genome/STAR_anno
GENCODE_GTF=/work/magroup/shared/4DN_project/ad/pipeline_RNA/annotation/gencode.v47.primary_assembly.annotation.gtf

~/Software/STAR_2.7.11a/STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ${out_STAR_folder} --genomeFastaFiles ${genome_fasta} --sjdbGTFfile ${GENCODE_GTF} --sjdbOverhang 100



