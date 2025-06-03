#!/bin/bash

############################
## test on demo config file
############################
python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/config_template.yaml 


############################
## run on real data for different batches for demultiplex, alignment, and bam2pair
############################
# 2024 06 batch
library_list=(1 2 3 4 5 6 7 8)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_06_24_lib-${library}.yaml --job_list demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_06_24_lib-${library}.log
done
# 2024 08 batch
library_list=(1 2 3 4 5 6 7 8)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_08_24_lib-${library}.yaml --job_list demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_08_24_lib-${library}.log
done
# 2024 09 batch A
library_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_09_24_A_lib-${library}.yaml --job_list demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_09_24_A_lib-${library}.log
done
# 2024 09 batch B
library_list=(1 2 3 4 5 6 7 8)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_09_24_B_lib-${library}.yaml --job_list demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_09_24_B_lib-${library}.log
done
# 2024 11
library_list=(9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_11_24_lib-${library}.yaml --job_list fastqc demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_11_24_lib-${library}.log
done
# 2024 12
library_list=(15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_12_24_lib-${library}.yaml --job_list fastqc demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_12_24_lib-${library}.log
do
library_list=(40 41 42 43 44 45)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_12_24_lib-${library}.yaml --job_list fastqc demultiplex alignment bam2pair 2>&1 | tee hic_pipeline_log/GAGE-seq_12_24_lib-${library}.log
done

############################
## merge pairs then dedup and perform on merged pairs
############################
library_list=(15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45)
for library in ${library_list[@]}
do
    python hic_pipeline_scripts/GAGE-seq_pipeline_v1.0.py --config hic_pipeline_config/GAGE-seq_merge_lib-${library}.yaml --job_list merge_pair dedup qc_pair 2>&1 | tee hic_pipeline_log/GAGE-seq_merge_lib-${library}.log
done



