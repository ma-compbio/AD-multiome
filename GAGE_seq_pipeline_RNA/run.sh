#!/bin/bash


## single library
for config_file in `ls pipeline_config/GAGE-seq_*.config | grep -v merge`
do
    python pipeline_script/build_sbatch_alignment.PSC.py --config ${config_file} --out_folder /ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/pipeline_job
done

## merge library
for config_file in `ls pipeline_config/*merge*.config` 
do
    python pipeline_script/build_sbatch_merge_bam.PSC.py --config ${config_file} --out_folder /ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/pipeline_job
done 

## dedup library
for config_file in `ls pipeline_config/*merge*.config`
do
    python pipeline_script/build_sbatch_dedup.PSC.py --config ${config_file} --out_folder /ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/pipeline_job
done

## first split the bam file into chromosome then dedup library

## feature count
for config_file in `ls pipeline_config/*merge*.config`
do
    python pipeline_script/build_sbatch_feature_count.PSC.py --config ${config_file} --out_folder /ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/pipeline_job
done
