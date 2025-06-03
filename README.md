

# Single-cell multiomics connects 3D genome and transcriptome alterations in Alzheimer’s disease

This repository contains the code used for the data processing and analysis.

The code is organized into the following directories:
- `GAGE_seq_pipeline_Hi-C`: Contains scripts for processing the raw GAGE-seq Hi-C, including demultiplexing, mapping, merging technical replicates, and deduplication
- `GAGE_seq_pipeline_RNA`: Contains scripts for processing the raw GAGE-seq RNA, including demultiplexing, mapping, deduplication, and quantification
- `analysis_scripts`: Contains scripts for the analysis of the processed data, including integrative analysis, gene expression analysis, 3D genome analysis, and visualization

--------------------------

# How to use the code

## GAGE-seq pipeline
The GAGE-seq pipeline scripts are designed to be run on a high-performance computing cluster. The scripts are written in Python and Bash. The required software and dependencies are listed in the `environment.yml` file, which can be used to create a conda environment.
The source of genome reference files and annotation files are provided in the `annotation` and `genome` folder with the corresponding scripts to download and prepare annoation files.
Note that there are some hard-coded paths in the scripts, so you may need to modify them to fit your local setup.
Script `run.sh` is the main entry point to run the GAGE-seq pipeline. It will run all the necessary steps for processing the raw data and generating the final output files.

## Analysis scripts
The scripts in the `analysis_scripts` directory are organized by analysis modules. The scripts are written in Python or R. Due to the privacy of the data, the scripts are not runable without proper sample annotation files. However, the scripts are well-documented and can be adapted to your own data. The main analysis modules are as follows:
- `3D_genome_feature_callin3D_genome_feature_callingg`: This module contains scripts for calling 3D genome features from the processed Hi-C data.
- `Gene_expression_analysis`: This module contains scripts for analyzing gene expression from the processed RNA data.
- `Integrative_analysis`: This module contains scripts for integrating the Hi-C, RNA, and snATAC data to identify correlations between 3D genome features and gene expression.
- `Hicformer_analysis`: This module contains scripts for the Hicformer model, including training and analysis.
- `Visualization`: This module contains scripts for creating visualizations of the results, including plots and figures for the manuscript.

