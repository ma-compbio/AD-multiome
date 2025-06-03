# Script for integrative analysis of multi-omics data

R notebook file `R_scripot/analysis_script.Rmd` contains the R script for integrative analysis of multi-omics data. This script is designed to perform the following tasks:
* Plotting the sample's clinical metadata
* Plotting the cell-type composition of the samples
* Integration with PFC427 scRNA-seq (Mathys et al. Cell 2023)
* Integration with snATAC-seq (Xiong et al. Cell 2023)
* Integration with Xenium data
* APA analysis of chromatin loops
* Saddle plot analysis of chromatin compartment
* Long-/short-range interaction analysis
* A/B compartment mingling analysis

Due to privacy concerns, some scripts are not directly runable. Instead, they are provided as templates for users to adapt to their own data. The scripts also requires specific result files from previous steps in the analysis pipeline, which are not included in this repository. For example, the loop APA analysis requires the output from the 3D_Genome_feature_calling step, and the saddle plot analysis requires the output from the 3D_Genome_feature_calling step as well.
