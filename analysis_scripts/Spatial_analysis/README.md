
# Scripts for Analyzing Xenium Data

## Overview
This directory contains scripts for analyzing Xenium spatial transcriptomics data.
Due to the privary requirement of the data, the scripts are not directly runnable without sample id information. These scripts are intended to be used as a reference for users who want to implement their own Xenium data analysis pipeline.

## Script Descriptions
1. **AD_analysis_Final_Celltypeproportion.R**: Extracts the proportions of major and sub cell types across samples..
2. **AD_analysis_Final_Expressionchange_compareXeniumGAGEseq.R**: Compares differential gene expression (log2 fold change between AD and non-AD) between Xenium and GAGE-seq data across representative major cell types.
3. **AD_analysis_Final_Neighboringcelltypeproportion.R**: Analyzes changes in neighboring cell type proportions between AD and non-AD across different brain regions.
4. **AD_analysis_Final_BivariateMoranI.py**: Calculates bivariate Moran’s I statistics between gene expression levels and imputed A/B compartment scores.
