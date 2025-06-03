
# Scripts for processing GAGE-seq Hi-C data and 3D genome feature calling

## Overview
This directory contains scripts for processing GAGE-seq Hi-C data and calling 3D genome features. The scripts are designed to work with the output from the GAGE-seq pipeline and generate various 3D genome features such as mcool file, long-/short-range summary per cell, etc.
Due to the privary requirement of the data, the scripts are not directly runnable without sample id information. These scripts are intended to be used as a reference for users who want to implement their own 3D genome feature calling pipeline.

## Overview of Pipeline
1. **Extract read pairs**: Extract read pairs from the GAGE-seq Hi-C data.
2. **Generate mcool file**: Convert the extracted read pairs into a multi-resolution cool file.
3. **Call 3D genome features**: Call various 3D genome features such as TADs, loops, and compartments from the mcool file.

### Extract read pairs
```bash
python src/cell_extract_hic_pair.py --pair_file hic_pipeline_results/hPFC-scHiC_merge_lib-43/complete/hPFC-scHiC_merge_lib-43_nodup_wo1k.pairs.gz --library_id hPFC-lib-43 --sample_anno annotation/sample_info.v2.tsv --output_folder processed_data/202501/raw_pair
```
Here, `--pair_file` is the path to the GAGE-seq Hi-C read pairs file, `--library_id` is the identifier for the library, `--sample_anno` is the path to the sample annotation file, and `--output_folder` is the directory where the processed read pairs will be saved. The script will create subdirectories for each sample (sample_id) and save the processed read pairs in the corresponding subdirectory. Each file will be named as cell_id.pairs.gz, where cell_id is the identifier for the cell.

### Generate mcool file
1. **Generate scool file from read pairs**:
```bash
python src/cell_scool_create.py --thread 80 --chrom_size annotation/hg38.chrom.sizes --resolution 10000 --input_pair_folder processed_data/202501/filtered_pair/${sample_id} --output_prefix processed_data/202501/filtered_pair_scool/${sample_id} --blacklist annotation/hg38_blacklist.bed
```

Here, `--thread` specifies the number of threads to use, `--chrom_size` is the path to the chromosome size file, `--resolution` is the resolution for the scool file, `--input_pair_folder` is the directory containing the processed read pairs, `--output_prefix` is the prefix for the output mcool file, and `--blacklist` is the path to the blacklist file. Blacklist region file is created using the following command:

```bash
cat ENCODE_DAC_blacklist.bed hg38_Gap.bed | sort -k1,1 -k2,2n | bedtools merge | grep -v 'chrUn\|alt\|random\|fix' > hg38_blacklist.bed
```
, where `ENCODE_DAC_blacklist.bed` is the ENCODE blacklist file and `hg38_Gap.bed` is the file containing gaps in the hg38 genome assembly.

2. **Convert scool file to mcool file**:
```bash
python src/cell_pseudo_bulk.scool.py --scool_file results/202501/filtered_pair_scool/${sample_id}.scool --cell_anno results/202501/cell_meta_table.tsv --genome_size annotation/hg38M.chrom.sizes --output_folder results/202501/cool_file/ --thread 12
```
Here, `--scool_file` is the path to the scool file generated in the previous step, `--cell_anno` is the path to the cell annotation file, `--genome_size` is the path to the genome size file, `--output_folder` is the directory where the mcool file will be saved, and `--thread` specifies the number of threads to use.

3. **Create psuedo bulk mcool file**:
```bash
python src/cell_pseudo_bulk.merge.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --thread 12
python src/cell_pseudo_bulk.merge.py --folder results/202501/cool_file/sub --cell_anno results/202501/cell_meta_table.tsv --cell_label sub --thread 12
## zoomify and balance the contact map
python src/cell_pseudo_bulk.zoomify.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --blacklist annotation/hg38_blacklist.bed --thread 12 --script_out log/pseudo_bulk_major_zoomify.202501.sh
python src/cell_pseudo_bulk.zoomify.py --folder results/202501/cool_file/sub --cell_anno results/202501/cell_meta_table.tsv --cell_label sub --blacklist annotation/hg38_blacklist.bed --thread 12 --script_out log/pseudo_bulk_sub_zoomify.202501.sh
## create the expected contact map
python src/cell_pseudo_bulk.expected_cis.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --script_out log/pseudo_bulk_major_expected_cis.202501.sh
python src/cell_pseudo_bulk.expected_cis.py --folder results/202501/cool_file/sub --cell_anno results/202501/cell_meta_table.tsv --cell_label sub --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --script_out log/pseudo_bulk_sub_expected_cis.202501.sh
```

Here, `--folder` is the directory containing the mcool files, `--cell_anno` is the path to the cell annotation file, `--cell_label` specifies whether the cells are major or sub, `--blacklist` is the path to the blacklist file, and `--script_out` is the path to save the generated script for running the command.
The above scripts will merge the mcool files from each sample into a single mcool file for each cell type (major and sub) and create the expected contact map for each cell type.

### Call 3D genome features

#### Call pseudo bulk features
We used the following script to perform saddle plot analysis, APA analysis, and long-/short-range summary analysis.

1. **Loop analysis**:
```bash
## get the APA results on loops
python src/cell_pseudo_bulk.loop_pileup.v2.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --loop_folder shared_results/result_GAGE-seq/pipeline_3D_feature/202501/loop --loop_processed_folder shared_results/result_GAGE-seq/pipeline_3D_feature/202501/loop_processed --output_folder results/202501/loop_pileup/major --script_out log/pseudo_bulk_major_loop_pileup.202504.sh
# parse all APA results to make a summary file
python src/cell_pseudo_bulk.loop_pileup_summary.py --folder results/202501/loop_pileup/major/AD --output results/202501/loop_pileup/major/AD_summary.tsv
python src/cell_pseudo_bulk.loop_pileup_summary.py --folder results/202501/loop_pileup/major/CT --output results/202501/loop_pileup/major/CT_summary.tsv
```
Here, `--folder` is the directory containing the mcool files, `--cell_anno` is the path to the cell annotation file, `--cell_label` specifies whether the cells are major or sub, `--view_file` is the path to the view file, `--loop_folder` is the directory containing the loop files, `--loop_processed_folder` is the directory containing the processed loop files, `--output_folder` is the directory where the results will be saved, and `--script_out` is the path to save the generated script for running the command.
Chromatin loops were called using the [scHiCluster](https://github.com/zhoujt1994/scHiCluster) pipeline for each major cell type. The loop files were then processed to generate the APA results for each cell type. The APA results were summarized and a differential analysis was performed between AD and CT cell types.

2. **Comparment analysis**:
```bash
# calculate the AB compartment (cis) using cooltools
python src/cell_pseudo_bulk.eigs_cis.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --phasing_track annotation/hg38M_100kb.GC_content.bed --script_out log/pseudo_bulk_major_eigs_cis.202501.sh 
python src/cell_pseudo_bulk.eigs_cis.py --folder results/202501/cool_file/sub --cell_anno results/202501/cell_meta_table.tsv --cell_label sub --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --phasing_track annotation/hg38M_100kb.GC_content.bed --script_out log/pseudo_bulk_sub_eigs_cis.202501.sh 
# get the saddle plots on AB score - cis
python src/cell_pseudo_bulk.saddle.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --script_out log/pseudo_bulk_major_saddle.202501.sh 
python src/cell_pseudo_bulk.saddle_summary.py --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major  --output results/202501/cool_file/major/saddle_summary.tsv  
python src/cell_pseudo_bulk.saddle.py --folder results/202501/cool_file/sub --cell_anno results/202501/cell_meta_table.tsv --cell_label sub --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --script_out log/pseudo_bulk_sub_saddle.202501.sh 
# get the saddle plots on AB score - trans
python src/cell_pseudo_bulk.saddle.py --use_trans --only_AD_and_CT --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --view_file results/202501/cool_file/view_hg38_10kb.tsv --thread 12 --script_out log/pseudo_bulk_major_saddle.202505.sh 
python src/cell_pseudo_bulk.saddle_summary.py --use_trans --folder results/202501/cool_file/major --cell_anno results/202501/cell_meta_table.tsv --cell_label major --output results/202501/cool_file/major/saddle_summary.trans.tsv  
#  subompartment vs distance
python src/cell_pseudo_bulk.compartment_distance.py --cool_folder results/202501/cool_file/sub --cell_anno results/202501/cell_meta_table.tsv --read_pair_folder results/202501/raw_pair --output results/202501/cool_file/sub/compartment_distance_summary.tsv --thread 12
```
Here, `--folder` is the directory containing the mcool files, `--cell_anno` is the path to the cell annotation file, `--cell_label` specifies whether the cells are major or sub, `--view_file` is the path to the view file, `--phasing_track` is the path to the phasing track file, and `--script_out` is the path to save the generated script for running the command. 
We used the GC content file `hg38M_100kb.GC_content.bed` to ensure that positive scores correspond to active chromatin. 

3. **Long-/short-range interaction analysis**:
```bash
# summarize the chromatin interaction into bins based on the distance between the two ends of the read pairs
python src/cell_stat_read_count.py --folder processed_data/202501/raw_pair/${sample_id} --patient_id ${sample_id} --read_count_cutoff 30000 --output processed_data/202501/stats/summary_read_pair_${sample_id} --folder_filtered processed_data/202501/filtered_pair/${sample_id}
# variantion of the above command for CTCF peaks
python src/cell_stat_read_count.CTCF.py --folder processed_data/202501/filtered_pair/${sample_id} --patient_id ${sample_id}  --ctcf_peak annotation/CTCF_recalibrated.final.strict.bed --output processed_data/202501/stats/summary_read_pair_${sample_id}
python src/cell_stat_read_count.CTCF.py --folder processed_data/202501/filtered_pair/${sample_id} --patient_id ${sample_id}  --ctcf_peak annotation/CTCF_recalibrated.final.bed --output processed_data/202501/stats/summary_read_pair_${sample_id}.loose
# variation of the above command for ATAC peaks
python src/cell_stat_read_count.ATAC.py --folder processed_data/202501/filtered_pair/${patient_id} --patient_id ${patient_id}  --cell_white_list processed_data/202501/cell_meta_table.tsv --cell_white_list_format cell_anno --ATAC_folder xiong2023/peak_folder --ATAC_format xiong --output processed_data/202501/stats/summary_read_pair_${patient_id}
```
Here, `--folder` is the directory containing the processed read pairs, `--patient_id` is the identifier for the patient, `--read_count_cutoff` is the minimum read count for filtering, `--output` is the path to save the summary file, and `--folder_filtered` is the directory where the filtered read pairs will be saved.

#### Call single-cell features
We used [Higashi](https://github.com/ma-compbio/Higashi) to impute scHi-C contact maps at 500 kb and 100 kb resolutions from the GAGE-seq data. Please refer to homepage of Higashi for the detailed usage. Based on the Higashi-imputed scHi-C contact maps, multiple 3D genome features were derived. Single-cell A/B compartment scores (at 100 kb and 500 kb resolutions) and single-cell insulation scores (at 100 kb resolution) were calculated using the scCompartment.py and scTAD.py from the Higashi software suite with default parameters. Additionally, [scGHOST](https://github.com/ma-compbio/scGHOST)  was used to calculate single-cell subcompartment annotations based on Higashi-imputed scHi-C contact maps and scA/B compartment scores at 500 kb resolution. 
