
# Hicformer: A transformer-based model to predict gene expression

This folder contains the code for the Hicformer model, which is a deep learning model to predict cell type-specific gene expression by integrating DNA sequence and 3D genome structure. Each input consists of a 409,600 bp genomic region and a corresponding Hi-C contact map at 1,024-bp resolution. The model was trained and evaluated using pseudo-bulk data aggregated at the cell type level.

## Input data
+ genome_regions_1024_200.bed: locations of genomic regions of interest. The first 400 rows are for the first gene; the next 400 rows are for the second gene; and so on. 
+ genome_regions_genes_1024_200.bed: locations of genes of interest
+ cells.tsv: metadata of cells
+ genes.tsv: metadata of genes
+ gene_id_conversion.tsv: mapping between Ensembl ID and gene name
+ cell_types.tsv: list of cell types
+ sequence_1024_200.tsv: DNA sequences of genes, one row per gene
+ expression_cov_1024_200_celltypebulk.pkl: a sparse matrix of cell-type-pseudo-bulk expression. Rows are cell types. Columns are expressions at individual 1024-bp loci. The first 400 columns are from the first gene, the following 400 columns are from the second gene, and so on
+ contact_1024_200_celltypebulk.h5: sparse matrices of cell-type-pseudo-bulk contact maps. The contact map for gene `gene_id` in cell type `cell_type_id` is stored under key `{cell_type_id}/{gene_id}` as a sparse COO matrix. The shape for each contact map is 400*400. The contact maps are symmetric and only the upper triangular part is stored
+ 1d-score-celltypebulk-10kb-*_1024_200_uint8.pkl: dense matrix of 1D scores, including A/B, insulation score (3 settings), and gene body score. Rows are cell types and columns are 1024-bp genomic regions

## How to run the model

+ Baseline  
  `python hcformer_pbulk_train.py --use_wandb --use_sweep --batch_size 64 --epochs 150 --gpu 0 1 2 3 4 5 6 7 --num 1 --gt_mode raw --dim 768`

+ HiCFormer 1d  
  `python hcformer_pbulk_train.py --use_wandb --use_sweep --batch_size 64 --epochs 150 --gpu 0 1 2 3 4 5 6 7 --num 1 --gt_mode raw --dim 768  --hic_1d`

+ HiCFormer 2d  
  `python hcformer_pbulk_train.py --use_wandb --use_sweep --batch_size 64 --epochs 150 --gpu 0 1 2 3 4 5 6 7 --num 1 --gt_mode raw --dim 768 --hic_2d`

+ HiCFormer 1d2d  
  `python hcformer_pbulk_train.py --use_wandb --use_sweep --batch_size 64 --epochs 150 --gpu 0 1 2 3 4 5 6 7 --num 1 --gt_mode raw --dim 768 --hic_1d --hic_2d`
