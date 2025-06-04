
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(ggrepel)
library(DESeq2)
library(Seurat)
library(pheatmap)
library(RColorBrewer)

load("D:/Project_AD/data/Xenium_Seuratanalysis.Rdata")
##CHoose the cells not in the data_spatial
domain_label <- read.csv("D:/Project_AD/Results/Results_Standard/Domain_label.csv")
rownames(domain_label)<-domain_label$cellid
domain_label$sample_id <- data_spatial$sample_id
data_spatial$domain_label <- domain_label$domain_label
table(data_spatial$domain_label,data_spatial$Cluster)
prop_table<-prop.table(table(data_spatial$domain_label,data_spatial$Cluster),1)
prop_table["WM",]
cell_ifinWM <- data_spatial$domain_label == "WM"
data_spatial$cell_ifinWM <- cell_ifinWM
data_spatial <- subset(data_spatial, cell_ifinWM == FALSE)
# ##Merge celltypes
# data_spatial$Cluster<-ifelse(data_spatial$Cluster %in% c("Exc L2/3 IT", "Exc L4 IT", "Exc L5 ET", "Exc L5 IT", "Exc L5/6 NP", "Exc L6 CT", "Exc L6 IT", "Exc L6 IT Car3", "Exc L6b"),"Exc",
#                              ifelse(data_spatial$Cluster %in% c("Inh Chandelier", "Inh Lamp5","Inh Pvalb","Inh Sncg","Inh Sst","Inh Vip"), "Inh",data_spatial$Cluster))
#####
#####
# data_spatial <- anndataR::read_h5ad("D:/Project_AD/data/processed_dataset_scvi.h5ad")
# data_spatial<-data_spatial$to_Seurat()
# data_spatial1 <- data_spatial
# load("D:/Project_AD/data/Xenium_Seuratanalysis.Rdata")
# data_spatial1$Cluster<- data_spatial$Cluster[names(data_spatial1$Cluster)]
# remove(data_spatial)
# data_spatial <-data_spatial1
# remove(data_spatial1)
# # load sample meta table
# data_sample <- read.table(file = "annotation/sample_info.v2.tsv", header = T, sep = '\t', stringsAsFactors = F)
# aggregate the counts based on donor-condition-celltype into  pseudobulk
# pseudo_obj_major <- AggregateExpression(data_spatial, assays = "RNA", return.seurat = T, group.by = c("major_cell_type", "condition", "sample_id"))
# pseudo_obj_major$cell_type <- pseudo_obj_major$major_cell_type
pseudo_obj_sub <- AggregateExpression(data_spatial, assays = "RNA", return.seurat = T, group.by = c("Cluster", "condition", "sample_id"))
pseudo_obj_sub$cell_type <- pseudo_obj_sub$Cluster
# get the cell count per sample
# pseudo_cell_count_major <- as.data.frame(table(obj[[]]$sample_id, obj[[]]$major_cell_type))
pseudo_cell_count_sub <- as.data.frame(table(data_spatial$sample_id, data_spatial$Cluster))
# colnames(pseudo_cell_count_major) <- c("projid", "cell_type", "cell_count")
colnames(pseudo_cell_count_sub) <- c("projid", "cell_type", "cell_count")

# ############################
# pseudo_obj<-pseudo_obj_sub
# pseudo_cell_count<-pseudo_cell_count_sub
# target_cell_type <- "Exc L5 ET"
# control_gene_set<-NULL
# is_gene_filter<-TRUE
# padj_cutoff<-0.05
# output_prefix<-"test"
############################
run_DESeq2_analysis <- function(pseudo_obj, pseudo_cell_count, target_cell_type, control_gene_set = c("TOP1", "PUM1", "PPIA", "CYC1", "UBE2D2", "RPL13", "TBP", "GAPDH", "ATP5F1B"),
                                is_gene_filter = TRUE,
                                padj_cutoff = 0.05, output_prefix = "test") {
  target_cell_type_show <- gsub(" ", "_", target_cell_type)
  target_cell_type_show <- gsub("/", "_", target_cell_type_show)
  ##
  # subset the cell type
  pseudo_obj_cell_type <- subset(pseudo_obj, cell_type == target_cell_type)
  # subse the cell count
  pseudo_cell_count_cell_type <- subset(pseudo_cell_count, cell_type == target_cell_type)
  # report the number of pseudobulk samples
  cat(paste0("Number of pseudobulk samples are ", ncol(pseudo_obj_cell_type)))
  # get sample list
  sample_list = colnames(pseudo_obj_cell_type)
  # get count matrix
  expr_count_matrix <- as.matrix(pseudo_obj_cell_type[['RNA']]$counts)
  # build column data
  df_sample_data <- data.frame(sample_id = sample_list)
  df_sample_data$projid <- gsub(".*_", "", df_sample_data$sample_id)
  # df_sample_data <- merge(df_sample_data, data_sample[, c("projid", "disease_group", "gender", "pmi", "age_group")], by = "projid", sort = F)
  # add the cell count information
  df_sample_data <- merge(df_sample_data, pseudo_cell_count_cell_type[, c("projid", "cell_count")], by = "projid", sort = F)
  # add rownames
  rownames(df_sample_data) <- df_sample_data$sample_id
  df_sample_data$sample_id <- NULL
  # convert character to factor
  df_sample_data$disease_group<- c("AD","AD","CT","CT")
  df_sample_data$disease_group <- factor(df_sample_data$disease_group, levels = c("CT", "AD"))
  # df_sample_data$gender <- factor(df_sample_data$gender, levels = c("male", "female"))
  # double check if the order the patients match exactly with colnames of the count matrix
  cat(paste0("Sample order match: ", all(rownames(df_sample_data) == colnames(expr_count_matrix))), '\n')
  # build DESeq2 object
  ##turn the expr_count_matrix to integer
  expr_count_matrix <- round(expr_count_matrix)
  dds <- DESeqDataSetFromMatrix(countData = expr_count_matrix,
                                colData = df_sample_data,
                                design = ~ disease_group)
  # filter genes 
  if (is_gene_filter) {
    smallestGroupSize <- 3
    keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
    cat("Do gene filtering, total gene left ", sum(keep), "\n")
    dds_filtered_by_gene <- dds[keep,]
  } else {
    dds_filtered_by_gene <- dds
  }
  # Transform counts for data visualization
  # rld <- rlog(dds_filtered_by_gene, blind=FALSE)
  # panel_sample_pca_rlog <- DESeq2::plotPCA(rld, ntop = min(5000, nrow(rld)), intgroup = "disease_group", returnData = F) +
  #   ggtitle("Sample PCA using rlog normalization") +
  #   scale_color_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
  #   theme_cowplot() +
  #   theme(aspect.ratio = 1)
  ##ggsave panel_sample_pca_rlog
  # ggsave(paste0("D:/Project_AD/Results/Results_pseudobulk/sample_pca_rlog_", target_cell_type_show, ".jpg"), plot = panel_sample_pca_rlog, width = 6, height = 6)
  # panel_sample_pca_rlog_cell_count <- DESeq2::plotPCA(rld, ntop = min(5000, nrow(rld)), intgroup = "cell_count", returnData = F) + 
  #   geom_point(color = "black", size = 3, pch = 21) +
  #   ggtitle("Sample PCA using normalized count with control gene set") +
  #   #scale_color_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
  #   scale_color_distiller(name="#cell", palette = "Reds", direction = 1) +
  #   theme_cowplot() +
  #   theme(aspect.ratio = 1)
  # ##ggsave
  # ggsave(paste0("D:/Project_AD/Results/Results_pseudobulk/sample_pca_rlog_cell_count_", target_cell_type_show, ".jpg"), plot = panel_sample_pca_rlog_cell_count, width = 6, height = 6)
  # Extract the rlog matrix from the object and compute pairwise correlation values
  # rld_mat <- assay(rld)
  # rld_cor <- cor(rld_mat)
  # plot heatmap
  #pheatmap(rld_cor, annotation = df_sample_data[, c("disease_group"), drop=F])
  # plot global gene expression
  # use rlog transform
  # raw_gene_count <- as.data.frame(rld_mat)
  gene_count_raw <- as.data.frame(log2(counts(dds_filtered_by_gene, normalized = FALSE) + 1))
  gene_count_raw$gene_id <- rownames(gene_count_raw)
  gene_count_raw <- gene_count_raw %>%
    pivot_longer(col = -gene_id, names_to = "sample_id", values_to = "expr_count") %>%
    mutate(projid = gsub(".*_", "", sample_id))
  gene_count_raw$sample_id <- NULL
  gene_count_raw <- merge(gene_count_raw, df_sample_data, by = "projid", sort = F)
  #
  gene_count_raw_mu <- gene_count_raw %>%
    group_by(disease_group) %>%
    summarise(mean_expr = mean(expr_count, na.rm = T))
  panel_gene_count_raw <- gene_count_raw %>%
    ggplot(aes(x = expr_count, fill = disease_group)) +
    geom_density(alpha = 0.5, linewidth = 0.5) +
    geom_vline(data = gene_count_raw_mu, aes(xintercept = mean_expr, color = disease_group), linetype = 2, linewidth = 1) +
    xlab("log2 UMI count +1") +
    ylab("Density") +
    scale_color_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
    scale_fill_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_cowplot()
  ggsave(paste0("D:/Project_AD/Results/Results_pseudobulk/gene_count_raw_", target_cell_type_show, ".jpg"), plot = panel_gene_count_raw, width = 6, height = 6)
  # use control gene sets to normalize 
  if (length(control_gene_set) > 0) {
    isControl <- rownames(dds_filtered_by_gene) %in% control_gene_set
    dds_filtered_by_gene <- estimateSizeFactors(dds_filtered_by_gene, controlGenes = isControl)
  }
  # Run DESeq2 differential expression analysis
  dds_filtered_by_gene <- DESeq(dds_filtered_by_gene)
  # Plot dispersion estimates
  #plotDispEsts(dds_filtered_by_gene)
  # Check the coefficients for the comparison 
  cat(paste0("check the contrast: ", resultsNames(dds_filtered_by_gene), '\n'))
  # Generate results object
  res <- results(dds_filtered_by_gene, contrast=c("disease_group","AD","CT"))
  # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
  # resLFC <- lfcShrink(dds_filtered_by_gene, 
  #                     coef = "disease_group_AD_vs_CT",
  #                     res=res,
  #                     type = "apeglm")
  ## sort the result by p-value
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene_name") %>%
    as_tibble() %>%
    arrange(padj)
  res_tbl_df<- data.frame(res_tbl)
  ##Save the res_tbl as tsv file
  write.table(res_tbl_df, file = paste0("D:/Project_AD/Results/Results_pseudobulk/res_tbl_", target_cell_type_show, ".tsv"), sep = '\t', quote = F, row.names = F)
  # Subset the significant results
  sig_res <- res_tbl %>%
    filter(padj < padj_cutoff) %>%
    arrange(padj)
  # MA plot
  data_MA_plot <- plotMA(res, ylim=c(-3,3), alpha = 0.5, returnData = TRUE)
  panel_MA_plot <- data_MA_plot %>%
    mutate(gene_group = ifelse(isDE, "DE", "non-DE")) %>%
    ggplot(aes(x = mean+1, y = lfc, color = gene_group)) +
    geom_point(alpha = 0.6, size = 0.4) +
    geom_hline(yintercept = 0.0, color = "black") +
    xlab("mean of normalized counts") +
    ylab("log fold change") +
    ggtitle("MA plot using normalized counts") +
    scale_color_manual(name = "", values = c("DE" = "blue", "non-DE" = "lightgrey")) +
    scale_x_continuous(trans = "log10") +
    theme_cowplot() +
    theme(aspect.ratio = 2.0)
  ## Extract normalized counts from dds object
  gene_count_normalized <- as.data.frame(log2(counts(dds_filtered_by_gene, normalized = TRUE) + 1))
  gene_count_normalized$gene_id <- rownames(gene_count_normalized)
  gene_count_normalized <- gene_count_normalized %>%
    pivot_longer(col = -gene_id, names_to = "sample_id", values_to = "expr_count") %>%
    mutate(projid = gsub(".*_", "", sample_id))
  gene_count_normalized$sample_id <- NULL
  gene_count_normalized <- merge(gene_count_normalized, df_sample_data, by = "projid", sort = F)
  if (target_cell_type == "Astro"){
    gene_set_astro_up_cell_death_oxygen <- c("APOE", "CST3", "PSAP", "RGCC", "CLU", "FAIM2", "GAPDH", "PRDX1", "MEF2C", "HGF", "HSP90AB1", "USP53", "EGLN3", "PTPRK", "PRDX6", "OXR1", "DDIT4", "SESN1", "THRB", "GSN", "SOD2", "HSPA1A", "PPARGC1A", "NR4A3", "SLC7A11", "ANO6", "HSPB1", "TFDP2", "FANCC", "BNIP3L", "PDE8A", "PRKN")
    gene_set_astro_up_lipid_storage <- c("AQP1", "OSBPL11", "LRAT", "PAM", "FOXP1", "MAOB", "C3", "ABCA1", "NR4A3", "FBXO32", "IRS2", "IRS2", "PADI2", "SPARC", "AKAP13", "PPARGC1A", "ACACB", "NTRK3", "FAM107A", "PLPP1", "TGFB2", "PLIN5", "FOXO1")
    gene_set_astro_down_angiogenesis <- c("PML", "F3", "NRP1", "SULF1", "GADD45A", "ANGPTL4", "VEGFA", "ACTG1", "RHOB", "ACTB", "ID1")
    gene_set_astro_up_cell_death_oxygen<-gene_set_astro_up_cell_death_oxygen[which(gene_set_astro_up_cell_death_oxygen %in% unique(gene_count_normalized$gene_id))]
    gene_set_astro_up_lipid_storage<-gene_set_astro_up_lipid_storage[which(gene_set_astro_up_lipid_storage %in% unique(gene_count_normalized$gene_id))]
    gene_set_astro_down_angiogenesis<-gene_set_astro_down_angiogenesis[which(gene_set_astro_down_angiogenesis %in% unique(gene_count_normalized$gene_id))]
    ##
    heat_colors <- rev(brewer.pal(11, "PuOr"))
    ##
    gene_count_normalized_part1 <- gene_count_normalized %>%
      filter(gene_id %in% gene_set_astro_up_cell_death_oxygen) %>%
      mutate(gene_group = "Astro_up_cell_death_oxygen")
    ##Form as projid by gene_id, value is cell_count
    gene_count_normalized_part1 <- gene_count_normalized_part1 %>%
      pivot_wider(names_from = gene_id, values_from = expr_count)
    gene_count_normalized_part1<-as.data.frame(gene_count_normalized_part1)
    gene_count_normalized_part1<-gene_count_normalized_part1[,c(1,5:ncol(gene_count_normalized_part1))]
    rownames(gene_count_normalized_part1)<-gene_count_normalized_part1$projid
    gene_count_normalized_part1<-gene_count_normalized_part1[,-1]
    gene_count_normalized_part1_heatmap <- pheatmap(t(gene_count_normalized_part1), 
                                                    color = heat_colors, 
                                                    # cluster_rows = TRUE, 
                                                    # cluster_cols = TRUE,
                                                    cluster_rows = FALSE,
                                                    cluster_cols = FALSE,
                                                    show_rownames = TRUE,
                                                    clustering_method = "average",
                                                    # annotation = df_sample_data[, c("disease_group", "gender")], 
                                                    border_color = NA, 
                                                    fontsize = 8, 
                                                    scale = "row", 
                                                    fontsize_row = 10, 
                                                    height = 20) 
    ##
    gene_count_normalized_part2 <- gene_count_normalized %>%
      filter(gene_id %in% gene_set_astro_up_lipid_storage) %>%
      mutate(gene_group = "Astro_up_lipid_storage")
    ##Form as projid by gene_id, value is cell_count
    gene_count_normalized_part2 <- gene_count_normalized_part2 %>%
      pivot_wider(names_from = gene_id, values_from = expr_count)
    gene_count_normalized_part2<-as.data.frame(gene_count_normalized_part2)
    gene_count_normalized_part2<-gene_count_normalized_part2[,c(1,5:ncol(gene_count_normalized_part2))]
    rownames(gene_count_normalized_part2)<-gene_count_normalized_part2$projid
    gene_count_normalized_part2<-gene_count_normalized_part2[,-1]
    gene_count_normalized_part2_heatmap <- pheatmap(t(gene_count_normalized_part2), 
                                                    color = heat_colors, 
                                                    # cluster_rows = TRUE, 
                                                    # cluster_cols = TRUE,
                                                    cluster_rows = FALSE,
                                                    cluster_cols = FALSE,
                                                    show_rownames = TRUE,
                                                    clustering_method = "average",
                                                    # annotation = df_sample_data[, c("disease_group
                                                    border_color = NA,
                                                    fontsize = 8,
                                                    scale = "row",
                                                    fontsize_row = 10,
                                                    height = 20)
    ##
    gene_count_normalized_part3 <- gene_count_normalized %>%
      filter(gene_id %in% gene_set_astro_down_angiogenesis) %>%
      mutate(gene_group = "Astro_down_angiogenesis")
    ##Form as projid by gene_id, value is cell_count
    gene_count_normalized_part3 <- gene_count_normalized_part3 %>%
      pivot_wider(names_from = gene_id, values_from = expr_count)
    gene_count_normalized_part3<-as.data.frame(gene_count_normalized_part3)
    gene_count_normalized_part3<-gene_count_normalized_part3[,c(1,5:ncol(gene_count_normalized_part3))]
    rownames(gene_count_normalized_part3)<-gene_count_normalized_part3$projid
    gene_count_normalized_part3<-gene_count_normalized_part3[,-1]
    gene_count_normalized_part3_heatmap <- pheatmap(t(gene_count_normalized_part3), 
                                                    color = heat_colors, 
                                                    # cluster_rows = TRUE, 
                                                    # cluster_cols = TRUE,
                                                    cluster_rows = FALSE,
                                                    cluster_cols = FALSE,
                                                    show_rownames = TRUE,
                                                    clustering_method = "average",
                                                    # annotation = df_sample_data[, c("disease_group
                                                    border_color = NA,
                                                    fontsize = 8,
                                                    scale = "row",
                                                    fontsize_row = 10,
                                                    height = 20)
    
  }
  # plot the density of normalized gene count
  gene_count_normalized_mu <- gene_count_normalized %>%
    group_by(disease_group) %>%
    summarise(mean_expr = mean(expr_count, na.rm = T))
  panel_gene_count_normalized <- gene_count_normalized %>%
    ggplot(aes(x = expr_count, fill = disease_group)) +
    geom_density(alpha = 0.5, linewidth = 0.5) +
    geom_vline(data = gene_count_normalized_mu, aes(xintercept = mean_expr, color = disease_group), linetype = 2, linewidth = 1) +
    xlab("log2 normalized count + 1") +
    ylab("Density") +
    scale_color_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
    scale_fill_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_cowplot()
  ## Extract normalized counts for significant genes only
  normalized_counts <- counts(dds_filtered_by_gene, normalized = TRUE)
  sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene_name, ]
  if (length(sig_res$gene_name) == 1) {
    sig_counts <- t(sig_counts)
    rownames(sig_counts) <- sig_res$gene_name
  }
  if (nrow(sig_counts)>=2){
    # plot the pca using normalized gene count
    se <- SummarizedExperiment(log2(counts(dds_filtered_by_gene, normalized=TRUE) + 1),
                               colData=colData(dds_filtered_by_gene))
    panel_sample_pca_normalized_count <- DESeq2::plotPCA(DESeqTransform(se), ntop = min(5000, nrow(normalized_counts)), intgroup = "disease_group", returnData = F) + 
      ggtitle("Sample PCA using gene set norm") +
      scale_color_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
      theme_cowplot() +
      theme(aspect.ratio = 1)
    panel_sample_pca_normalized_count_cell_count <- DESeq2::plotPCA(DESeqTransform(se), ntop = min(5000, nrow(normalized_counts)), intgroup = "cell_count", returnData = F) + 
      geom_point(color = "black", size = 3, pch = 21) +
      ggtitle("Sample PCA using gene set norm") +
      #scale_color_manual(name = "", values = c("AD" = "#f1a340", "CT" = "#998ec3")) +
      scale_color_distiller(name="#cell", palette = "Reds", direction = 1) +
      theme_cowplot() +
      theme(aspect.ratio = 1)
    ## Set a color-blind friendly palette
    heat_colors <- rev(brewer.pal(11, "PuOr"))
    ## Run pheatmap using the metadata data frame for the annotation
    panel_sig_gene_heatmap <- pheatmap(sig_counts, 
                                       color = heat_colors, 
                                       cluster_rows = TRUE, 
                                       cluster_cols = TRUE,
                                       show_rownames = TRUE,
                                       clustering_method = "average",
                                       # annotation = df_sample_data[, c("disease_group", "gender")], 
                                       border_color = NA, 
                                       fontsize = 8, 
                                       scale = "row", 
                                       fontsize_row = 10, 
                                       height = 20) 
    ##ggsave
    ggsave(paste0("D:/Project_AD/Results/Results_pseudobulk/sample_pca_normalized_count_", target_cell_type_show, ".jpg"), plot = panel_sample_pca_normalized_count, width = 6, height = 6)
    
  }
  # plot the figures
  panel_gene_density <- plot_grid(panel_gene_count_raw, panel_gene_count_normalized, align = "v", axis = 'lr', ncol = 1)
  ggsave(paste0("D:/Project_AD/Results/Results_pseudobulk/gene_density_", target_cell_type_show, ".jpg"), plot = panel_gene_density, width = 6, height = 6)
  # return(res_tbl)
  return(list(gene_count_normalized,res_tbl_df))
}
############################
gene_count_normalized_list <- list()
res_tbl_df_list <- list()
celltype_unique <- unique(data_spatial$Cluster)
for (celltyoe_cur in celltype_unique) {
  cat(paste0("Start to run DESeq2 analysis for cell type ", celltyoe_cur, '\n'))
  run_DESeq2_analysis_obj <- run_DESeq2_analysis(pseudo_obj_sub, pseudo_cell_count_sub,
                                                 celltyoe_cur, control_gene_set = NULL, is_gene_filter = TRUE, padj_cutoff = 0.05, output_prefix = "test")
  gene_count_normalized_list[[celltyoe_cur]] <- run_DESeq2_analysis_obj[[1]]
  res_tbl_df_list[[celltyoe_cur]] <- run_DESeq2_analysis_obj[[2]]
}

##Load the result of GAGEseq
file_GAGEseq <- "D:/Project_AD/Results/Results_pseudobulk/GAGEseq/"
if_have_major <-FALSE
res_tbl_df_list_GAGEseq <- list()
celltype_unique <- unique(data_spatial$Cluster)
for (celltyoe_cur in celltype_unique) {
  celltyoe_cur_show <- gsub(" ", "_", celltyoe_cur)
  celltyoe_cur_show <- gsub("/", "_", celltyoe_cur_show)
  ##load the tsv file
  if (celltyoe_cur %in% c("Exc","Inh")){
    file_cur<-paste0(file_GAGEseq, "major_noControlGeneSet.", celltyoe_cur_show, ".tsv") 
    if_have_major <-TRUE
  }else{
    file_cur<-paste0(file_GAGEseq, "sub_noControlGeneSet.", celltyoe_cur_show, ".tsv")
  }
  ##Chech whether the file exists
  if (!file.exists(file_cur)){
    cat(paste0("The file ", file_cur, " does not exist, skip this cell type\n"))
    next
  }else{
    data_cur <- read.table(file =file_cur, header = T, sep = '\t', stringsAsFactors = F)
    res_tbl_df_list_GAGEseq[[celltyoe_cur]] <- data_cur
  }
}
##
# file_GAGEseq <- "D:/Project_AD/Results/Results_pseudobulk/GAGEseq_female/"
# 
# res_tbl_df_list_GAGEseq <- list()
# celltype_unique <- unique(data_spatial$Cluster)
# for (celltyoe_cur in celltype_unique) {
#   celltyoe_cur_show <- gsub(" ", " ", celltyoe_cur)
#   celltyoe_cur_show <- gsub("/", " ", celltyoe_cur_show)
#   if (celltyoe_cur == "Exc L2/3 IT"){
#     celltyoe_cur_show <- "Exc L2-3 IT"
#   }
#   if (celltyoe_cur == "Exc L5/6 IT"){
#     celltyoe_cur_show <- "Exc L5-6 IT"
#   }
#   if (celltyoe_cur == "Exc L5/6 NP"){
#     celltyoe_cur_show <- "Exc L5-6 NP"
#   }
#   ##load the tsv file
#   file_cur<-paste0(file_GAGEseq, celltyoe_cur_show,"-AD-vs-CT-MAST", ".csv")
#   ##Chech whether the file exists
#   if (!file.exists(file_cur)){
#     cat(paste0("The file ", file_cur, " does not exist, skip this cell type\n"))
#     next
#   }else{
#     data_cur <- read.table(file =file_cur, header = T, sep = ',', stringsAsFactors = F)
#     names(data_cur)[1]<-"gene_name"
#     names(data_cur)[3]<-"log2FoldChange"
#     names(data_cur)[2]<-"pvalue"
#     names(data_cur)[6]<-"padj"
#     res_tbl_df_list_GAGEseq[[celltyoe_cur]] <- data_cur 
#   }
# }
##
celltype_common <- intersect( names(res_tbl_df_list), names(res_tbl_df_list_GAGEseq))

res_tbl_df_twotech_merge <- data.frame(gene = character(), log2FoldChange_Xenium = numeric(), log2FoldChange_GAGEseq = numeric(), cell_type = character())
corr_list <- list()
for (celltype_cur in celltype_common) {
  res_tbl_df_cur <- res_tbl_df_list[[celltype_cur]]
  res_tbl_df_cur_GAGEseq <- res_tbl_df_list_GAGEseq[[celltype_cur]]
  rownames(res_tbl_df_cur)<-res_tbl_df_cur$gene_name
  rownames(res_tbl_df_cur_GAGEseq)<-res_tbl_df_cur_GAGEseq$gene_name
  commom_gene_cur <- intersect(res_tbl_df_cur$gene_name, res_tbl_df_cur_GAGEseq$gene_name)
  res_tbl_df_cur<-res_tbl_df_cur[commom_gene_cur,]
  res_tbl_df_cur_GAGEseq<-res_tbl_df_cur_GAGEseq[commom_gene_cur,]
  sig_gene_cur <- union(res_tbl_df_cur$gene_name[res_tbl_df_cur$pvalue<0.05],res_tbl_df_cur_GAGEseq$gene_name[res_tbl_df_cur_GAGEseq$pvalue<0.05])
  ##Merge the two data frame
  res_tbl_df_twotech <- data.frame(gene = commom_gene_cur, log2FoldChange_Xenium = res_tbl_df_cur$log2FoldChange,
                                   log2FoldChange_GAGEseq = res_tbl_df_cur_GAGEseq$log2FoldChange,
                                   cell_type = celltype_cur)
  res_tbl_df_twotech<-res_tbl_df_twotech[res_tbl_df_twotech$gene %in% sig_gene_cur,]
  res_tbl_df_twotech_merge <- rbind(res_tbl_df_twotech_merge, res_tbl_df_twotech)
  corr_obj <- cor.test(res_tbl_df_twotech_merge$log2FoldChange_GAGEseq, res_tbl_df_twotech_merge$log2FoldChange_Xenium)
  corr_list[[celltype_cur]] <- as.numeric(corr_obj$estimate)
}
corr_vec <- unlist(corr_list)
corr_vec
summary(corr_vec)

ncol_use <- 5
if (if_have_major){
  ncol_use <- 4
}
library(ggplot2)
res_tbl_df_twotech_merge$gene_type <- ifelse(sign(res_tbl_df_twotech_merge$log2FoldChange_GAGEseq) != sign(res_tbl_df_twotech_merge$log2FoldChange_Xenium), "opposite", "same")
ggplot(data = res_tbl_df_twotech_merge,aes(x = log2FoldChange_Xenium, y = log2FoldChange_GAGEseq)) +
  geom_point(aes(color = gene_type)) +
  facet_wrap(~cell_type, scales = "free",ncol = ncol_use) + ##Add regression line, line's color is black, also add the confidence interval
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  xlab("log2FoldChange_Xenium") +
  ylab("log2FoldChange_GAGEseq") +
  # ggtitle("Comparison of log2FoldChange between Xenium and GAGEseq") +
  theme_cowplot()+
  scale_color_manual(values = c("opposite" = "#CBE0F1", "same" = "#164893"))



library(ggplot2)
library(ggpmisc)
library(cowplot)
library(ggpubr)
ncol_use <- 7
if (if_have_major){
  ncol_use <- 4
}
res_tbl_df_twotech_merge <- res_tbl_df_twotech_merge %>%
  mutate(quadrant = case_when(
    log2FoldChange_Xenium >= 0 & log2FoldChange_GAGEseq >= 0 ~ "Q1",
    log2FoldChange_Xenium <  0 & log2FoldChange_GAGEseq >= 0 ~ "Q2",
    log2FoldChange_Xenium <  0 & log2FoldChange_GAGEseq <  0 ~ "Q3",
    log2FoldChange_Xenium >= 0 & log2FoldChange_GAGEseq <  0 ~ "Q4"
  ))

quad_counts <- res_tbl_df_twotech_merge %>%
  group_by(cell_type, quadrant) %>%
  summarise(n = n(), .groups = "drop")

quad_labels <- quad_counts %>%
  mutate(x = case_when(
    quadrant == "Q1" ~ 2,
    quadrant == "Q2" ~ -2,
    quadrant == "Q3" ~ -2,
    quadrant == "Q4" ~ 2
  ),
  y = case_when(
    quadrant == "Q1" ~ 2.5,
    quadrant == "Q2" ~ 2.5,
    quadrant == "Q3" ~ -2.5,
    quadrant == "Q4" ~ -2.5
  ),
  label = paste0("n = ", n))

ggplot(data = res_tbl_df_twotech_merge, 
       aes(x = log2FoldChange_Xenium, y = log2FoldChange_GAGEseq)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#A2A2A2", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2, color = "#A2A2A2", linewidth = 1) +
  geom_point(aes(color = gene_type), size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", fullrange = FALSE) +
  geom_text(data = quad_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6) +
  facet_wrap(~cell_type, scales = "free", ncol = ncol_use) +
  # facet_grid(cell_type ~ ., scales = "free", switch = "y") +
  xlim(-3, 3) +
  ylim(-3, 3) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-2, 2, 1)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-2, 2, 1)) +
  xlab("log2 fold change \n (Xenium)") +
  ylab("log2 fold change \n (GAGE-seq)") +
  theme_cowplot() +
  scale_color_manual(values = c("opposite" = "#CBE0F1", "same" = "#164893")) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"),  # 白色背景+白边
    strip.text = element_text(size = 16, face = "bold", color = "black"),  # 字体设置
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) 
  )


# ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq.jpg", width = 16, height = 7,
#        dpi = 600)
# ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM.jpg", width = 20, height = 9,
#        dpi = 600)
# ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM.pdf", width = 20, height = 9,
#        dpi = 600)
# ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM_onlyfemale.jpg", width = 16, height = 7,
#        dpi = 600)
ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM_major.jpg", width = 12, height = 7,
       dpi = 600)
ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM_major.pdf", width = 12, height = 7,
       dpi = 600)



library(ggplot2)
library(ggpmisc)
library(cowplot)
library(ggpubr)
ncol_use <- 7
if (if_have_major){
  ncol_use <- 4
}
res_tbl_df_twotech_merge<-res_tbl_df_twotech_merge[(abs(res_tbl_df_twotech_merge$log2FoldChange_Xenium) > 0.25) & (abs(res_tbl_df_twotech_merge$log2FoldChange_GAGEseq) > 0.25),]
res_tbl_df_twotech_merge <- res_tbl_df_twotech_merge %>%
  mutate(quadrant = case_when(
    log2FoldChange_Xenium >= 0 & log2FoldChange_GAGEseq >= 0 ~ "Q1",
    log2FoldChange_Xenium <  0 & log2FoldChange_GAGEseq >= 0 ~ "Q2",
    log2FoldChange_Xenium <  0 & log2FoldChange_GAGEseq <  0 ~ "Q3",
    log2FoldChange_Xenium >= 0 & log2FoldChange_GAGEseq <  0 ~ "Q4"
  ))

quad_counts <- res_tbl_df_twotech_merge %>%
  group_by(cell_type, quadrant) %>%
  summarise(n = n(), .groups = "drop")

quad_labels <- quad_counts %>%
  mutate(x = case_when(
    quadrant == "Q1" ~ 1.5,
    quadrant == "Q2" ~ -1.5,
    quadrant == "Q3" ~ -1.5,
    quadrant == "Q4" ~ 1.5
  ),
  y = case_when(
    quadrant == "Q1" ~ 1.5,
    quadrant == "Q2" ~ 1.5,
    quadrant == "Q3" ~ -1.5,
    quadrant == "Q4" ~ -1.5
  ),
  label = paste0("n = ", n))

ggplot(data = res_tbl_df_twotech_merge, 
       aes(x = log2FoldChange_Xenium, y = log2FoldChange_GAGEseq)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#A2A2A2", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "#A2A2A2", linewidth = 1.5) +
  geom_point(aes(color = gene_type), size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", fullrange = TRUE) +
  geom_text(data = quad_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6) +
  facet_wrap(~cell_type, scales = "free", ncol = ncol_use) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  xlab("log2 FC (Xenium)") +
  ylab("log2 FC (GAGE-seq)") +
  theme_cowplot() +
  scale_color_manual(values = c("opposite" = "#CBE0F1", "same" = "#164893")) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 16)
  )


# ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq.jpg", width = 16, height = 7,
#        dpi = 600)
ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM.jpg", width = 20, height = 9,
       dpi = 600)
ggsave("D:/Project_AD/Results/Results_pseudobulk/Comparison_log2FoldChange_Xenium_GAGEseq_rmWM.pdf", width = 20, height = 9,
       dpi = 600)