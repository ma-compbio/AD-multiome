library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)
library(ggplot2)
library(dplyr)
library(Seurat)
library(pheatmap)

path_result <- "D:/Project_AD/Results/Results_final/"
##Data preparation
######################################################################################################################################################
load("D:/Project_AD/data/Xenium_Seuratanalysis.Rdata")
data_spatial
##Obtain the cell class removing the low expression cell type
Cellclass_rmlow<-setdiff(unique(data_spatial$Cluster),c("Oligo-Lowexpr","Astro-Lowexpr","Micro-Lowexpr"))
#
##
Idents(data_spatial) <- "Cluster"
data_spatial_rmlow<-subset(data_spatial, idents = Cellclass_rmlow)

##Load the domain label
domain_label <- read.csv("D:/Project_AD/Results/Results_Standard/Domain_label.csv")
rownames(domain_label)<-domain_label$cellid
domain_label$sample_id <- data_spatial$sample_id
data_spatial_rmlow$domain_label<-domain_label[names(data_spatial_rmlow$orig.ident),]$domain_label
######################################################################################################################################################


##Cell type and subtype proportion visualization
######################################################################################################################################################
##Cell type proportion
##Calculate the cell type proportion
cluster_proportions <- data_spatial_rmlow@meta.data %>%
  group_by(sample_id, Cluster) %>%
  summarise(count = n()) %>%
  group_by(sample_id) %>%
  mutate(percentage = count / sum(count) * 100)
cluster_proportions$sample_id<-factor(cluster_proportions$sample_id,levels = c("AD-1","AD-2","CT-1","CT-2"))

##Show the cell type proportion
ggplot(cluster_proportions, aes(x = sample_id, y = percentage, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Sample ID", y = "% Cell Type", fill = "Cell Type") +
  scale_fill_manual(values = c("Astro" = "#665D48", "Endo" = "#8E6C62", "Exc L2/3 IT" = "#B4D335", "Exc L4 IT" = "#61C4B2",
                               "Exc L5 ET" = "#0C5B79", "Exc L5 IT" = "#4EB2AD", "Exc L5/6 NP" = "#3C9E64", "Exc L6 CT" = "#2D8DB9",
                               "Exc L6 IT" = "#A19935", "Exc L6 IT Car3" = "#5352A3", "Exc L6b" = "#6F4A9E",
                               "Inh Chandelier" = "#E94B9B", "Inh Lamp5" = "#DA808D", "Inh PAX6" = "#702B8B", "Inh Pvalb" = "#D93037",
                               "Inh Sncg" = "#B57BB5", "Inh Sst" = "#F8981D", "Inh Vip" = "#9A62A8",
                               "Micro" = "#95AE96", "Oligo" = "#53776C", "OPC" = "#384A46", "VLMC" = "#6B7257",
                               "Exc-Lowexpr" = "#B3CDE3", "Oligo-Lowexpr" = "#FBB4AE")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(path_result,"cluster_proportions.jpg",sep = ""), width = 6, height = 6, dpi = 300)
ggsave(paste(path_result,"cluster_proportions.pdf",sep = ""), width = 6, height = 6, dpi = 300)

##

data_spatial_rmlow$domain_label_ifGM <- ifelse(data_spatial_rmlow$domain_label == "WM", "No", "Yes")
data_spatial_rmlow_copy <- data_spatial_rmlow
# Idents(data_spatial_rmlow_copy) <-"Cluster"
data_spatial_rmlow_copy <- data_spatial_rmlow_copy[,data_spatial_rmlow_copy$Cluster != "Unk"]
Idents(data_spatial_rmlow_copy) <-"domain_label_ifGM"
data_spatial_rmlowM_GM<- data_spatial_rmlow_copy[,data_spatial_rmlow_copy$domain_label_ifGM == "Yes"]
cluster_proportions_GM <- data_spatial_rmlowM_GM@meta.data %>%
  group_by(sample_id, Cluster) %>%
  summarise(count = n()) %>%
  group_by(sample_id) %>%
  mutate(percentage = count / sum(count) * 100)
cluster_proportions_GM$sample_id<-factor(cluster_proportions_GM$sample_id,levels = c("AD-1","AD-2","CT-1","CT-2"))
##
ggplot(cluster_proportions_GM, aes(x = sample_id, y = percentage, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Sample ID", y = "% Cell Type", fill = "Cell Type") +
  scale_fill_manual(values = c("Astro" = "#665D48", "Endo" = "#8E6C62", "Exc L2/3 IT" = "#B4D335", "Exc L4 IT" = "#61C4B2",
                               "Exc L5 ET" = "#0C5B79", "Exc L5 IT" = "#4EB2AD", "Exc L5/6 NP" = "#3C9E64", "Exc L6 CT" = "#2D8DB9",
                               "Exc L6 IT" = "#A19935", "Exc L6 IT Car3" = "#5352A3", "Exc L6b" = "#6F4A9E",
                               "Inh Chandelier" = "#E94B9B", "Inh Lamp5" = "#DA808D", "Inh PAX6" = "#702B8B", "Inh Pvalb" = "#D93037",
                               "Inh Sncg" = "#B57BB5", "Inh Sst" = "#F8981D", "Inh Vip" = "#9A62A8",
                               "Micro" = "#95AE96", "Oligo" = "#53776C", "OPC" = "#384A46", "VLMC" = "#6B7257",
                               "Exc-Lowexpr" = "#B3CDE3", "Oligo-Lowexpr" = "#FBB4AE")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("% Cell Subtype")+ ##Change the name of legend title
  labs(fill = "Cell Subtype")
ggsave(paste(path_result,"cluster_proportions_GM.jpg",sep = ""), width = 6, height = 6, dpi = 300)
ggsave(paste(path_result,"cluster_proportions_GM.pdf",sep = ""), width = 6, height = 6, dpi = 300)

CV_df<-data.frame(type = character(),
                  CV = numeric())
for(type in unique(cluster_proportions_GM$Cluster)){
  cluster_proportions_GM_cur <-cluster_proportions_GM[cluster_proportions_GM$Cluster == type,]
  CV_value <- sd(cluster_proportions_GM_cur$percentage)/mean(cluster_proportions_GM_cur$percentage)
  CV_df <- rbind(CV_df, data.frame(type = type, CV = CV_value))
}
CV_df<-CV_df[CV_df$type!= "Unk",]

data_spatial_rmlowM_GM$Major_Cluster <- ifelse(data_spatial_rmlowM_GM$Cluster %in% c("Exc L2/3 IT","Exc L4 IT","Exc L5 ET","Exc L5 IT","Exc L5/6 NP","Exc L6 CT","Exc L6 IT","Exc L6 IT Car3","Exc L6b"),
                                               "Exc",ifelse(data_spatial_rmlowM_GM$Cluster %in% c("Inh Chandelier","Inh Lamp5","Inh Pvalb","Inh Sncg","Inh Sst","Inh Vip"),"Inh",data_spatial_rmlowM_GM$Cluster))
cluster_proportions_GM <- data_spatial_rmlowM_GM@meta.data %>%
  group_by(sample_id, Major_Cluster) %>%
  summarise(count = n()) %>%
  group_by(sample_id) %>%
  mutate(percentage = count / sum(count) * 100)
cluster_proportions_GM$sample_id<-factor(cluster_proportions_GM$sample_id,levels = c("AD-1","AD-2","CT-1","CT-2"))
##
ggplot(cluster_proportions_GM, aes(x = sample_id, y = percentage, fill = Major_Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Sample ID", y = "% Cell Type", fill = "Cell Type") +
  scale_fill_manual(values = c("Astro" = "#665D48", "Endo" = "#8E6C62", "Exc" = "#B4D335",
                               "Inh" = "#E94B9B", 
                               "Micro" = "#95AE96", "Oligo" = "#53776C", "OPC" = "#384A46", "VLMC" = "#6B7257",
                               "Exc-Lowexpr" = "#B3CDE3", "Oligo-Lowexpr" = "#FBB4AE")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("% Major Cell Type")+ ##Change the name of legend title
  labs(fill = "Major Cell Type")
CV_df<-data.frame(type = character(),
                  CV = numeric())
for(type in unique(cluster_proportions_GM$Major_Cluster)){
  cluster_proportions_GM_cur <-cluster_proportions_GM[cluster_proportions_GM$Major_Cluster == type,]
  CV_value <- sd(cluster_proportions_GM_cur$percentage)/mean(cluster_proportions_GM_cur$percentage)
  CV_df <- rbind(CV_df, data.frame(type = type, CV = CV_value))
}
CV_df<-CV_df[CV_df$type!= "Unk",]
ggsave(paste(path_result,"Major_Cluster_proportions_GM.jpg",sep = ""), width = 4.6, height = 6, dpi = 300)
ggsave(paste(path_result,"Major_Cluster_proportions_GM.pdf",sep = ""), width = 4.6, height = 6, dpi = 300)

######################################################################################################################################################