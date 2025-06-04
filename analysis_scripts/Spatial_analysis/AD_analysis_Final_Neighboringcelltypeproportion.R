library(Seurat)
library(ggplot2)
path_result <- "D:/Project_AD/Results/Results_final/"
##Data preparation
######################################################################################################################################################
load("D:/Project_AD/data/Xenium_Seuratanalysis.Rdata")
data_spatial

##Load the domain label
domain_label <- read.csv("D:/Project_AD/Results/Results_Standard/Domain_label.csv")
rownames(domain_label)<-domain_label$cellid
domain_label$sample_id <- data_spatial$sample_id
data_spatial$domain_label<-domain_label$domain_label

celltype_list <- unique(data_spatial$Cluster)
ntp_data<-read.csv("D:/Project_AD/data/neighboringcelltype_proportion_merge.csv",header = T,check.names = F)
ntp_data<-ntp_data[,!colnames(ntp_data) %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr")]
ntp_data[,3:ncol(ntp_data)]<-ntp_data[,3:ncol(ntp_data)]/rowSums(ntp_data[,3:ncol(ntp_data)])
mean_ntp<-apply(X = ntp_data[,3:ncol(ntp_data)],MARGIN = 2,mean,na.rm = TRUE)
data_spatial$domain_label_merge <- ifelse(data_spatial$domain_label %in% c("L1","L2","L3","L4","L5","L6"),"GM","WM")
######################################################################################################################################################

##Calculate the cell type proportion
######################################################################################################################################################
ctp_list <- list()
for (sample_id in unique(data_spatial$sample_id)){
  print(sample_id)
  data_spatial_cur<-data_spatial[,data_spatial$sample_id == sample_id]
  Idents(data_spatial_cur) <- "domain_label"
  ##
  table_use <- table(data_spatial_cur$Cluster)
  table_use<-table_use[!(names(table_use) %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr","Unk"))]
  table_use_prop <- prop.table(table_use)
  ##
  ctp_list[[as.character(sample_id)]] <- table_use_prop
}
######################################################################################################################################################

##Calculate the cell type proportion of each region
######################################################################################################################################################
ctp_list_region <- list()
for (domain_cur in sort(unique(data_spatial$domain_label))){
  print(domain_cur)
  data_spatial_cur<-data_spatial[,data_spatial$domain_label == domain_cur]
  ##
  table_use <- table(data_spatial_cur$Cluster)
  table_use<-table_use[!(names(table_use) %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr","Unk"))]
  table_use_prop <- prop.table(table_use)
  ##
  ctp_list_region[[domain_cur]] <- table_use_prop > 0.01
}
ctp_list_region$L5["Exc L5 ET"] <- TRUE
ctp_list_region$L1["Inh Chandelier"] <- TRUE
ctp_list_region$L2["Inh Chandelier"] <- TRUE
ctp_list_region$L3["Inh Chandelier"] <- TRUE
ctp_list_region$L4["Inh Chandelier"] <- TRUE
ctp_list_region$L5["Inh Chandelier"] <- TRUE
ctp_list_region$L6["Inh Chandelier"] <- TRUE
ctp_list_region$WM["Inh Chandelier"] <- TRUE

ctp_list_region_merge <- list()
for (domain_cur in sort(unique(data_spatial$domain_label_merge))){
  print(domain_cur)
  data_spatial_cur<-data_spatial[,data_spatial$domain_label_merge == domain_cur]
  ##
  table_use <- table(data_spatial_cur$Cluster)
  table_use<-table_use[!(names(table_use) %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr","Unk"))]
  table_use_prop <- prop.table(table_use)
  ##
  ctp_list_region_merge[[domain_cur]] <- table_use_prop > 0.01
}
######################################################################################################################################################

##Calcualte the neighboring cell type proportion and difference
##Regions
mean_ntp_diff_list <- list()
mean_ntp_reldiff_list <- list()
mean_ntp_AD_list <- list()
mean_ntp_CT_list <- list()
mean_ntp_test_list <- list()
for (celltype_cur in celltype_list){
  print(celltype_cur)
  data_spatial_cur<-subset(data_spatial,Cluster == celltype_cur)
  mean_ntp_diff_list0<-list()
  mean_ntp_reldiff_list0<-list()
  mean_ntp_AD_list0<-list()
  mean_ntp_CT_list0<-list()
  mean_ntp_test_list0 <-list()
  Idents(data_spatial_cur) <- "domain_label"
  for (domain_choose in unique(data_spatial_cur$domain_label)){
    if (sum(data_spatial_cur$domain_label == domain_choose)/length(data_spatial_cur$domain_label) > 0.1){
      # print(domain_choose)
      data_spatial_cur_sub <- subset(data_spatial_cur, idents = domain_choose)
      cellid_AD<-names(data_spatial_cur_sub$orig.ident)[which(data_spatial_cur_sub$condition == "AD")]
      cellid_CT<-names(data_spatial_cur_sub$orig.ident)[which(data_spatial_cur_sub$condition == "CT")]
      sampleid_AD<-(data_spatial_cur_sub$sample_id)[which(data_spatial_cur_sub$condition == "AD")]
      sampleid_CT<-(data_spatial_cur_sub$sample_id)[which(data_spatial_cur_sub$condition == "CT")]
      ntp_data_AD<-ntp_data[which(ntp_data[,2] %in% cellid_AD),]
      ntp_data_CT<-ntp_data[which(ntp_data[,2] %in% cellid_CT),]
      ##
      for (celltype_cur_temp in unique(c(names(ctp_list[['AD-1']]),names(ctp_list[['AD-2']]),names(ctp_list[['CT-1']]),names(ctp_list[['CT-2']])))){
        ntp_data_AD[sampleid_AD == AD-1,celltype_cur_temp] <- ntp_data_AD[sampleid_AD == AD-1,celltype_cur_temp]/as.numeric(ctp_list[["AD-1"]][celltype_cur_temp])
        ntp_data_AD[sampleid_AD == AD-2,celltype_cur_temp] <- ntp_data_AD[sampleid_AD == AD-2,celltype_cur_temp]/as.numeric(ctp_list[["AD-2"]][celltype_cur_temp])
        ntp_data_CT[sampleid_CT == CT-1,celltype_cur_temp] <- ntp_data_CT[sampleid_CT == CT-1,celltype_cur_temp]/as.numeric(ctp_list[["CT-1"]][celltype_cur_temp])
        ntp_data_CT[sampleid_CT == CT-2,celltype_cur_temp] <- ntp_data_CT[sampleid_CT == CT-2,celltype_cur_temp]/as.numeric(ctp_list[["CT-2"]][celltype_cur_temp])
      }
      #
      ntp_AD_mean<-colMeans(ntp_data_AD[,3:ncol(ntp_data_AD)])
      ntp_CT_mean<-colMeans(ntp_data_CT[,3:ncol(ntp_data_CT)])
      ##wilcox test
      pvalue_log10_vec <-c()
      for (neighborcelltype in colnames(ntp_data_AD)[3:ncol(ntp_data_AD)]){
        ntp_AD_mean_cur<-ntp_data_AD[,neighborcelltype]
        ntp_CT_mean_cur<-ntp_data_CT[,neighborcelltype]
        ##test
        wilcox_test_result<-wilcox.test(ntp_AD_mean_cur,ntp_CT_mean_cur,
                                        alternative = "two.sided",paired = FALSE,exact = FALSE)
        pvalue_log10<- -log10(wilcox_test_result$p.value)
        pvalue_log10_vec<-c(pvalue_log10_vec,pvalue_log10)
      }
      names(pvalue_log10_vec)<-colnames(ntp_data_AD)[3:ncol(ntp_data_AD)]
      mean_ntp_test_list0[[domain_choose]]<-pvalue_log10_vec
      ##
      # mean_ntp_diff_list0[[domain_choose]]<-ntp_AD_mean - ntp_CT_mean
      mean_ntp_diff_list0[[domain_choose]]<- log10((ntp_AD_mean + 1e-3)/(ntp_CT_mean + 1e-3))
      reldiff<-(ntp_AD_mean - ntp_CT_mean)/ntp_CT_mean
      mean_ntp_reldiff_list0[[domain_choose]]<-reldiff
      mean_ntp_AD_list0[[domain_choose]]<-ntp_AD_mean
      mean_ntp_CT_list0[[domain_choose]]<-ntp_CT_mean 
    }
  }
  mean_ntp_diff_list[[celltype_cur]]<-mean_ntp_diff_list0
  mean_ntp_reldiff_list[[celltype_cur]]<-mean_ntp_reldiff_list0
  mean_ntp_AD_list[[celltype_cur]]<-mean_ntp_AD_list0
  mean_ntp_CT_list[[celltype_cur]]<-mean_ntp_CT_list0
  mean_ntp_test_list[[celltype_cur]]<-mean_ntp_test_list0
}

##Regions (Merged: GM, WM)
mean_ntp_diff_list_merge <- list()
mean_ntp_reldiff_list_merge <- list()
mean_ntp_AD_list_merge <- list()
mean_ntp_CT_list_merge <- list()
mean_ntp_test_list_merge <- list()
for (celltype_cur in celltype_list){
  print(celltype_cur)
  data_spatial_cur<-subset(data_spatial,Cluster == celltype_cur)
  mean_ntp_diff_list0<-list()
  mean_ntp_reldiff_list0<-list()
  mean_ntp_AD_list0<-list()
  mean_ntp_CT_list0<-list()
  mean_ntp_test_list0 <-list()
  Idents(data_spatial_cur) <- "domain_label_merge"
  for (domain_choose in unique(data_spatial_cur$domain_label_merge)){
    if (sum(data_spatial_cur$domain_label_merge == domain_choose)/length(data_spatial_cur$domain_label_merge) > 0.1){
      # print(domain_choose)
      data_spatial_cur_sub <- subset(data_spatial_cur, idents = domain_choose)
      cellid_AD<-names(data_spatial_cur_sub$orig.ident)[which(data_spatial_cur_sub$condition == "AD")]
      cellid_CT<-names(data_spatial_cur_sub$orig.ident)[which(data_spatial_cur_sub$condition == "CT")]
      sampleid_AD<-(data_spatial_cur_sub$sample_id)[which(data_spatial_cur_sub$condition == "AD")]
      sampleid_CT<-(data_spatial_cur_sub$sample_id)[which(data_spatial_cur_sub$condition == "CT")]
      ntp_data_AD<-ntp_data[which(ntp_data[,2] %in% cellid_AD),]
      ntp_data_CT<-ntp_data[which(ntp_data[,2] %in% cellid_CT),]
      ##
      for (celltype_cur_temp in names(ctp_list[["AD-1"]])){
        ntp_data_AD[sampleid_AD == AD-1,celltype_cur_temp] <- ntp_data_AD[sampleid_AD == AD-1,celltype_cur_temp]/as.numeric(ctp_list[["AD-1"]][celltype_cur_temp])
        ntp_data_AD[sampleid_AD == AD-2,celltype_cur_temp] <- ntp_data_AD[sampleid_AD == AD-2,celltype_cur_temp]/as.numeric(ctp_list[["AD-2"]][celltype_cur_temp])
        ntp_data_CT[sampleid_CT == CT-1,celltype_cur_temp] <- ntp_data_CT[sampleid_CT == CT-1,celltype_cur_temp]/as.numeric(ctp_list[["CT-1"]][celltype_cur_temp])
        ntp_data_CT[sampleid_CT == CT-2,celltype_cur_temp] <- ntp_data_CT[sampleid_CT == CT-2,celltype_cur_temp]/as.numeric(ctp_list[["CT-2"]][celltype_cur_temp])
      }
      #
      ntp_AD_mean<-colMeans(ntp_data_AD[,3:ncol(ntp_data_AD)])
      ntp_CT_mean<-colMeans(ntp_data_CT[,3:ncol(ntp_data_CT)])
      ##wilcox test
      pvalue_log10_vec <-c()
      for (neighborcelltype in colnames(ntp_data_AD)[3:ncol(ntp_data_AD)]){
        ntp_AD_mean_cur<-ntp_data_AD[,neighborcelltype]
        ntp_CT_mean_cur<-ntp_data_CT[,neighborcelltype]
        ##test
        wilcox_test_result<-wilcox.test(ntp_AD_mean_cur,ntp_CT_mean_cur,
                                        alternative = "two.sided",paired = FALSE,exact = FALSE)
        pvalue_log10<- -log10(wilcox_test_result$p.value)
        pvalue_log10_vec<-c(pvalue_log10_vec,pvalue_log10)
      }
      names(pvalue_log10_vec)<-colnames(ntp_data_AD)[3:ncol(ntp_data_AD)]
      mean_ntp_test_list0[[domain_choose]]<-pvalue_log10_vec
      ##
      # mean_ntp_diff_list0[[domain_choose]]<-ntp_AD_mean - ntp_CT_mean
      mean_ntp_diff_list0[[domain_choose]]<- log10((ntp_AD_mean + 1e-3)/(ntp_CT_mean + 1e-3))
      reldiff<-(ntp_AD_mean - ntp_CT_mean)/ntp_CT_mean
      mean_ntp_reldiff_list0[[domain_choose]]<-reldiff
      mean_ntp_AD_list0[[domain_choose]]<-ntp_AD_mean
      mean_ntp_CT_list0[[domain_choose]]<-ntp_CT_mean 
    }
  }
  mean_ntp_diff_list_merge[[celltype_cur]]<-mean_ntp_diff_list0
  mean_ntp_reldiff_list_merge[[celltype_cur]]<-mean_ntp_reldiff_list0
  mean_ntp_AD_list_merge[[celltype_cur]]<-mean_ntp_AD_list0
  mean_ntp_CT_list_merge[[celltype_cur]]<-mean_ntp_CT_list0
  mean_ntp_test_list_merge[[celltype_cur]]<-mean_ntp_test_list0
}

##
mean_ntp_diff_df<-data.frame("celltype"=character(),neighcelltype = character(),"domain"=character(),"mean_ntp_diff"=numeric(),"mean_npt_AD" = numeric(),"mean_npt_CT" = numeric())
mean_ntp_reldiff_df<-data.frame("celltype"=character(),neighcelltype = character(),"domain"=character(),"mean_ntp_reldiff"=numeric())
mean_ntp_test_list_df <-data.frame("celltype"=character(),neighcelltype = character(),"domain"=character(),"mean_ntp_test"=numeric())
for (celltype_cur in celltype_list){
  mean_ntp_diff_list0<-mean_ntp_diff_list[[celltype_cur]]
  mean_ntp_reldiff_list0<-mean_ntp_reldiff_list[[celltype_cur]]
  for (domain_cur in names(mean_ntp_diff_list0)){
    mean_ntp_diff_df<-rbind(mean_ntp_diff_df,
                            data.frame("celltype"=celltype_cur,
                                       "neighcelltype"=names(mean_ntp_diff_list0[[domain_cur]]),
                                       "domain"=domain_cur,
                                       "mean_ntp_diff"=mean_ntp_diff_list0[[domain_cur]],
                                       "mean_npt_AD"=mean_ntp_AD_list[[celltype_cur]][[domain_cur]],
                                       "mean_npt_CT"=mean_ntp_CT_list[[celltype_cur]][[domain_cur]]
                            ))
    mean_ntp_reldiff_df<-rbind(mean_ntp_reldiff_df,data.frame("celltype"=celltype_cur,"neighcelltype"=names(mean_ntp_reldiff_list0[[domain_cur]]),"domain"=domain_cur,"mean_ntp_reldiff"=mean_ntp_reldiff_list0[[domain_cur]]))
    mean_ntp_test_list_df <- rbind(mean_ntp_test_list_df,
                                   data.frame("celltype"=celltype_cur,
                                              "neighcelltype"=names(mean_ntp_reldiff_list0[[domain_cur]]),
                                              "domain"=domain_cur,
                                              "mean_ntp_test"=mean_ntp_test_list[[celltype_cur]][[domain_cur]]))
    
  }
}
mean_ntp_diff_merge<-merge(mean_ntp_diff_df,mean_ntp_reldiff_df,by = c("celltype","neighcelltype","domain"))
mean_ntp_diff_merge$mean_npt_conditions<-(mean_ntp_diff_merge$mean_npt_AD + mean_ntp_diff_merge$mean_npt_CT)/2
mean_ntp_diff_merge<-merge(mean_ntp_diff_merge,mean_ntp_test_list_df,by = c("celltype","neighcelltype","domain"))
mean_ntp_diff_merge$mean_ntp_test[is.infinite(mean_ntp_diff_merge$mean_ntp_test)]<-5
mean_ntp_diff_merge$mean_ntp_test[mean_ntp_diff_merge$mean_ntp_test < -log10(0.05)]<-NA
mean_ntp_diff_merge$mean_ntp_test[mean_ntp_diff_merge$mean_ntp_test > 5] <- 5
##
mean_ntp_diff_df_merge<-data.frame("celltype"=character(),neighcelltype = character(),"domain"=character(),"mean_ntp_diff"=numeric(),"mean_npt_AD" = numeric(),"mean_npt_CT" = numeric())
mean_ntp_reldiff_df_merge<-data.frame("celltype"=character(),neighcelltype = character(),"domain"=character(),"mean_ntp_reldiff"=numeric())
mean_ntp_test_list_df_merge <-data.frame("celltype"=character(),neighcelltype = character(),"domain"=character(),"mean_ntp_test"=numeric())
for (celltype_cur in celltype_list){
  mean_ntp_diff_list0<-mean_ntp_diff_list_merge[[celltype_cur]]
  mean_ntp_reldiff_list0<-mean_ntp_reldiff_list_merge[[celltype_cur]]
  for (domain_cur in names(mean_ntp_diff_list0)){
    mean_ntp_diff_df_merge<-rbind(mean_ntp_diff_df_merge,
                                  data.frame("celltype"=celltype_cur,
                                             "neighcelltype"=names(mean_ntp_diff_list0[[domain_cur]]),
                                             "domain"=domain_cur,
                                             "mean_ntp_diff"=mean_ntp_diff_list0[[domain_cur]],
                                             "mean_npt_AD"=mean_ntp_AD_list_merge[[celltype_cur]][[domain_cur]],
                                             "mean_npt_CT"=mean_ntp_CT_list_merge[[celltype_cur]][[domain_cur]]
                                  ))
    mean_ntp_reldiff_df_merge<-rbind(mean_ntp_reldiff_df_merge,data.frame("celltype"=celltype_cur,"neighcelltype"=names(mean_ntp_reldiff_list0[[domain_cur]]),"domain"=domain_cur,"mean_ntp_reldiff"=mean_ntp_reldiff_list0[[domain_cur]]))
    mean_ntp_test_list_df_merge <- rbind(mean_ntp_test_list_df_merge,
                                         data.frame("celltype"=celltype_cur,
                                                    "neighcelltype"=names(mean_ntp_reldiff_list0[[domain_cur]]),
                                                    "domain"=domain_cur,
                                                    "mean_ntp_test"=mean_ntp_test_list_merge[[celltype_cur]][[domain_cur]]))
    
  }
}
mean_ntp_diff_merge2<-merge(mean_ntp_diff_df_merge,mean_ntp_reldiff_df_merge,by = c("celltype","neighcelltype","domain"))
mean_ntp_diff_merge2$mean_npt_conditions<-(mean_ntp_diff_merge2$mean_npt_AD + mean_ntp_diff_merge2$mean_npt_CT)/2
mean_ntp_diff_merge2<-merge(mean_ntp_diff_merge2,mean_ntp_test_list_df_merge,by = c("celltype","neighcelltype","domain"))
mean_ntp_diff_merge2$mean_ntp_test[is.infinite(mean_ntp_diff_merge2$mean_ntp_test)]<-5
mean_ntp_diff_merge2$mean_ntp_test[mean_ntp_diff_merge2$mean_ntp_test < -log10(0.05)]<-NA
mean_ntp_diff_merge2$mean_ntp_test[mean_ntp_diff_merge2$mean_ntp_test > 5] <- 5
##

mean_ntp_diff_merge$mean_ntp_reldiff_abs <- abs(mean_ntp_diff_merge$mean_ntp_reldiff)
mean_ntp_diff_merge$mean_ntp_diff_abs <- abs(mean_ntp_diff_merge$mean_ntp_diff)
mean_ntp_diff_merge$mean_ntp_reldiff_abs_log<-log1p(mean_ntp_diff_merge$mean_ntp_reldiff_abs)
mean_ntp_diff_merge$mean_ntp_reldiff_log<-log1p(mean_ntp_diff_merge$mean_ntp_reldiff_abs) * sign(mean_ntp_diff_merge$mean_ntp_reldiff)
# mean_ntp_diff_merge$mean_ntp_reldiff_log<-sqrt(mean_ntp_diff_merge$mean_ntp_reldiff_abs) * sign(mean_ntp_diff_merge$mean_ntp_reldiff)
mean_ntp_diff_merge$mean_ntp_diff_abs[is.infinite(mean_ntp_diff_merge$mean_ntp_reldiff_log)]<-NaN
mean_ntp_diff_merge$mean_ntp_reldiff_log[is.infinite(mean_ntp_diff_merge$mean_ntp_reldiff_log)]<-NaN
# mean_ntp_diff_merge<-mean_ntp_diff_merge[!is.na(mean_ntp_diff_merge$mean_ntp_reldiff_log),]
celltype_order<-sort(unique(mean_ntp_diff_merge$celltype))
mean_ntp_diff_merge$celltype<-factor(mean_ntp_diff_merge$celltype,
                                     levels = celltype_order)
mean_ntp_diff_merge$neighcelltype<-factor(mean_ntp_diff_merge$neighcelltype,
                                          levels = rev(celltype_order))
mean_ntp_diff_merge<-mean_ntp_diff_merge[!(mean_ntp_diff_merge$celltype %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr")) &
                                           !(mean_ntp_diff_merge$neighcelltype %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr")),]
mean_ntp_diff_merge$mean_ntp_diff_change<-mean_ntp_diff_merge$mean_ntp_diff
# mean_ntp_diff_merge$mean_ntp_diff_change[mean_ntp_diff_merge$mean_ntp_diff_change > 0.2]<-0.2
# mean_ntp_diff_merge$mean_ntp_diff_change[mean_ntp_diff_merge$mean_ntp_diff_change < -0.2]<--0.2
# mean_ntp_diff_merge$mean_ntp_diff_change[mean_ntp_diff_merge$mean_ntp_diff_change > 0.1]<-0.1
# mean_ntp_diff_merge$mean_ntp_diff_change[mean_ntp_diff_merge$mean_ntp_diff_change < -0.1]<--0.1
mean_ntp_diff_merge$mean_ntp_reldiff_change<-mean_ntp_diff_merge$mean_ntp_reldiff
# mean_ntp_diff_merge$mean_ntp_reldiff_change[mean_ntp_diff_merge$mean_ntp_reldiff_change > 1]<-1
# mean_ntp_diff_merge$mean_ntp_reldiff_change[mean_ntp_diff_merge$mean_ntp_reldiff_change < -1]<--1
# mean_ntp_diff_merge$mean_ntp_reldiff_change<-abs(mean_ntp_diff_merge$mean_ntp_reldiff_change)
mean_ntp_diff_merge$mean_ntp_reldiff_change_abs <- abs(mean_ntp_diff_merge$mean_ntp_reldiff_change)
# mean_ntp_diff_merge$mean_ntp_diff<-sqrt(abs(mean_ntp_diff_merge$mean_ntp_diff)) * sign(mean_ntp_diff_merge$mean_ntp_diff)
mean_ntp_diff_merge$mean_npt_conditions_norm<-mean_ntp_diff_merge$mean_npt_conditions/mean_ntp[as.character(mean_ntp_diff_merge$neighcelltype)]
mean_ntp_diff_merge$mean_npt_conditions_norm<-log1p(mean_ntp_diff_merge$mean_npt_conditions_norm)
mean_ntp_diff_merge<- mean_ntp_diff_merge[!(is.na(mean_ntp_diff_merge$mean_ntp_diff)),]
mean_ntp_diff_merge <- mean_ntp_diff_merge[!((mean_ntp_diff_merge$celltype == "Unk") | (mean_ntp_diff_merge$neighcelltype == "Unk")),]
##filter by ctp_list_region
# mean_ntp_diff_merge
mean_ntp_diff_merge_useindex_all <- c()
for (domain_cur in unique(data_spatial$domain_label)){
  ctp_list_region_cur <- ctp_list_region[[domain_cur]]
  ctp_list_region_cur_true <- names(ctp_list_region_cur[which(ctp_list_region_cur)])
  mean_ntp_diff_merge_useindex <-  which((mean_ntp_diff_merge$celltype %in% ctp_list_region_cur_true) & (mean_ntp_diff_merge$neighcelltype %in% ctp_list_region_cur_true) & (mean_ntp_diff_merge$domain == domain_cur))
  mean_ntp_diff_merge_useindex_all <- c(mean_ntp_diff_merge_useindex_all,mean_ntp_diff_merge_useindex)
}
mean_ntp_diff_merge <- mean_ntp_diff_merge[mean_ntp_diff_merge_useindex_all,]
##
ggplot(data = mean_ntp_diff_merge,aes(x = domain, y = neighcelltype))+
  geom_point(aes(size = mean_npt_conditions,color = mean_ntp_diff_change))+
  # geom_point(aes(size = mean_npt_conditions_norm,color = mean_npt_conditions_norm))+
  # geom_point(aes(size = mean_ntp_reldiff_change_abs,color = mean_ntp_reldiff_change))+
  theme_bw()+
  facet_grid(~celltype)+
  # facet_wrap(~celltype,scales = "free",nrow = 2)+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 60,hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  scale_colour_gradient2(low = "#702B8B",
                         mid = "white",
                         high = "#FB8661",
                         midpoint = 0,
                         # limits = c(-0.2,0.2))+ ##increase the size of facet
                         limits = c(-0.1,0.1))+ ##increase the size of facet
  theme(strip.text.x = element_text(size = 12))+
  ylab("Neighboring cell type proportion")+
  xlab("Domain")
ggsave(paste(path_result,"mean_ntp_diff_merge_facet.jpg",sep = ""),width = 30,height = 5,dpi = 300)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet.pdf",sep = ""),width = 35,height = 5,dpi = 1000)


ggplot(data = mean_ntp_diff_merge,aes(x = domain, y = neighcelltype))+
  geom_point(aes(size = mean_ntp_test,color = mean_ntp_diff_change))+
  # geom_point(aes(size = mean_npt_conditions_norm,color = mean_npt_conditions_norm))+
  # geom_point(aes(size = mean_ntp_reldiff_change_abs,color = mean_ntp_reldiff_change))+
  theme_bw()+
  facet_grid(~celltype)+
  # facet_wrap(~celltype,scales = "free",nrow = 2)+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 60,hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  scale_colour_gradient2(low = "#702B8B",
                         mid = "white",
                         high = "#FB8661",
                         midpoint = 0,
                         # limits = c(-0.2,0.2))+ ##increase the size of facet
                         limits = c(-0.1,0.1))+ ##increase the size of facet
  theme(strip.text.x = element_text(size = 12))+
  ylab("Neighboring cell type proportion")+
  xlab("Domain")

ggsave(paste(path_result,"mean_ntp_diff_merge_facet_region.jpg",sep = ""),width = 30,height = 5,dpi = 300)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_region.pdf",sep = ""),width = 35,height = 5,dpi = 1000)
##

##
##Save the mean_ntp_diff_merge to csv
write.csv(mean_ntp_diff_merge, paste(path_result,"mean_ntp_diff_merge.csv",sep = ""),row.names = F)
##
ggplot(data = mean_ntp_diff_merge,aes(x = celltype, y = neighcelltype))+
  geom_point(aes(size = mean_npt_conditions,color = mean_ntp_diff_change))+
  # geom_point(aes(size = mean_npt_conditions_norm,color = mean_npt_conditions_norm))+
  # geom_point(aes(size = mean_ntp_reldiff_change_abs,color = mean_ntp_reldiff_change))+
  theme_bw()+
  facet_grid(~domain)+
  # facet_wrap(~celltype,scales = "free",nrow = 2)+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 60,hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+ ##increase the text size of the facet
  scale_colour_gradient2(low = "#702B8B",
                         mid = "white",
                         high = "#FB8661",
                         midpoint = 0,
                         # limits = c(-0.2,0.2))+ ##increase the size of facet
                         limits = c(-0.1,0.1))+ ##increase the size of facet
  theme(strip.text.x = element_text(size = 12))+
  ylab("Neighboring cell type")+
  xlab("Cell type")
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2.jpg",sep = ""),width = 30,height = 5,dpi = 300)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2.pdf",sep = ""),width = 35,height = 10,dpi = 1000)

ggplot(data = mean_ntp_diff_merge,aes(x = celltype, y = neighcelltype))+
  geom_point(aes(size = mean_ntp_test,color = mean_ntp_diff))+
  # geom_point(aes(size = mean_npt_conditions_norm,color = mean_npt_conditions_norm))+
  # geom_point(aes(size = mean_ntp_reldiff_change_abs,color = mean_ntp_reldiff_change))+
  theme_bw()+
  facet_grid(~domain)+
  # facet_wrap(~celltype,scales = "free",nrow = 2)+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 60,hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+ ##increase the text size of the facet
  scale_colour_gradient2(low = "#702B8B",
                         mid = "white",
                         high = "#FB8661",
                         midpoint = 0)+
  # limits = c(-0.2,0.2))+ ##increase the size of facet
  # limits = c(-0.1,0.1))+ ##increase the size of facet
  theme(strip.text.x = element_text(size = 12))+
  ylab("Neighboring cell type")+
  xlab("Cell type")
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2_1.jpg",sep = ""),width = 30,height = 5,dpi = 300)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2_1.pdf",sep = ""),width = 35,height = 10,dpi = 1000)

write.csv(mean_ntp_diff_merge,paste(path_result,"mean_ntp_diff_region.csv",sep = ""),row.names = F)

ggplot(data = mean_ntp_diff_merge,aes(x = celltype, y = neighcelltype))+
  geom_point(aes(size = mean_ntp_test,color = mean_ntp_diff))+
  # geom_point(aes(size = mean_npt_conditions_norm,color = mean_npt_conditions_norm))+
  # geom_point(aes(size = mean_ntp_reldiff_change_abs,color = mean_ntp_reldiff_change))+
  theme_bw()+
  # facet_grid(~domain,scales = "free_x",space = "fixed")+
  facet_wrap(~domain, scales = "free_x", nrow = 1)+
  # facet_wrap(~celltype,scales = "free",nrow = 2)+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 60,hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 20))+ ##increase the text size of the facet
  scale_colour_gradient2(low = "#702B8B",
                         mid = "white",
                         high = "#FB8661",
                         midpoint = 0)+
  # limits = c(-0.2,0.2))+ ##increase the size of facet
  # limits = c(-0.1,0.1))+ ##increase the size of facet
  theme(strip.text.x = element_text(size = 12))+
  ylab("Neighboring cell type")+
  xlab("Cell type")

write.csv(mean_ntp_diff_merge,paste(path_result,"mean_ntp_diff_region.csv",sep = ""),row.names = F)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2_1_1.jpg",sep = ""),width = 30,height = 6,dpi = 300)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2_1_1.jpg",sep = ""),width = 30,height = 6,dpi = 1000)




mean_ntp_diff_merge2$mean_ntp_reldiff_abs <- abs(mean_ntp_diff_merge2$mean_ntp_reldiff)
mean_ntp_diff_merge2$mean_ntp_diff_abs <- abs(mean_ntp_diff_merge2$mean_ntp_diff)
mean_ntp_diff_merge2$mean_ntp_reldiff_abs_log<-log1p(mean_ntp_diff_merge2$mean_ntp_reldiff_abs)
mean_ntp_diff_merge2$mean_ntp_reldiff_log<-log1p(mean_ntp_diff_merge2$mean_ntp_reldiff_abs) * sign(mean_ntp_diff_merge2$mean_ntp_reldiff)
# mean_ntp_diff_merge2$mean_ntp_reldiff_log<-sqrt(mean_ntp_diff_merge2$mean_ntp_reldiff_abs) * sign(mean_ntp_diff_merge2$mean_ntp_reldiff)
mean_ntp_diff_merge2$mean_ntp_diff_abs[is.infinite(mean_ntp_diff_merge2$mean_ntp_reldiff_log)]<-NaN
mean_ntp_diff_merge2$mean_ntp_reldiff_log[is.infinite(mean_ntp_diff_merge2$mean_ntp_reldiff_log)]<-NaN
# mean_ntp_diff_merge2<-mean_ntp_diff_merge2[!is.na(mean_ntp_diff_merge2$mean_ntp_reldiff_log),]
celltype_order<-sort(unique(mean_ntp_diff_merge2$celltype))
mean_ntp_diff_merge2$celltype<-factor(mean_ntp_diff_merge2$celltype,
                                      levels = celltype_order)
mean_ntp_diff_merge2$neighcelltype<-factor(mean_ntp_diff_merge2$neighcelltype,
                                           levels = rev(celltype_order))
mean_ntp_diff_merge2<-mean_ntp_diff_merge2[!(mean_ntp_diff_merge2$celltype %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr")) &
                                             !(mean_ntp_diff_merge2$neighcelltype %in% c("Astro-Lowexpr","Oligo-Lowexpr","Micro-Lowexpr")),]
mean_ntp_diff_merge2$mean_ntp_diff_change<-mean_ntp_diff_merge2$mean_ntp_diff
# mean_ntp_diff_merge2$mean_ntp_diff_change[mean_ntp_diff_merge2$mean_ntp_diff_change > 0.2]<-0.2
# mean_ntp_diff_merge2$mean_ntp_diff_change[mean_ntp_diff_merge2$mean_ntp_diff_change < -0.2]<--0.2
# mean_ntp_diff_merge2$mean_ntp_diff_change[mean_ntp_diff_merge2$mean_ntp_diff_change > 0.1]<-0.1
# mean_ntp_diff_merge2$mean_ntp_diff_change[mean_ntp_diff_merge2$mean_ntp_diff_change < -0.1]<--0.1
mean_ntp_diff_merge2$mean_ntp_reldiff_change<-mean_ntp_diff_merge2$mean_ntp_reldiff
# mean_ntp_diff_merge2$mean_ntp_reldiff_change[mean_ntp_diff_merge2$mean_ntp_reldiff_change > 1]<-1
# mean_ntp_diff_merge2$mean_ntp_reldiff_change[mean_ntp_diff_merge2$mean_ntp_reldiff_change < -1]<--1
# mean_ntp_diff_merge2$mean_ntp_reldiff_change<-abs(mean_ntp_diff_merge2$mean_ntp_reldiff_change)
mean_ntp_diff_merge2$mean_ntp_reldiff_change_abs <- abs(mean_ntp_diff_merge2$mean_ntp_reldiff_change)
# mean_ntp_diff_merge2$mean_ntp_diff<-sqrt(abs(mean_ntp_diff_merge2$mean_ntp_diff)) * sign(mean_ntp_diff_merge2$mean_ntp_diff)
mean_ntp_diff_merge2$mean_npt_conditions_norm<-mean_ntp_diff_merge2$mean_npt_conditions/mean_ntp[as.character(mean_ntp_diff_merge2$neighcelltype)]
mean_ntp_diff_merge2$mean_npt_conditions_norm<-log1p(mean_ntp_diff_merge2$mean_npt_conditions_norm)
mean_ntp_diff_merge2<- mean_ntp_diff_merge2[!(is.na(mean_ntp_diff_merge2$mean_ntp_diff)),]
mean_ntp_diff_merge2 <- mean_ntp_diff_merge2[!((mean_ntp_diff_merge2$celltype == "Unk") | (mean_ntp_diff_merge2$neighcelltype == "Unk")),]
##filter by ctp_list_region
# mean_ntp_diff_merge2
mean_ntp_diff_merge2_useindex_all <- c()
for (domain_cur in unique(data_spatial$domain_label_merge)){
  ctp_list_region_cur <- ctp_list_region_merge[[domain_cur]]
  ctp_list_region_cur_true <- names(ctp_list_region_cur[which(ctp_list_region_cur)])
  mean_ntp_diff_merge2_useindex <-  which((mean_ntp_diff_merge2$celltype %in% ctp_list_region_cur_true) & (mean_ntp_diff_merge2$neighcelltype %in% ctp_list_region_cur_true) & (mean_ntp_diff_merge2$domain == domain_cur))
  mean_ntp_diff_merge2_useindex_all <- c(mean_ntp_diff_merge2_useindex_all,mean_ntp_diff_merge2_useindex)
}
mean_ntp_diff_merge2 <- mean_ntp_diff_merge2[mean_ntp_diff_merge2_useindex_all,]

ggplot(data = mean_ntp_diff_merge2,aes(x = celltype, y = neighcelltype))+
  geom_point(aes(size = mean_ntp_test,color = mean_ntp_diff))+
  # geom_point(aes(size = mean_npt_conditions_norm,color = mean_npt_conditions_norm))+
  # geom_point(aes(size = mean_ntp_reldiff_change_abs,color = mean_ntp_reldiff_change))+
  theme_bw()+
  facet_grid(~domain)+
  # facet_wrap(~celltype,scales = "free",nrow = 2)+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 60,hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+ ##increase the text size of the facet
  scale_colour_gradient2(low = "#702B8B",
                         mid = "white",
                         high = "#FB8661",
                         midpoint = 0)+
  # limits = c(-0.2,0.2))+ ##increase the size of facet
  # limits = c(-0.1,0.1))+ ##increase the size of facet
  theme(strip.text.x = element_text(size = 12))+
  ylab("Neighboring cell subtype")+
  xlab("Cell subtype")
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2_2.jpg",sep = ""),width = 15,height = 5,dpi = 300)
ggsave(paste(path_result,"mean_ntp_diff_merge_facet_V2_2.pdf",sep = ""),width = 15,height = 10,dpi = 1000)
write.csv(mean_ntp_diff_merge,paste(path_result,"mean_ntp_diff_region.csv",sep = ""),row.names = F)
