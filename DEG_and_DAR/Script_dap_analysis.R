setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference")
library(Seurat)
library(Signac)

tIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/CIO_and_HIO/tHIO_vs_tCIO_epi_ATAC/redo_after_adding_fragment_file/redo_dap_after_refining_cell_type_annotation/Res_tIO_atac_with_unified_cell_type_annotation.rds")
SCpubr::do_DimPlot(tIO_atac, group.by = "High_resolution_cell_type")
saveRDS(tIO_atac, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_ATAC_seurat_object.rds")

# transfer the tIO_rna regional score to the tIO_atac cells
HIO_CIO_meta_cell <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/CIO_and_HIO/tHIO_vs_tCIO_epi_ATAC/redo_after_adding_fragment_file/Res_CIO_HIO_merged_meta_cell_info.rds")
saveRDS(HIO_CIO_meta_cell, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_HIO_CIO_meta_cell.rds")

mapped_region_prob_mat <- tIO_epi@meta.data[,grep("Intestine_region_L2:prob", colnames(tIO_epi@meta.data))]
aa <- apply(mapped_region_prob_mat, 2, as.numeric)
rownames(aa) <- rownames(mapped_region_prob_mat)
mapped_region_prob_mat <- aa
mat <- HIO_CIO_meta_cell[which(HIO_CIO_meta_cell$RNA_cell_id2%in%colnames(tIO_epi) & HIO_CIO_meta_cell$ATAC_cell%in%colnames(tIO_atac)),c("RNA_cell_id2","ATAC_cell")]
atac_cell_score <- t(sapply(unique(mat$ATAC_cell), function(atac_cell){
  rna_cells <- mat$RNA_cell_id2[which(mat$ATAC_cell==atac_cell)]
  if(length(rna_cells)>1){
    value <- colMeans(mapped_region_prob_mat[rna_cells,])
    value <- value/sum(value)
  }else{
    value <- mapped_region_prob_mat[rna_cells,]
  }
  return(value)
}))

saveRDS(atac_cell_score, file="Dat_atac_cell_regional_score_from_matched_RNA_cells.rds")

all_atac_cell_score <- matrix(NA, nrow=ncol(tIO_atac), ncol=3)
rownames(all_atac_cell_score) <- colnames(tIO_atac)
colnames(all_atac_cell_score) <- colnames(atac_cell_score)
all_atac_cell_score[rownames(atac_cell_score),] <- atac_cell_score
#seu_obj_list <- SplitObject(tIO_atac, split.by = "Sample.name")
for(x in names(seu_obj_list)){
  seu_obj <- seu_obj_list[[x]]
  #det_rates <- rowMeans(seu_obj@assays$peaks_species@counts > 0) 
  #top <- names(which(det_rates>0.05 & det_rates<0.9))
  #seu_obj <- FindTopFeatures(seu_obj, min.cutoff = ncol(seu_obj) * 0.05)
  #seu_obj <- RunSVD(seu_obj)
  #seu_obj <- RunUMAP(seu_obj, reduction="lsi", dims = 2:20)
  #seu_obj <- FindNeighbors(seu_obj, reduction="lsi", dims=2:20) %>% FindClusters(resolution=2)
  #seu_obj_list[[x]] <- seu_obj
  
  nn_mat <- seu_obj@graphs$peaks_species_nn
  que_cells <- setdiff(rownames(nn_mat), rownames(atac_cell_score))
  # take the mean score of the 20 NN as the score 
  que_score_mat <- t(sapply(que_cells, function(i){
    vec <- as.vector(as.matrix(nn_mat[i,]))
    nn_cells <- colnames(nn_mat)[which(vec==1)]
    ref_cells <- intersect(rownames(atac_cell_score),nn_cells)
    if(length(ref_cells)<1){
      cl <- seu_obj@meta.data[i,"peaks_species_snn_res.2"]
      
      ref_cells <- intersect(rownames(atac_cell_score),
                             colnames(seu_obj)[which(seu_obj$peaks_species_snn_res.2==cl)])
      if(length(ref_cells)<1){
        return(colMeans(atac_cell_score))
      }else if(length(ref_cells)==1){
        return(atac_cell_score[ref_cells,])
      }else{
        return(colMeans(atac_cell_score[ref_cells,]))
      }
    }else if(length(ref_cells)>1){
      return(colMeans(atac_cell_score[ref_cells,]))
    }else if(length(ref_cells)==1){
      return(atac_cell_score[ref_cells,])
    }
  }))
  all_atac_cell_score[que_cells,] <- que_score_mat
}
colnames(all_atac_cell_score) <- sub(pattern = "Intestine_region_L2:prob.",
                                     replacement = "",
                                     colnames(all_atac_cell_score))
saveRDS(all_atac_cell_score, file="Dat_all_atac_cell_regional_score.rds")

tIO_atac@meta.data[,colnames(all_atac_cell_score)] <- all_atac_cell_score
tIO_atac@meta.data[,"Intestine_region_L2:pred.id"] <- colnames(all_atac_cell_score)[apply(all_atac_cell_score,1,which.max)]
saveRDS(tIO_atac, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_ATAC_seurat_object.rds")

n1 <- sapply(sort(unique(tIO_atac$`Intestine_region_L2:pred.id`)), function(tissue){
  sapply(sort(unique(tIO_atac$Sample.name)), function(sample){
    sum(tIO_atac$Sample.name==sample & tIO_atac$`Intestine_region_L2:pred.id`==tissue)
  })
})

tissue.cols <- setNames(c("#d9d9d9", "#696969","#303030"), c("Colon", "Ileum", "Prox_SI"))
row_num=1
col_num=4
pdf("Plot_barplot_tIO_epi_scATAC_inferred_regional_identity.pdf", height=5*row_num, width=5*col_num)
par(mfrow=c(row_num,col_num), mar=c(12,5,5,5))
for(x in unique(tIO_atac$Sample.name)){
  barplot(res[[ct]], main=ct, las=2, ylab="Cell number", col = tissue.cols[rownames(res[[ct]])])
  legend("topleft", legend=names(tissue.cols), fill=tissue.cols, bty="n")
}
dev.off()

tIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_ATAC_seurat_object.rds")
tIO_epi <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
p1 <- SCpubr::do_ViolinPlot(tIO_atac, features = "Intestine_region_L2:prob.Prox_SI", group.by = "Sample.name", ylab="Prob. Prox_SI")+ggtitle(label="ATAC")
p2 <- SCpubr::do_ViolinPlot(tIO_atac, features = "Intestine_region_L2:prob.Colon", group.by = "Sample.name", ylab="Prob. Colon")+ggtitle(label="ATAC")
p3 <- SCpubr::do_ViolinPlot(tIO_epi, features = "Intestine_region_L2:prob.Prox_SI", group.by = "Sample.name", ylab="Prob. Prox_SI")+ggtitle(label="RNA")
p4 <- SCpubr::do_ViolinPlot(tIO_epi, features = "Intestine_region_L2:prob.Colon", group.by = "Sample.name", ylab="Prob. Colon")+ggtitle(label="RNA")
pdf("Plot_violinPlot_tIO_epi_RNA_and_ATAC_regional_score_distribution.pdf")
p3+p4+p1+p2+patchwork::plot_layout(ncol=2)
dev.off()



# identify species DAP
tIO_atac$High_resolution_cell_type_per_species <- paste(tIO_atac$High_resolution_cell_type, tIO_atac$Species, sep="@")
sort(table(tIO_atac$High_resolution_cell_type_per_species))
saveRDS(tIO_atac, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_ATAC_seurat_object.rds")

selected_cell_type <- setdiff(sort(unique(tIO_atac$High_resolution_cell_type)), "Putative_secretory_progenitors") # putative secretory progenitors are only identified in human, and on the integrated UMAP embedding, putative secretory cells do not collocalize with human cells
## get cell type per species detection rate
source("~/Work/commonScript/Script_other_functions.R")
cell_type_det_rate <- getExpressedProp(seu.obj=tIO_atac, feature.to.calc = "High_resolution_cell_type_per_species", colname.prefix=NULL, assay.type="peaks_species")
saveRDS(cell_type_det_rate, file="Dat_tIO_epi_ATAC_hg38_coor_high_resolution_cell_type_per_species_detection_rate.rds")

## run test per cell type
library(doParallel)
ddp_by_cell_type <- list()
for(x in selected_cell_type){
  print(paste(x, "start"))
  
  # get detected peaks
  rate_mat <- cell_type_det_rate[,paste(x,c("Human","Chimp"),sep="@")]
  expressed_peaks <- unique(unlist(lapply(seq(ncol(rate_mat)), function(j){
    vec <- rate_mat[,j]
    peaks1 <- rownames(cell_type_det_rate)[which(vec>0.01)]
    num <- min(c(150000, length(peaks1)))
    peaks1[order(rate_mat[peaks1,j], decreasing = T)[1:num]]
  })))
  
  
  idx <- which(tIO_atac$High_resolution_cell_type==x)
  combined_count <- as.matrix(tIO_atac@assays$peaks_species@counts[expressed_peaks,idx])
  count_cutoff <- 0
  bi_mat <- combined_count>count_cutoff
  species_vec <- as.factor(tIO_atac$Species[idx])
  umi_vec <- tIO_atac$nCount_peaks_species[idx]
  colon_score_vec <- tIO_atac$`Intestine_region_L2:prob.Colon`[idx]
  proxSI_score_vec <- tIO_atac$`Intestine_region_L2:prob.Prox_SI`[idx]
  registerDoParallel(20)
  species_dd_test_res <- foreach(k=seq(nrow(bi_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.numeric(as.vector(bi_mat[k,]))
    m0 <- glm(e ~ umi_vec+colon_score_vec+proxSI_score_vec, family = binomial)
    m1 <- glm(e ~ umi_vec+colon_score_vec+proxSI_score_vec+species_vec, family = binomial)
    a0 <- anova(m0)
    a1 <- anova(m1)
    p_anova <- anova(m1,m0, test = "Chisq")$Pr[2]
    p_resi <- pf((a0[nrow(a0),"Resid. Dev"]/a0[nrow(a0),"Resid. Df"]) / (a1[nrow(a1),"Resid. Dev"]/a1[nrow(a1),"Resid. Df"]), 
                 df1 = a0[nrow(a0),"Resid. Df"], df2 = a1[nrow(a1),"Resid. Df"], lower.tail = F)
    coef <- coef(m1)[length(coef(m1))]
    return(c(p_anova, p_resi ,coef))
  }
  stopImplicitCluster()
  rownames(species_dd_test_res) <- expressed_peaks
  colnames(species_dd_test_res) <- c("p_ANOVA", "p_Resi", "Coef")
  resi_p_adj <- p.adjust(species_dd_test_res[,2], method = "BH")
  anova_p_adj <- p.adjust(species_dd_test_res[,1], method = "BH")
  species_dd_test_res <- as.data.frame(species_dd_test_res)
  
  expressed_prop <- sapply(sort(unique(species_vec)), function(x){
    idx <- which(species_vec==x)
    prop_vec <- rowMeans(bi_mat[,idx])
    return(prop_vec)
  })
  colnames(expressed_prop) <- sort(unique(species_vec))
  expressed_prop <- as.data.frame(expressed_prop)
  expressed_prop$Prop_diff <- expressed_prop$Human-expressed_prop$Chimp
  res_mat <- data.frame("Species_test_anova_P"=species_dd_test_res$p_ANOVA,
                        "Species_test_resi_P"=species_dd_test_res$p_Resi,
                        "Corrected_anova_P"=anova_p_adj,
                        "Corrected_resi_P"=resi_p_adj,
                        "Species_coef"=species_dd_test_res$Coef,
                        "Human_expressed_prop"=expressed_prop$Human,
                        "Chimp_expressed_prop"=expressed_prop$Chimp,
                        "Human_Chimp_prop_diff"=expressed_prop$Prop_diff)
  
  X <- tIO_atac@assays$peaks_species@data[expressed_peaks,idx]
  y <- tIO_atac@meta.data[idx, "Species"]
  da <- presto::wilcoxauc(X=X, y=y)
  res_mat$Species_test_wilcox_P <- da$pval[which(da$group=="Human")]
  res_mat$Corrected_wilcox_P <- p.adjust(res_mat$Species_test_wilcox_P, method="BH")
  res_mat$Human_avgExpr <- da$avgExpr[which(da$group=="Human")]
  res_mat$Chimp_avgExpr <- da$avgExpr[which(da$group=="Chimp")]
  res_mat$Human_Chimp_logFC <- da$logFC[which(da$group=="Human")]
  
  human_high_features <- rownames(res_mat)[which(res_mat$Species_test_resi_P<0.3 
                                                 & res_mat$Human_Chimp_prop_diff> 0.05
                                                 & res_mat$Human_expressed_prop>0.1)]
  chimp_high_features <- rownames(res_mat)[which(res_mat$Species_test_resi_P<0.3 
                                                 & res_mat$Human_Chimp_prop_diff< -0.05 
                                                 & res_mat$Chimp_expressed_prop>0.1)]
  
  res_mat$Human_high <- rownames(res_mat)%in%human_high_features
  res_mat$Chimp_high <- rownames(res_mat)%in%chimp_high_features
  print(paste("Human_high_feature =", length(human_high_features)))
  print(paste("Chimp_high_feature =", length(chimp_high_features)))  
  saveRDS(res_mat, file=paste0("Res_tHIO_vs_tCIO_",x,"_differential_detection_test_res_no_wilcoxauc_filtering-2.rds"))
  
  ddp_by_cell_type[[x]] <- res_mat
}
saveRDS(ddp_by_cell_type, file="Res_DDP_by_cell_type.rds")


# merge the EEC subtypes, using EEC subtype label as covariates
x <- "EEC"
print(paste(x, "start"))

# get detected peaks
tIO_atac$Cell_type_per_species <- paste(tIO_atac$Cell_type, tIO_atac$Species, sep="@")
cell_type_det_rate <- getExpressedProp(seu.obj=tIO_atac, feature.to.calc = "Cell_type_per_species", colname.prefix=NULL, assay.type="peaks_species")
saveRDS(cell_type_det_rate, file="Dat_tIO_epi_ATAC_hg38_coor_cell_type_per_species_detection_rate.rds")

cell_type_det_rate <- readRDS("Dat_tIO_epi_ATAC_hg38_coor_high_resolution_cell_type_per_species_detection_rate.rds")
rate_mat <- cell_type_det_rate[,grep(x,colnames(cell_type_det_rate))]
expressed_peaks <- unique(unlist(lapply(seq(ncol(rate_mat)), function(j){
  vec <- rate_mat[,j]
  peaks1 <- rownames(cell_type_det_rate)[which(vec>0.01)]
  num <- min(c(150000, length(peaks1)))
  peaks1[order(rate_mat[peaks1,j], decreasing = T)[1:num]]
})))


idx <- which(tIO_atac$Cell_type==x)
combined_count <- as.matrix(tIO_atac@assays$peaks_species@counts[expressed_peaks,idx])
count_cutoff <- 0
bi_mat <- combined_count>count_cutoff
species_vec <- as.factor(tIO_atac$Species[idx])
subtype_vec <- as.factor(tIO_atac$High_resolution_cell_type[idx])
umi_vec <- tIO_atac$nCount_peaks_species[idx]
colon_score_vec <- tIO_atac$`Intestine_region_L2:prob.Colon`[idx]
proxSI_score_vec <- tIO_atac$`Intestine_region_L2:prob.Prox_SI`[idx]
library(doParallel)
registerDoParallel(10)
species_dd_test_res <- foreach(k=seq(nrow(bi_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.numeric(as.vector(bi_mat[k,]))
  m0 <- glm(e ~ umi_vec+colon_score_vec+proxSI_score_vec+subtype_vec, family = binomial)
  m1 <- glm(e ~ umi_vec+colon_score_vec+proxSI_score_vec+subtype_vec+species_vec, family = binomial)
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0, test = "Chisq")$Pr[2]
  p_resi <- pf((a0[nrow(a0),"Resid. Dev"]/a0[nrow(a0),"Resid. Df"]) / (a1[nrow(a1),"Resid. Dev"]/a1[nrow(a1),"Resid. Df"]), 
               df1 = a0[nrow(a0),"Resid. Df"], df2 = a1[nrow(a1),"Resid. Df"], lower.tail = F)
  coef <- coef(m1)[length(coef(m1))]
  return(c(p_anova, p_resi ,coef))
}
saveRDS(species_dd_test_res, file="tmp.rds")
stopImplicitCluster()
rownames(species_dd_test_res) <- expressed_peaks
colnames(species_dd_test_res) <- c("p_ANOVA", "p_Resi", "Coef")
resi_p_adj <- p.adjust(species_dd_test_res[,2], method = "BH")
anova_p_adj <- p.adjust(species_dd_test_res[,1], method = "BH")
species_dd_test_res <- as.data.frame(species_dd_test_res)

expressed_prop <- sapply(sort(unique(species_vec)), function(x){
  idx <- which(species_vec==x)
  prop_vec <- rowMeans(bi_mat[,idx])
  return(prop_vec)
})
colnames(expressed_prop) <- sort(unique(species_vec))
expressed_prop <- as.data.frame(expressed_prop)
expressed_prop$Prop_diff <- expressed_prop$Human-expressed_prop$Chimp
res_mat <- data.frame("Species_test_anova_P"=species_dd_test_res$p_ANOVA,
                      "Species_test_resi_P"=species_dd_test_res$p_Resi,
                      "Corrected_anova_P"=anova_p_adj,
                      "Corrected_resi_P"=resi_p_adj,
                      "Species_coef"=species_dd_test_res$Coef,
                      "Human_expressed_prop"=expressed_prop$Human,
                      "Chimp_expressed_prop"=expressed_prop$Chimp,
                      "Human_Chimp_prop_diff"=expressed_prop$Prop_diff)

X <- tIO_atac@assays$peaks_species@data[expressed_peaks,idx]
y <- tIO_atac@meta.data[idx, "Species"]
da <- presto::wilcoxauc(X=X, y=y)
res_mat$Species_test_wilcox_P <- da$pval[which(da$group=="Human")]
res_mat$Corrected_wilcox_P <- p.adjust(res_mat$Species_test_wilcox_P, method="BH")
res_mat$Human_avgExpr <- da$avgExpr[which(da$group=="Human")]
res_mat$Chimp_avgExpr <- da$avgExpr[which(da$group=="Chimp")]
res_mat$Human_Chimp_logFC <- da$logFC[which(da$group=="Human")]

human_high_features <- rownames(res_mat)[which(res_mat$Species_test_resi_P<0.3 
                                               & res_mat$Human_Chimp_prop_diff> 0.05
                                               & res_mat$Human_expressed_prop>0.1)]
chimp_high_features <- rownames(res_mat)[which(res_mat$Species_test_resi_P<0.3 
                                               & res_mat$Human_Chimp_prop_diff< -0.05 
                                               & res_mat$Chimp_expressed_prop>0.1)]

res_mat$Human_high <- rownames(res_mat)%in%human_high_features
res_mat$Chimp_high <- rownames(res_mat)%in%chimp_high_features
print(paste("Human_high_feature =", length(human_high_features)))
print(paste("Chimp_high_feature =", length(chimp_high_features)))  
saveRDS(res_mat, file="Res_tHIO_vs_tCIO_EEC_with_subtype_identity_as_covariate_differential_detection_test_res_no_wilcoxauc_filtering-2.rds")

ddp_by_cell_type[["EEC_no_subtype_covar"]] <- ddp_by_cell_type[["EEC"]]
ddp_by_cell_type[["EEC"]] <- res_mat
saveRDS(ddp_by_cell_type, file="Res_DDP_by_cell_type_2.rds")

p_cutoff = 0.3
prop_diff_cutoff = 0.05
prop_cutoff = 0.1
selected_cell_type <- setdiff(sort(unique(tIO_atac$Cell_type)), "Putative_secretory_progenitors")
dap_num <- t(sapply(selected_cell_type, function(x){
  res_mat <- ddp_by_cell_type[[x]]
  human_high_features <- rownames(res_mat)[which(res_mat$Species_test_resi_P < p_cutoff
                                                 & res_mat$Human_Chimp_prop_diff > prop_diff_cutoff
                                                 & res_mat$Human_expressed_prop > prop_cutoff)]
  chimp_high_features <- rownames(res_mat)[which(res_mat$Species_test_resi_P < p_cutoff 
                                                 & res_mat$Human_Chimp_prop_diff < -prop_diff_cutoff 
                                                 & res_mat$Chimp_expressed_prop > prop_cutoff)]
  return(c(length(human_high_features), length(chimp_high_features)))
}))
colnames(dap_num) <- c("Human_high", "Chimp_high")
dap_num

input <- t(dap_num)
colnames(input) <- sub(pattern="_", replacement = " ", colnames(input))
input[2,] <- -input[2,]

pdf("Plot_human_chimp_DAP_number_per_cell_type.pdf")
par(mar=c(10,10,5,5))
barplot(input[1,], ylim=c(-max(abs(input))-1000, max(input)+1000), las=2, ylab="DAP number", col=epi.ct.cols[colnames(input)])
barplot(input[2,], add=T, density= 20, las=2, col=epi.ct.cols[colnames(input)])
legend("topleft",legend = "Human high", bty="n")
legend("bottomleft",legend = "Chimp high", bty="n")
dev.off()


# run GREAT analysis on DAP per cell type
# read species DA peak list
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/rGREAT")
ddp_by_cell_type <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/Res_DDP_by_cell_type_2.rds")
orth_peaks <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/Res_combined_peaks_in_panTro6_hg38_coor.rds")
# run GREAT analysis
library(rGREAT)
great_res_list <- list()
for(x in sort(unique(tIO_atac$Cell_type))){
  print(x)
  mat <- ddp_by_cell_type[[x]]
  rownames(mat) <- orth_peaks$combined_peaks_hg38_coor[match(rownames(mat), orth_peaks$peak_ID)]
  human_peaks <- rownames(mat)[which(mat$Human_high)]
  gr <- StringToGRanges(human_peaks)
  job = submitGreatJob(gr, species = "hg38")
  tbl = getEnrichmentTables(job, download_by = "tsv")
  great_res_list[[x]][["Human"]] <- tbl
  
  chimp_peaks <- rownames(mat)[which(mat$Chimp_high)]
  gr <- StringToGRanges(chimp_peaks)
  job = submitGreatJob(gr, species = "hg38")
  tbl = getEnrichmentTables(job, download_by = "tsv")
  great_res_list[[x]][["Chimp"]] <- tbl
}
saveRDS(great_res_list, file="Res_online_rGREAT_res_per_cell_type_DA_peaks.rds")

enriched_term_list <- list()
for(ct in names(great_res_list)){
  for(sp in c("Human","Chimp")){
    id <- paste(ct,sp,sep=":")
    res <- great_res_list[[ct]][[sp]][["GO Biological Process"]]
    enriched_terms <- res$Desc[which(res$BinomFdrQ<0.05 & res$RegionFoldEnrich>2)]
    enriched_term_list[[sp]][[ct]] <- enriched_terms
  }
}
enrichment_by_species <- list()
for(sp in names(enriched_term_list)){
  all_terms <- unique(unlist(enriched_term_list[[sp]]))
  enrichment_mat <- matrix(NA, nrow=length(all_terms), ncol=length(enriched_term_list[[sp]]))
  rownames(enrichment_mat) <- all_terms
  colnames(enrichment_mat) <- names(enriched_term_list[[sp]])
  for(ct in colnames(enrichment_mat)){
    terms <- enriched_term_list[[sp]][[ct]]
    enrichment_mat[terms,ct] <- TRUE
  }
  enrichment_by_species[[sp]] <- enrichment_mat
}
mat <- enrichment_by_species[["Human"]]
freq <- rowSums(!is.na(enrichment_by_species[["Human"]]))
terms <- rownames(mat)[which(mat[,"Enterocyte"])]
saveRDS(enrichment_by_species, file="Res_enriched_GO_BP_terms_by_species_and_cell_type.rds")
saveRDS(enrichment_by_species, file="/projects/site/pred/ihb-intestine-evo/used_object/GO/Res_tIO_atac_enriched_GO_BP_terms_by_species_and_cell_type.rds")


# combine multi-layer information
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation")
## detected? differentially detected?
ddp_by_cell_type <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/Res_DDP_by_cell_type_2.rds")
saveRDS(ddp_by_cell_type, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_DDP_by_cell_type_2.rds")
# defined as epithelial cell type enriched?
ct_marker_peak_res <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/cell_type_markers/Res_tIO_epi_cell_type_marker_peaks_per_species_dap_res.rds")
human_mat <- ct_marker_peak_res$Human$sig_res
chimp_mat <- ct_marker_peak_res$Chimp$sig_res
tIO_peak_id <- rownames(tIO_atac)
peak_anno <- matrix(F, nrow=nrow(tIO_atac), ncol=3*length(ddp_by_cell_type)+2*length(unique(human_mat$group)))
rownames(peak_anno) <- rownames(tIO_atac)
colnames(peak_anno) <- c(paste(rep(c("Expressed","Human_high","Chimp_high"), each=length(ddp_by_cell_type)),
                             rep(names(ddp_by_cell_type), 3),
                             sep=":"),
                         paste(rep(c("Human","Chimp"), each=length(unique(human_mat$group))),
                               rep(sort(unique(human_mat$group)),2),
                               sep=":"))
colname_vec1 <- c(paste(rep(c("Expressed","Human_high","Chimp_high"), each=length(ddp_by_cell_type)),
        rep(names(ddp_by_cell_type), 3),
        sep=":"),
  paste(rep(c("Human","Chimp"), each=length(unique(human_mat$group))),
        rep(sort(unique(human_mat$group)),2),
        sep=":"))

for(x in names(ddp_by_cell_type)){
  mat <- ddp_by_cell_type[[x]]
  peak_anno[,paste("Expressed",x,sep=":")] <- tIO_peak_id %in% rownames(mat)
  peak_anno[,paste("Human_high",x,sep=":")] <- tIO_peak_id %in% rownames(mat)[which(mat$Human_high)]
  peak_anno[,paste("Chimp_high",x,sep=":")] <- tIO_peak_id %in% rownames(mat)[which(mat$Chimp_high)]
}
for(x in sort(unique(human_mat$group))){
  peak_anno[,paste("Human",x,sep=":")] <- tIO_peak_id %in% human_mat$feature[which(human_mat$group==x)]
}
for(x in sort(unique(chimp_mat$group))){
  peak_anno[,paste("Chimp",x,sep=":")] <- tIO_peak_id %in% chimp_mat$feature[which(chimp_mat$group==x)]
}

orth_peaks <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/Res_combined_peaks_in_panTro6_hg38_coor.rds")
rownames(peak_anno) <- orth_peaks$combined_peaks_hg38_coor[match(rownames(peak_anno), orth_peaks$peak_ID)]

# overlap with evolution signatures? 
evo_sig_overlap_mat <- readRDS("/projects/site/pred/ihb-intestine-evo/evo_signature/GO_enrichment/Res_tIO_all_peaks_overlap_with_evolution_signatures_per_paper.rds")
peak_anno <- data.frame(peak_anno, evo_sig_overlap_mat, stringsAsFactors = F)
colnames(peak_anno)[(ncol(peak_anno)-(ncol(evo_sig_overlap_mat)-1)):ncol(peak_anno)] <- colnames(evo_sig_overlap_mat)
saveRDS(peak_anno, file="Res_tIO_merged_peak_annotation-2.rds")


# ChIP-seeker-based peak annotation
anno <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/human_peak_annotation/Data_frame_peak_gene_pairs.rds")
rownames(anno) <- anno$name
vec <- setNames(rep(NA, nrow(peak_anno)), rownames(peak_anno))
vec[rownames(anno)] <- anno$annotation_category
vec2 <- setNames(rep(NA, nrow(peak_anno)), rownames(peak_anno))
vec2[rownames(anno)] <- anno$gene_symbol
peak_anno$Closest_gene <- vec2
peak_anno$annotation_category <- vec

# get human-chimp DEGs
hc_degs <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff/Res_human_chimp_DEGs_with_input_gene_list.rds")
saveRDS(hc_degs, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_DEGs_with_input_gene_list.rds")
# get lineage-specific genes
hcm_degs <- readRDS("/projects/site/pred/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/pairwise_DEGs/Res_human_chimp_mouse_DEG_mat.rds")
saveRDS(hcm_degs, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_mouse_DEG_mat.rds")
# get human-mouse DEGs
hm_degs <- readRDS("/projects/site/pred/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/pairwise_DEGs/Res_fetal_human_and_mouse_ANCOVA_DEG_per_cell_type_res.rds")
saveRDS(hm_degs, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_fetal_human_and_mouse_ANCOVA_DEG_per_cell_type_res.rds")

shared_genes <- intersect(rownames(hc_degs), rownames(hm_degs))
length(shared_genes)

colnames(hc_degs) <- paste("H_and_C", colnames(hc_degs), sep=":")
colnames(hm_degs) <- paste("H_and_M", colnames(hm_degs), sep=":")
colnames(hcm_degs) <- paste("H_C_M", colnames(hcm_degs), sep=":")
combined_degs <- cbind(hc_degs[shared_genes,], hm_degs[shared_genes,], hcm_degs[shared_genes,])
saveRDS(combined_degs, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_mouse_pairwise_DEG_combined_mat.rds")

for(x in colnames(hc_degs)){
  degs <- rownames(hc_degs)[which(hc_degs[,x]%in%c("Human_high","Chimp_high"))]
  peak_anno[[x]] <- peak_anno$Closest_gene%in%degs
}
for(x in colnames(hm_degs)){
  degs <- rownames(hm_degs)[which(hm_degs[,x]%in%c("Human_high","Mouse_high"))]
  peak_anno[[x]] <- peak_anno$Closest_gene%in%degs
}
for(x in colnames(hcm_degs)){
  degs <- rownames(hcm_degs)[which(hcm_degs[,x]%in%c("human_specific_high","human_specific_low"))]
  x <- sub(pattern="H_C_M", replacement = "Human_specific", x)
  peak_anno[[x]] <- peak_anno$Closest_gene%in%degs
}
for(x in colnames(hcm_degs)){
  degs <- rownames(hcm_degs)[which(hcm_degs[,x]%in%c("primate_specific_high","primate_specific_low"))]
  x <- sub(pattern="H_C_M", replacement = "Primate_specific", x)
  peak_anno[[x]] <- peak_anno$Closest_gene%in%degs
}
for(x in colnames(hcm_degs)){
  degs <- rownames(hcm_degs)[which(hcm_degs[,x]=="conserved")]
  x <- sub(pattern="H_C_M", replacement = "Conserved", x)
  peak_anno[[x]] <- peak_anno$Closest_gene%in%degs
}
saveRDS(peak_anno, file="Res_tIO_merged_peak_annotation-2.rds")
saveRDS(peak_anno, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation.rds")

tIO_glm <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/Pando_on_tIO_ct_markers_and_DEGs/Res_tIO_epi_Pando_glm_model_grn.rds")
df <- tIO_glm@assays$peaks_species_human@links

coef <- coef(tIO_glm)
grn <- NetworkModules(tIO_glm)
library(Signac)
region_gr <- StringToGRanges(coef$region)
peak_gr <- StringToGRanges(rownames(peak_anno))
overlap_reg <- findOverlaps(query = region_gr, subject = peak_gr)
coef$peak <- rownames(peak_anno)[overlap_reg@to]
pando_pairs <- unique(paste(coef$peak, coef$target, sep=":")[which(coef$pval<0.05)])
cor_pairs <- unique(paste(df$peak, df$gene, sep=":"))
anno_pairs <- paste(rownames(peak_anno), peak_anno$Closest_gene, sep=":")
peak_anno$Linked_by_cor <- anno_pairs%in%cor_pairs
peak_anno$Linked_by_Pando <- anno_pairs%in%pando_pairs
colnames(peak_anno)[seq(length(colname_vec1))] <- colname_vec1
colnames(peak_anno)[c(1:37,54:83)] <- paste(colnames(peak_anno)[c(1:37,54:83)], "fetal", sep="@")
saveRDS(peak_anno, file="Res_tIO_merged_peak_annotation-2.rds")
saveRDS(peak_anno, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation.rds")

# evolutionary node assignment
node_res <- readRDS("/projects/site/pred/ihb-intestine-evo/regulome_and_evolution/christi.node.assignments.rds")
tHIO_peak_gr <- StringToGRanges(node_res$name)
overlap_node <- findOverlaps(query = tHIO_peak_gr, subject = peak_gr)
freq1 <- table(overlap_node@from)
idx <- which(overlap_node@from %in% names(freq1)[which(freq1==1)])
used_overlap_node <- overlap_node[idx]

node_res$combined_peaks <- NA
node_res$combined_peaks[used_overlap_node@from] <- rownames(peak_anno)[used_overlap_node@to]
mat <- as.matrix(unique(node_res[,c("combined_peaks","node.label")]))
idx <- which(rowSums(is.na(mat))==0)
mat <- mat[idx,]
freq <- table(mat[,1])
aa <- names(freq)[which(freq==1)]
unique_mat <- mat[which(mat[,1] %in% aa),]
peak_anno$node.label <- NA
peak_anno[unique_mat[,1], "node.label"] <- unique_mat[,2]
saveRDS(peak_anno, file="Res_tIO_merged_peak_annotation-2.rds")
saveRDS(peak_anno, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation.rds")

write.table(peak_anno, file="Table_tIO_merged_peak_annotation.txt",sep="\t",quote=F)
write.table(combined_degs, file="Table_fetal_human_chimp_mouse_DEGs.txt",sep="\t",quote=F)
genes <- c("SLC5A12","SLC26A3","ETV4","FEV","MYC")


# compare the GREAT enrichment for DARs of enterocytes and stem cells
# get DAR enriched terms that are specific to enterocyte, stem cells and shared
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/DE_vs_DA/GO")

GO_anno <- read.csv("/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v109/Ensembl_v109_GO.csv")
idx <- which(GO_anno$GO.domain=="biological_process")
GO_anno <- GO_anno[idx,]
GO_term_names <- unique(GO_anno$GO.term.name)
# load the GO enrichment for DARs classifying the species differential direction
great_res_list <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/rGREAT/Res_online_rGREAT_res_per_cell_type_DA_peaks.rds")



pval_list <- list()
or_list <- list()


all_top_terms <- c()
for(ct in names(great_res_list)){
  for(sp in c("Human", "Chimp")){
    res <- great_res_list[[ct]][[sp]][["GO Biological Process"]]
    all_top_terms <- union(all_top_terms, res$Desc)
    
    pval_vec <- setNames(res$BinomFdrQ, res$Desc)
    or_vec <-  setNames(res$RegionFoldEnrich, res$Desc)
    
    if(length(pval_list[[ct]])<1){
      pval_list[[ct]] <- pval_vec
      or_list[[ct]] <- or_vec
      
    }else{
      union_terms <- union(names(pval_list[[ct]]), names(pval_vec))
      shared_terms <- intersect(names(pval_list[[ct]]), names(pval_vec))
      newly_added_terms <- setdiff(names(pval_vec), names(pval_list[[ct]]))
      
      pval_vec_combined <- setNames(rep(NA, length(union_terms)), union_terms)
      pval_vec_combined[names(pval_list[[ct]])] <- pval_list[[ct]]
      pval_vec_combined[newly_added_terms] <- pval_vec[newly_added_terms]
      pval_vec_combined[shared_terms] <- apply(cbind(pval_list[[ct]][shared_terms], pval_vec[shared_terms]),1,min)
      pval_list[[ct]] <- pval_vec_combined
      
      or_vec_combined <- setNames(rep(NA, length(union_terms)), union_terms)
      or_vec_combined[names(or_list[[ct]])] <- or_list[[ct]]
      or_vec_combined[newly_added_terms] <- or_vec[newly_added_terms]
      or_vec_combined[shared_terms] <- apply(cbind(or_list[[ct]][shared_terms], or_vec[shared_terms]),1,max)
      or_list[[ct]] <- or_vec_combined
      
    }
  }
}

pval_dar <- matrix(NA, nrow=length(all_top_terms), ncol=length(pval_list))
rownames(pval_dar) <- all_top_terms
colnames(pval_dar) <- names(pval_list)

or_dar <- matrix(NA, nrow=length(all_top_terms), ncol=length(pval_list))
rownames(or_dar) <- all_top_terms
colnames(or_dar) <- names(pval_list)

for(ct in names(pval_list)){
  pval_dar[names(pval_list[[ct]]), ct] <- pval_list[[ct]]
  or_dar[names(pval_list[[ct]]), ct] <- or_list[[ct]]
}
saveRDS(pval_dar, file="Res_tIO_per_cell_type_DAR_pval_mat.rds")
saveRDS(or_dar, file="Res_tIO_per_cell_type_DAR_odds_ratio_mat.rds")


# fill the NA in pval matrix with 1, and in odds ratio matrix with 0
pval_dar <- apply(pval_dar, 2, function(vec){
  vec[is.na(vec)] <- 1
  return(vec)
})
or_dar <- apply(or_dar, 2, function(vec){
  vec[is.na(vec)] <- 0
  return(vec)
})
saveRDS(pval_dar, file="Res_tIO_per_cell_type_DAR_pval_mat_filled_with_1.rds")
saveRDS(or_dar, file="Res_tIO_per_cell_type_DAR_odds_ratio_mat_filled_with_0.rds")

# get enriched GO terms
sig_idx <- (pval_dar<0.05) * (or_dar>2)

# exclude the GO terms that are enriched in at least one cell type but not all cell types
cell_type_sig_idx <- rowSums(sig_idx)>0 & rowSums(sig_idx)<ncol(sig_idx)
input <- sig_idx[cell_type_sig_idx,]
dar_enriched_GO_terms <- list(
  "all_top_involved"=sig_idx,
  "cell_type_biased_sig"=input
)
saveRDS(dar_enriched_GO_terms, file="Res_DAR_enriche_GO_terms.rds")


list <- lapply(seq(ncol(input)), function(j){
  rownames(input)[which(input[,j]>0)]
})
el <- cbind(unlist(list),
            rep(colnames(input), sapply(list, length)))
saveRDS(el, file="Res_DAR_enriche_GO_terms_igraph_input_edge_list.rds")

library(igraph)
igraph_obj <- graph_from_edgelist(el)
layout_nicely(igraph_obj) -> layout
rownames(layout) <- V(igraph_obj)$name
saveRDS(layout, file="Res_DAR_enriche_GO_terms_igraph_layout.rds")


pdf("Plot_igraph_DAR_enriched_GO_terms.pdf")
plot(layout, pch=16, col="#d0d0d0", xlab="", ylab="", yaxt="n", xaxt="n", yaxt="n", bty="n")
for(ct in colnames(sig_idx)){
  terms <- el[el[,2]==ct,1]
  for(x in terms){
    lines(layout[c(x,ct), 1], layout[c(x,ct), 2], lwd=0.5, col="#d0d0d0")
  }
}
points(layout[colnames(sig_idx),], pch=16, cex=3)
text(layout[colnames(sig_idx),]+2, labels =colnames(sig_idx))
dev.off()

terms <- grep("lactate", rownames(input), value = T)
points(layout[terms,], pch=16, cex=1.5)
text(layout[terms,], labels = terms)

# determine the species differential direction
all_enriched_GO_terms <- c()
all_id <- c()
enriched_GO_res <- list()
for(ct in names(great_res_list)){
  for(sp in c("Human", "Chimp")){
    res <- great_res_list[[ct]][[sp]][["GO Biological Process"]]
    idx <- which(res$BinomFdrQ<0.05 & res$RegionFoldEnrich>2)
    terms <- res$Desc[idx]
    pval <- setNames(res$BinomFdrQ[idx], terms)
    or <- setNames(res$RegionFoldEnrich[idx], terms)
    id <- paste(ct,sp,sep="@")
    enriched_GO_res[[id]][["pval"]] <- pval
    enriched_GO_res[[id]][["or"]] <- or
    all_enriched_GO_terms <- union(all_enriched_GO_terms, terms)
    all_id <- c(all_id, id)
  }
}
DAR_enriched_GO <- matrix(0, nrow=length(all_enriched_GO_terms), ncol=length(all_id))
rownames(DAR_enriched_GO) <- all_enriched_GO_terms
colnames(DAR_enriched_GO) <- all_id
for(id in names(enriched_GO_res)){
  DAR_enriched_GO[names(enriched_GO_res[[id]][["pval"]]),id] <- 1
}

DAR_enriched_GO_pval <- matrix(0, nrow=length(all_enriched_GO_terms), ncol=length(all_id))
rownames(DAR_enriched_GO_pval) <- all_enriched_GO_terms
colnames(DAR_enriched_GO_pval) <- all_id
DAR_enriched_GO_or <- matrix(0, nrow=length(all_enriched_GO_terms), ncol=length(all_id))
rownames(DAR_enriched_GO_or) <- all_enriched_GO_terms
colnames(DAR_enriched_GO_or) <- all_id

for(id in names(enriched_GO_res)){
  DAR_enriched_GO_pval[names(enriched_GO_res[[id]][["pval"]]),id] <- enriched_GO_res[[id]][["pval"]]
  DAR_enriched_GO_or[names(enriched_GO_res[[id]][["or"]]),id] <- enriched_GO_res[[id]][["or"]]
}

DAR_enriched_GO_res_list <- list(
  "direction"=DAR_enriched_GO,
  "pval"=DAR_enriched_GO_pval,
  "or"=DAR_enriched_GO_or
)
saveRDS(DAR_enriched_GO_res_list, file="Res_per_cell_type_per_species_DAR_enriched_GO_res_list.rds")

DAR_GO_direction <- matrix(NA, nrow=nrow(DAR_enriched_GO), ncol=length(great_res_list))
rownames(DAR_GO_direction) <- rownames(DAR_enriched_GO)
colnames(DAR_GO_direction) <- names(great_res_list)
DAR_GO_pval <- matrix(NA, nrow=nrow(DAR_enriched_GO), ncol=length(great_res_list))
rownames(DAR_GO_pval) <- rownames(DAR_enriched_GO)
colnames(DAR_GO_pval) <- names(great_res_list)
DAR_GO_or <- matrix(NA, nrow=nrow(DAR_enriched_GO), ncol=length(great_res_list))
rownames(DAR_GO_or) <- rownames(DAR_enriched_GO)
colnames(DAR_GO_or) <- names(great_res_list)

for(ct in names(great_res_list)){
  # determine the direction
  mat <- DAR_enriched_GO[,paste(ct,c("Human","Chimp"),sep="@")]
  idx0 <- which(mat[,1]>0 & mat[,2]>0)
  idx1 <- which(mat[,1]>0 & mat[,2]==0)
  idx2 <- which(mat[,1]==0 & mat[,2]>0)
  DAR_GO_direction[idx0,ct] <- 0
  DAR_GO_direction[idx1,ct] <- 1
  DAR_GO_direction[idx2,ct] <- -1
  
  # determine pval
  mat <- -log10(DAR_enriched_GO_pval[,paste(ct,c("Human","Chimp"),sep="@")])
  if(length(idx0)>1){
    DAR_GO_pval[idx0,ct] <- rowMeans(mat[idx0,])
  }else{
    DAR_GO_pval[idx0,ct] <- mean(mat[idx0,])
  }
  DAR_GO_pval[idx1,ct] <- mat[idx1,1]
  DAR_GO_pval[idx2,ct] <- mat[idx2,2]
  
  # determine odds ratio
  mat <- DAR_enriched_GO_or[,paste(ct,c("Human","Chimp"),sep="@")]
  if(length(idx0)>1){
    DAR_GO_or[idx0,ct] <- rowMeans(mat[idx0,])
  }else{
    DAR_GO_or[idx0,ct] <- mean(mat[idx0,])
  }
  DAR_GO_or[idx1,ct] <- mat[idx1,1]
  DAR_GO_or[idx2,ct] <- mat[idx2,2]
  
}

DAR_enriched_GO_res_list <- list(
  "direction"=DAR_GO_direction,
  "pval"=DAR_GO_pval,
  "or"=DAR_GO_or
)
saveRDS(DAR_enriched_GO_res_list, file="Res_per_cell_type_DAR_enriched_GO_res_list.rds")



DAR_GO_direction_across_cell_type <- rowSums(DAR_GO_direction, na.rm=T)
DAR_GO_direction_id_across_cell_type <- ifelse(DAR_GO_direction_across_cell_type>0, "Human_high", ifelse(DAR_GO_direction_across_cell_type<0, "Chimp_high", "Bidirection"))
DAR_GO_pval_across_cell_type <- rowMeans(DAR_GO_pval, na.rm=T)
DAR_GO_or_across_cell_type <- rowMeans(DAR_GO_or, na.rm=T)


DAR_enriched_GO_res_df <- data.frame(
  "direction"=DAR_GO_direction_across_cell_type,
  "group"=DAR_GO_direction_id_across_cell_type,
  "pval"=DAR_GO_pval_across_cell_type,
  "or"=DAR_GO_or_across_cell_type,
  stringsAsFactors = F
)
saveRDS(DAR_enriched_GO_res_df, file="Res_across_cell_type_DAR_enriched_GO_res_df.rds")


# exclude the GO terms that are enriched in at least one cell type but not all cell types
cell_type_sig_idx <- rowSums(!is.na(DAR_GO_direction))<ncol(DAR_GO_direction)
input <- DAR_GO_direction[cell_type_sig_idx,]

list <- lapply(seq(ncol(input)), function(j){
  rownames(input)[which(!is.na(input[,j]))]
})
el <- cbind(unlist(list),
            rep(colnames(input), sapply(list, length)))
saveRDS(el, file="Res_DAR_enriche_GO_terms_igraph_input_edge_list.rds")

direction_per_cell_type_vec <- unlist(lapply(seq(ncol(input)), function(ct){
  input[which(!is.na(input[,ct])),ct]
}))
input <- DAR_GO_pval[cell_type_sig_idx,]
pval_per_cell_type_vec <- unlist(lapply(seq(ncol(input)), function(ct){
  input[which(!is.na(input[,ct])),ct]
}))
input <- DAR_GO_or[cell_type_sig_idx,]
or_per_cell_type_vec <- unlist(lapply(seq(ncol(input)), function(ct){
  input[which(!is.na(input[,ct])),ct]
}))
edge_df <- data.frame(
  "from_name"=el[,1],
  "to_name"=el[,2],
  "pair_name"=paste(el[,1], el[,2], sep=":"),
  "direction"=direction_per_cell_type_vec,
  "pval"=pval_per_cell_type_vec,
  "or"=or_per_cell_type_vec,
  stringsAsFactors = F
)
saveRDS(edge_df, file="Dat_igraph_edge_df.rds")



library(igraph)
igraph_obj <- graph_from_edgelist(el)
layout_nicely(igraph_obj) -> layout
rownames(layout) <- V(igraph_obj)$name
saveRDS(layout, file="Res_DAR_enriche_GO_terms_igraph_layout.rds")

pval_vec <- setNames(rep(5, nrow(layout)), rownames(layout))
or_vec <- setNames(rep(100, nrow(layout)), rownames(layout))
group_vec <- setNames(rep("Cell_type_label", nrow(layout)), rownames(layout))
terms <- intersect(rownames(layout), rownames(DAR_enriched_GO_res_df))
pval_vec[terms] <- DAR_enriched_GO_res_df[terms, "pval"]
or_vec[terms] <- DAR_enriched_GO_res_df[terms, "or"]
group_vec[terms] <- DAR_enriched_GO_res_df[terms, "group"]

node_df <- data.frame(
  "name"=rownames(layout),
  "X"=layout[,1],
  "Y"=layout[,2],
  "size_pval"=pval_vec,
  "size_or"=or_vec,
  "group"=group_vec,
  stringsAsFactors = F
)

saveRDS(node_df, file="Dat_igraph_node_df.rds")

group_cols <- setNames(c("#d0d0d0", "#303030", "#0571b0","#ca0020"),
                       c("Bidirection", "Cell_type_label", "Chimp_high", "Human_high"))



# focus on terms that are enriched in no more than 2 cell types
cell_type_specific_idx <- rowSums(!is.na(DAR_GO_direction))<3
# select top 2 enriched GO terms per cell type per species
input_metric <- DAR_enriched_GO_or[cell_type_specific_idx,]
top_terms <- apply(input_metric, 2, function(vec){
  rownames(input_metric)[order(vec, decreasing = T)[1:2]]
})
top_terms <- unique(as.vector(top_terms))
cell_types <- 
to_label <- c(top_terms, unique(el[,2]))


coor <- cbind(df$X, df$Y)
rownames(coor) <- df$name

plot(coor, pch=16, col="#d0d0d0", xlab="", ylab="", yaxt="n", xaxt="n", yaxt="n", bty="n", type="n")
for(ct in unique(el[,2])){
  terms <- el[el[,2]==ct,1]
  for(x in terms){
    lines(coor[c(x,ct), 1], coor[c(x,ct), 2], lwd=0.3, col="#d0d0d0")
  }
}
points(coor, pch=21, cex=df$size*0.5, bg=group_cols[df$group], col="#303030", lwd=0.3)
points(coor[to_label,], cex=df[to_label,"size"]*0.5, lwd=1)
text(coor[to_label,], labels = to_label, col=group_cols[df[to_label,"group"]])
legend("bottomright", pch=16, legend=names(group_cols), col=group_cols, bty="n")
dev.off()



library(ggraph)
library(tidygraph)
ggraph_obj <- as_tbl_graph(igraph_obj)

ggraph_obj <- ggraph_obj %>% 
  tidygraph::activate(nodes) %>%
  left_join(node_df, by = c("name" = "name")) %>% 
  tidygraph::activate(edges) %>%
  mutate(from_name = (.N()$name[from])) %>%
  mutate(to_name = (.N()$name[to])) %>%
  mutate(pair_name = paste(from_name, to_name, sep=":")) %>%
  left_join(edge_df, by=c("pair_name"="pair_name")) %>%
  tidygraph::activate(nodes)
saveRDS(ggraph_obj, file="Res_ggraph_object_of_DAR_enriched_GO_terms.rds")



p1 <- ggraph(ggraph_obj, x=X, y=Y) +
  geom_edge_diagonal(aes(alpha=or, color=factor(direction)),
                     width=0.5, arrow = arrow(length = unit(1,"mm"))) +
  scale_edge_alpha_continuous(range=c(0.1,0.8), guide = "none") +
  scale_edge_color_manual(values = c('-1'='#7FB3D5', '1'='#EC7063')) +
  geom_node_point(aes(size = size_or, fill = group),
                  shape = 21, color='darkgrey') +
  scale_fill_manual(values = group_cols) +
  scale_size_continuous(range = c(1,5), trans = "sqrt") +
  geom_node_point(aes(size = size_or, filter = name%in%top_terms),
                  shape = 1, color='#303030') +
  geom_node_label(aes(label=name, filter = group=="Cell_type_label"),
                 size=3, repel=T, max.overlaps = 13) +
  geom_node_text(aes(label=name, filter = name %in% top_terms),
                  min.segment.length = unit(0, 'lines'), size=3, repel=T, max.overlaps = 99999) +
  theme_void()

pdf("Plot_ggraph_DAR_enriched_GO_terms.pdf")
p1
dev.off()


