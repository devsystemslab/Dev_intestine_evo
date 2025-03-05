setwd("/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet/update_cell_type_annotation")
# get dn/ds ratio of epithelial cell types
# get dN/dS score per cell based on all genes and weighted by normalized expression
fetal_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_fetal_tissue_all_cell_class_unified_cell_type_annotation.rds")

dnds_df <- readRDS("/pstore/data/ihb-intestine-evo/evo_signature/DN_DS_ratio/data/Dat_human_to_10sp_one2one_orthologs_dnds_ratio.rds")
shared_genes <- intersect(rownames(dnds_df), rownames(fetal_rna))

for(sp in c("primate","mammal")){
  ave_dnds <- setNames(dnds_df[[sp]], rownames(dnds_df))
  weighted_expr <- fetal_rna@assays$RNA@data[shared_genes,]*ave_dnds[shared_genes]
  sum_weighted_expr <- colSums(weighted_expr, na.rm=T)
  sum_log_expr <- colSums(fetal_rna@assays$RNA@data[shared_genes,])
  expr_weighted_ave_dnds_score <- sum_weighted_expr/sum_log_expr
  fetal_rna[[paste0(sp, "_average_dnds_score")]] <- expr_weighted_ave_dnds_score 
}

dnds_ratio_list <- sapply(unique(fetal_rna$Unified_cell_type), function(x){
  fetal_rna$primate_average_dnds_score[which(fetal_rna$Unified_cell_type==x)]
})
names(dnds_ratio_list) <- unique(fetal_rna$Unified_cell_type)

median_vec <- sapply(unique(fetal_rna$Unified_cell_type), function(x){
  median(fetal_rna$primate_average_dnds_score[which(fetal_rna$Unified_cell_type==x)])
})
names(median_vec) <- unique(fetal_rna$Unified_cell_type)

cell_type_order <- names(rev(sort(median_vec)))
pval <- sapply(1:(length(cell_type_order)-1), function(i){
  ct1 <- cell_type_order[i]
  ct2 <- cell_type_order[(i+1)]
  wilcox.test(dnds_ratio_list[[ct1]], dnds_ratio_list[[ct2]], alternative="greater")$p.value
})
fetal_rna@active.ident  <-  factor(x = fetal_rna$Unified_cell_type, 
                                  levels = names(rev(sort(median_vec))))
p1 <- SCpubr::do_FeaturePlot(fetal_rna, 
                         reduction = "cca_umap",
                         features="primate_average_dnds_score",
                         pt.size = 4)+NoLegend()+NoAxes()
p1
png("Plot_UMAP_dnds_score_all_cell_class.png", height=2000, width=2000)
p1
dev.off()

p1 <- SCpubr::do_ViolinPlot(fetal_rna, 
features="primate_average_dnds_score", 
colors.use=cl_cols) 
pdf("Plot_violin_plot_dnds_score_all_cell_type.pdf", height=7, width=15)
p1
dev.off()
saveRDS(fetal_rna, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_fetal_tissue_all_cell_class_unified_cell_type_annotation.rds")

# get dN/dS ratio of expressed genes without weighting by expression
fetal_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_fetal_tissue_all_cell_class_unified_cell_type_annotation.rds")
dnds_df <- readRDS("/pstore/data/ihb-intestine-evo/evo_signature/DN_DS_ratio/data/Dat_human_to_10sp_one2one_orthologs_dnds_ratio.rds")
shared_genes <- intersect(rownames(dnds_df), rownames(fetal_rna))

sp <- "primate"
ave_dnds <- setNames(dnds_df[[sp]], rownames(dnds_df))
expressed_gene_dnds <- (fetal_rna@assays$RNA@data[shared_genes,]>0)*1*ave_dnds[shared_genes]
fetal_rna$expressed_gene_mean_dnds <- apply(expressed_gene_dnds, 2, function(vec){
  mean(vec[which(vec>0)])
})  
saveRDS(fetal_rna, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_fetal_tissue_all_cell_class_unified_cell_type_annotation_with_dnds.rds")

fetal_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_fetal_tissue_all_cell_class_unified_cell_type_annotation_with_dnds.rds")
median_vec <- sapply(unique(fetal_rna$Unified_cell_type), function(x){
  median(fetal_rna$expressed_gene_mean_dnds[which(fetal_rna$Unified_cell_type==x)])
})
names(median_vec) <- unique(fetal_rna$Unified_cell_type)

dnds_ratio_list <- sapply(unique(fetal_rna$Unified_cell_type), function(x){
  fetal_rna$expressed_gene_mean_dnds[which(fetal_rna$Unified_cell_type==x)]
})
names(dnds_ratio_list) <- unique(fetal_rna$Unified_cell_type)


fetal_rna@active.ident  <-  factor(x = fetal_rna$Unified_cell_type, 
                                  levels = names(rev(sort(median_vec))))

# define colors
mat <- unique(fetal_rna@meta.data[,c("Major_cell_type", "Unified_cell_type")])
mat <- mat[order(mat[,1]),]
epi_cols <- c("#fee5d9","#fcae91","#fb6a4a","#de2d26","#a50f15")
mes_cols <- c("#edf8e9","#bae4b3","#74c476","#31a354","#006d2c")
imm_cl_cols <- setNames(c("#deebf7", "#9ecae1","#3182bd"), mat[which(mat[,1]=="Immune"),2])
neu_cl_cols <- setNames(c("#fde0dd","#fa9fb5","#c51b8a"), mat[which(mat[,1]=="Neural"),2])
epi_cl_cols <- setNames(colorRampPalette(epi_cols)(sum(mat[,1]=="Epithelial")), mat[which(mat[,1]=="Epithelial"),2])
mes_cl_cols <- setNames(colorRampPalette(mes_cols)(sum(mat[,1]=="Mesenchymal")), mat[which(mat[,1]=="Mesenchymal"),2])
endo_cols <- setNames("#8E44AD", "Endothelial cell")
cl_cols <- c(imm_cl_cols, neu_cl_cols, epi_cl_cols, mes_cl_cols, endo_cols)
saveRDS(cl_cols, file="Dat_fetal_tissue_all_cell_class_color.rds")
cell_cols <- cl_cols[fetal_rna$Unified_cell_type]
fetal_rna$Cell_cols <- cell_cols

p1 <- SCpubr::do_ViolinPlot(fetal_rna, 
features="expressed_gene_mean_dnds", 
colors.use=cl_cols) 
pdf("Plot_violin_plot_dnds_score_all_cell_type_mean_dnds_no_expr_weighting.pdf", height=7, width=15)
p1
dev.off()

# get cell type markers and show the distribution of dN/dS score of cell type markers
deg_res <- presto::wilcoxauc(fetal_rna, group_by="Unified_cell_type")
deg_res$pct_diff <- deg_res$pct_in - deg_res$pct_out
sig_res <- deg_res[which(deg_res$padj<0.05 & deg_res$logFC>0.1 & deg_res$pct_in>10 & deg_res$pct_diff>10),]
saveRDS(sig_res, file="Res_fetal_tissue_all_cell_class_unified_cell_type_marker_genes.rds")
saveRDS(sig_res, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Res_fetal_tissue_all_cell_cl
ass_unified_cell_type_marker_genes.rds")

marker_dnds <- lapply(sort(unique(sig_res$group)), function(ct){
  genes <- sig_res$feature[which(sig_res$group==ct)]
  shared_genes <- intersect(names(ave_dnds)[!is.nan(ave_dnds)], genes)
  ave_dnds[shared_genes]
})
names(marker_dnds) <- sort(unique(sig_res$group))

median_vec <- sapply(names(marker_dnds), function(ct){median(marker_dnds[[ct]])})
idx <- order(median_vec, decreasing=T)
ordered_marker_dnds <- marker_dnds[idx]
pdf("Plot_boxplot_marker_gene_dnds.pdf", width=15, height=7)
boxplot(ordered_marker_dnds, las=2, col=cl_cols[names(ordered_marker_dnds)], outline=F)
dev.off()

library(ggpubr)
df <- sig_res[,c("feature", "group")]
df$dnds <- ave_dnds[df$feature]
df <- df[which(!is.na(df$dnds)),]
df$group <- factor(x = df$group, 
                  levels = names(rev(sort(median_vec))))
df$log_dnds <- log(df$dnds)
p1 <- ggviolin(df, x = "group", y = "dnds", fill = "group",
         palette = cl_cols,
         add = "boxplot", add.params = list(fill = "white"))+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
         geom_hline(yintercept=mean(median_vec), linetype='dotted', col = '#696969')
pdf("Plot_violinPlot_marker_gene_dnds.pdf", width=15, height=7)
p1
dev.off()



library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(Seurat)
# load major cell type markers
folder <- "/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet/update_cell_type_annotation/expression_specificity/epithelial_focused"
# load the gene AA conservation and cell type expression specificity data
df_tau <- readRDS("/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet/update_cell_type_annotation/Res_gene_dnds_and_tau_score.rds")

# get the genes that are expressed in enterocyte and enriched in either in epithelium in general or enterocyte
#gene_list <- list()
#coef_mat <- c()
#for(ct in c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell", "Stem_cell")){
#  expressed_genes <- rownames(cell_type_rate)[which(cell_type_rate[,ct]>0.1)]
#  genes <- intersect(
#    expressed_genes,
#    union(cell_class_res$feature[which(cell_class_res$group=="Epithelial")], epi_res$feature[which(epi_res$group==ct)])
#  )
#  gene_list[[ct]] <- genes
#  df_marker <- df_tau[which(df_tau$feature %in% genes),]
#  df_marker$log1p_conservation <- log1p(exp(df_marker$conservation))
#  m0 <- lm(log1p_conservation ~ specificity, data=df_marker)
#  coef_mat <- rbind(coef(m0), coef_mat)
#}
#saveRDS(gene_list, file="Res_expressed_pan_epi_and_cell_type_specific_marker_of_each_major_epi_cell_type.rds")
gene_list <- readRDS(paste0(folder, "/Res_expressed_pan_epi_and_cell_type_specific_marker_of_each_major_epi_cell_type.rds"))

fetal_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_fetal_tissue_all_cell_class_unified_cell_type_annotation_with_dnds.rds")
# identify epithelial cell type markers
## exclude early enterocyte from DEG identification
cells <- colnames(fetal_rna)[which(fetal_rna$Unified_cell_type != "Early_enterocyte" & fetal_rna$Major_cell_type=="Epithelial")]
epi_subset <- subset(fetal_rna, cells=cells)
saveRDS(epi_subset, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_tissue_epi_excluding_early_enterocyte.rds")
deg_res <- presto::wilcoxauc(epi_subset, group_by="Unified_cell_type")

# add expr. fold change info to the df_tau dataframe
for(ct in names(gene_list)){
  deg_res_sub <- deg_res[which(deg_res$group==ct),]
  df_tau[[paste0(ct, "_logFC")]] <- deg_res_sub[match(rownames(df_tau), deg_res_sub$feature), "logFC"]
}
df_tau$Most_enriched_cell_type <- apply(df_tau[,paste0(names(gene_list), "_logFC")], 1, function(vec){
  sorted_value <- sort(as.numeric(as.matrix(vec)), decreasing=T)
  if(sorted_value[1] - sorted_value[2] > 0.4){
    return(names(gene_list)[which.max(vec)])
  }else{
    return("Undetermined")
  }
})
saveRDS(df_tau, file="/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet/update_cell_type_annotation/Res_gene_dnds_tau_score_ct_logFC.rds")
df_tau <- readRDS("/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet/update_cell_type_annotation/Res_gene_dnds_tau_score_ct_logFC.rds")

other_cell_type_specific_genes <- lapply(names(gene_list), function(ct){
  genes <- gene_list[[ct]]
  df_marker <- df_tau[which(df_tau$feature %in% genes),]
  dnds_ratio_high_cutoff <- quantile(df_marker$conservation, 0.8)
  dnds_ratio_low_cutoff <- quantile(df_marker$conservation, 0.2)
  tau_high_cutoff <- quantile(df_marker[[paste0(ct, "_logFC")]], 0.8)
  tau_low_cutoff <- quantile(df_marker[[paste0(ct, "_logFC")]], 0.2)

  # genes relaxed of negative selection
  df_fast_genes <- df_marker[which(df_marker$conservation > dnds_ratio_high_cutoff),]
  n_fast_genes <- ifelse(nrow(df_fast_genes)>15, 15, nrow(df_fast_genes))
  fast_genes <- df_fast_genes$feature[order(df_fast_genes$conservation, decreasing=T)[1:n_fast_genes]]

  # genes under strong negative selection
  df_slow_genes <- df_marker[which(df_marker$conservation < dnds_ratio_low_cutoff),]
  n_slow_genes <- ifelse(nrow(df_slow_genes)>15, 15, nrow(df_slow_genes))
  slow_genes <- df_slow_genes$feature[order(df_slow_genes$conservation, decreasing=F)[1:n_slow_genes]]

  other_cell_type_specific_genes <- df_marker$feature[which(df_marker$specificity > tau_high_cutoff & df_marker$Most_enriched_cell_type!=ct & df_marker$Most_enriched_cell_type!="Undetermined")]
  gene_to_exclude_from_highlight <- intersect(c(fast_genes, slow_genes), other_cell_type_specific_genes)
  
})
names(other_cell_type_specific_genes) <- names(gene_list)
saveRDS(other_cell_type_specific_genes, file="Res_other_cell_type_specific_genes.rds")

plot_list <- list()
ct_order <- c("Enterocyte", "Stem_cell", "EEC", "Goblet cell", "BEST4+ cell")
for(ct in ct_order){
  genes <- gene_list[[ct]]
  other_cell_type_specific_genes <- df_tau$feature[which(df_tau$Most_enriched_cell_type!=ct & df_tau$Most_enriched_cell_type!="Undetermined")]
  genes <- setdiff(genes, other_cell_type_specific_genes)
  df_marker <- df_tau[which(df_tau$feature %in% genes),]
  dnds_ratio_high_cutoff <- quantile(df_marker$conservation, 0.8)
  dnds_ratio_low_cutoff <- quantile(df_marker$conservation, 0.2)
  tau_high_cutoff <- quantile(df_marker$specificity, 0.8)
  tau_low_cutoff <- quantile(df_marker$specificity, 0.2)
  
  # genes relaxed of negative selection
  df_fast_genes <- df_marker[which(df_marker$conservation > dnds_ratio_high_cutoff),]
  n_fast_genes <- ifelse(nrow(df_fast_genes)>15, 15, nrow(df_fast_genes))
  fast_genes <- df_fast_genes$feature[order(df_fast_genes$conservation, decreasing=T)[1:n_fast_genes]]
  
  # genes under strong negative selection
  df_slow_genes <- df_marker[which(df_marker$conservation < dnds_ratio_low_cutoff),]
  n_slow_genes <- ifelse(nrow(df_slow_genes)>15, 15, nrow(df_slow_genes))
  slow_genes <- df_slow_genes$feature[order(df_slow_genes$conservation, decreasing=F)[1:n_slow_genes]]
  
  highlight_genes <- c(fast_genes, slow_genes)
  p1 <- ggplot(df_marker, aes(x=conservation, y=specificity)) +
    stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(level)), colour = "#303030",
                    bins = 10)+
    scale_fill_distiller(palette = "Reds", direction = 1) +
    geom_vline(xintercept =  dnds_ratio_high_cutoff, linetype="dashed", color = "#303030") +
    geom_vline(xintercept = dnds_ratio_low_cutoff, linetype="dashed", color = "#303030") +
    #geom_hline(yintercept =  tau_high_cutoff, linetype="dashed", color = "#696969") +
    #geom_hline(yintercept = tau_low_cutoff, linetype="dashed", color = "#696969") +
    #geom_text(data=subset(df_marker, feature%in%slow_genes),
    #          aes(label=feature))+
    coord_cartesian(clip = "off") +
    geom_point(data = subset(df_marker, feature%in%highlight_genes))+
    ggrepel::geom_text_repel(data = subset(df_marker, feature%in%highlight_genes),
                             aes(label=feature),
                             max.overlaps=999,
                             #box.padding = 1,
                             show.legend = FALSE) + #this removes the 'a' from the legend
    theme_bw() +
    labs(x="log(dN/dS)", y="Tau", color="", title=paste(ct, "and pan-epithelial"))
  
  plot_list[[ct]] <- p1
  
}
p <- plot_grid(plotlist = plot_list, align = "hv",ncol = 1)  
pdf("Plot_scatterPlot_between_tau_and_dNdS_removing_other_cell_type_specific_genes_from_input_laptop.pdf", height=5*length(plot_list), width=7)
p
dev.off()

coef_mat <- c()
for(ct in c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell", "Stem_cell")){
  genes <- gene_list[[ct]]
  other_cell_type_specific_genes <- df_tau$feature[which(df_tau$Most_enriched_cell_type!=ct & df_tau$Most_enriched_cell_type!="Undetermined")]
  genes <- setdiff(genes, other_cell_type_specific_genes)
  gene_list[[ct]] <- genes
  df_marker <- df_tau[which(df_tau$feature %in% genes),]
  m0 <- lm(specificity ~ log1p_dnds, data=df_marker)
  pred <- predict(m0, newdata = data.frame(log1p_dnds=1))
  coef_mat <- rbind(c(coef(m0), pred), coef_mat)
}
rownames(coef_mat) <- c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell", "Stem_cell")
saveRDS(coef_mat, file="Res_epi_coef_mat_log1p_dnds_and_tau_excluding_other_cell_type_specific_marker.rds")

coor <- df_marker[,c("log1p_dnds","specificity")]

colors <- setNames(
  c("#fba486", "#e2352b", "#f2583f", "#fcc5af", "#c5201e"),
  c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell", "Stem_cell"))


pdf("Plot_slope_dnds_and_tau_excluding_other_cell_type_specific_marker.pdf")
plot(coor, pch=16, type="n", ylim=c(0.4, 0.7),  bty="n", xlab="log(dN/dS + 1)", ylab="Tau")
for(i in seq(nrow(coef_mat))){
  abline(a=coef_mat[i,1], b=coef_mat[i,2], col=colors[rownames(coef_mat)[i]], lwd=2)
}
legend("topleft", legend=rownames(coef_mat)[order(coef_mat[,3], decreasing = T)], fill=colors[rownames(coef_mat)[order(coef_mat[,3], decreasing = T)]], bty="n")
dev.off()



plot_list <- list()
ct_order <- c("Enterocyte", "Stem_cell", "EEC", "Goblet cell", "BEST4+ cell")
for(ct in ct_order){
  genes <- gene_list[[ct]]
  other_cell_type_specific_genes <- df_tau$feature[which(df_tau$Most_enriched_cell_type!=ct & df_tau$Most_enriched_cell_type!="Undetermined")]
  genes <- setdiff(genes, other_cell_type_specific_genes)
  gene_list[[ct]] <- genes
  deg_res_sub <- deg_res[which(deg_res$group==ct),]
  df_tau[[paste0(ct, "_logFC")]] <- deg_res_sub[match(rownames(df_tau), deg_res_sub$feature), "logFC"]
  df_marker <- df_tau[which(df_tau$feature %in% genes),]
  
  df_marker$enrichment <- df_marker[[paste0(ct, "_logFC")]]
  dnds_ratio_high_cutoff <- quantile(df_marker$conservation, 0.8)
  dnds_ratio_low_cutoff <- quantile(df_marker$conservation, 0.2)
  fc_high_cutoff <- quantile(df_marker$enrichment, 0.8)
  fc_low_cutoff <- quantile(df_marker$enrichment, 0.2)
  
  # genes relaxed of negative selection
  df_fast_genes <- df_marker[which(df_marker$conservation > dnds_ratio_high_cutoff),]
  n_fast_genes <- ifelse(nrow(df_fast_genes)>15, 15, nrow(df_fast_genes))
  fast_genes <- df_fast_genes$feature[order(df_fast_genes$conservation, decreasing=T)[1:n_fast_genes]]
  
  # genes under strong negative selection
  df_slow_genes <- df_marker[which(df_marker$conservation < dnds_ratio_low_cutoff),]
  n_slow_genes <- ifelse(nrow(df_slow_genes)>15, 15, nrow(df_slow_genes))
  slow_genes <- df_slow_genes$feature[order(df_slow_genes$conservation, decreasing=F)[1:n_slow_genes]]
  
  highlight_genes <- c(fast_genes, slow_genes)
  p1 <- ggplot(df_marker, aes(x=conservation, y=enrichment)) +
    stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(level)), colour = "#303030",
                    bins = 10)+
    scale_fill_distiller(palette = "Reds", direction = 1) +
    geom_vline(xintercept =  dnds_ratio_high_cutoff, linetype="dashed", color = "#303030") +
    geom_vline(xintercept = dnds_ratio_low_cutoff, linetype="dashed", color = "#303030") +
    #geom_hline(yintercept =  tau_high_cutoff, linetype="dashed", color = "#696969") +
    #geom_hline(yintercept = tau_low_cutoff, linetype="dashed", color = "#696969") +
    #geom_text(data=subset(df_marker, feature%in%slow_genes),
    #          aes(label=feature))+
    coord_cartesian(clip = "off") +
    geom_point(data = subset(df_marker, feature%in%highlight_genes))+
    ggrepel::geom_text_repel(data = subset(df_marker, feature%in%highlight_genes),
                             aes(label=feature),
                             max.overlaps=999,
                             #box.padding = 1,
                             show.legend = FALSE) + #this removes the 'a' from the legend
    theme_bw() +
    labs(x="log(dN/dS)", y="Enrichment", color="", title=ct)
  
  plot_list[[ct]] <- p1
  
}
saveRDS(gene_list, file=paste0(folder, "/Res_expressed_pan_epi_and_cell_type_specific_marker_of_each_major_epi_cell_type_excluding_specific_markers_of_other_cell_types.rds"))
p <- plot_grid(plotlist = plot_list, align = "hv",ncol = 1)  
pdf("Plot_scatterPlot_between_fc_and_dNdS_laptop.pdf", height=5*length(plot_list), width=7)
p
dev.off()


# redo the GO enrichment for fast enterocyte markers
source("/home/yuq22/ihb-intestine-evo/common_script/Script_statiscal_test.R")
source("/home/yuq22/ihb-intestine-evo/colors/colors.R")
gene_list <- readRDS(paste0(folder, "/Res_expressed_pan_epi_and_cell_type_specific_marker_of_each_major_epi_cell_type_excluding_specific_markers_of_other_cell_types.rds"))
genes <- gene_list[["Enterocyte"]]
df_marker <- df_tau[which(df_tau$feature %in% genes),]
dnds_ratio_high_cutoff <- quantile(df_marker$conservation, 0.8)
dnds_ratio_low_cutoff <- quantile(df_marker$conservation, 0.2)
tau_high_cutoff <- quantile(df_marker$specificity, 0.8)
tau_low_cutoff <- quantile(df_marker$specificity, 0.2)

fast_specific_genes <- df_marker$feature[which(df_marker$conservation>dnds_ratio_high_cutoff & df_marker$specificity>tau_high_cutoff)]
fast_pan_genes <- df_marker$feature[which(df_marker$conservation>dnds_ratio_high_cutoff & df_marker$specificity<tau_low_cutoff)]
slow_specific_genes <- df_marker$feature[which(df_marker$conservation<dnds_ratio_low_cutoff & df_marker$specificity>tau_high_cutoff)]
slow_pan_genes <- df_marker$feature[which(df_marker$conservation<dnds_ratio_low_cutoff & df_marker$specificity<tau_low_cutoff)]

gene_sets <- list(
  "fast_specific_genes"=fast_specific_genes,
  "fast_pan_genes"=fast_pan_genes,
  "slow_specific_genes"=slow_specific_genes,
  "slow_pan_genes"=slow_pan_genes
)
gene_folder <- "/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet/update_cell_type_annotation/expression_specificity/epithelial_focused"
saveRDS(gene_sets, file=paste0(gene_folder, "/Res_enterocyte_evo_and_cell_type_marker_group_genes_excluding_specific_markers_of_other_cell_types.rds"))


# get enterocyte markers and epithelial markers that are not significantly depleted in enterocyte
expressed_genes <- readRDS(paste0(gene_folder,"/Res_enterocyte_expressed_genes.rds"))
dnds_df <- readRDS("/pstore/data/ihb-intestine-evo/evo_signature/DN_DS_ratio/data/Dat_human_to_10sp_one2one_orthologs_dnds_ratio.rds")

go_res <- list()
for(gene_group in names(gene_sets)){
  go_res[[gene_group]] <- GO_enrichment_test(deg=gene_sets[[gene_group]], GO_path="/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v109/Ensembl_v109_GO.csv", GO_domain="biological_process", expressed_genes=intersect(expressed_genes, rownames(dnds_df)))

}
saveRDS(go_res, file="Res_GO_biological_process_enrichment_res_for_Enterocyte_expressed_genes_excluding_specific_markers_of_other_cell_types.rds")

# get enriched terms for each gene category
terms <- lapply(names(go_res), function(x){
    df <- go_res[[x]]
    rownames(df)[which(df$"BH_corrected_P"<0.1)]
  })
names(terms) <- names(go_res)
saveRDS(terms, file="Res_ent_exressed_gene_groups_enriched_terms_excl_other_markers.rds")

logp_mat <- sapply(names(go_res), function(x){
    -log10(go_res[[x]][, "Hypogeometric_test_nominal_P"])

})
rownames(logp_mat) <- rownames(go_res[[1]])
saveRDS(logp_mat, file="Res_GO_biological_process_enrichment_pval_for_Enterocyte_expressed_genes_excl_other_marker.rds")

selected_groups <- c("fast_specific_genes", "slow_pan_genes")
data <- logp_mat[unlist(terms[selected_groups]), selected_groups]
data <- data[order(data[,1]-data[,2], decreasing=T),]
#data[,2] <- -data[,2]
input <- sqrt(data)

pdf("Plot_heatmap_Enterocyte_fast_specific_and_slow_pan_genes_enriched_GOs_excl_other_markers.pdf")
gplots::heatmap.2(input,trace="none",main="", key=TRUE, keysize=0.5, density.info="none",dendrogram="none", scale="none", col=grey_scale.heatmap,Rowv=FALSE, cexRow=1, Colv=FALSE, margins = c(8, 26))
dev.off()
