library(Seurat)
library(gplots)
library(ggplot2)
library(splines)
library(doParallel)
library(simspec)
library(dplyr)
library(destiny)
library(patchwork)
source("~/Work/commonScript/Script_functions.R")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte")
# read day 17.5 fetal mouse proximal SI stem cell to enterocyte data 
mouse_se <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/refine_fetal_mouse_annotation/pt/Res_d17.5_fetal_mouse_stem_cell_to_enterocyte.rds")
cells <- colnames(mouse_se)[which(mouse_se$Phase=="G2M" & mouse_se$Cell_type=="Stem_cell")]
mouse_se@meta.data[cells, "Cell_type"] <- "Transit_amplifying_stem_cells"
saveRDS(mouse_se, file="/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/refine_fetal_mouse_annotation/pt/Res_d17.5_fetal_mouse_stem_cell_to_enterocyte.rds")
DimPlot(mouse_se, group.by = "Cell_type")
# read fetal human duodenum stem cell to enterocyte data 
human_se <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/Res_integrated_human_stemCell2Enterocyte.rds")
## extract PCW 11, 12 14 cells that map the mouse cells best
cells <- colnames(human_se)[which(human_se$Age.week%in%c(11,12,14))]
human_mid_se <- subset(human_se, cells=cells)

p1 <- DimPlot(human_se, reduction = "umap_css", group.by = "Age.week")
p2 <- DimPlot(human_mid_se, reduction = "umap_css", group.by = "Age.week")
p1+p2

human_mid_se <- FindVariableFeatures(object = human_mid_se, selection.method = "vst", nfeatures = 3000)
human_mid_se <- ScaleData(object = human_mid_se, verbose = T)
human_mid_se <- RunPCA(object = human_mid_se, features = VariableFeatures(human_mid_se), verbose = F, npcs = 50)
usefulPCs <- 1:20
human_mid_se <- FindNeighbors(object = human_mid_se, dims = usefulPCs)
human_mid_se <- FindClusters(object = human_mid_se, resolution = 1)
human_mid_se <- RunUMAP(object = human_mid_se, dims = usefulPCs)
# CSS integration
human_mid_se <- cluster_sim_spectrum(human_mid_se, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
human_mid_se <- RunUMAP(human_mid_se, reduction = "css", dims = 1:ncol(human_mid_se@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
human_mid_se <- FindNeighbors(object = human_mid_se, reduction = "css", dims = 1:ncol(human_mid_se@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
DimPlot(human_mid_se, reduction = "umap_css", label=T, group.by = "Cell_type")+NoLegend()
input <- human_mid_se@reductions$css@cell.embeddings
dm <- DiffusionMap(input, k=20)
vec <- rank(dm$DC1)
aa <- median(vec[which(human_mid_se$Cell_type=="Stem_cell")]) > median(vec[which(human_mid_se$Cell_type=="Enterocyte")])
if(aa){
  human_mid_se$Pt <- rank(-dm$DC1)
}else{
  human_mid_se$Pt <- rank(dm$DC1)
}
FeaturePlot(human_mid_se, reduction = "umap_css", features = "Pt")
saveRDS(human_mid_se, file="Res_PCW11-14_fetal_human_stemCell2Enterocyte_integrated_with_CSS.rds")

pt_by_age <- lapply(c("80d", "85d", "101d"), function(age){
  human_mid_se$Pt[which(human_mid_se$Age==age)]
})
boxplot(pt_by_age)

expr_mat <- as.matrix(human_mid_se@assays$RNA@data)
pt_vec <- human_mid_se$Pt
g2m_score_vec <- human_mid_se$G2M.Score
s_score_vec <- human_mid_se$S.Score
age_vec <- as.numeric(human_mid_se$Age.week)
registerDoParallel(50)
res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.vector(expr_mat[k,])
  m0 <- lm(e ~ g2m_score_vec + s_score_vec + age_vec)
  m1 <- lm(e ~ g2m_score_vec + s_score_vec + age_vec + ns(pt_vec, df=6))
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0)$Pr[2]
  p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], df1=a0["Residuals", "Df"], df2=a1["Residuals","Df"], lower.tail = F)
  coef <- coef(m1)[4]
  return(c(p_anova, p_resi ,coef))
}
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_Resi", "Coef")
df <- data.frame(res)
saveRDS(df, file="Res_PCW11-14_fetal_human_stem_cell_to_enterocyte_Pt_test_cc_score_and_sample_age_as_covariates.rds")

aa <- getExprByPt(pt.vec=human_mid_se$Pt,
                  expr.mat = human_mid_se@assays$RNA@data,
                  mode = "fix.bin.num",
                  bin.num=62,
                  return.idx = T)
saveRDS(aa, file="Res_PCW11-14_fetal_human_stemCell2Enterocyte_expr_pseudotime_bin_expr.rds")

pt_expr <- aa$expr.mat
expr_sd <- apply(pt_expr, 1, sd)
expressed_df <- df[names(expr_sd)[which(expr_sd>0)],]
expressed_df$p_ANOVA_adj <- p.adjust(expressed_df$p_ANOVA, method="BH")
saveRDS(expressed_df, file="Res_PCW11-14_fetal_human_stem_cell_to_enterocyte_Pt_test_cc_score_as_covariates_for_expressed_genes_with_BH_adjustment.rds")

plot(expressed_df$Coef, -log10(expressed_df$p_ANOVA), pch=16)
pt_genes <- rownames(expressed_df)[which(expressed_df$p_ANOVA_adj<0.05 & abs(expressed_df$Coef)>0.01)]
writeLines(pt_genes, con="List_PCW11-14_fetal_human_Pt_genes.csv")
print(paste("Identify", length(pt_genes), "pseudotime-dependent genes"))

# load mouse pseudotime bin average expression and mouse differentiaion pseudotime dependent genes
mouse_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/refine_fetal_mouse_annotation/pt"
mouse_aa <- readRDS(file.path(mouse_dir, "Res_fetal_mouse_d17.5_prox_stemCell2Enterocyte_expr_pseudotime_bin_expr.rds"))
mouse_pt_expr <- mouse_aa$expr.mat
mouse_pt_genes <- readLines(con=file.path(mouse_dir, "List_d17.5_proxSI_fetal_mouse_stemCell2Enterocyte_pt_dependent_gene.csv"))

human_mouse_orth <- readRDS("~/Work/Annotation/Ensembl/Human/Dat_human_mouse_one2one_symbol_only.rds")
idx1 <- which(human_mouse_orth[,"Human_symbol"] %in% pt_genes & human_mouse_orth[,"Mouse_symbol"] %in% rownames(mouse_pt_expr))
idx2 <- which(human_mouse_orth[,"Human_symbol"] %in% rownames(pt_expr) & human_mouse_orth[,"Mouse_symbol"] %in% mouse_pt_genes)
idx <- union(idx1, idx2)
idx <- intersect(idx1, idx2)
length(idx)
both_pt_orth <- human_mouse_orth[idx,]
saveRDS(both_pt_orth, file="Res_stemCell2Enterocyte_pt_genes_both_human_and_mouse.rds")
human_expr <- pt_expr[both_pt_orth[,"Human_symbol"],]
mouse_expr <- mouse_pt_expr[both_pt_orth[,"Mouse_symbol"],]
rownames(mouse_expr) <- both_pt_orth[,"Human_symbol"]
dim(human_expr)
dim(mouse_expr)
cor_mat <- cor(human_expr, mouse_expr, method="spearman")
heatmap.2(cor_mat, trace="none",
          Rowv=FALSE, Colv=FALSE, dendrogram = "none", col = darkBlue2Red.heatmap)

source("~/Work/commonScript/pt_alignment.r")
expr_ref <- as.matrix(human_mid_se@assays$RNA@data[both_pt_orth[,"Human_symbol"],])
pt_ref <- human_mid_se$Pt/ncol(human_mid_se)
expr_query <- as.matrix(mouse_se@assays$RNA@data[both_pt_orth[,"Mouse_symbol"],])
rownames(expr_query) <- both_pt_orth[,"Human_symbol"]
pt_query <- mouse_se$Pt/ncol(mouse_se)
align_res <- align_pt_traj_with_ref(expr_ref, pt_ref, expr_query, pt_query, # data input
                                    num_breaks_ref = 20, num_breaks_query = 20, dist_method = "spearman", degree = 1, # dtw alignment parameters
                                    ref_name = "human", query_name = "mouse", mode = "global", rev = FALSE, # dtw alignment parameters
                                    nknots_cobs = 20, degree_cobs = 2)
saveRDS(align_res, file="Res_human_mouse_Pt_alignment_global_mode.rds")
heatmap.2(align_res$alignment$diff_mat, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)
aligned_pt_query <- predict(align_res$model, pt_query)
mouse_se$aligned_scToEnt_pt <- aligned_pt_query[,"fit"]
human_mid_se$aligned_scToEnt_pt <- human_mid_se$Pt/ncol(human_mid_se)
saveRDS(human_mid_se, file="Res_human_mid_se_with_Pt_alignment.rds")
saveRDS(mouse_se, file="Res_mouse_se_with_Pt_alignment.rds")

num_breaks <- 20
aligned_pt <- human_mid_se$aligned_scToEnt_pt
expr_mat <- as.matrix(human_mid_se@assays$RNA@data)
human_aligned_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(human_aligned_expr, file="Dat_human_se_aligned_expr.rds")

aligned_pt <- mouse_se$aligned_scToEnt_pt
expr_mat <- as.matrix(mouse_se@assays$RNA@data)
mouse_aligned_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(mouse_aligned_expr, file="Dat_mouse_se_aligned_expr.rds")

# get mouse expression before alignment
pt <- mouse_se$Pt/ncol(mouse_se)
expr_mat <- as.matrix(mouse_se@assays$RNA@data)
mouse_expr_before_alignment <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(mouse_expr_before_alignment, file="Res_d17.5_fetal_mouse_stemCell2Enterocyte_pt_bin_expr_before_alignment.rds")

mouse_expr_before <- mouse_expr_before_alignment[both_pt_orth[,"Mouse_symbol"],]
mouse_expr_after <- mouse_aligned_expr[both_pt_orth[,"Mouse_symbol"],]
rownames(mouse_expr_before) <- rownames(mouse_expr_after) <- both_pt_orth[,"Human_symbol"]
colnames(mouse_expr_before) <- colnames(mouse_expr_after) <- paste("M", seq(num_breaks), sep="_")
human_expr <- human_aligned_expr[both_pt_orth[,"Human_symbol"],]
colnames(human_expr) <- paste("H", seq(num_breaks), sep="_")
cor_before <- cor(human_expr, mouse_expr_before, method="spearman")
heatmap.2(cor_before, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)
cor_after <- cor(human_expr, mouse_expr_after, method="spearman")
heatmap.2(cor_after, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)

plot(mouse_se$Pt/ncol(mouse_se), mouse_se$aligned_scToEnt_pt, pch=16) # basically no change at all
g1 <- c("MKI67", "OLFM4", "LGR5","APOA4", "FABP2", "GSTA1", "GSTA2", "AADAC", "REEP6","DPP4")
mat <- both_pt_orth[which(both_pt_orth[,"Human_symbol"] %in% g1),]
par(mfrow=c(2,3))
for(i in seq(nrow(mat))){
  human_vec <- human_expr[mat[i,"Human_symbol"],]
  mouse_vec <- mouse_expr_before[mat[i,"Human_symbol"],]
  mouse_vec_after <- mouse_expr_after[mat[i,"Human_symbol"],]
  ymax <- max(c(human_vec, mouse_vec, mouse_vec_after))
  ymin <- min(c(human_vec, mouse_vec, mouse_vec_after))
  plot(seq(ncol(mouse_expr_before)), human_vec, type = "l", lwd=2,
       ylim=c(ymin, ymax), main=mat[i,"Human_symbol"], col="red",
       xlab="Pt", ylab="Normed. expr", bty="n")
  lines(seq(ncol(mouse_expr_before)), mouse_vec, lwd=2)
  lines(seq(ncol(mouse_expr_before)), mouse_vec_after, col="blue", lty=2)
}
cols <- c("#d73027","#fc8d59","#fee08b","#1a9850")
ct <- c("Stem_cell", "Transit_amplifying_stem_cells", "Early_enterocyte", "Enterocyte")
gCols <- setNames(cols, ct)
df1 <- data.frame("Cell_index"=seq(ncol(human_mid_se)),
                  "Pt"=human_mid_se$aligned_scToEnt_pt,
                  "Cell_type"=human_mid_se$Cell_type,
                  stringsAsFactors = F)
p1 <- ggplot(df1, aes(x=Cell_index, y=Pt, color=Cell_type)) +
  geom_point()+
  scale_color_manual(values = gCols)+
  theme_minimal()+
  NoLegend()
df2 <- data.frame("Cell_index"=seq(ncol(mouse_se)),
                  "Pt_unaligned"=mouse_se$Pt/ncol(mouse_se),
                  "Pt_aligned"=mouse_se$aligned_scToEnt_pt,
                  "Cell_type"=mouse_se$Cell_type,
                  stringsAsFactors = F)
p2 <- ggplot(df2, aes(x=Cell_index, y=Pt_unaligned, color=Cell_type)) +
  ylab("Pt")+
  geom_point()+
  scale_color_manual(values = gCols)+
  theme_minimal()+
  NoLegend()
p3 <- ggplot(df2, aes(x=Cell_index, y=Pt_aligned, color=Cell_type)) +
  ylab("Pt")+
  geom_point()+
  scale_color_manual(values = gCols)+
  theme_minimal()
p1+p2+p3

par(mfrow=c(1,2), xpd=TRUE)
plotFeature2(coor=Embeddings(human_mid_se, reduction = "umap_css"),
             values = human_mid_se$Cell_type,
             point.order = "random",
             main="Human",
             gCols = gCols,
             add.label = T)
plotFeature2(coor=Embeddings(mouse_se, reduction = "umap"),
             values = mouse_se$Cell_type,
             point.order = "random",
             main="Mouse",
             gCols = gCols,
             add.label = T)

par(mfrow=c(2,4))
for(g in c("Olfm4", "Lgr5", "Mki67", "Cdk1")){
  plotFeature2(coor=Embeddings(mouse_se, reduction = "umap"),
               values = mouse_se@assays$RNA@data[g,],
               point.order = "sorted",
               main=g)
}

orth_mat <- human_mouse_orth[which(human_mouse_orth[,"Human_symbol"]%in%rownames(human_mid_se) & human_mouse_orth[,"Mouse_symbol"]%in%rownames(mouse_se)),]
human_expr <- as.matrix(human_mid_se@assays$RNA@data[orth_mat[,"Human_symbol"],])
mouse_expr <- as.matrix(mouse_se@assays$RNA@data[orth_mat[,"Mouse_symbol"],])
expr_mat <- cbind(human_expr, mouse_expr)
pt_vec <- c(human_mid_se$aligned_scToEnt_pt, mouse_se$aligned_scToEnt_pt)
g2m_score_vec <- c(human_mid_se$G2M.Score, mouse_se$G2M.Score)
s_score_vec <- c(human_mid_se$S.Score, mouse_se$S.Score)
age_vec <- as.factor(c(human_mid_se$Age, mouse_se$Age))
species_vec <- as.factor(rep(c("Human", "Mouse"), c(ncol(human_mid_se), ncol(mouse_se))))
registerDoParallel(50)
res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.vector(expr_mat[k,])
  m0 <- lm(e ~ g2m_score_vec + s_score_vec + age_vec + ns(pt_vec, df=6))
  m1 <- lm(e ~ g2m_score_vec + s_score_vec + age_vec + ns(pt_vec, df=6) + species_vec)
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0)$Pr[2]
  p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], df1=a0["Residuals", "Df"], df2=a1["Residuals","Df"], lower.tail = F)
  coef <- coef(m1)[4]
  return(c(p_anova, p_resi ,coef))
}
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_Resi", "Coef")
df <- data.frame(res)

saveRDS(df, file="Res_PCW11-14_fetal_human_stem_cell_to_enterocyte_versus_fetal_mouse.rds")

# run Pt alignment between human samples to see the correlation heatmap - would it be two blocks or diagonally maximal
human_mid_se_list <- list()
pt_idx <- matrix(F, nrow=nrow(human_mid_se), ncol=3)
rownames(pt_idx) <- rownames(human_mid_se)
colnames(pt_idx) <- unique(human_mid_se$Age)
main_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/stem_cell_to_enterocyte_by_sample"
for(x in unique(human_mid_se$Age)){
  pt_gene <- readLines(file.path(main_dir, x, "List_Pt_genes.csv"))
  human_mid_se_list[[x]] <- readRDS(file.path(main_dir, x, "Res_stem_cell_to_enterocyte.rds"))
  pt_idx[pt_gene,x] <- TRUE
}
pairs <- combn(unique(human_mid_se$Age), 2)
dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte"
css_pt_genes <- readLines(file.path(dir,"List_PCW11-14_fetal_human_Pt_genes.csv"))
main_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/alignment_between_human_samples"
for(j in seq(ncol(pairs))){
  age_1 <- pairs[1,j]
  age_2 <- pairs[2,j]
  sub_dir <- paste(age_1, age_2, "css_integrated", sep="_")
  wd <- file.path(main_dir, sub_dir)
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)
  ref_seu <- subset(human_mid_se, cells=colnames(human_mid_se)[human_mid_se$Age==age_1])
  que_seu <- subset(human_mid_se, cells=colnames(human_mid_se)[human_mid_se$Age==age_2])
  features <- css_pt_genes
  expr_ref <- as.matrix(ref_seu@assays$RNA@data[features,])
  pt_ref <- rank(ref_seu$Pt)/ncol(ref_seu)
  expr_query <- as.matrix(que_seu@assays$RNA@data[features,])
  pt_query <- rank(que_seu$Pt)/ncol(que_seu)
  align_res <- align_pt_traj_with_ref(expr_ref, pt_ref, expr_query, pt_query, # data input
                                      num_breaks_ref = 20, num_breaks_query = 20, dist_method = "spearman", degree = 1, # dtw alignment parameters
                                      ref_name = age_1, query_name = age_2, mode = "global", rev = FALSE, # dtw alignment parameters
                                      nknots_cobs = 20, degree_cobs = 2)
  saveRDS(align_res, file=paste0("Res_fetal_human_",age_1, "_", age_2, "_Pt_alignment_global_mode.rds"))
  pdf(paste0("Plot_heatmap_fetal_human_",age_1, "_", age_2, "_Pt_alignment.pdf"))
  heatmap.2(align_res$alignment$diff_mat, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)
  dev.off()
}

# human and mouse
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css")
human_count <- as.matrix(human_mid_se@assays$RNA@counts)
human_data <- as.matrix(human_mid_se@assays$RNA@data)
human_meta <- human_mid_se@meta.data
mouse_count <- as.matrix(mouse_se@assays$RNA@counts)
mouse_data <- as.matrix(mouse_se@assays$RNA@data)
mouse_meta <- mouse_se@meta.data
## only take genes with one2one orthologs between human and mouse
## do not redo expression normalization
expressed_orth <- human_mouse_orth[human_mouse_orth[,"Human_symbol"]%in%rownames(human_count) & human_mouse_orth[,"Mouse_symbol"]%in%rownames(mouse_count),]
human_count <- human_count[expressed_orth[,"Human_symbol"],]
human_data <- human_data[expressed_orth[,"Human_symbol"],]
mouse_count <- mouse_count[expressed_orth[,"Mouse_symbol"],]
mouse_data <- mouse_data[expressed_orth[,"Mouse_symbol"],]
rownames(mouse_count) <- rownames(mouse_data) <- rownames(human_count)
combined_count <- cbind(human_count, mouse_count)
combined_data <- cbind(human_data, mouse_data)
att <- intersect(colnames(human_meta), colnames(mouse_meta))[c(1:13,16,17,18,19)]
combined_meta <- rbind(human_meta[,att], mouse_meta[,att])
combined_obj <- CreateSeuratObject(counts = combined_count,
                                   meta.data=combined_meta)
combined_obj@assays$RNA@data <- combined_data
## use the Pt defined in both human and mouse as features
VariableFeatures(combined_obj) <- both_pt_orth[,"Human_symbol"]
combined_obj <- ScaleData(object = combined_obj, verbose = T)
combined_obj <- RunPCA(object = combined_obj, features = VariableFeatures(combined_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined_obj <- FindNeighbors(object = combined_obj, dims = usefulPCs)
combined_obj <- FindClusters(object = combined_obj, resolution = 1)
combined_obj <- RunUMAP(object = combined_obj, dims = usefulPCs)
## calculate similarity between each cell to each de novo cluster defined in each sample separately
## central scaled the cluster similarity spectrum in each sample then combine across samples
hvg <- VariableFeatures(combined_obj)
que_expr <- as.matrix(combined_obj@assays$RNA@data[hvg,])
cor_list <- lapply(unique(combined_obj$Age), function(x){
  ref_obj <- subset(combined_obj, cells=colnames(combined_obj)[combined_obj$Age==x]) %>% 
    FindNeighbors(dims = usefulPCs) %>% 
    FindClusters(resolution = 1)
  ref_expr <- getAveExpr(seu.obj=ref_obj, feature.to.calc = "RNA_snn_res.1", colname.prefix = x, genes = hvg)
  cor_mat <- t(scale(t(cor(que_expr, ref_expr, method="spearman"))))
  return(cor_mat)
})
css_mat <- do.call('cbind', cor_list)
combined_obj[['css']] <- CreateDimReducObject(embeddings = css_mat, key = "css_", assay = "RNA")
combined_obj <- RunUMAP(combined_obj, reduction = "css", dims = 1:ncol(combined_obj@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined_obj <- FindNeighbors(object = combined_obj, reduction = "css", dims = 1:ncol(combined_obj@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
saveRDS(combined_obj, file="Res_PCW11-14_fetal_human_and_d17.5_fetal_mouse_stemCell2Enterocyte_integrated_with_CSS.rds")


## Construct Pt on CSS space and compare the transcriptome similarity on pseudotime bin after alignment
cor_mat <- c()
name_vec <- c()
### approach 1: Mouse to human combined, reconstruct Pt in the two species together
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_mouse_integrated")
#### step 1: get pseudotime
input <- css_mat
dm <- DiffusionMap(input, k=20)
vec <- rank(dm$DC1)
aa <- median(vec[which(combined_obj$Cell_type=="Stem_cell")]) > median(vec[which(combined_obj$Cell_type=="Enterocyte")])
if(aa){
  combined_obj$Pt <- rank(-dm$DC1)
}else{
  combined_obj$Pt <- rank(dm$DC1)
}
#FeaturePlot(combined_obj, reduction = "umap", features = "Pt")
#### step 2: run Pt alignment
ref_seu <- subset(combined_obj, cells=colnames(combined_obj)[which(combined_obj$Species=="Mouse")])
que_seu <- subset(combined_obj, cells=colnames(combined_obj)[which(combined_obj$Species=="Human")])
expr_ref <- as.matrix(ref_seu@assays$RNA@data[features,])
pt_ref <- rank(ref_seu$Pt)/ncol(ref_seu)
expr_query <- as.matrix(que_seu@assays$RNA@data[features,])
pt_query <- rank(que_seu$Pt)/ncol(que_seu)
align_res <- align_pt_traj_with_ref(expr_ref, pt_ref, expr_query, pt_query, # data input
                                    num_breaks_ref = 20, num_breaks_query = 20, dist_method = "spearman", degree = 1, # dtw alignment parameters
                                    ref_name = "Human", query_name = "Mouse", mode = "global", rev = FALSE, # dtw alignment parameters
                                    nknots_cobs = 20, degree_cobs = 2)
saveRDS(align_res, file="Res_fetal_human_mouse_Pt_alignment_global_mode.rds")
pdf("Plot_heatmap_fetal_human_mouse_Pt_alignment.pdf")
heatmap.2(align_res$alignment$diff_mat, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)
dev.off()
#### step 3: calculate Pt bin transcriptome similarity before and after alignment
expr_mat <- as.matrix(que_seu@assays$RNA@data)
##### query expression after alignment
aligned_pt_query <- predict(align_res$model, pt_query)
aligned_pt <- aligned_pt_query[,"fit"]
plot(aligned_pt_query[,"fit"], pt_query, pch=16)
cor(aligned_pt_query[,"fit"], pt_query) # no change before and after alignment

que_aligned_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(que_aligned_expr, file="Dat_human_aligned_expr.rds")
##### query expression before alignment
aligned_pt <- pt_query
que_unaligned_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(que_unaligned_expr, file="Dat_human_unaligned_expr.rds")
##### reference expression
expr_mat <- as.matrix(ref_seu@assays$RNA@data)
aligned_pt <- pt_ref
ref_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(ref_expr, file="Dat_mouse_expr.rds")

cor_before <- sapply(seq(num_breaks), function(j){
  cor(ref_expr[features,j], que_unaligned_expr[features,j], method="spearman")
})
cor_after <- sapply(seq(num_breaks), function(j){
  cor(ref_expr[features,j], que_aligned_expr[features,j], method="spearman")
})
cor_mat <- rbind(cor_mat, rbind(cor_before, cor_after))
name_vec <- c(name_vec, c("Integrated_before", "Integrated_after"))

### approach 2: Mouse to each human time point, reconstruct Pt in each sample separately
seu_obj_list <- SplitObject(combined_obj, split.by = "Age")
for(x in names(seu_obj_list)){
  seu_obj <- seu_obj_list[[x]]
  input <- seu_obj@reductions$css@cell.embeddings
  dm <- DiffusionMap(input, k=20)
  vec <- rank(dm$DC1)
  aa <- median(vec[which(seu_obj$Cell_type=="Stem_cell")]) > median(vec[which(seu_obj$Cell_type=="Enterocyte")])
  if(aa){
    seu_obj$Pt_by_age <- rank(-dm$DC1)
  }else{
    seu_obj$Pt_by_age <- rank(dm$DC1)
  }
  seu_obj_list[[x]] <- seu_obj
}

main_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/align_by_age"
setwd(main_dir)
ref_seu <- seu_obj_list[["17.5d"]]
features <- VariableFeatures(combined_obj)
expr_ref <- as.matrix(ref_seu@assays$RNA@data[features,])
pt_ref <- rank(ref_seu$Pt_by_age)/ncol(ref_seu)

num_breaks <- 20
aligned_pt <- pt_ref
expr_mat <- as.matrix(ref_seu@assays$RNA@data)
ref_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(ref_expr, file="Dat_mouse_expr.rds")

for(age in unique(combined_obj$Age[which(combined_obj$Species=="Human")])){
  print(paste(age, "start"))
  wd <- file.path(main_dir, age)
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)
  que_seu <- seu_obj_list[[age]]
  expr_query <- as.matrix(que_seu@assays$RNA@data[features,])
  pt_query <- rank(que_seu$Pt_by_age)/ncol(que_seu)
  align_res <- align_pt_traj_with_ref(expr_ref, pt_ref, expr_query, pt_query, # data input
                                      num_breaks_ref = 20, num_breaks_query = 20, dist_method = "spearman", degree = 1, # dtw alignment parameters
                                      ref_name = "Human", query_name = "Mouse", mode = "global", rev = FALSE, # dtw alignment parameters
                                      nknots_cobs = 20, degree_cobs = 2)
  saveRDS(align_res, file=paste0("Res_",age,"_fetal_human_mouse_Pt_alignment_global_mode.rds"))
  pdf(paste0("Plot_heatmap_",age,"_fetal_human_mouse_Pt_alignment.pdf"))
  heatmap.2(align_res$alignment$diff_mat, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)
  dev.off()
  
  aligned_pt_query <- predict(align_res$model, pt_query)
  aligned_pt <- aligned_pt_query[,"fit"]
  pt_pcc <- cor(aligned_pt, pt_query)
  print(paste("Pt PCC before and after alignment =", pt_pcc))
  
  expr_mat <- as.matrix(que_seu@assays$RNA@data)
  que_aligned_expr <- sapply(1:num_breaks, function(i){
    idx <- which(ceiling(aligned_pt * num_breaks) == i)
    if (i == 1)
      idx <- c(idx, which(aligned_pt == 0))
    if (length(idx) > 1){
      return(rowMeans(expr_mat[,idx]))
    } else if (length(idx) == 1){
      return(expr_mat[,idx])
    } else{
      return(rnorm(nrow(expr_mat)))
    }
  })
  saveRDS(que_aligned_expr, file="Dat_human_aligned_expr.rds")
  
  cor_vec <- sapply(seq(num_breaks), function(j){
    cor(ref_expr[features,j], que_aligned_expr[features,j], method="spearman")
  })
  cor_mat <- rbind(cor_mat, cor_vec)
  name_vec <- c(name_vec, age)
}

### approach 3: Mouse to human combined, reconstruct Pt in each species separately
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined")
seu_obj_list <- SplitObject(combined_obj, split.by = "Species")
for(x in names(seu_obj_list)){
  seu_obj <- seu_obj_list[[x]]
  input <- seu_obj@reductions$css@cell.embeddings
  dm <- DiffusionMap(input, k=20)
  vec <- rank(dm$DC1)
  aa <- median(vec[which(seu_obj$Cell_type=="Stem_cell")]) > median(vec[which(seu_obj$Cell_type=="Enterocyte")])
  if(aa){
    seu_obj$Pt_by_species <- rank(-dm$DC1)
  }else{
    seu_obj$Pt_by_species <- rank(dm$DC1)
  }
  seu_obj_list[[x]] <- seu_obj
}
ref_seu <- seu_obj_list[["Mouse"]]
que_seu <- seu_obj_list[["Human"]]
expr_ref <- as.matrix(ref_seu@assays$RNA@data[features,])
pt_ref <- rank(ref_seu$Pt_by_species)/ncol(ref_seu)
expr_query <- as.matrix(que_seu@assays$RNA@data[features,])
pt_query <- rank(que_seu$Pt_by_species)/ncol(que_seu)
align_res <- align_pt_traj_with_ref(expr_ref, pt_ref, expr_query, pt_query, # data input
                                    num_breaks_ref = 20, num_breaks_query = 20, dist_method = "spearman", degree = 1, # dtw alignment parameters
                                    ref_name = "Human", query_name = "Mouse", mode = "global", rev = FALSE, # dtw alignment parameters
                                    nknots_cobs = 20, degree_cobs = 2)
saveRDS(align_res, file="Res_fetal_human_mouse_Pt_alignment_global_mode.rds")
pdf("Plot_heatmap_fetal_human_mouse_Pt_alignment.pdf")
heatmap.2(align_res$alignment$diff_mat, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=darkBlue2Orange)
dev.off()

aligned_pt_query <- predict(align_res$model, pt_query)
aligned_pt <- aligned_pt_query[,"fit"]
pt_pcc <- cor(aligned_pt, pt_query)
print(paste("Pt PCC before and after alignment =", pt_pcc))
expr_mat <- as.matrix(que_seu@assays$RNA@data)
que_aligned_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(que_aligned_expr, file="Dat_human_aligned_expr.rds")

num_breaks <- 20
aligned_pt <- pt_ref
expr_mat <- as.matrix(ref_seu@assays$RNA@data)
ref_expr <- sapply(1:num_breaks, function(i){
  idx <- which(ceiling(aligned_pt * num_breaks) == i)
  if (i == 1)
    idx <- c(idx, which(aligned_pt == 0))
  if (length(idx) > 1){
    return(rowMeans(expr_mat[,idx]))
  } else if (length(idx) == 1){
    return(expr_mat[,idx])
  } else{
    return(rnorm(nrow(expr_mat)))
  }
})
saveRDS(ref_expr, file="Dat_mouse_expr.rds")

cor_vec <- sapply(seq(num_breaks), function(j){
  cor(ref_expr[features,j], que_aligned_expr[features,j], method="spearman")
})
cor_mat <- rbind(cor_mat, cor_vec)
name_vec <- c(name_vec, "Human_combined")

df <- data.frame("Method"=rep(rownames(cor_mat), ncol(cor_mat)),
                 "Pt"=rep(seq(ncol(cor_mat)), each=nrow(cor_mat)),
                 "SCC"=as.vector(cor_mat),
                 stringsAsFactors = F)

gCols <- setNames(c(c("#fcc5c0","#7a0177"),c("#c6dbef","#6baed6", "#08519c"), "#006837"), rownames(cor_mat))
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css")
pdf("Plot_transcriptome_similarity.pdf", height=7, width=10)
ggplot(df, aes(x=Pt, y=SCC, color=Method))+
  geom_point(size=2)+
  geom_line(lwd=1.2, lty=2)+
  scale_color_manual(values=gCols)
dev.off()

pt_ref <- rank(ref_seu$Pt_by_species)/ncol(ref_seu)
pt_query <- rank(que_seu$Pt_by_species)/ncol(que_seu)
combined_obj@meta.data[colnames(ref_seu), "Pt_by_species"] <- pt_ref
combined_obj@meta.data[colnames(que_seu), "Pt_by_species"] <- pt_query
saveRDS(combined_obj, file="Res_PCW11-14_fetal_human_and_d17.5_fetal_mouse_stemCell2Enterocyte_integrated_with_CSS.rds")

df <- data.frame("Age"=combined_obj$Age,
                 "Pt"=combined_obj$Pt,
                 "Pt_by_species"=combined_obj$Pt_by_species,
                 stringsAsFactors = F)

p1 <- ggplot(df, aes(x=Age, y=Pt))+
  geom_boxplot()+
  ylab("Pt")+
  ggtitle(labe="Integrated")+
  scale_x_discrete(limits=c("17.5d", "101d", "85d", "80d"))
p2 <- ggplot(df, aes(x=Age, y=Pt_by_species))+
  geom_boxplot()+
  ylab("Pt")+
  ggtitle(label="Human_combined")+
  scale_x_discrete(limits=c("17.5d", "101d", "85d", "80d"))

pdf("Plot_boxplot_Pt_distribution_by_sample.pdf", height=5, width=10)
p1+p2
dev.off()

cols <- c("#d73027","#fc8d59","#fee08b","#1a9850")
ct <- c("Stem_cell", "Transit_amplifying_stem_cells", "Early_enterocyte", "Enterocyte")
gCols <- setNames(cols, ct)
seu_obj_list <- SplitObject(combined_obj, split.by = "Species")
que_seu <- seu_obj_list[["Human"]]
ref_seu <- seu_obj_list[["Mouse"]]
df1 <- data.frame("Cell_index"=seq(ncol(que_seu)),
                  "Pt"=que_seu$Pt_by_species,
                  "Cell_type"=que_seu$Cell_type,
                  stringsAsFactors = F)
p1 <- ggplot(df1, aes(x=Cell_index, y=Pt, color=Cell_type)) +
  geom_point(size=0.5)+
  scale_color_manual(values = gCols)+
  theme_minimal()+
  ggtitle(label="Human")+
  NoLegend()
df2 <- data.frame("Cell_index"=seq(ncol(ref_seu)),
                  "Pt"=ref_seu$Pt_by_species,
                  "Cell_type"=ref_seu$Cell_type,
                  stringsAsFactors = F)
p2 <- ggplot(df2, aes(x=Cell_index, y=Pt, color=Cell_type)) +
  ylab("Pt")+
  geom_point(size=0.5)+
  scale_color_manual(values = gCols)+
  ggtitle(label="Mouse")+
  theme_minimal()
pdf("Plot_cell_type_distribution_along_pt.pdf")
p1+p2
dev.off()


cells <- colnames(combined_obj)[combined_obj$Species=="Human"]
combined_obj@meta.data[cells, "Pt_normalized_by_species"] <-  rank(combined_obj@meta.data[cells, "Pt"])/length(cells)
cells <- colnames(combined_obj)[combined_obj$Species=="Mouse"]
combined_obj@meta.data[cells, "Pt_normalized_by_species"] <-  rank(combined_obj@meta.data[cells, "Pt"])/length(cells)
saveRDS(combined_obj, file="Res_PCW11-14_fetal_human_and_d17.5_fetal_mouse_stemCell2Enterocyte_integrated_with_CSS.rds")

png("Plot_UMAP_CSS_Pt.png", height=2000*2, width=2000*3)
par(mfrow=c(2,3), mar=c(5,5,10,5))
plotFeature2(coor=Embeddings(combined_obj, reduction = "umap_css"),
             values = combined_obj$Species,
             point.order = "random",
             main="Species",
             cex=5,
             lwd=0.8,
             cex.main=10,
             add.legend = T,
             legend.cex = 10)
plotFeature2(coor=Embeddings(combined_obj, reduction = "umap_css"),
             values = combined_obj$Cell_type,
             point.order = "random",
             main="Cell_type",
             cex=5,
             lwd=0.8,
             cex.main=10,
             add.legend = T,
             legend.cex = 7)
plotFeature2(coor=Embeddings(combined_obj, reduction = "umap_css"),
             values = combined_obj$G2M.Score,
             point.order = "random",
             main="G2M score",
             cex=5,
             lwd=0.8,
             cex.main=10)

plotFeature2(coor=Embeddings(combined_obj, reduction = "umap_css"),
             values = combined_obj$Pt,
             point.order = "random",
             main="Pt",
             cex=5,
             lwd=0.8,
             cex.main=10)
plotFeature2(coor=Embeddings(combined_obj, reduction = "umap_css"),
             values = combined_obj$Pt_normalized_by_species,
             point.order = "random",
             main="Pt_normed",
             cex=5,
             lwd=0.8,
             cex.main=10)
plotFeature2(coor=Embeddings(combined_obj, reduction = "umap_css"),
             values = combined_obj$Pt_by_species,
             point.order = "random",
             main="Pt_human_combined",
             cex=5,
             lwd=0.8,
             cex.main=10)
dev.off()

genes <- c("AADAC", "REEP6", "FABP2", "DPP4", "APOA4", "OLFM4", "LGR5", "ASCL2", "MKI67", "CDK1")
plotFeature(seu.obj=combined_obj,
            dr="umap_css",
            genes.to.plot = genes,
            col.num = 6)

cols <- c("#d73027","#fc8d59","#fee08b","#1a9850")
ct <- c("Stem_cell", "Transit_amplifying_stem_cells", "Early_enterocyte", "Enterocyte")
gCols <- setNames(cols, ct)
df1 <- data.frame("Cell_index"=seq(ncol(human_mid_se)),
                  "Pt"=human_mid_se$aligned_scToEnt_pt,
                  "Cell_type"=human_mid_se$Cell_type,
                  stringsAsFactors = F)
p1 <- ggplot(df1, aes(x=Cell_index, y=Pt, color=Cell_type)) +
  geom_point()+
  scale_color_manual(values = gCols)+
  theme_minimal()+
  NoLegend()
df2 <- data.frame("Cell_index"=seq(ncol(mouse_se)),
                  "Pt_unaligned"=mouse_se$Pt/ncol(mouse_se),
                  "Pt_aligned"=mouse_se$aligned_scToEnt_pt,
                  "Cell_type"=mouse_se$Cell_type,
                  stringsAsFactors = F)
p2 <- ggplot(df2, aes(x=Cell_index, y=Pt_unaligned, color=Cell_type)) +
  ylab("Pt")+
  geom_point()+
  scale_color_manual(values = gCols)+
  theme_minimal()+
  NoLegend()

# integrated Pt 
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_mouse_integrated")
mouse_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_mouse_integrated/Dat_mouse_expr.rds")
human_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_mouse_integrated/Dat_human_unaligned_expr.rds")
pdf("Plot_selected_genes_across_pt.pdf", height=3*2, width=3*5)
par(mfrow=c(2,5))
for(g in genes){
  human_vec <- human_expr[g,]
  mouse_vec <- mouse_expr[g,]
  ymax <- max(c(human_vec, mouse_vec))
  ymin <- min(c(human_vec, mouse_vec))
  plot(seq(ncol(mouse_expr)), human_vec, type = "l", lwd=2,
       ylim=c(ymin, ymax), main=g, col="dark blue",
       xlab="Pt", ylab="Normed. expr", bty="n")
  lines(seq(ncol(mouse_expr)), mouse_vec, lwd=2, col="dark red")
}
dev.off()
cor_vec_integrated <- sapply(features, function(g){
  cor(mouse_expr[g,], human_expr[g,], method="spearman")
})


# human combined Pt
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined")
mouse_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_mouse_expr.rds")
human_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_human_aligned_expr.rds")
pdf("Plot_selected_genes_across_pt.pdf", height=3*2, width=3*5)
par(mfrow=c(2,5))
for(g in genes){
  human_vec <- human_expr[g,]
  mouse_vec <- mouse_expr[g,]
  ymax <- max(c(human_vec, mouse_vec))
  ymin <- min(c(human_vec, mouse_vec))
  plot(seq(ncol(mouse_expr)), human_vec, type = "l", lwd=2,
       ylim=c(ymin, ymax), main=g, col="dark blue",
       xlab="Pt", ylab="Normed. expr", bty="n")
  lines(seq(ncol(mouse_expr)), mouse_vec, lwd=2, col="dark red")
}
dev.off()

cor_vec_human_combined <- sapply(features, function(g){
  cor(mouse_expr[g,], human_expr[g,], method="spearman")
})
wilcox.test(cor_vec_integrated, cor_vec_human_combined, paired = T)
plot(density(cor_vec_integrated, from=-1, to=1, bw=0.1))
lines(density(cor_vec_human_combined, from=-1, to=1, bw=0.1), col="blue")
# Temporal expression patterns of feature genes along Pt is more similar between species when we use Pt_human_combined (i.e. Pt constructed in each species separately), rather than Pt
# so I choose the Pt_human_combined for downstream analysis

# approach 1: identify genes with differential expression difference along Pt 
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/DEG")
mouse_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_mouse_expr.rds")
human_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_human_aligned_expr.rds")
colnames(mouse_expr) <- paste0("M", seq(ncol(mouse_expr)))
colnames(human_expr) <- paste0("H", seq(ncol(human_expr)))
expr_mat <- cbind(human_expr, mouse_expr)
cor_mat <- cor(expr_mat, method="spearman")
hc <- hclust(as.dist(1-cor_mat), method="ward.D2")
plot(hc, hang=-1)

expr_mat <- human_expr - mouse_expr
genes <- setdiff(rownames(human_expr), features)
cor_vec_human_combined_non_feature <- sapply(genes, function(g){
  cor(human_expr[g,], mouse_expr[g,], method="spearman")
})
boxplot(list("aligned_features"=cor_vec_human_combined, 
             "non_aligned_features"=cor_vec_human_combined_non_feature))
pt_vec <- seq(ncol(expr_mat))
registerDoParallel(50)
p_anova <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'c')%dopar%{
  e <- as.vector(expr_mat[k,])
  res <- anova(lm(e~ns(pt_vec, df=3)))
  p <- res$`Pr(>F)`[1]
  return(p)
}
p_BH <- p.adjust(p_anova, method="BH")
res <- data.frame("p_ANOVA"=p_anova,
                  "p_BH"=p_BH,
                  stringsAsFactors = F)
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_BH")
saveRDS(res, file="Res_PCW11-14_fetal_human_vs_d17.5_fetal_mouse_stem_cell_to_enterocyte_pval.rds")

deg <- rownames(res)[which(res$p_BH<0.001)]
writeLines(deg, con = "List_human_mouse_deg_with_variable_expression_level_difference_along_pt.csv")
deg_pt_expr <- expr_mat[deg,]
saveRDS(deg_pt_expr, file="Dat_human_mouse_deg_expresssion_level_difference_across_pt_bin.rds")
hc <- hclust(as.dist(1-cor(t(deg_pt_expr))))
plot(hc, hang=-1, cex=0.01)
saveRDS(hc, file="Res_human_mouse_DEG_module_hc.rds")

deg_module_id <- cutree(hc, h=1.5)
deg_module <- lapply(sort(unique(deg_module_id)),
                     function(i){names(deg_module_id[which(deg_module_id==i)])})
deg_module_expr_diff <- plotClusterExprProfile(expr=deg_pt_expr, 
                                               time.vec=seq(ncol(deg_pt_expr)), 
                                               group.vec=rep("Human", ncol(deg_pt_expr)), 
                                               cluster.vec=deg_module_id, 
                                               group.cols="#31a354", 
                                               return.value=T, 
                                               to.plot=T, 
                                               plot.name="Plot_human_vs_mouse_stemCell2Enterocyte_pt_dependent_deg_module_average_expr_diff_profile.pdf", 
                                               add.legend=F, 
                                               legend.pos="topleft", 
                                               cex.legend=2, 
                                               col.num=3, 
                                               border.do.smooth=T, 
                                               mean.do.smooth=T, 
                                               df=8,
                                               ylab="Expr. diff.")
saveRDS(deg_module_expr_diff, file="Res_human_vs_mouse_stemCell2Enterocyte_pt_dependent_gene_module_expr_diff.rds")

normed_diff <- t(apply(deg_pt_expr, 1, function(vec){
  (vec-min(vec))/(max(vec)-min(vec))
}))
normed_diff_sum <- colSums(normed_diff)
boxplot(normed_diff)

# plot expression profile in human and mouse separately 
two_sp_expr <- cbind(human_expr, mouse_expr)
two_sp_expr <- two_sp_expr[deg,]
deg_module_expr <- plotClusterExprProfile(expr=two_sp_expr, 
                                          time.vec=rep(seq(ncol(deg_pt_expr)),2), 
                                          group.vec=rep(c("Human","Mouse"), each=ncol(deg_pt_expr)), 
                                          cluster.vec=deg_module_id, 
                                          group.cols=c("#8966A9","#DC3838"), 
                                          return.value=T, 
                                          to.plot=T, 
                                          plot.name="Plot_human_vs_mouse_stemCell2Enterocyte_pt_dependent_deg_module_average_expr_profile.pdf", 
                                          add.legend=T, 
                                          legend.pos="topright", 
                                          cex.legend=2, 
                                          col.num=3, 
                                          border.do.smooth=T, 
                                          mean.do.smooth=T, 
                                          df=8,
                                          ylab="Relative expr.")
saveRDS(deg_module_expr, file="Res_human_vs_mouse_stemCell2Enterocyte_pt_dependent_gene_module_expr.rds")

## GO enrichment analysis with DAVID
load("/home/yuq/Work/Annotation/Ensembl/Human/v93/ensembl.v93.hg38.RData")
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/DEG/GO")
for(i in seq(length(deg_module_expr$g.list))){
  gene_symbol <- deg_module_expr$g.list[[i]]
  gene_id <- unique(ensembl.v93.hg38$gene_id[which(ensembl.v93.hg38$gene_name%in%gene_symbol)])
  writeLines(gene_id, con = paste0("List_expr_diff_module_",i,"_ensemblID.csv"))
}
gene_id <- unique(ensembl.v93.hg38$gene_id[which(ensembl.v93.hg38$gene_name%in%rownames(human_expr))])
writeLines(gene_id, con="List_expressed_human_mouse_one2one_orthologs.csv")
gene_id <- unique(ensembl.v93.hg38$gene_id[which(ensembl.v93.hg38$gene_name%in%rownames(two_sp_expr))])
writeLines(gene_id, con="List_human_mouse_DEG.csv")

# additional approaches for trajectory alignment
## approach 1: dtw - done, and the result is no change before and after alignment
# use additional approaches to align Pt bin which is defiend in each species separately.
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/align_pt_bin")
## approach 2: Seurat label transfer - transfer the pseudotime bin index on mouse trajectory
seu_obj_list <- SplitObject(combined_obj, split.by = "Species")
num_breaks = 20
idx_vec <- setNames(rep(NA, ncol(combined_obj)),
                    colnames(combined_obj))
for(sp in names(seu_obj_list)){
  seu_obj <- seu_obj_list[[sp]]
  pt_bin_idx <- rep(NA, ncol(seu_obj))
  aligned_pt <- seu_obj$Pt_by_species
  for(i in 1:num_breaks){
    idx <- which(ceiling(aligned_pt * num_breaks) == i)
    if (i == 1)
      idx <- c(idx, which(aligned_pt == 0))
    pt_bin_idx[idx] <- i
  }
  idx_vec[colnames(seu_obj)] <- pt_bin_idx
}
combined_obj@meta.data[names(idx_vec), "Pt_by_species_index"] <- idx_vec
DimPlot(combined_obj, group.by = "Pt_by_species_index", reduction = "umap_css")
saveRDS(combined_obj, file="Res_combined_obj_with_Pt_by_species_bin_index.rds")

seu_obj_list <- SplitObject(combined_obj, split.by = "Species")
for(sp in names(seu_obj_list)){
  seu_obj_list[[sp]] <- FindVariableFeatures(seu_obj_list[[sp]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
seu_query <- seu_obj_list[["Mouse"]]
seu_ref <- seu_obj_list[["Human"]]
seu_anchors <- FindTransferAnchors(reference = seu_ref, query = seu_query, 
                                   dims = 1:30)
predictions <- TransferData(anchorset = seu_anchors, refdata = paste0("H",seu_ref$Pt_by_species_index), 
                            dims = 1:30)
seu_query$Pred_Pt_bin <- predictions$predicted.id
DimPlot(seu_query, group.by = "Pred_Pt_bin", reduction = "umap_css")

### try rPCA for human mouse cell co-embeddings - merge human cells of different samples but human and mouse cell keep separated
#library(future)
#library(future.apply)
#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 20000 * 1024^2)
#seu_obj_list <- SplitObject(combined_obj, split.by = "orig.ident")
#bm280k.list <- future_lapply(X = seu_obj_list, FUN = function(x) {
#  x <- NormalizeData(x, verbose = FALSE)
#  x <- FindVariableFeatures(x, verbose = FALSE)
#})
#features <- SelectIntegrationFeatures(object.list = bm280k.list)
#bm280k.list <- future_lapply(X = bm280k.list, FUN = function(x) {
#  x <- ScaleData(x, features = features, verbose = FALSE)
#  x <- RunPCA(x, features = features, verbose = FALSE)
#})
#anchors <- FindIntegrationAnchors(object.list = bm280k.list, reference = c(1, 2), reduction = "rpca", 
#                                  dims = 1:50)
#bm280k.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
#bm280k.integrated <- ScaleData(bm280k.integrated, verbose = FALSE)
#bm280k.integrated <- RunPCA(bm280k.integrated, verbose = FALSE)
#bm280k.integrated <- RunUMAP(bm280k.integrated, dims = 1:50)
#DimPlot(bm280k.integrated, group.by = "orig.ident")


## approach 5: smoothed cell-to-cluster similarity-based matching
### reference option 1 - pt bin expr 
mouse_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_mouse_expr.rds")
human_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_human_aligned_expr.rds")
colnames(human_expr) <- paste0("H", seq(ncol(human_expr)))
colnames(mouse_expr) <- paste0("M", seq(ncol(mouse_expr)))
cl2cl_cor <- cor(mouse_expr[VariableFeatures(combined_obj),], human_expr[VariableFeatures(combined_obj),], method="spearman")
max_id <- setNames(colnames(cl2cl_cor)[apply(cl2cl_cor, 1, which.max)], rownames(cl2cl_cor)) 
max_id
ref_expr <- human_expr[VariableFeatures(combined_obj),]

### reference option 2 - cluster expr
#### use cluster defined in fetal human integrated data under resolution = 1
human_mid_se$integrated_RNA_snn_res.1 <- human_se@meta.data[colnames(human_mid_se), "RNA_snn_res.1"]
p1 <- DimPlot(human_mid_se, reduction = "umap_css", group.by = "integrated_RNA_snn_res.1", label=T)
ref_expr <- getAveExpr(seu.obj=human_mid_se,
                       #feature.to.calc = "RNA_snn_res.0.3",
                       feature.to.calc = "integrated_RNA_snn_res.1",
                       colname.prefix = "H",
                       genes = VariableFeatures(combined_obj))
#saveRDS(ref_expr, file="Dat_ref_expr_human_mid_se_res0.3_cluster_expr_for_human_mouse_pt_genes.rds")
saveRDS(ref_expr, file="Dat_ref_expr_human_mid_se_integrated_res1_cluster_expr_for_human_mouse_pt_genes.rds")

que_expr <- as.matrix(seu_query@assays$RNA@data[VariableFeatures(combined_obj),])
score_mat <- cor(que_expr, ref_expr, method="spearman")
id <- colnames(score_mat)[apply(score_mat, 1, which.max)]
seu_query$singleCell2Cl <- id 
p5 <- DimPlot(seu_query, group.by = "singleCell2Cl", reduction = "umap_css", label = T)

seu_query <- ScaleData(object = seu_query, verbose = T)
seu_query <- RunPCA(object = seu_query, features = VariableFeatures(seu_query), verbose = F, npcs = 20)
seu_query <- FindNeighbors(object = seu_query)
nn_mat <- seu_query@graphs$RNA_nn
smoothed_score_mat <- t(sapply(seq(nrow(nn_mat)), function(i){
  vec <- as.vector(as.matrix(nn_mat[i,]))
  nn_cells <- colnames(nn_mat)[which(vec==1)]
  nn_score <- colMeans(rbind(score_mat[nn_cells,], score_mat[i,]))
}))
rownames(smoothed_score_mat) <- rownames(nn_mat)
id <- colnames(smoothed_score_mat)[apply(smoothed_score_mat, 1, which.max)]
#seu_query$cell2Pt_bin_Pred_Pt_bin <- id # use pt bin expr as reference
seu_query$smoothed_cell2Cl <- id     # use cluster expr as reference


p1 <- DimPlot(seu_ref, group.by = "Pt_by_species_index", reduction = "umap_css", label=T)+NoLegend()
p2 <- DimPlot(seu_query, group.by = "cell2Pt_bin_Pred_Pt_bin", reduction = "umap_css", label = T)+NoLegend()
#p3 <- DimPlot(human_mid_se, group.by = "RNA_snn_res.0.3", reduction = "umap_css", label=T)+NoLegend()
p3 <- DimPlot(human_mid_se, group.by = "integrated_RNA_snn_res.1", reduction = "umap_css", label=T)+NoLegend()
p4 <- DimPlot(seu_query, group.by = "smoothed_cell2Cl", reduction = "umap_css", label = T)+NoLegend()
p3+p4+p5+plot_layout(ncol=2)
p1+p2+p3+p4+plot_layout(ncol=2)
FeaturePlot(human_mid_se, reduction = "umap_css", features = c("MKI67", "CDK1", "OLFM4", "LGR5"), order = T)
p7+p3
## approach 3: CSS-based MCMF, compare the Pt_species for matched human and mouse cells
library(reticulate)
source_python("~/Work/commonScript/matching.py")
css_mat <- combined_obj@reductions$css@cell.embeddings
ref_cells <- colnames(combined_obj)[which(combined_obj$Species=="Human")]
que_cells <- colnames(combined_obj)[which(combined_obj$Species=="Mouse")]
cost_graph <- get_cost_knn_graph(
  source = css_mat[ref_cells,],
  target = css_mat[que_cells,],
  knn_k = 5,
  knn_n_jobs = 10,
  null_cost_percentile = 99,
  capacity_method = "uniform"
)
grp_idx <- setNames(data.frame(do.call(cbind, mcmf(cost_graph))+1), c("Human","Mouse"))
grp_cell <- data.frame("Human"=ref_cells[grp_idx[,"Human"]], 
                       "Mouse"=que_cells[grp_idx[,"Mouse"]])
idx <- which(!is.na(grp_cell[,2]))
grp_cell <- grp_cell[idx,]
saveRDS(grp_cell, file="Res_MCMF_grp_cell.rds")

### approach 4: Get kNN (k=20) human fetal cells for mouse cell on the CSS space
ref_mat <- combined_obj@reductions$css@cell.embeddings[combined_obj$Species=="Human",]
que_mat <- combined_obj@reductions$css@cell.embeddings[combined_obj$Species=="Mouse",]
knn <- RANN::nn2(ref_mat, que_mat, k = 20)$nn.idx
dim(knn)
#ref_idx <- human_mid_se@meta.data[rownames(ref_mat), "RNA_snn_res.0.3"]
ref_idx <- paste0("T",combined_obj@meta.data[rownames(ref_mat), "Pt_by_species_index"])
nn_idx <- matrix(ref_idx[as.vector(knn)], nrow=nrow(knn))
pred_id <- apply(nn_idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
#seu_query@meta.data[rownames(que_mat), "kNN_pred_id"] <- pred_id
seu_query@meta.data[rownames(que_mat), "kNN_pred_Pt"] <- pred_id
DimPlot(seu_query, group.by = "kNN_pred_Pt", reduction = "umap_css", label=T)
seu_query$Pt_by_species_index <- paste0("T", seu_query$Pt_by_species_index)
n1 <- sapply(sort(unique(seu_query$Pt_by_species_index)), function(real_idx){
  sapply(sort(unique(seu_query$kNN_pred_Pt)), function(pred_idx){
    sum(seu_query$Pt_by_species_index==real_idx & seu_query$kNN_pred_Pt==pred_idx)
  })
})

num_mat <- sapply(sort(unique(human_mid_se$RNA_snn_res.0.3)), function(i){
  n1 <- sum(human_mid_se$RNA_snn_res.0.3==i)
  n2 <- sum(seu_query$kNN_pred_id==i)
  return(c(n1,n2))
})
rownames(num_mat) <- c("Human","Mouse")
colnames(num_mat) <- paste0("H",sort(unique(human_mid_se$RNA_snn_res.0.3)))
prop_mat <- num_mat/rowSums(num_mat)
df <- data.frame("Species"=rep(rownames(prop_mat), ncol(prop_mat)),
                 "Cluster"=rep(colnames(prop_mat), each=nrow(prop_mat)),
                 "Proportion"=as.vector(prop_mat),
                 stringsAsFactors = F)
p1 <- ggplot(df, aes(x=Species,y=Proportion,fill=Cluster))+
  geom_bar(stat="identity", width=0.5)
p1

pdf("Plot_Pt_correspondence.pdf")
plot(x=combined_obj@meta.data[grp_cell[,"Human"], "Pt_by_species"],
     y=combined_obj@meta.data[grp_cell[,"Mouse"], "Pt_by_species"],
     pch=16,
     xlab="Human",
     ylab="Mouse",
     bty="n")
abline(a=0, b=1, lty=2)
dev.off()

n1 <- sapply(sort(unique(combined_obj$Pt_by_species_index)), function(idx_h){
  sapply(sort(unique(combined_obj$Pt_by_species_index)), function(idx_m){
    cell_h <- colnames(combined_obj)[which(combined_obj$Pt_by_species_index==idx_h & combined_obj$Species=="Human")]
    cell_m <- colnames(combined_obj)[which(combined_obj$Pt_by_species_index==idx_m & combined_obj$Species=="Mouse")]
    sum(grp_cell[,"Human"]%in%cell_h & grp_cell[,"Mouse"]%in%cell_m)
  })
})
colnames(n1) <- paste0("H",sort(unique(combined_obj$Pt_by_species_index)))
rownames(n1) <- paste0("M",sort(unique(combined_obj$Pt_by_species_index)))


n2 <- sapply(sort(unique(combined_obj$Cell_type)), function(idx_h){
  sapply(sort(unique(combined_obj$Cell_type)), function(idx_m){
    cell_h <- colnames(combined_obj)[which(combined_obj$Cell_type==idx_h & combined_obj$Species=="Human")]
    cell_m <- colnames(combined_obj)[which(combined_obj$Cell_type==idx_m & combined_obj$Species=="Mouse")]
    sum(grp_cell[,"Human"]%in%cell_h & grp_cell[,"Mouse"]%in%cell_m)
  })
})
colnames(n2) <- paste0("H-",sort(unique(combined_obj$Cell_type)))
rownames(n2) <- paste0("M-",sort(unique(combined_obj$Cell_type)))


id <- human_mid_se@meta.data[grp_cell[,"Human"],"RNA_snn_res.0.3"]
seu_1 <- subset(human_mid_se, cells=unique(grp_cell[,"Human"]))
seu_2 <- subset(seu_query, cells=unique(grp_cell[,"Mouse"]))
seu_2@meta.data[grp_cell[,"Mouse"], "MCMF_pred_id"] <- id 
DimPlot(seu_2, reduction = "umap_css", group.by = "MCMF_pred_id", label=T)
p5 <- DimPlot(seu_1, group.by = "RNA_snn_res.0.3", reduction = "umap_css", label=T)+NoLegend()
p6 <- DimPlot(seu_2, group.by = "MCMF_pred_id", reduction = "umap_css", label = T)+NoLegend()
saveRDS(seu_query, file="Res_mouse_with_pred_human_ident.rds")

num_mat <- sapply(sort(unique(human_mid_se$RNA_snn_res.0.3)), function(i){
  n1 <- sum(human_mid_se$RNA_snn_res.0.3==i)
  n2 <- sum(seu_2$MCMF_pred_id==i)
  return(c(n1,n2))
})
rownames(num_mat) <- c("Human","Mouse")
colnames(num_mat) <- paste0("H",sort(unique(human_mid_se$RNA_snn_res.0.3)))
prop_mat <- num_mat/rowSums(num_mat)
df <- data.frame("Species"=rep(rownames(prop_mat), ncol(prop_mat)),
                 "Cluster"=rep(colnames(prop_mat), each=nrow(prop_mat)),
                 "Proportion"=as.vector(prop_mat),
                 stringsAsFactors = F)
p1 <- ggplot(df, aes(x=Species,y=Proportion,fill=Cluster))+
  geom_bar(stat="identity", width=0.5)



# align human-mouse cell populations on the stem-cell-to-enterocyte trajectory
aligned_ct <- setNames(rep(NA, ncol(combined_obj)), colnames(combined_obj))
aligned_ct[colnames(human_mid_se)] <- paste("H", human_mid_se$integrated_RNA_snn_res.1, sep="_")
aligned_ct[colnames(seu_query)] <- seu_query$singleCell2Cl
combined_obj$Aligned_cl <- aligned_ct
p8 <- DimPlot(combined_obj, reduction = "umap_css", group.by = "Species")
saveRDS(combined_obj, file="Res_PCW11-14_fetal_human_and_d17.5_fetal_mouse_stemCell2Enterocyte_integrated_with_CSS_with_ref_cl_defined_in_integrated_fetal_human.rds")
p9 <- DimPlot(combined_obj, reduction = "umap_css", group.by = "Aligned_cl")
p1+p3+p5+p2+p4+p6+p8+p7+p9+plot_layout(ncol=3)

p1 <- DimPlot(seu_2, reduction = "umap_css", group.by = "Cell_type")+NoLegend()
id <- human_mid_se@meta.data[grp_cell[,"Human"],"Cell_type"]
seu_2@meta.data[grp_cell[,"Mouse"], "MCMF_pred_ct"] <- id 
p2 <- DimPlot(seu_2, reduction = "umap_css", group.by = "MCMF_pred_ct")
p1+p2
file="/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/refine_fetal_mouse_annotation/pt/Res_d17.5_fetal_mouse_stem_cell_to_enterocyte.rds"
pca_mouse_se <- readRDS(file)
coor <- Embeddings(pca_mouse_se, reduction = "umap")
DimPlot(pca_mouse_se, group.by = "Cell_type")
seu_2[['pca_umap']] <- CreateDimReducObject(embeddings=coor[colnames(seu_2),],
                                            key = "PCAUMAP_")
p3 <- DimPlot(seu_2, reduction = "pca_umap", group.by = "MCMF_pred_ct")+NoLegend()
p4 <- DimPlot(seu_2, reduction = "pca_umap", group.by = "Cell_type")
p3+p4
pca_mouse_se$Aligned_cl <- combined_obj@meta.data[colnames(pca_mouse_se), "Aligned_cl"]
p1 <- FeaturePlot(pca_mouse_se, features = "Mki67", order = T)
p2 <- FeaturePlot(human_mid_se, features = "MKI67", order = T, reduction = "umap_css")
p3 <- DimPlot(pca_mouse_se, group.by = "Aligned_cl", label=T)+NoLegend()
p4 <- DimPlot(human_mid_se, reduction = "umap_css", group.by = 'RNA_snn_res.0.3', label=T)+NoLegend()
p5 <- DimPlot(pca_mouse_se, group.by = "Cell_type")
p6 <- DimPlot(human_mid_se, group.by = "Cell_type", reduction = "umap_css")
p1+p3+p5+p2+p4+p6+plot_layout(ncol=3)

n1 <- sapply(sort(unique(human_se$RNA_snn_res.1)), function(i){
  sapply(paste0(sort(as.numeric(sub("d", "", unique(human_se$Age)))),"d"), function(x){
    sum(human_se$Age==x & human_se$RNA_snn_res.1==i)
  })
})
vec <- sapply(paste("H",sort(unique(human_se$RNA_snn_res.1)),sep="_"), function(i){
  sum(combined_obj$Aligned_cl==i & combined_obj$Species=="Mouse")
})
colnames(n1) <- paste("H",sort(unique(human_se$RNA_snn_res.1)),sep="_")
n1 <- rbind(n1, vec)
rownames(n1)[nrow(n1)] <- "17.5d"
p1 <- n1/rowSums(n1)

df <- data.frame("Age"=rep(rownames(p1), ncol(p1)),
                 "Cluster"=rep(colnames(p1), each=nrow(p1)),
                 "Proportion"=as.vector(p1),
                 stringsAsFactors = F)
p1 <- ggplot(df, aes(x=Age,y=Proportion,fill=Cluster))+
  geom_bar(stat="identity", width=0.5)+
  scale_x_discrete(limits=rownames(n1))
p2 <- DimPlot(human_se, reduction = "umap_css", label=T,group.by = 'RNA_snn_res.1')
p1+p2

# approach 2 of DEG identification: identify DEG between species, including the ones with global expression level difference 
# --> then further identify those with Pt-dependent expression level difference from those DEG
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/DEG_control_ct")
expr_mat <- as.matrix(combined_obj@assays$RNA@data)
sp_vec <- as.factor(combined_obj$Species)
pt_vec <- combined_obj$Pt_by_species_index
g2m_score_vec <- combined_obj$G2M.Score
s_score_vec <- combined_obj$S.Score
registerDoParallel(50)
res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.vector(expr_mat[k,])
  m0 <- lm(e ~ g2m_score_vec+s_score_vec+ns(pt_vec, df=3))
  m1 <- lm(e ~ g2m_score_vec+s_score_vec+ns(pt_vec, df=3)+sp_vec)
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0)$Pr[2]
  p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], 
               df1=a0["Residuals", "Df"], 
               df2=a1["Residuals","Df"], 
               lower.tail = F)
  coef <- coef(m1)[7]
  return(c(p_anova, p_resi ,coef))
}
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_Resi", "Coef")
df <- data.frame(res)
df$p_ANOVA_BH <- p.adjust(df$p_ANOVA, method="BH")
saveRDS(df, file="Res_human_mouse_SE_with_Pt_bin_index_and_cc_score_as_covariates.rds")
stopImplicitCluster()

plot(df$Coef, -log10(df$p_ANOVA), pch=16)

deg <- rownames(df)[which(df$p_ANOVA_BH<0.01)]
writeLines(deg, con = "List_overall_human_mouse_s2e_DEG.csv")

# focus on genes with Pt-dependent expression pattern in either human or mouse
## human
cells <- which(combined_obj$Species=="Human")
expr_mat <- as.matrix(combined_obj@assays$RNA@data[,cells])
pt_vec <- combined_obj$Pt_by_species[cells]
g2m_score_vec <- combined_obj$G2M.Score[cells]
s_score_vec <- combined_obj$S.Score[cells]
registerDoParallel(50)
res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.vector(expr_mat[k,])
  m0 <- lm(e ~ g2m_score_vec+s_score_vec)
  m1 <- lm(e ~ g2m_score_vec+s_score_vec+ns(pt_vec, df=6))
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0)$Pr[2]
  p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], 
               df1=a0["Residuals", "Df"], 
               df2=a1["Residuals","Df"], 
               lower.tail = F)
  return(c(p_anova, p_resi))
}
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_Resi")
df <- data.frame(res)
df$p_ANOVA_BH <- p.adjust(df$p_ANOVA, method="BH")
saveRDS(df, file="Res_human_Pt-dependent_gene_with_cc_score_as_covariates.rds")
stopImplicitCluster()
human_pt <- df

## mouse
cells <- which(combined_obj$Species=="Mouse")
expr_mat <- as.matrix(combined_obj@assays$RNA@data[,cells])
pt_vec <- combined_obj$Pt_by_species[cells]
g2m_score_vec <- combined_obj$G2M.Score[cells]
s_score_vec <- combined_obj$S.Score[cells]
registerDoParallel(50)
res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.vector(expr_mat[k,])
  m0 <- lm(e ~ g2m_score_vec+s_score_vec)
  m1 <- lm(e ~ g2m_score_vec+s_score_vec+ns(pt_vec, df=6))
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0)$Pr[2]
  p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], 
               df1=a0["Residuals", "Df"], 
               df2=a1["Residuals","Df"], 
               lower.tail = F)
  return(c(p_anova, p_resi))
}
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_Resi")
mouse_pt <- data.frame(res)
mouse_pt$p_ANOVA_BH <- p.adjust(mouse_pt$p_ANOVA, method="BH")
saveRDS(mouse_pt, file="Res_mouse_Pt-dependent_gene_with_cc_score_as_covariates.rds")
stopImplicitCluster()
human_pt_genes <- rownames(human_pt)[which(human_pt$p_ANOVA_BH<0.01)]
mouse_pt_genes <- rownames(mouse_pt)[which(mouse_pt$p_ANOVA_BH<0.01)]
writeLines(human_pt_genes, con="List_human_pt_genes.csv")
writeLines(mouse_pt_genes, con="List_mouse_pt_genes.csv")
union_pt_genes <- union(mouse_pt_genes, human_pt_genes)
length(union_pt_genes)
deg_pt_genes <- intersect(deg, union_pt_genes)
length(deg_pt_genes)
writeLines(deg_pt_genes, con="List_human_mouse_DE_union_Pt-dependent_genes.csv")

mouse_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_mouse_expr.rds")
human_expr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/human_combined/Dat_human_aligned_expr.rds")
colnames(mouse_expr) <- paste0("M", seq(ncol(mouse_expr)))
colnames(human_expr) <- paste0("H", seq(ncol(human_expr)))
expr_mat <- cbind(human_expr, mouse_expr)
deg_expr <- expr_mat[deg_pt_genes,]
dim(deg_expr)
hc <- hclust(as.dist(1-cor(t(deg_expr))))
plot(hc, hang=-1, cex=0.01)
saveRDS(hc, file="Res_human_mouse_DEG_module_hc.rds")

deg_module_id <- cutree(hc, h=1.5)
deg_module <- lapply(sort(unique(deg_module_id)),
                     function(i){names(deg_module_id[which(deg_module_id==i)])})
deg_module_expr <- plotClusterExprProfile(expr=deg_expr, 
                                          time.vec=c(seq(ncol(human_expr)),seq(ncol(mouse_expr))), 
                                          group.vec=rep(c("Human","Mouse"), each=ncol(human_expr)), 
                                          cluster.vec=deg_module_id, 
                                          group.cols=c("#DC3838","#8966A9"), 
                                          return.value=T, 
                                          to.plot=T, 
                                          plot.name="Plot_human_vs_mouse_stemCell2Enterocyte_pt_dependent_deg_module_average_expr_profile.pdf", 
                                          add.legend=T, 
                                          legend.pos="topright", 
                                          cex.legend=2, 
                                          col.num=3, 
                                          border.do.smooth=T, 
                                          mean.do.smooth=T, 
                                          df=8,
                                          ylab="Relative expr")
saveRDS(deg_module_expr, file="Res_human_vs_mouse_stemCell2Enterocyte_pt_dependent_gene_module_expr.rds")

# further distinguish genes with Pt-dependent expression level difference
expr_mat <- human_expr[deg_pt_genes,] - mouse_expr[deg_pt_genes,]
dim(expr_mat)
pt_vec <- seq(ncol(expr_mat))
registerDoParallel(50)
p_anova <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'c')%dopar%{
  e <- as.vector(expr_mat[k,])
  res <- anova(lm(e~ns(pt_vec, df=3)))
  p <- res$`Pr(>F)`[1]
  return(p)
}
p_BH <- p.adjust(p_anova, method="BH")
res <- data.frame("p_ANOVA"=p_anova,
                  "p_BH"=p_BH,
                  stringsAsFactors = F)
rownames(res) <- rownames(expr_mat)
colnames(res) <- c("p_ANOVA", "p_BH")
saveRDS(res, file="Res_human_mouse_expr_dif_age_test_pval.rds")

expr_dif_pt_genes <- rownames(res)[which(res$p_BH<0.001)]
length(expr_dif_pt_genes)
writeLines(expr_dif_pt_genes, con = "List_human_mouse_deg_with_variable_expression_level_difference_along_pt.csv")

# fold-change of human-mouse DEGs in different time points  
#
