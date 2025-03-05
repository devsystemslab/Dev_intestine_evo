library(Seurat)
library(destiny)
library(dplyr)
library(doParallel)
library(splines)
library(gplots)
library(patchwork)
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/stem_cell_to_enterocyte_by_sample")
human_epi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/Res_fetal_8-19_PCW_hg38_epi_with_CSS_integration.rds")
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi")
human_se <- subset(human_epi, cells=colnames(human_epi)[human_epi$Cell_type %in% c("Stem_cell", "Transit_amplifying_stem_cells", "Early_enterocyte", "Enterocyte")])
saveRDS(human_se, file="Res_integrated_human_stemCell2Enterocyte.rds")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/CSS_based_Pt")
human_se <- FindVariableFeatures(object = human_se, selection.method = "vst", nfeatures = 3000)
human_se <- ScaleData(object = human_se, verbose = T)
human_se <- RunPCA(object = human_se, features = VariableFeatures(human_se), verbose = F, npcs = 50)
usefulPCs <- 1:20
human_se <- FindNeighbors(object = human_se, dims = usefulPCs)
human_se <- FindClusters(object = human_se, resolution = 2)
human_se <- RunUMAP(object = human_se, dims = usefulPCs)
# run CSS integration
human_se <- cluster_sim_spectrum(human_se, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
human_se <- RunUMAP(human_se, reduction = "css", dims = 1:ncol(human_se@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
human_se <- FindNeighbors(object = human_se, reduction = "css", dims = 1:ncol(human_se@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 0.3)
input <- human_se@reductions$css@cell.embeddings
dm <- DiffusionMap(input, k=20)
vec <- rank(dm$DC1)
aa <- median(vec[which(human_se$Cell_type=="Stem_cell")]) > median(vec[which(human_se$Cell_type=="Enterocyte")])
if(aa){
  human_se$Pt <- rank(-dm$DC1)
}else{
  human_se$Pt <- rank(dm$DC1)
}
FeaturePlot(human_se, reduction = "umap_css", features = "Pt")
saveRDS(human_se, file="Res_fetal_human_stemCell2Enterocyte_integrated_with_CSS.rds")


DimPlot(human_se, group.by = "Cell_type", reduction = "umap_css")
human_se_list <- SplitObject(human_se, split.by = "Age")
main_dir <- getwd()

start_bin_idx <- 3
for(x in names(human_se_list)){
  
  print(paste(x, "start"))
  sub_dir <- x
  if(!file.exists(file.path(main_dir, sub_dir))){
    dir.create(file.path(main_dir, sub_dir))
  }
  setwd(file.path(main_dir, sub_dir))
  seu_obj <- human_se_list[[x]]
  
  # preprocessing
  print("Seurat preprocessing")
  seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
  seu_obj <- ScaleData(object = seu_obj, verbose = T)
  seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
  usefulPCs <- 1:20
  seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
  seu_obj <- FindClusters(object = seu_obj, resolution = 2)
  seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
  
  # get pseudotime
  print("Get pseudotime")
  input <- seu_obj@reductions$pca@cell.embeddings
  dm <- DiffusionMap(input, k=20)
  vec <- rank(dm$DC1)
  aa <- median(vec[which(seu_obj$Cell_type=="Stem_cell")]) > median(vec[which(seu_obj$Cell_type=="Enterocyte")])
  if(aa){
    seu_obj$Pt <- rank(-dm$DC1)
  }else{
    seu_obj$Pt <- rank(dm$DC1)
  }
  saveRDS(seu_obj, file="Res_stem_cell_to_enterocyte.rds")
  human_se_list[[x]] <- seu_obj
  
  # identify stem cell-to-enterocyte pseudotime dependent genes, using cell cycle phase score as covariates
  print("Get pseudotime-dependent genes")
  expr_mat <- as.matrix(seu_obj@assays$RNA@data)
  pt_vec <- seu_obj$Pt
  g2m_score_vec <- seu_obj$G2M.Score
  s_score_vec <- seu_obj$S.Score
  registerDoParallel(50)
  res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.vector(expr_mat[k,])
    m0 <- lm(e ~ g2m_score_vec + s_score_vec )
    m1 <- lm(e ~ g2m_score_vec + s_score_vec + ns(pt_vec, df=6))
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
  saveRDS(df, file="Res_stem_cell_to_enterocyte_Pt_test_cc_score_as_covariates.rds")
  
  aa <- getExprByPt(pt.vec=seu_obj$Pt,
                    expr.mat = seu_obj@assays$RNA@data,
                    mode = "fix.cell.num",
                    cell.num.per.bin=20,
                    return.idx = T)
  saveRDS(aa, file="Res_stemCell2Enterocyte_expr_pseudotime_bin_expr.rds")
  
  pt_expr <- aa$expr.mat
  expr_sd <- apply(pt_expr, 1, sd)
  expressed_df <- df[names(expr_sd)[which(expr_sd>0)],]
  expressed_df$p_ANOVA_adj <- p.adjust(expressed_df$p_ANOVA, method="BH")
  saveRDS(expressed_df, file="Res_stem_cell_to_enterocyte_Pt_test_cc_score_as_covariates_for_expressed_genes_with_BH_adjustment.rds")
  
  pt_genes <- rownames(expressed_df)[which(expressed_df$p_ANOVA_adj<0.05 & abs(expressed_df$Coef)>0.1)]
  writeLines(pt_genes, con="List_Pt_genes.csv")
  print(paste("Identify", length(pt_genes), "pseudotime-dependent genes"))
  
  # identify genes with significant turning point in spline smoothed expression profiles
  print("Get expression pattern turning point")
  expr_mat <- pt_expr[pt_genes,]
  pt_vec <- seq(ncol(expr_mat))
  res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.vector(expr_mat[k,])
    m0 <- lm(e ~ pt_vec)
    a0 <- anova(m0)
    
    p_vec <- sapply(3:(ncol(expr_mat)-3), function(j){
      m1_1 <- lm(e[1:j] ~ pt_vec[1:j]) 
      m1_2 <- lm(e[j:ncol(expr_mat)] ~ pt_vec[j:ncol(expr_mat)])
      a1_1 <- anova(m1_1)
      a1_2 <- anova(m1_2)
      resi_sum_sq_1 <- a1_1$`Sum Sq`[2]+a1_2$`Sum Sq`[2]
      resi_df_1 <- a1_1$Df[2]+a1_2$Df[2]
      resi_mean_sq_1 <- resi_sum_sq_1/resi_df_1
      p <- pf(a0$`Mean Sq`[2]/resi_mean_sq_1, 
              df1=a0$Df[2], 
              df2=resi_df_1,
              lower.tail = F) 
      return(p)
    })
    return(p_vec)
    
  }
  stopImplicitCluster()
  rownames(res) <- rownames(expr_mat)
  colnames(res) <- paste0("Pt-bin_", start_bin_idx:(ncol(expr_mat)-start_bin_idx))
  saveRDS(res, file="Res_stem_cell_to_enterocyte_turning_point_pval.rds")
  padj <- t(apply(res, 1, p.adjust, method="bonferroni"))
  saveRDS(padj, file="Res_stem_cell_to_enterocyte_turning_point_padj.rds")
  min_pval_idx <- colnames(res)[apply(res, 1, which.min)]
  freq <- table(min_pval_idx)
  a <- setNames(rep(0, ncol(padj)), colnames(padj))
  a[names(freq)] <- freq
  
  pdf("Plot_most_significant_turning_point_distribution.pdf")
  plot(seq(length(a))+(start_bin_idx-1), a, pch=16, bty="n", xlab="Pt bin", ylab="Turning point number")
  lines(smooth.spline(seq(length(a))+(start_bin_idx-1), a, df=5))
  dev.off()
  
}


# compare pseudotime across samples
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/stem_cell_to_enterocyte_by_sample/integrate")
main_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/stem_cell_to_enterocyte_by_sample"

for(age in names(human_se_list)){
  seu_obj <- human_se_list[[age]]
  aa <- readRDS(file.path(main_dir, age, "Res_stemCell2Enterocyte_expr_pseudotime_bin_expr.rds"))
  seu_obj$Pt_bin_idx <- aa$cell.idx.vec
  human_se_list[[age]] <- seu_obj
}
saveRDS(human_se_list, file="Res_human_se_list_by_sample.rds")

# get pt-genes defined in more than 1 sample
pt_idx_mat <- matrix(F, nrow=nrow(human_se), ncol=length(human_se_list))
rownames(pt_idx_mat) <- rownames(human_se)
colnames(pt_idx_mat) <- names(human_se_list)
for(age in names(human_se_list)){
  pt_genes <- readLines(file.path(main_dir, age, "List_Pt_genes.csv"))
  pt_idx_mat[pt_genes, age] <- T
}
union_pt_genes <- rownames(pt_idx_mat)[which(rowSums(pt_idx_mat)>1)]
expr_list <- lapply(names(human_se_list), function(age){
  aa <- readRDS(file.path(main_dir, age, "Res_stemCell2Enterocyte_expr_pseudotime_bin_expr.rds"))
  expr_mat <- aa$expr.mat
  colnames(expr_mat) <- paste(age, seq(ncol(expr_mat)), sep="_")
  return(expr_mat)
})
expr_mat <- do.call('cbind', expr_list)
saveRDS(expr_mat, file="Res_pt_bin_expr_mat_by_sample.rds")
pt_expr_mat <- expr_mat[union_pt_genes,]
saveRDS(pt_expr_mat, file="Res_pt_bin_expr_mat_by_sample_for_human_union_pt_genes.rds")

pt_hc <- hclust(as.dist(1-cor(pt_expr_mat)), method="ward.D2")
pdf("Plot_pt_bin_hc.pdf")
plot(pt_hc, hang=-1, cex=0.2)
dev.off()
g <- cutree(pt_hc, 3)
saveRDS(g, file="Res_pt_bin_group.rds")

res_list <- list()
for(age in names(human_se_list)){
  res <- readRDS(file.path(main_dir, age, "Res_stem_cell_to_enterocyte_turning_point_pval.rds"))
  min_pval_idx <- colnames(res)[apply(res, 1, which.min)]
  freq <- table(min_pval_idx)
  a <- setNames(rep(0, ncol(res)), colnames(res))
  a[names(freq)] <- freq
  
  g_age <- g[grep(age, names(g))]
  stage_2 <- sub(paste0(age,"_"), "", names(g_age)[which(g_age==2)])
  stage_2_start <- paste0("Pt-bin_", min(as.numeric(stage_2)))
  stage_2_end <- paste0("Pt-bin_", max(as.numeric(stage_2)))
  res_list[[age]] <- list("a"=a,
                          "stage_2_start"=stage_2_start,
                          "stage_2_end"=stage_2_end)
  seu_obj <- human_se_list[[age]]
  seu_obj$Stage_idx <- NA
  for(i in sort(unique(g_age))){
    pt_idx <- sub(paste0(age,"_"), "", names(g_age)[g_age==i])
    seu_obj$Stage_idx[which(seu_obj$Pt_bin_idx%in%pt_idx)] <- i
  }
  human_se_list[[age]] <- seu_obj
}
saveRDS(human_se_list, file="Res_human_se_list_by_sample.rds")

col_num=4
row_num=2
png("Plot_UMAP_se_stages_by_sample.png", height=1000*row_num, width=1000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,10,5))
for(age in names(human_se_list)){
  seu_obj <- human_se_list[[age]]
  plotFeature2(coor=Embeddings(seu_obj, reduction = "umap"),
               values = seu_obj$Stage_idx,
               point.order = "random",
               main=age,
               cex=4,
               cex.main=8,
               add.label = T,
               label.cex = 5)
}
dev.off()

png("Plot_UMAP_se_cell_type_by_sample.png", height=1000*row_num, width=1000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,10,5))
for(age in names(human_se_list)){
  seu_obj <- human_se_list[[age]]
  plotFeature2(coor=Embeddings(seu_obj, reduction = "umap"),
               values = seu_obj$Cell_type,
               point.order = "random",
               main=age,
               cex=4,
               cex.main=8,
               add.label = T,
               label.cex = 5)
}
dev.off()


stage_idx_vec <- setNames(rep(NA, ncol(human_se)), colnames(human_se))
for(age in names(human_se_list)){
  seu_obj <- human_se_list[[age]]
  stage_idx_vec[colnames(seu_obj)] <- seu_obj$Stage_idx
}
human_se$Stage_idx <- stage_idx_vec
DimPlot(human_se, reduction = "umap_css", group.by = "Stage_idx", split.by = "Age")


pdf("Plot_most_significant_turning_point_distribution.pdf", height=5*2, width=5*4)
par(mfrow=c(2,4))
for(age in names(res_list)){
  a <- res_list[[age]]$a
  stage_2_start <- res_list[[age]]$stage_2_start
  stage_2_end <- res_list[[age]]$stage_2_end
  plot(seq(length(a)), a, pch=16, bty="n", xlab="Pt bin", ylab="Turning point number", main=age)
  lines(smooth.spline(seq(length(a)), a, df=5))
  abline(v=which(names(a)==stage_2_start))
  abline(v=which(names(a)==stage_2_end))
}
dev.off()


