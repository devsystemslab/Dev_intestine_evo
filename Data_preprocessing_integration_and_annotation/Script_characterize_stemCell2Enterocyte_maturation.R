library(simspec)
library(Seurat)
library(presto)
library(dplyr)
source("~/Work/commonScript/Script_functions.R")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age")
human_se <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/Res_integrated_human_stemCell2Enterocyte.rds")
DimPlot(human_se, reduction = "umap_css", group.by = "Cell_type", label=T)

ct_expr <- getAveExpr(seu.obj=human_se,
                      feature.to.calc = "Cell_type",
                      colname.prefix = "Human")
saveRDS(ct_expr, file="Dat_human_fetal_cell_type_average_expr.rds")

human_se_list_by_ct <- SplitObject(human_se, split.by = "Cell_type")
main_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age"
for(cell_type in names(human_se_list_by_ct)){
  print(paste(cell_type, "start"))
  sub_dir <- cell_type
  if(!file.exists(file.path(main_dir, sub_dir))){
    dir.create(file.path(main_dir, sub_dir))
  }
  setwd(file.path(main_dir, sub_dir))
  seu_obj <- human_se_list_by_ct[[cell_type]]
  
  # preprocessing
  print("Seurat preprocessing")
  seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
  seu_obj <- ScaleData(object = seu_obj, verbose = T)
  seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
  usefulPCs <- 1:20
  seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
  seu_obj <- FindClusters(object = seu_obj, resolution = 2)
  seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
  
  expr_mat <- as.matrix(seu_obj@assays$RNA@data)
  age_vec <- as.numeric(seu_obj$Age.week)
  diff_vec <- seu_obj$Pt # stem cell to enterocyte differentiation score
  g2m_score_vec <- seu_obj$G2M.Score
  s_score_vec <- seu_obj$S.Score
  registerDoParallel(50)
  res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.vector(expr_mat[k,])
    m0 <- lm(e ~ g2m_score_vec + s_score_vec + ns(diff_vec, df=6))
    m1 <- lm(e ~ g2m_score_vec + s_score_vec + ns(diff_vec, df=6) + age_vec)
    a0 <- anova(m0)
    a1 <- anova(m1)
    p_anova <- anova(m1,m0)$Pr[2]
    p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], 
                 df1=a0["Residuals", "Df"], 
                 df2=a1["Residuals","Df"], 
                 lower.tail = F)
    coef <- coef(m1)[4]
    return(c(p_anova, p_resi ,coef))
  }
  rownames(res) <- rownames(expr_mat)
  colnames(res) <- c("p_ANOVA", "p_Resi", "Coef")
  df <- data.frame(res)
  saveRDS(df, file=paste0("Res_", cell_type, "_age_test_cc_and_differentiation_score_as_covariates.rds"))
  stopImplicitCluster()
  
  age_expr <- getAveExpr(seu.obj=seu_obj,
                         feature.to.calc = "Age.week",
                         colname.prefix = "PCW",
                         size.cutoff = 3,
                         specified.order = sort(as.numeric(unique(seu_obj$Age.week))))
  saveRDS(age_expr, file=paste0("Res_", cell_type, "_average_expr_by_age.rds"))
  
  expr_sd <- apply(age_expr, 1, sd)
  expressed_df <- df[names(expr_sd)[which(expr_sd>0)],]
  expressed_df$p_ANOVA_adj <- p.adjust(expressed_df$p_ANOVA, method="BH")
  saveRDS(expressed_df, file=paste0("Res_", cell_type, "_age_test_cc_score_as_covariates_for_expressed_genes_with_BH_adjustment.rds"))
  
  age_genes <- rownames(expressed_df)[which(expressed_df$p_ANOVA_adj<0.05)]
  cor_vec <- cor(t(age_expr[age_genes,]), as.numeric(sub("PCW_","",colnames(age_expr))), method="spearman")[,1]
  age_up_genes <- names(cor_vec)[order(cor_vec, decreasing = T)[1:200]]
  summary(cor_vec[age_up_genes])
  age_down_genes <- names(cor_vec)[order(cor_vec)[1:200]]
  summary(cor_vec[age_down_genes])
  age_genes <- list("up"=age_up_genes,
                    "down"=age_down_genes)
  saveRDS(age_genes, file=paste0("Res_", cell_type, "_age_dependent_genes.rds"))
  
  seu_obj$Age_up_score <- colSums(expr_mat[age_up_genes,])
  seu_obj$Age_down_score <- colSums(expr_mat[age_down_genes,])
  human_se_list_by_ct[[cell_type]] <- seu_obj
  saveRDS(seu_obj, file=paste0("Res_fetal_human_", cell_type, "_seu_obj.rds"))
  up_score <- lapply(sort(as.numeric(unique(seu_obj$Age.week))), function(age){
    seu_obj$Age_up_score[which(seu_obj$Age.week==age)]
  })
  names(up_score) <- paste0("PCW", sort(as.numeric(unique(seu_obj$Age.week))))
  down_score <- lapply(sort(as.numeric(unique(seu_obj$Age.week))), function(age){
    seu_obj$Age_down_score[which(seu_obj$Age.week==age)]
  })
  names(down_score) <- paste0("PCW", sort(as.numeric(unique(seu_obj$Age.week))))
  pdf(paste0("Plot_boxplot_",cell_type, "_age_score.pdf"), height = 5, width=10)
  par(mfrow=c(1,2))
  boxplot(up_score, main="Up_score", las=2)
  boxplot(down_score, main="Down_score", las=2)
  dev.off()
  
}

age_up_gene_list <- lapply(sort(unique(human_se$Cell_type)), function(cell_type){
  path <- file.path(main_dir, cell_type, paste0("Res_", cell_type, "_age_dependent_genes.rds"))
  age_genes <- readRDS(path)
  age_up_genes <- age_genes$up
})
names(age_up_gene_list) <- sort(unique(human_se$Cell_type))
n1 <- sapply(seq(length(age_up_gene_list)), function(i){
  sapply(seq(length(age_up_gene_list)), function(j){
    sum(age_up_gene_list[[i]] %in% age_up_gene_list[[j]])
  })
})
rownames(n1) <- colnames(n1) <- names(age_up_gene_list)
idx <- matrix(F, nrow=nrow(human_se), ncol=length(human_se_list_by_ct))
rownames(idx) <- rownames(human_se)
colnames(idx) <- names(human_se_list_by_ct)
for(x in names(human_se_list_by_ct)){
  genes <- age_up_gene_list[[x]]
  idx[genes,x] <- TRUE
}
union_up_genes <- rownames(idx)[which(rowSums(idx)>1)]

age_down_gene_list <- lapply(sort(unique(human_se$Cell_type)), function(cell_type){
  path <- file.path(main_dir, cell_type, paste0("Res_", cell_type, "_age_dependent_genes.rds"))
  age_genes <- readRDS(path)
  age_down_genes <- age_genes$down
})
names(age_down_gene_list) <- sort(unique(human_se$Cell_type))
idx_down <- matrix(F, nrow=nrow(human_se), ncol=length(human_se_list_by_ct))
rownames(idx_down) <- rownames(human_se)
colnames(idx_down) <- names(human_se_list_by_ct)
for(x in names(human_se_list_by_ct)){
  genes <- age_down_gene_list[[x]]
  idx_down[genes,x] <- TRUE
}
union_down_genes <- rownames(idx)[which(rowSums(idx_down)>1)]

human_se$Age_up_score <- colSums(human_se@assays$RNA@data[union_up_genes,])
human_se$Age_down_score <- colSums(human_se@assays$RNA@data[union_down_genes,])

up_score <- lapply(sort(as.numeric(unique(human_se$Age.week))), function(age){
  human_se$Age_up_score[which(human_se$Age.week==age)]
})
names(up_score) <- paste0("PCW", sort(as.numeric(unique(human_se$Age.week))))
down_score <- lapply(sort(as.numeric(unique(human_se$Age.week))), function(age){
  human_se$Age_down_score[which(human_se$Age.week==age)]
})
names(down_score) <- paste0("PCW", sort(as.numeric(unique(human_se$Age.week))))
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age")
pdf("Plot_age_up_and_down_score_for_fetal_human_stemCell2Enterocyte.pdf", height=5, width=10)
par(mfrow=c(1,2))
boxplot(up_score, main="Up_score", las=2)
boxplot(down_score, main="Down_score", las=2)
dev.off()

union_age_genes <- c(union_up_genes, union_down_genes)
age_genes <- list("up"=union_up_genes,
                  "down"=union_down_genes,
                  "combined"=union_age_genes)
saveRDS(age_genes, file="Res_union_age_genes.rds")
n1_list <- list()
for(cell_type in sort(unique(human_se$Cell_type))){
  print(cell_type)
  path <- file.path(main_dir, cell_type, paste0("Res_", cell_type, "_average_expr_by_age.rds"))
  age_expr <- readRDS(path) 
  ref_expr <- age_expr[union_age_genes,]
  que_expr <- as.matrix(human_se@assays$RNA@data[union_age_genes, which(human_se$Cell_type==cell_type)])
  cor_mat <- cor(que_expr, ref_expr, method="spearman")
  pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                   rownames(cor_mat))
  real <- paste("PCW", human_se@meta.data[names(pred), "Age.week"], sep="_")
  n1 <- sapply(colnames(age_expr), function(pred_age){
    sapply(colnames(age_expr), function(real_age){
      sum(pred==pred_age & real==real_age)
    })
  })
  n1_list[[cell_type]] <- n1
}


ref_expr_list <- lapply(sort(unique(human_se$Cell_type)), function(cell_type){
  path <- file.path(main_dir, cell_type, paste0("Res_", cell_type, "_average_expr_by_age.rds"))
  age_expr <- readRDS(path) 
  ref_expr <- age_expr[union_age_genes,]
  colnames(ref_expr) <- paste(cell_type, colnames(ref_expr), sep=":")
  return(ref_expr)
})
combined_ref_expr <- do.call('cbind', ref_expr_list)
saveRDS(combined_ref_expr, file="Dat_fetal_human_cell_type_combined_age_estimation_ref_expr.rds")
que_expr <- as.matrix(human_se@assays$RNA@data[union_age_genes, ])
cor_mat <- cor(que_expr, combined_ref_expr, method="spearman")
pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
df <- data.frame("real_cell_type"=human_se@meta.data[rownames(mat), "Cell_type"],
                 "pred_cell_type"=mat[,1],
                 "real_sample_age"=paste("PCW", human_se@meta.data[rownames(mat), "Age.week"], sep="_"),
                 "pred_sample_age"=mat[,2],
                 stringsAsFactors = F)
age_n1 <- sapply(colnames(age_expr), function(pred_age){
  sapply(colnames(age_expr), function(real_age){
    sum(df$pred_sample_age==pred_age & df$real_sample_age==real_age)
  })
})
ct_n1 <- sapply(sort(unique(human_se$Cell_type)), function(pred_ct){
  sapply(sort(unique(human_se$Cell_type)), function(real_ct){
    sum(df$pred_cell_type==pred_ct & df$real_cell_type==real_ct)
  })
})

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age/to_mouse")
human_mouse_orth <- readRDS("~/Work/Annotation/Ensembl/Human/Dat_human_mouse_one2one_symbol_only.rds")
mouse_se <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/refine_fetal_mouse_annotation/pt/Res_d17.5_fetal_mouse_stem_cell_to_enterocyte.rds")

age_genes_orth <- human_mouse_orth[which(human_mouse_orth[,"Human_symbol"]%in%union_age_genes & human_mouse_orth[,"Mouse_symbol"]%in%rownames(mouse_se)),]
dim(age_genes_orth)
ref_expr <- combined_ref_expr[age_genes_orth[,"Human_symbol"],]
saveRDS(ref_expr, file="Dat_fetal_human_cell_type_combined_age_estimation_ref_expr_only_include_human_mouse_one2one_expressed_orthologs.rds")

que_expr <- as.matrix(human_se@assays$RNA@data[age_genes_orth[,"Human_symbol"], ])
cor_mat <- cor(que_expr, ref_expr, method="spearman")
pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
df <- data.frame("real_cell_type"=human_se@meta.data[rownames(mat), "Cell_type"],
                 "pred_cell_type"=mat[,1],
                 "real_sample_age"=paste("PCW", human_se@meta.data[rownames(mat), "Age.week"], sep="_"),
                 "pred_sample_age"=mat[,2],
                 stringsAsFactors = F)
age_n1 <- sapply(colnames(age_expr), function(pred_age){
  sapply(colnames(age_expr), function(real_age){
    sum(df$pred_sample_age==pred_age & df$real_sample_age==real_age)
  })
})
ct_n1 <- sapply(sort(unique(human_se$Cell_type)), function(pred_ct){
  sapply(sort(unique(human_se$Cell_type)), function(real_ct){
    sum(df$pred_cell_type==pred_ct & df$real_cell_type==real_ct)
  })
})


que_expr <- as.matrix(mouse_se@assays$RNA@data[age_genes_orth[,"Mouse_symbol"], ])
cor_mat <- cor(que_expr, ref_expr, method="spearman")
pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
df <- data.frame("real_cell_type"=mouse_se@meta.data[rownames(mat), "Cell_type"],
                 "pred_cell_type"=mat[,1],
                 "pred_sample_age"=mat[,2],
                 stringsAsFactors = F)
mouse_se$Pred_human_age <- df$pred_sample_age
mouse_se$Pred_human_cell_type <- df$pred_cell_type
saveRDS(mouse_se, file="Res_d17.5_fetal_mouse_stemCell2Enterocyte_with_predicted_human_age.rds")

freq <- table(mouse_se$Pred_human_age)
res <- data.frame(Pred_human_age=names(freq),
                 Cell_num=as.numeric(freq))

p1 <- DimPlot(mouse_se, group.by="Pred_human_age")
p2 <-ggplot(data=res, aes(x=Pred_human_age, y=Cell_num)) +
  geom_bar(stat="identity", width=0.5, fill='steelblue')+
  theme_minimal()+
  coord_flip() +
  scale_x_discrete(limits=colnames(age_expr))
p3 <- DimPlot(mouse_se, group.by="Cell_type")
p4 <- DimPlot(mouse_se, group.by="Pred_human_cell_type")
p1+p2+p3+p4

# apply to mouse cells of different ages
mouse_epi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/Res_mouse_SI_epi_with_coarse_cell_type_annotation.rds")
p1 <- DimPlot(mouse_epi, reduction = "umap_css", group.by = "Inferred_SI_region_identity")
p2 <- DimPlot(mouse_epi, reduction = "umap_css", group.by = "Coarse_cell_type", label=T)
p1+p2
cells <- colnames(mouse_epi)[which(mouse_epi$Coarse_cell_type %in% c("Enterocyte", "Enterocyte_precursor", "Stem_cell"))]
mouse_multi_se <- subset(mouse_epi, cells=cells)
que_expr <- as.matrix(mouse_multi_se@assays$RNA@data[age_genes_orth[,"Mouse_symbol"], ])
cor_mat <- cor(que_expr, ref_expr, method="spearman")
pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
df <- data.frame("real_cell_type"=mouse_multi_se@meta.data[rownames(mat), "Coarse_cell_type"],
                 "pred_cell_type"=mat[,1],
                 "real_sample_age"=mouse_multi_se@meta.data[rownames(mat), "Age"],
                 "pred_sample_age"=mat[,2],
                 stringsAsFactors = F)
mouse_multi_se$Pred_human_age <- df$pred_sample_age
mouse_multi_se$Pred_human_cell_type <- df$pred_cell_type
saveRDS(mouse_multi_se, file="Res_fetal_mouse_stemCell2Enterocyte_with_predicted_human_age.rds")

freq <- table(mouse_multi_se$Pred_human_age)
res <- data.frame(Pred_human_age=names(freq),
                  Cell_num=as.numeric(freq))

real_age <- paste0(sort(as.numeric(sub("d", "", unique(mouse_multi_se$Age)))), "d")
pred_age <- paste0("PCW_", sort(as.numeric(sub("PCW_", "", unique(mouse_multi_se$Pred_human_age)))))
n1 <- sapply(real_age, function(x){
  sapply(pred_age, function(y){
    sum(mouse_multi_se$Age==x & mouse_multi_se$Pred_human_age==y)
  })
})
p1 <- t(t(n1)/colSums(n1))
gCols <- setNames(colorRampPalette(c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))(nrow(n1)),
                  rownames(p1))
barplot(p1, col=gCols[rownames(p1)], border = NA, ylab="Pred_human_age_prop", xlab="Mouse_sample_age")

df <- data.frame("Mouse_sample_age"=rep(colnames(p1), each=nrow(p1)),
                 "Pred_human_age"=rep(rownames(p1), ncol(p1)),
                 "Proportion"=as.vector(p1))


p0 <- DimPlot(mouse_multi_se, group.by="Age", reduction = "umap_css")
p1 <- DimPlot(mouse_multi_se, group.by="Pred_human_age", reduction = "umap_css")
p2 <-ggplot(data=res, aes(x=Pred_human_age, y=Cell_num)) +
  geom_bar(stat="identity", width=0.5, fill='steelblue')+
  theme_minimal()+
  coord_flip() +
  scale_x_discrete(limits=colnames(age_expr))
p3 <- DimPlot(mouse_multi_se, group.by="Coarse_cell_type", reduction = "umap_css")
p4 <- DimPlot(mouse_multi_se, group.by="Pred_human_cell_type", reduction = "umap_css")
p5 <- DimPlot(mouse_multi_se, group.by="Inferred_SI_region_identity", reduction = "umap_css")
p6 <- ggplot(data=df, aes(x=Mouse_sample_age, y=Proportion, fill=Pred_human_age)) +
  geom_bar(stat="identity", width = 0.8)+
  scale_fill_manual(values=gCols)+
  theme_minimal()
p0+p1+p3+p4+p5+p6+plot_layout(ncol=3)

# to tCIOs
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age/to_tCIO")
ref_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age"
combined_ref_expr <- readRDS(file.path(ref_dir,"Dat_fetal_human_cell_type_combined_age_estimation_ref_expr.rds"))
tCIO_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/signac_unified_peak_list_based/RNA_ATAC_integration"
epi_clusters <- list(
  "C7-tCIO-W12"=c(2,3,4,12),
  "C2-tCIO-W12"=c(0,1,3,4,6,8,9,10,12)
)
seurat_rna_list <- list()
p_list <- list()
pred_age_list <- list()
for(dataset in c("C2-tCIO-W12", "C7-tCIO-W12")){
  seurat_rna <- readRDS(file.path(tCIO_dir, dataset, paste0("Res_",dataset,"_seurat_rna.rds")))
  seurat_rna <- subset(seurat_rna, cells=colnames(seurat_rna)[seurat_rna$RNA_snn_res.1%in%epi_clusters[[dataset]]])
  genes <- intersect(rownames(combined_ref_expr), rownames(seurat_rna))
  ref_expr <- combined_ref_expr[genes,]
  que_expr <- as.matrix(seurat_rna@assays$RNA@data[genes,])
  cor_mat <- cor(que_expr, ref_expr, method="spearman")
  pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                   rownames(cor_mat))
  mat <- do.call('rbind', strsplit(pred, ":"))
  seurat_rna$Pred_human_age <- mat[,2]
  p_list[[dataset]] <- DimPlot(seurat_rna, group.by = "Pred_human_age")
  seurat_rna_list[[dataset]] <- seurat_rna
  pred_age_list[[dataset]] <- table(seurat_rna$Pred_human_age)
}
saveRDS(seurat_rna_list, file="Res_tCIO_epi_subset_with_pred_human_age.rds")

n1 <- matrix(NA, nrow=length(pred_age_list), ncol=length(pred_age_list[[1]]))
rownames(n1) <- names(pred_age_list)
colnames(n1) <- paste("PCW", sort(as.numeric(sub("PCW_", "", names(pred_age_list[[1]])))), sep="_")
for(x in names(pred_age_list)){
  n1[x, names(pred_age_list[[x]])] <- pred_age_list[[x]]
}
p1 <- n1/rowSums(n1)
df <- data.frame("Age"=rep(colnames(p1), each=nrow(p1)),
                 "Sample"=rep(rownames(p1), ncol(p1)),
                 "Prop"=as.vector(p1),
                 stringsAsFactors = F)
pdf("Plot_pred_human_age_prop.pdf")
ggplot(df, aes(x=Age, y=Prop, fill=Sample))+
  geom_bar(stat="identity",position=position_dodge())+
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))
dev.off()

# to H9-tHIO

# to iPSC-tHIO
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age/to_tHIO/iPSC_tHIO")
tHIO <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/suspension_HIO_derived_tHIO/RNA_human_chimp_consensus/Epithelial/remove_mesenchymal_doublets_and_foregut_cells/Res_suspention-HIO-derived_tHIO_epithelium_scRNA-seq.rds")
seurat_rna <- subset(tHIO, cells=colnames(tHIO)[which(tHIO$Epithelial_cell_type %in% c("Enterocytes", "Enterocyte precursors", "Stem cells"))])
genes <- intersect(rownames(combined_ref_expr), intersect(rownames(seurat_rna), rownames(tCIO_rna)))
ref_expr <- combined_ref_expr[genes,]
dim(ref_expr)
que_expr <- as.matrix(seurat_rna@assays$RNA@data[genes,])
dim(que_expr)
cor_mat <- cor(que_expr, ref_expr, method="spearman")
pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
seurat_rna$Pred_human_age <- mat[,2]
table(seurat_rna$Pred_human_age)
saveRDS(seurat_rna, file="Res_suspention_HIO_derived_tHIO_scRNA-seq_with_pred_human_age.rds")
DimPlot(seurat_rna, group.by = "Pred_human_age")
freq <- table(seurat_rna$Pred_human_age)
df <- data.frame("Age"=names(freq),
                 "Number"=as.numeric(freq),
                 stringsAsFactors = F)
pdf("Plot_pred_human_age_prop.pdf")
ggplot(df, aes(x=Age, y=Number))+
  geom_bar(stat="identity",position=position_dodge())+
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))
dev.off()

#res_DE <- expressed_df
#res_DE <- res_DE %>%
#  mutate(DE = p_ANOVA_adj<0.01 & abs(Coef)>3) %>%
#  mutate(DEG = ifelse(DE, rownames(res_DE), NA))
#sum(res_DE$DE)
#
#res_DE[which(!is.na(res_DE$DEG)),] -> mat
#g_up <- rownames(mat)[which(mat$Age_change>0)]
#g_down <- rownames(mat)[which(mat$Age_change<0)]
#plotFeature(seu.obj=human_se, 
#            genes.to.plot = g_down,
#            plot.name = "Plot_UMAP_human_se_enterocyte_age_down_genes.png")
#
##library(ggrepel)
#pdf("Plot_scater_plot_age_change_and_enterocyte_enrichment.pdf")
#ggplot(res_DE, aes(x = Enterocyte_enrichment, y = Age_change, col=Age_dependent, label=DEG)) +
#  geom_hex(bins=30) +
#  scale_fill_continuous(type = "viridis") +
#  geom_text_repel() +
#  geom_hline(yintercept=c(1.5, -1.5), col="#303030", linetype="dotted") +
#  scale_color_manual(values=c("#909090", "red")) +
#  theme_minimal()
#dev.off()
#
#
## Basic scatterplot
#ggplot(data, aes(x=x, y=y) ) +
#  geom_point()
#
## Hexbin chart with default option
#ggplot(data, aes(x=x, y=y) ) +
#  geom_hex() +
#  theme_bw()
#
## Bin size control + color palette
#ggplot(data, aes(x=x, y=y) ) +
#  geom_hex(bins = 30) +
#  scale_fill_continuous(type = "viridis") +
#  theme_bw()
# 
#
#
## Area + contour
#ggplot(data, aes(x=x, y=y) ) +
#  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
  


