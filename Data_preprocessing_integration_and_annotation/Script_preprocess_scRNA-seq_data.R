setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Human_tHIO_hg38/scRNA-seq")
library(Seurat)
library(simspec)
source("~/Work/commonScript/Script_functions.R")

# step 1. read scRNA-seq data of tHIO samples
# read paths to tHIO scRNA-seq count matrix
library(readxl)
readout <- "scRNA-seq_respective_genome"
species <- "Human"
tissue <- "transplanted"

sample_info <- read_excel("~/Work/Endoderm/intestine_evolution/sample/HIO_CIO_metadata.xlsx", sheet=readout)
additional_sample_to_include <- c(
  "H9-NRG1-W11"
)
selected_sample_info <- sample_info[which(sample_info$Species==species & sample_info$Tissue %in% tissue | sample_info$`Sample name`%in%additional_sample_to_include),]
dim(selected_sample_info)


# read data
seu_obj_list <- list()
for(i in seq(nrow(selected_sample_info))){
  aa <- selected_sample_info$`Sample name`[i]
  print(paste("Reading", aa))
  vec <- as.vector(as.matrix(selected_sample_info[i,3:11]))
  info.mat <- cbind(names(selected_sample_info)[3:11], vec)
  path <- paste0(selected_sample_info[i,12],"/outs/filtered_feature_bc_matrix/")
  seu_obj <- prepareSeuratObject2(rawData = Read10X(path), 
                                  namePrefix = paste0("HR-",i), 
                                  additional.sample.info.mat=info.mat, seu.version=3, 
                                  mito.cutoff=Inf, min.cell.cutoff=1, 
                                  min.gene.cutoff=500, gene.num.low=500, gene.num.high=Inf)
  seu_obj_list[[i]] <- seu_obj
}

combined <- merge(x=seu_obj_list[[1]], y=seu_obj_list[-1])
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size=0 )
dim(combined)
combined_sub <- subset(combined, subset = percent.mito < 0.25 & nFeature_RNA > 1000 & nCount_RNA < 50000)
VlnPlot(combined_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size=0 )
table(combined_sub$orig.ident)
bk <- combined
combined <- combined_sub
dim(combined)
rm(combined_sub)

# remove confounding genes
dir <- "~/Work/Annotation/confound_genes/" 
all.files <- list.files(dir, pattern=".txt")
confound.genes <- lapply(seq(length(all.files)), function(i){
  file <- paste0(dir, all.files[i])
  g <- readLines(file)
  return(g)
})
names(confound.genes) <- c("Cell_cycle", "Experiment_induced", "HB", "MT", "Ribosome", "Sex_chr") 
all.confound.genes <- unique(unlist(confound.genes))

genes.to.remove <- unique(c(confound.genes[["MT"]], confound.genes[["Ribosome"]], confound.genes[["Sex_chr"]]))
counts <- combined@assays$RNA@counts
counts <- counts[setdiff(rownames(counts), genes.to.remove),]
meta <- combined@meta.data[,4:13]
combined <- CreateSeuratObject(counts=counts, meta=meta)

combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 1e4)
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)
saveRDS(combined, file="Res_tHIO_hg38_no_batch_correction.rds")
DimPlot(combined, reduction = "umap", group.by = "Sample.name")

# CSS integration
combined <- cluster_sim_spectrum(combined, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]

DimPlot(combined, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1") + DimPlot(combined, reduction = "umap_css", group.by = "Sample.name")
saveRDS(combined, file="Res_tHIO_hg38_with_CSS_integration.rds")
# MNN integration
library(SeuratWrappers)
sf.list <- SplitObject(combined,split.by="orig.ident")
for (i in 1:length(sf.list)) {
  sf.list[[i]] <- NormalizeData(sf.list[[i]],assay="RNA",normalization.method="LogNormalize")
  sf.list[[i]] <- FindVariableFeatures(sf.list[[i]], selection.method="vst", nfeatures=2000)
}

sff <- RunFastMNN(object.list = sf.list)
combined[['mnn']] <- CreateDimReducObject(Embeddings(sff, "mnn")[colnames(combined),], key="MNN_")
combined <- RunUMAP(combined, reduction="mnn", dims = 1:20, reduction.name = "umap_mnn", reduction.key = "UMAPMNN_")
combined <- FindNeighbors(object = combined, reduction = "mnn", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_MNN_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(combined, file="Res_tHIO_hg38_with_CSS_and_MNN_integration.rds")
p1 <- DimPlot(combined, reduction = "umap_mnn", group.by = "Sample.name")+NoLegend()
p2 <- DimPlot(combined, reduction = "umap_css", group.by = "Sample.name")
p3 <- DimPlot(combined, reduction = "umap_mnn", group.by = "RNA_MNN_snn_res.1", label=T)+NoLegend()
p4 <- DimPlot(combined, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", label=T)+NoLegend()
p1+p2+p3+p4+plot_layout(ncol=2)
# plot cell class marker expression
g1 <- c("EPCAM","COL1A2","CDH5","PTPRC","ASCL1")
plotFeature.batch(seu.obj=combined, dr="umap_css", genes.to.plot=g1, 
                  col.num = 5, plot.name = "Plot_UMAP_CSS_major_cell_type_marker_expr.png",
                  nCols = blue.cols, cex=2)
plotFeature.batch(seu.obj=combined, dr="umap_mnn", genes.to.plot=g1, 
                  col.num = 5, plot.name = "Plot_UMAP_MNN_major_cell_type_marker_expr.png",
                  nCols = blue.cols, cex=2)


g1 <- c("CDX2", "PDX1", "ONECUT2", "LGALS4", "SOX2", "SATB2", "CLDN18", "GATA4", "GATA6", "HOXA10", "OSR2")
plotFeature.batch(seu.obj=combined, dr="umap_css", genes.to.plot=g1, 
                  col.num = 6, plot.name = "Plot_UMAP_CSS_regional_marker_expr.png",
                  nCols = blue.cols, cex=2)
plotFeature.batch(seu.obj=combined, dr="umap_mnn", genes.to.plot=g1, 
                  col.num = 6, plot.name = "Plot_UMAP_MNN_regional_marker_expr.png",
                  nCols = blue.cols, cex=2)

# classify clusters into cell classes
known.markers <- read.table("~/Work/Annotation/cellTypeMarker/Intestine/Table_major_cell_type_markers_from_literature_search.txt", head=T, sep="\t",stringsAsFactors=F)  
cell.type.order <- c("Neuronal", "Epithelial", "Mesenchymal", "Endothelial", "Immune", "Erythroid", "Hepatocyte")
pan.cm <- unique(unlist(lapply(cell.type.order, function(ct){
  g1 <- known.markers$Gene[which(known.markers$Used_pan_cell_type_markers==ct)]
  intersect(g1, rownames(combined))
}))) 

# cell class annotation
cell_class <- list(
  "Neural" = 28,
  "Epithelial" = c(4,6,7,8,13,14,15,18,20,21,22),
  "Epithelial_and_mesenchymal_doublet"=c(9,24),
  "Mesenchymal" = c(0,1,2,3,5,10,11,16,17,19,23,25,27),
  "Endothelial_and_mesenchymal" = 26,
  "Immune_and_others" = 12
)
cl_order <- unlist(cell_class)
seu_obj_expr <- getAveExpr(seu.obj=combined, feature.to.calc = "RNA_MNN_snn_res.1", specified.order = cl_order, genes=pan.cm)
seu_obj_prop <- getExpressedProp(seu.obj = combined, feature.to.calc = "RNA_MNN_snn_res.1", specified.order = cl_order, genes=pan.cm)
getDotPlot(cl.expr.mat=seu_obj_expr, cl.prop.mat=seu_obj_prop, gene.reorder=FALSE, specified.order=NULL, cl.reorder.by.hc=FALSE, 
           genes=pan.cm, colors=c("#d9d9d9", "#252525"), point.size.factor=4.5, 
           plot.name="Plot_dot_plot_cell_class_marker_RNA.pdf", 
           plot.height=20, plot.width=30, plot.margin=c(8,12,5,8), 
           plot.cex=2, max.diag=FALSE, col.space.factor=0.5, row.space.factor=0.7)


cl_vec <- combined$RNA_MNN_snn_res.1
class_vec <- rep(NA, length(cl_vec))
for(x in names(cell_class)){
  cl_x <- cell_class[[x]]
  class_vec[cl_vec %in% cl_x] <- x
}
combined$Cell_class <- class_vec
DimPlot(combined, group.by = "Cell_class", reduction = "umap_mnn")
saveRDS(combined, file="Res_tHIO_hg38_scRNA-seq.rds")

# epithelial cell type or cell state markers
## extract epithelial cells and perform subclustering
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Human_tHIO_hg38/scRNA-seq/epi")
epi <- subset(combined, Cell_class=="Epithelial")
epi <- FindVariableFeatures(object = epi, selection.method = "vst", nfeatures = 3000)
epi <- ScaleData(object = epi, verbose = T)
epi <- RunPCA(object = epi, features = VariableFeatures(epi), verbose = F, npcs = 50)
usefulPCs <- 1:20
epi <- FindNeighbors(object = epi, dims = usefulPCs)
epi <- FindClusters(object = epi, resolution = 2)
epi <- RunUMAP(object = epi, dims = usefulPCs)
DimPlot(epi, label=T)+DimPlot(epi, group.by = "Sample.name")
saveRDS(epi, file="Res_tHIO_epithelium_scRNA-seq.rds")
epi$RNA_PCA_snn_res.2 <- epi$RNA_snn_res.1
# CSS integration
epi <- cluster_sim_spectrum(epi, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
epi <- RunUMAP(epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
epi <- FindNeighbors(object = epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
epi[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- epi[[paste0("RNA_snn_res.", 0.2*5)]]

# MNN integration
seu_list <- SplitObject(epi,split.by="orig.ident")
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- NormalizeData(seu_list[[i]],assay="RNA",normalization.method="LogNormalize")
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method="vst", nfeatures=2000)
}

seu_obj <- RunFastMNN(object.list = seu_list)
epi[['mnn']] <- CreateDimReducObject(Embeddings(seu_obj, "mnn")[colnames(epi),], key="MNN_")
epi <- RunUMAP(epi, reduction="mnn", dims = 1:20, reduction.name = "umap_mnn", reduction.key = "UMAPMNN_")
epi <- FindNeighbors(object = epi, reduction = "mnn", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
epi[[paste0("RNA_MNN_snn_res.", 0.2*5)]] <- epi[[paste0("RNA_snn_res.", 0.2*5)]]
DimPlot(epi, reduction = "umap_css", group.by = "Sample.name")+DimPlot(epi, reduction = "umap_mnn", group.by = "Sample.name")
saveRDS(epi, file="Res_tHIO_epi_hg38_with_CSS_and_MNN_integration.rds")


g1 <- c("EPCAM","COL1A2","CDH5","PTPRC","ASCL1")
plotFeature.batch(seu.obj=epi, dr="umap", genes.to.plot=g1, 
                  col.num = 5, plot.name = "Plot_UMAP_major_cell_type_marker_expr.png",
                  nCols = blue.cols, cex=2)

g1 <- c("CDX2", "PDX1", "ONECUT2", "LGALS4", "SOX2", "SATB2", "CLDN18", "GATA4", "GATA6","HOXA10", "OSR2","PTF1A", "PROX1", "NKX6-1", "ONECUT1")
plotFeature.batch(seu.obj=epi, dr="umap_css", genes.to.plot=g1, 
                  col.num = 6, plot.name = "Plot_UMAP_CSS_epithelium_regional_marker_expr.png",
                  nCols = blue.cols, cex=2)



# Epithelial cell type annotation
p1 <- DimPlot(epi, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", label=T)+NoLegend()
p2 <- DimPlot(epi, reduction = "umap_css", group.by = "Sample.name")
p1+p2
g1 <- c("CDX2","PDX1","ONECUT2", "SATB2","SOX2", "FOXJ1", "TP63", "CLDN18", "HOXA10", "OSR2","LGALS4", "SHH","LGR5", "ASCL2", "OLFM4", "MUC2", "SPINK4", "SPDEF", "MKI67","CDK1",
        "CHGA", "ISL1", "BEST4", "SPIB", "GP2", "CA7", "MUC1","DPP4", "APOA4", "FABP2", "SI","DEFA5", "DEFA6","LYZ", "RGS13")
plotFeature.batch(seu.obj = epi, dr="umap_css", genes.to.plot = g1, nCols = blue.cols, cex=3,
                  plot.name = "Plot_UMAP_CSS_tHIO_epithelium_epi_cell_type_marker_expression.png", 
                  col.num = 8)

plotFeature2(coor = Embeddings(epi, reduction = "umap"),
             values = epi$nCount_RNA,
             point.order = 'sorted',
             cex=1)

# remove non-intestinal cell clusters
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Human_tHIO_hg38/scRNA-seq/epi/intestinal")
epi <- subset(epi, cells=colnames(epi)[which(!epi$RNA_CSS_snn_res.1 %in% c(9,12))])
epi <- FindVariableFeatures(object = epi, selection.method = "vst", nfeatures = 3000)
epi <- ScaleData(object = epi, verbose = T)
epi <- RunPCA(object = epi, features = VariableFeatures(epi), verbose = F, npcs = 50)
usefulPCs <- 1:20
epi <- FindNeighbors(object = epi, dims = usefulPCs)
epi <- FindClusters(object = epi, resolution = 2)
epi <- RunUMAP(object = epi, dims = usefulPCs)
DimPlot(epi, label=T)+DimPlot(epi, group.by = "Sample.name")
epi$RNA_PCA_snn_res.2 <- epi$RNA_snn_res.2
# CSS integration
epi <- cluster_sim_spectrum(epi, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
epi <- RunUMAP(epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
epi <- FindNeighbors(object = epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
epi[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- epi[[paste0("RNA_snn_res.", 0.2*5)]]

# MNN integration
seu_list <- SplitObject(epi,split.by="orig.ident")
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method="vst", nfeatures=2000)
}

seu_obj <- RunFastMNN(object.list = seu_list)
epi[['mnn']] <- CreateDimReducObject(Embeddings(seu_obj, "mnn")[colnames(epi),], key="MNN_")
epi <- RunUMAP(epi, reduction="mnn", dims = 1:20, reduction.name = "umap_mnn", reduction.key = "UMAPMNN_")
epi <- FindNeighbors(object = epi, reduction = "mnn", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
epi[[paste0("RNA_MNN_snn_res.", 0.2*5)]] <- epi[[paste0("RNA_snn_res.", 0.2*5)]]
DimPlot(epi, reduction = "umap_css", group.by = "Sample.name")+DimPlot(epi, reduction = "umap_mnn", group.by = "Sample.name")
saveRDS(epi, file="Res_tHIO_intestinal_epi_hg38_with_CSS_and_MNN_integration.rds")

# cell type annotation
cell_anno <- list(
  "EEC_ISL1+"=10,
  "EEC_ISL1-"=14,
  "Goblet_and_Paneth"=5,
  "BEST4+_enterocyte"=15
)

cl_vec <- epi$RNA_CSS_snn_res.1
ct_vec <- rep("Stem-cell-to-enterocyte", ncol(epi))
for(x in names(cell_anno)){
  cl_x <- cell_anno[[x]]
  ct_vec[which(cl_vec %in% cl_x)] <- x
}
epi$Cell_type <- ct_vec
DimPlot(epi, reduction = "umap_css", group.by = "Cell_type", label=T)
ref_expr <- colSums(epi@assays$RNA@data[c("DEFA5","DEFA6"), which(epi$Cell_type=="Goblet_and_Paneth")])
que_expr <- as.matrix(epi@assays$RNA@data[,which(epi$Cell_type=="Goblet_and_Paneth")])
cor_vec <- cor(t(que_expr), ref_expr)[,1]
genes <- names(cor_vec)[order(cor_vec, decreasing = T)[1:10]]
score <- colSums(epi@assays$RNA@data[genes,which(epi$Cell_type=="Goblet_and_Paneth")])
plot(seq(length(score)), sort(score))
cutoff=quantile(score, 0.975)
abline(h=cutoff)
paneth_cells <- names(score)[which(score>cutoff)]
goblet_cells <- setdiff(colnames(epi)[which(epi$Cell_type=="Goblet_and_Paneth")],
                        paneth_cells)
epi@meta.data[paneth_cells, "Cell_type"] <- "Paneth_cell"
epi@meta.data[goblet_cells, "Cell_type"] <- "Goblet_cell"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
s2e_features <- list(
  "Stem_cell"=c("LGR5","OLFM4","ASCL2"),
  "Early_enterocyte"=c("GSTA1", "GSTA2", "AADAC", "REEP6"),
  "Enterocyte"=c("FABP2","APOA4")
)
s2e <- subset(epi, cells=colnames(epi)[which(epi$Cell_type=="Stem-cell-to-enterocyte")])
DimPlot(s2e, reduction = "umap_css")
score_mat <- sapply(names(s2e_features), function(x){
  scale(colSums(s2e@assays$RNA@data[s2e_features[[x]],]))
})
rownames(score_mat) <- colnames(s2e)
epi <- FindNeighbors(object = epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), force.recalc = T) 

nn_mat <- epi@graphs$RNA_nn[colnames(s2e), colnames(s2e)] 
smoothed_score_mat <- t(sapply(seq(nrow(nn_mat)), function(i){
  vec <- as.vector(as.matrix(nn_mat[i,]))
  nn_cells <- colnames(nn_mat)[which(vec==1)]
  nn_score <- colMeans(rbind(score_mat[nn_cells,], score_mat[i,]))
}))
rownames(smoothed_score_mat) <- rownames(nn_mat)
id <- colnames(smoothed_score_mat)[apply(smoothed_score_mat, 1, which.max)]
s2e$s2E_smoothed_ident_scaled_expr_based <- id
DimPlot(s2e, reduction = "umap_css", group.by = "s2E_smoothed_ident_scaled_expr_based", shuffle = T, label = T)+NoLegend()
epi@meta.data[colnames(s2e), "Cell_type"] <- s2e$s2E_smoothed_ident_scaled_expr_based
DimPlot(epi, reduction = "umap_css", group.by = "Cell_type", label=T)+NoLegend()
saveRDS(epi, file="Res_tHIO_intestinal_epi_hg38_with_CSS_and_MNN_integration.rds")

png("Plot_UMAP_CSS_tHIO_intestinal_epithelial.png", height=2000, width=2000*3)
par(mar=c(5,10,15,10), mfrow=c(1,3))
plotFeature2(Embeddings(epi, reduction = "umap_css"),
             values = epi$Cell_type,
             point.order = "random",
             add.label = T,
             label.cex = 7,
             cex=3,
             lwd=0.8,
             main="Cell type",
             cex.main=10)
plotFeature2(Embeddings(epi, reduction = "umap_css"),
             values = epi$RNA_CSS_snn_res.1,
             point.order = "random",
             add.label = T,
             label.cex = 7,
             cex=3,
             lwd=0.8,
             main="Cell cluster",
             cex.main=10)
plotFeature2(Embeddings(epi, reduction = "umap_css"),
             values = epi$Sample.name,
             point.order = "random",
             add.legend = T,
             legend.cex = 8,
             cex=3,
             lwd=0.8,
             main="Sample.name",
             cex.main=10)

dev.off()

# match tHIO ages to fetal human ages
ref_dir <-"/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/se_across_age"
combined_ref_expr <- readRDS(file.path(ref_dir,"Dat_fetal_human_cell_type_combined_age_estimation_ref_expr.rds"))
genes <- intersect(rownames(combined_ref_expr), rownames(epi))
length(genes)
s2e <- subset(epi, cells=colnames(epi)[which(epi$Cell_type %in% c("Stem_cell", "Early_enterocyte", "Enterocyte"))])
DimPlot(s2e, reduction = "umap_css")
que_expr <- as.matrix(s2e@assays$RNA@data[genes, ])
ref_expr <- combined_ref_expr[genes,]
cor_mat <- cor(que_expr, ref_expr, method="spearman")

pred <- setNames(rep(NA, ncol(s2e)), colnames(s2e))
s2e_list <- SplitObject(s2e, split.by = "Sample.name")
for(x in names(s2e_list)){
  seu_obj <- s2e_list[[x]]
  seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
  seu_obj <- ScaleData(object = seu_obj, verbose = T)
  seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
  usefulPCs <- 1:20
  seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
  nn_mat <- seu_obj@graphs$RNA_nn
  smoothed_cor_mat <- t(sapply(seq(nrow(nn_mat)), function(i){
    vec <- as.vector(as.matrix(nn_mat[i,]))
    nn_cells <- colnames(nn_mat)[which(vec==1)]
    nn_cor <- colMeans(rbind(cor_mat[nn_cells,], cor_mat[rownames(nn_mat)[i],]))
  }))
  rownames(smoothed_cor_mat) <- rownames(nn_mat)
  id <- colnames(smoothed_cor_mat)[apply(smoothed_cor_mat, 1, which.max)]
  pred[rownames(nn_mat)] <- id
}

#pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
#                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
df <- data.frame("real_cell_type"=s2e@meta.data[rownames(mat), "Cell_type"],
                 "pred_cell_type"=mat[,1],
                 "real_sample_age"=s2e@meta.data[rownames(mat), "Age"],
                 "pred_sample_age"=mat[,2],
                 stringsAsFactors = F)
s2e$smoothed_Pred_human_age <- df$pred_sample_age
s2e$smoothed_Pred_human_cell_type <- df$pred_cell_type

pred <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)],
                 rownames(cor_mat))
mat <- do.call('rbind', strsplit(pred, ":"))
df <- data.frame("real_cell_type"=s2e@meta.data[rownames(mat), "Cell_type"],
                 "pred_cell_type"=mat[,1],
                 "real_sample_age"=s2e@meta.data[rownames(mat), "Age"],
                 "pred_sample_age"=mat[,2],
                 stringsAsFactors = F)
s2e$Pred_human_age <- df$pred_sample_age
s2e$Pred_human_cell_type <- df$pred_cell_type

#p1 <- DimPlot(s2e, reduction = "umap_css", group.by = "Age")
#p2 <- DimPlot(s2e, reduction = "umap_css", group.by = "Pred_human_age")
#p1+p2
saveRDS(s2e, file="Res_tHIO_stemCell2Enterocyte_with_predicted_fetal_human_age.rds")

real_age <- sort(unique(s2e$Sample.name))
pred_age <- paste0("PCW_", sort(as.numeric(sub("PCW_", "", unique(s2e$smoothed_Pred_human_age)))))
n1 <- sapply(real_age, function(x){
  sapply(pred_age, function(y){
    sum(s2e$Sample.name==x & s2e$smoothed_Pred_human_age==y)
  })
})
p1 <- t(t(n1)/colSums(n1))

df <- data.frame("tHIO_sample_name"=rep(colnames(p1), each=nrow(p1)),
                 "Pred_human_age"=rep(rownames(p1), ncol(p1)),
                 "Proportion"=as.vector(p1),
                 stringsAsFactors = F)

pdf("Plot_barplot_tHIO_pred_fetal_human_age.pdf", height=5, width=7)
ggplot(df, aes(x=Pred_human_age, y=Proportion, fill=tHIO_sample_name))+
  geom_bar(stat="identity",position=position_dodge())+
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))
dev.off()

DimPlot(s2e, reduction="umap_css", group.by="smoothed_Pred_human_age", split.by = "Sample.name", pt.size=2)
plotFeature2(coor=Embeddings(s2e, reduction = "umap_css"),
             values = s2e$smoothed_Pred_human_age,
             emphasize = which(s2e$Sample.name=="H9-tHIO-EGF-mesentery-W8"),
             add.legend = T)

# combine stem-cell-to-enterocyte trajectory of tHIOs, fetal human and fetal mouse primary
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Human_tHIO_hg38/scRNA-seq/epi/intestinal/with_fetal_human_and_mouse")
tHIO_s2e <- s2e
mouse_se <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/refine_fetal_mouse_annotation/pt/Res_d17.5_fetal_mouse_stem_cell_to_enterocyte.rds")
DimPlot(mouse_se, group.by = "Cell_type")
# read fetal human duodenum stem cell to enterocyte data 
human_se <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/Res_integrated_human_stemCell2Enterocyte.rds")
DimPlot(human_se, group.by = "Cell_type", reduction = "umap_css")
human_se <- subset(human_se, cells=colnames(human_se)[!human_se$Age%in%c("59d","72d")])
human <- merge(human_se, tHIO_s2e)

human_count <- as.matrix(human@assays$RNA@counts)
human_data <- as.matrix(human@assays$RNA@data)
human_meta <- human@meta.data
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
combined_obj <- FindVariableFeatures(combined_obj)
# load human_mid_se and mouse_se integrated data
folder = "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css"
file="Res_PCW11-14_fetal_human_and_d17.5_fetal_mouse_stemCell2Enterocyte_integrated_with_CSS.rds"
aa <- readRDS(file.path(folder, file))
DimPlot(aa, reduction = "umap_css")
VariableFeatures(combined_obj) <- VariableFeatures(aa)
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
cor_list <- lapply(unique(combined_obj$Sample.name), function(x){
  ref_obj <- subset(combined_obj, cells=colnames(combined_obj)[combined_obj$Sample.name==x]) %>% 
    FindNeighbors(dims = usefulPCs) %>% 
    FindClusters(resolution = 0.6)
  ref_expr <- getAveExpr(seu.obj=ref_obj, feature.to.calc = "RNA_snn_res.1", colname.prefix = x, genes = hvg)
  cor_mat <- t(scale(t(cor(que_expr, ref_expr, method="spearman"))))
  return(cor_mat)
})
css_mat <- do.call('cbind', cor_list)
combined_obj[['css']] <- CreateDimReducObject(embeddings = css_mat, key = "css_", assay = "RNA")
combined_obj <- RunUMAP(combined_obj, reduction = "css", dims = 1:ncol(combined_obj@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined_obj <- FindNeighbors(object = combined_obj, reduction = "css", dims = 1:ncol(combined_obj@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
combined_obj[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined_obj[[paste0("RNA_snn_res.", 0.2*5)]]
DimPlot(combined_obj, reduction = "umap_css", group.by = "Cell_type")

### reconstruct Pt in each sample separately
seu_obj_list <- SplitObject(combined_obj, split.by = "Sample.name")
expr_list <- list()
for(x in names(seu_obj_list)){
  seu_obj <- seu_obj_list[[x]]
  input <- seu_obj@reductions$css@cell.embeddings
  dm <- DiffusionMap(input, k=20)
  vec <- rank(dm$DC1)
  aa <- median(vec[which(seu_obj$Cell_type=="Stem_cell")]) > median(vec[which(seu_obj$Cell_type=="Enterocyte")])
  if(aa){
    seu_obj$Pt_by_sample <- rank(-dm$DC1)/ncol(seu_obj)
  }else{
    seu_obj$Pt_by_sample <- rank(dm$DC1)/ncol(seu_obj)
  }
  seu_obj_list[[x]] <- seu_obj
  
  num_breaks <- 20
  pt <- seu_obj$Pt_by_sample
  expr_mat <- as.matrix(seu_obj@assays$RNA@data)
  pt_bin_expr <- sapply(1:num_breaks, function(i){
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
  expr_list[[x]] <- pt_bin_expr
}
combined_expr_mat <- do.call('cbind', expr_list)
rownames(combined_expr_mat) <- rownames(combined_obj)
saveRDS(expr_list, file="Res_pt_by_sample_bin_average_expr.rds")
saveRDS(seu_obj_list, file="Res_seurat_object_by_sample.rds")
saveRDS(combined_obj, file="Res_fetal_human_without_two_most_early_sample_E17.7_mouse_and_tHIO_s2e.rds")
folder = "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/compare_fetal_human_and_mouse_duodenum/stem_cell_to_enterocyte/human_mouse_css/DEG_control_ct"
file="Res_human_vs_mouse_stemCell2Enterocyte_pt_dependent_gene_module_expr.rds"
deg_module_expr <- readRDS(file.path(folder, file)) 

group.cols <- setNames(c("#feebe2","#fcc5c0","#fa9fb5","#f768a1","#c51b8a","#7a0177",
                         "#c7e9b4","#7fcdbb","#2c7fb8","#253494",
                         "#F7DC6F"),
                       names(expr_list))
group.vec <- rep(names(expr_list), each=num_breaks)
time.vec <- rep(seq(num_breaks), length(expr_list))
g.list <- deg_module_expr$g.list
g.size = deg_module_expr$g.size
genes <- intersect(unlist(g.list), rownames(combined_expr_mat))
length(genes)
scale.expr = t(scale(t(combined_expr_mat[genes,])))

mean.list <- list()
for(group in unique(group.vec)){
  print(paste(group, "start"))
  e <- scale.expr[,which(group.vec==group)]
  t <- time.vec[which(group.vec==group)]
  mean.mat <- sapply(seq(length(g.list)), function(i){
    apply(e[intersect(g.list[[i]],rownames(e)),], 2, mean)
  })
  mean.list[[group]] <- mean.mat
}
mean.combined.mat <- do.call("rbind", mean.list)

plot.name <- "Plot_fetal_human_mouse_DEG_cluster_average_expr_including_all_fetal_human_and_tHIO.pdf"
col.num=4
row.num=ceiling(length(g.list)/col.num)
pdf(plot.name, height=5*row.num, width=5*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,5,5))
for(i in seq(length(g.list))){
  mean.vec <- mean.combined.mat[,i]
  ymax <- max(mean.vec)
  ymin <- min(mean.vec)
  plot(c(min(time.vec), max(time.vec)), c(ymax, ymin), type="n", xlab="Pt bin",ylab="Relative expr", main=paste("Cluster",i,"(",g.size[i],")"), bty="n", cex.axis=2, cex.lab=2, cex.main=3)
  for(group in unique(group.vec)){
    g.idx <- which(group.vec==group)
    g.time <- time.vec[g.idx]
    g.mean <- mean.vec[g.idx]
    lines(g.time, smooth.spline(g.time, g.mean, df=6)$y, col=group.cols[group], lwd=3)
  }
  if(i==8){
    legend("topright", legend=names(group.cols), text.col=group.cols, bty="n", cex=1.3)
  }
  
}
dev.off()
