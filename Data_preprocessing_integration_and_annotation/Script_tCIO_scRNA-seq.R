setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO")
library(Seurat)
library(simspec)
library(harmony)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(presto)
library(dplyr)
library(destiny)
library(doParallel)
library(splines)
source("~/Work/commonScript/Script_functions.R")

# step 1. read scRNA-seq data of tCIO samples
# read paths to tCIO scRNA-seq count matrix
library(readxl)
readout <- "scRNA-seq_respective_genome"
species <- "Chimp"
sample <- c("C7-tCIO-W12","C2-tCIO-W12")

sample_info <- read_excel("~/Work/Endoderm/intestine_evolution/sample/HIO_CIO_metadata.xlsx", sheet=readout)
selected_sample_info <- sample_info[which(sample_info$`Sample name` %in% sample),]
selected_sample_info

# read data
seu_obj_list <- list()
for(i in seq(nrow(selected_sample_info))){
  aa <- selected_sample_info$`Sample name`[i]
  print(paste("Reading", aa))
  vec <- as.vector(as.matrix(selected_sample_info[i,3:11]))
  info.mat <- cbind(names(selected_sample_info)[3:11], vec)
  path <- paste0(selected_sample_info[i,12],"/outs/filtered_feature_bc_matrix/")
  seu_obj <- prepareSeuratObject2(rawData = Read10X(path), 
                                  namePrefix = paste0("CR-",i), 
                                  species=species, 
                                  mito.gene.list = "~/Work/Annotation/panTro6/List_chimp_mitochondrial_genes.txt",
                                  additional.sample.info.mat=info.mat, seu.version=3, 
                                  mito.cutoff=Inf, min.cell.cutoff=1, 
                                  min.gene.cutoff=500, gene.num.low=500, gene.num.high=Inf)
  seu_obj_list[[i]] <- seu_obj
}

combined <- merge(x=seu_obj_list[[1]], y=seu_obj_list[-1])
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size=0 )
dim(combined)
combined_sub <- subset(combined, subset = percent.mito < 0.2 & nFeature_RNA > 1000 & nCount_RNA < 50000)
VlnPlot(combined_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size=0 )
bk <- combined
combined <- combined_sub
dim(combined)
rm(combined_sub)

# chimp does not have annotation of ribosomal genes, so we remove the chimp orthologs of human ribosomal genes
human_ribo_genes <- readLines("~/Work/Annotation/confound_genes/List_ribosome_protein_genes.txt")
human_chimp_orth <- read.table("~/Work/Annotation/Ensembl/Human/v107/Ensembl_v107_biomart_human_chimp_ortholog.txt",
                               sep="\t",
                               stringsAsFactors = F,
                               head=T)
chimp_ribo_orth <- setdiff(human_chimp_orth$Chimpanzee.gene.name[which(human_chimp_orth$Gene.name%in%human_ribo_genes)], "")
writeLines(chimp_ribo_orth, con = "~/Work/Annotation/panTro6/confound_genes/List_EnsV107_chimp_ribosomal_orthologs_of_human_ones.txt")

dir <- "~/Work/Annotation/panTro6/confound_genes/" 
all.files <- list.files(dir, pattern=".txt")
confound.genes <- lapply(seq(length(all.files)), function(i){
  file <- paste0(dir, all.files[i])
  g <- readLines(file)
  return(g)
})
names(confound.genes) <- c("MT", "Sex_chr", "Ribosomal") 

genes.to.remove <- unique(c(confound.genes[["MT"]], confound.genes[["Sex_chr"]], confound.genes[["Ribosomal"]]))
counts <- combined@assays$RNA@counts
counts <- counts[setdiff(rownames(counts), genes.to.remove),]
meta <- combined@meta.data[,4:13]
seu_obj <- CreateSeuratObject(counts=counts, meta=meta)
seu_obj <- NormalizeData(object = seu_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
saveRDS(seu_obj, file="Res_C2_and_C7_tCIO_panTro6_scRNA-seq_merged_without_integration.rds")
DimPlot(seu_obj, group.by = "Sample.name")
g1 <- c("CDX2", "LGALS4", "PDX1", "ONECUT2", "SATB2","SOX2","CLDN18",
        "GATA4","GATA6","OSR2","HOXA10",
        "CDH1", "EPCAM", "KRT8", "KRT18","FXYD3","PERP",
        "COL1A2", "VIM","DCN","MFAP4","COL3A1","COL6A2",
        "PTPRC", "CORO1A","CD53","CD37","LCP1","LAPTM5",
        "PHOX2B", "PRPH", "ASCL1","HAND2","TBX3",
        "CDH5", "PECAM1","FLT1","CLDN5","EGFL7","ESAM","ECSCR")
length(g1)
plotFeature.batch(seu.obj = seu_obj, dr="umap", genes.to.plot = g1, nCols = blue.cols, 
                  plot.name = "Plot_UMAP_C2_and_C7-tCIO_cell_class_marker_expression.png", col.num = 7)

group_anno <- list(
  "Epithelial" = c(1,2,3,8,9,11,13,17,18,19),
  "Mesenchymal" = c(0,4,5,6,7,10,12,15,16),
  "Immune" = 14
)

cl_vec <- seu_obj$RNA_snn_res.1
group_vec <- rep(NA, length(cl_vec))
for(ct in names(group_anno)){
  cl <- group_anno[[ct]]
  group_vec[which(cl_vec%in%cl)] <- ct
}

seu_obj$Cell_class <- group_vec
DimPlot(seu_obj, group.by = "Cell_class", label=T)
saveRDS(seu_obj, file="Res_C2_and_C7_tCIO_panTro6_scRNA-seq_merged_without_integration.rds")

# harmony and CSS integration
## CSS integration
seu_obj_css <- cluster_sim_spectrum(seu_obj, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
seu_obj_css <- RunUMAP(seu_obj_css, reduction = "css", dims = 1:ncol(seu_obj_css@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
seu_obj_css <- FindNeighbors(object = seu_obj_css, reduction = "css", dims = 1:ncol(seu_obj_css@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
saveRDS(seu_obj_css, file="Res_C2_and_C7_tCIO_panTro6_scRNA-seq_with_CSS_integration.rds")

# harmony integration
library(harmony)
seu_obj_harmony <- RunHarmony(seu_obj,
                              group.by.vars = "Sample.name",
                              reduction="pca",
                              dims.use = 1:20,
                              assay.use = "RNA",
                              project.dim = F)
seu_obj_harmony <- RunUMAP(seu_obj_harmony,
                           reduction = "harmony",
                           dims = 1:20,
                           reduction.key = "UMAPHARMONY_",
                           reduction.name="umap_harmony")
seu_obj_harmony <- FindNeighbors(object = seu_obj_harmony, reduction = "harmony", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
saveRDS(seu_obj_harmony, file="Res_C2_and_C7_tCIO_panTro6_scRNA-seq_with_harmony.rds")
p0 <- DimPlot(seu_obj, reduction = "umap", shuffle = T, group.by = "Sample.name")+ggtitle(label = "Non-integrated")
p0_2 <- DimPlot(seu_obj, reduction = "umap", shuffle = T, group.by = "Cell_class")+ggtitle(label = "Non-integrated")
p1 <- DimPlot(seu_obj_css, reduction = "umap_css", shuffle = T, group.by = "Sample.name")+ggtitle(label="CSS")
p2 <- DimPlot(seu_obj_harmony, reduction = "umap_harmony", shuffle = T, group.by = "Sample.name")+ggtitle(label="Harmony")
p0+p0_2+ p1+p2


# extract epithelial cells
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi")
tCIO_rna <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/Res_C2_and_C7_tCIO_panTro6_scRNA-seq_with_harmony.rds")
epi <- subset(tCIO_rna, subset = Cell_class=="Epithelial")
DimPlot(epi)
epi <- FindVariableFeatures(object = epi, selection.method = "vst", nfeatures = 3000)
epi <- ScaleData(object = epi, verbose = T)
epi <- RunPCA(object = epi, features = VariableFeatures(epi), verbose = F, npcs = 50)
usefulPCs <- 1:20
epi <- FindNeighbors(object = epi, dims = usefulPCs)
epi <- FindClusters(object = epi, resolution = 1)
epi <- RunUMAP(object = epi, dims = usefulPCs)
saveRDS(epi, file="Res_C2_and_C7_tCIO_epi_panTro6_without_integration.rds")
g1 <- c("EPCAM","COL1A2","CDH5","PTPRC","ASCL1")
plotFeature.batch(seu.obj=epi, dr="umap", genes.to.plot=g1, 
                  col.num = 5, plot.name = "Plot_UMAP_major_cell_type_marker_expr.png",
                  nCols = blue.cols, cex=2)

g1 <- c("CDX2", "CDH1", "EPCAM", "COL1A2", "DCN","PTPRC", "CORO1A","PHOX2B", "CDH5", "PECAM1", "SATB2", "ONECUT2",
        "SOX2", "FOXJ1", "TP63", "CLDN18", "LGALS4","LGR5", "ASCL2", "OLFM4", "MUC2", "SPINK4", "SPDEF",
        "CHGA", "ISL1", "BEST4", "SPIB", "GP2", "DPP4", "APOA4", "FABP2", "SI","DEFA5", "LYZ", "RGS13", "MKI67","PCNA")
length(g1)
plotFeature.batch(seu.obj = epi, dr="umap", genes.to.plot = g1, nCols = blue.cols, 
                  plot.name = "Plot_UMAP_C7-tCIO_epi_cell_class_and_epi_cell_type_marker_expression.png", col.num = 5)

#C13 is a doublet cluster with mesenchymal cells
epi_2 <- subset(epi, subset = RNA_snn_res.1!=13)
DimPlot(epi_2, label=T)
epi_2 <- FindVariableFeatures(object = epi_2, selection.method = "vst", nfeatures = 3000)
epi_2 <- ScaleData(object = epi_2, verbose = T)
epi_2 <- RunPCA(object = epi_2, features = VariableFeatures(epi_2), verbose = F, npcs = 50)
usefulPCs <- 1:20
epi_2 <- FindNeighbors(object = epi_2, dims = usefulPCs)
epi_2 <- FindClusters(object = epi_2, resolution = 1)
epi_2 <- RunUMAP(object = epi_2, dims = usefulPCs)
saveRDS(epi_2, file="Res_C2_C7_tCIO_panTro6_epi_2_removing_mesenchyme_doublet_cluster.rds")
g1 <- c("EPCAM","COL1A2","CDH5","PTPRC","ASCL1")
plotFeature.batch(seu.obj=epi_2, dr="umap", genes.to.plot=g1, 
                  col.num = 5, plot.name = "Plot_UMAP_epi_2_major_cell_type_marker_expr.png",
                  nCols = blue.cols, cex=2)

# Epithelial cell type annotation
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi")
tCIO_epi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/Res_C2_C7_tCIO_panTro6_epi_2_removing_mesenchyme_doublet_cluster.rds")

# MNN, harmony and CSS integration
## CSS integration
tCIO_epi <- cluster_sim_spectrum(tCIO_epi, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
tCIO_epi <- RunUMAP(tCIO_epi, reduction = "css", dims = 1:ncol(tCIO_epi@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
tCIO_epi <- FindNeighbors(object = tCIO_epi, reduction = "css", dims = 1:ncol(tCIO_epi@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
tCIO_epi[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- tCIO_epi[[paste0("RNA_snn_res.", 0.2*5)]]

# harmony integration
tCIO_epi <- RunHarmony(tCIO_epi,
                       group.by.vars = "Sample.name",
                       reduction="pca",
                       dims.use = 1:20,
                       assay.use = "RNA",
                       project.dim = F)
tCIO_epi <- RunUMAP(tCIO_epi,
                    reduction = "harmony",
                    dims = 1:20,
                    reduction.key = "UMAPHARMONY_",
                    reduction.name="umap_harmony")
tCIO_epi <- FindNeighbors(object = tCIO_epi, reduction = "harmony", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
tCIO_epi[[paste0("RNA_Harmony_snn_res.", 0.2*5)]] <- tCIO_epi[[paste0("RNA_snn_res.", 0.2*5)]]

# MNN integration
seu_list <- SplitObject(tCIO_epi,split.by="orig.ident")
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method="vst", nfeatures=2000)
}

seu_obj <- RunFastMNN(object.list = seu_list)
tCIO_epi[['mnn']] <- CreateDimReducObject(Embeddings(seu_obj, "mnn")[colnames(tCIO_epi),], key="MNN_")
tCIO_epi <- RunUMAP(tCIO_epi, reduction="mnn", dims = 1:20, reduction.name = "umap_mnn", reduction.key = "UMAPMNN_")
tCIO_epi <- FindNeighbors(object = tCIO_epi, reduction = "mnn", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
tCIO_epi[[paste0("RNA_MNN_snn_res.", 0.2*5)]] <- tCIO_epi[[paste0("RNA_snn_res.", 0.2*5)]]
p1 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "Sample.name")+NoLegend()
p2 <- DimPlot(tCIO_epi, reduction = "umap_mnn", group.by = "Sample.name")+ theme(legend.position="bottom")
p3 <- DimPlot(tCIO_epi, reduction = "umap_harmony", group.by = "Sample.name")+NoLegend()
p1+p2+p3+plot_spacer()+plot_layout(ncol=2)
saveRDS(tCIO_epi, file="Res_C2_and_C7_tCIO_epi_panTro6_scRNA-seq_with_CSS_harmony_MNN.rds")

g1 <- c("CDX2","PDX1","ONECUT2", "SATB2","SOX2", "FOXJ1", "TP63", "CLDN18", "HOXA10", "OSR2","LGALS4", "SHH","LGR5", "ASCL2", "OLFM4", "MUC2", "SPINK4", "SPDEF", "MKI67","CDK1",
        "CHGA", "ISL1", "BEST4", "SPIB", "GP2", "CA7", "MUC1","DPP4", "APOA4", "FABP2", "SI","DEFA5", "DEFA6","LYZ", "RGS13")
plotFeature.batch(seu.obj = tCIO_epi, dr="umap_css", genes.to.plot = g1, nCols = blue.cols, cex=3,
                  plot.name = "Plot_UMAP_CSS_tCIO_epi_cell_type_marker_expression.png", 
                  col.num = 8)

de_res <- wilcoxauc(tCIO_epi, group_by = "RNA_CSS_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[which(de_res$padj<0.01 & de_res$auc>0.6 & de_res$logFC>0.1 & de_res$pct_diff>20 & de_res$pct_in>20),]
saveRDS(sig_res, file="Res_RNA_CSS_res.1_cluster_markers.rds")
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=auc)
genes <- top_res$feature[which(top_res$group==4)]
plotFeature.batch(seu.obj = tCIO_epi, dr="umap_css", genes.to.plot = genes, nCols = blue.cols, cex=3,
                  plot.name = "Plot_UMAP_CSS_tCIO_epi_C4_marker_expression.png", 
                  col.num = 5)

# integrate with differentiated HIOs
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_differentiated_HIO")
HIO_epi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Human_tHIO_hg38/scRNA-seq/with_NRG1_HIO/epi/Res_differentiated_HIO_epithelium_scRNA-seq.rds")
combined <- merge(x=HIO_epi, y=tCIO_epi)
hvg <- intersect(VariableFeatures(HIO_epi), VariableFeatures(tCIO_epi))
length(hvg)
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
VariableFeatures(combined) <- hvg
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)

# MNN, harmony and CSS integration
## CSS integration
combined <- cluster_sim_spectrum(combined, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]

# harmony integration
combined <- RunHarmony(combined,
                       group.by.vars = "Sample.name",
                       reduction="pca",
                       dims.use = 1:20,
                       assay.use = "RNA",
                       project.dim = F)
combined <- RunUMAP(combined,
                    reduction = "harmony",
                    dims = 1:20,
                    reduction.key = "UMAPHARMONY_",
                    reduction.name="umap_harmony")
combined <- FindNeighbors(object = combined, reduction = "harmony", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_Harmony_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]

# MNN integration
seu_list <- SplitObject(combined,split.by="orig.ident")
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method="vst", nfeatures=2000)
}

seu_obj <- RunFastMNN(object.list = seu_list)
combined[['mnn']] <- CreateDimReducObject(Embeddings(seu_obj, "mnn")[colnames(combined),], key="MNN_")
combined <- RunUMAP(combined, reduction="mnn", dims = 1:20, reduction.name = "umap_mnn", reduction.key = "UMAPMNN_")
combined <- FindNeighbors(object = combined, reduction = "mnn", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_MNN_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]

# run rPCA integration
# split the dataset into a list of seurat objects
seu_obj_list <- SplitObject(combined, split.by = "Sample.name")
# normalize and identify variable features for each dataset independently
seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seu_obj_list)
seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = seu_obj_list, anchor.features = features, reduction = "rpca",
                                  k.anchor = 20)
combined.rpca <- IntegrateData(anchorset = anchors)
DefaultAssay(combined.rpca) <- "integrated"
combined.rpca <- ScaleData(combined.rpca, verbose = FALSE)
combined.rpca <- RunPCA(combined.rpca, npcs = 30, verbose = FALSE)
combined.rpca <- RunUMAP(combined.rpca, reduction = "pca", dims = 1:30)
combined.rpca <- FindNeighbors(combined.rpca, reduction = "pca", dims = 1:30)
combined.rpca <- FindClusters(combined.rpca, resolution = 0.5)
combined[['umap_rpca']] <- CreateDimReducObject(embeddings = Embeddings(combined.rpca, reduction = "umap"),
                                                key="UMAPRPCA_")
p1 <- DimPlot(combined, reduction = "umap_css", group.by = "Sample.name")+NoLegend()
p2 <- DimPlot(combined, reduction = "umap_harmony", group.by = "Sample.name")+theme(legend.position = "bottom")
p3 <- DimPlot(combined, reduction = "umap_mnn", group.by = "Sample.name")+NoLegend()
p4 <- DimPlot(combined, reduction = "umap_rpca", group.by = "Sample.name")+NoLegend()
p1+p2+p3+p4+plot_layout(ncol=2)
saveRDS(combined, file="Res_tCIO_and_differentiated_HIO_integrated_with_CSS_MNN_Harmony_rPCA.rds")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO")
tHIO_epi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Human_tHIO_hg38/scRNA-seq/epi/intestinal/Res_tHIO_intestinal_epi_hg38_with_CSS_and_MNN_integration.rds")
combined <- merge(x=tHIO_epi, y=tCIO_epi)
hvg <- intersect(VariableFeatures(tHIO_epi), VariableFeatures(tCIO_epi))
length(hvg)
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
VariableFeatures(combined) <- hvg
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)

# MNN, harmony and CSS integration
## CSS integration
combined <- cluster_sim_spectrum(combined, label_tag = "Sample.name", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]
DimPlot(combined, reduction = "umap_css", group.by = "Individual.or.stem.cell.line")
# harmony integration
combined <- RunHarmony(combined,
                       group.by.vars = "Sample.name",
                       reduction="pca",
                       dims.use = 1:20,
                       assay.use = "RNA",
                       project.dim = F)
combined <- RunUMAP(combined,
                    reduction = "harmony",
                    dims = 1:20,
                    reduction.key = "UMAPHARMONY_",
                    reduction.name="umap_harmony")
combined <- FindNeighbors(object = combined, reduction = "harmony", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_Harmony_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]

# MNN integration
seu_list <- SplitObject(combined,split.by="orig.ident")
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method="vst", nfeatures=2000)
}

seu_obj <- RunFastMNN(object.list = seu_list)
combined[['mnn']] <- CreateDimReducObject(Embeddings(seu_obj, "mnn")[colnames(combined),], key="MNN_")
combined <- RunUMAP(combined, reduction="mnn", dims = 1:20, reduction.name = "umap_mnn", reduction.key = "UMAPMNN_")
combined <- FindNeighbors(object = combined, reduction = "mnn", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_MNN_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]

# run rPCA integration
# split the dataset into a list of seurat objects
seu_obj_list <- SplitObject(combined, split.by = "Sample.name")
# normalize and identify variable features for each dataset independently
seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seu_obj_list)
seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = seu_obj_list, anchor.features = features, reduction = "rpca",
                                  k.anchor = 20)
combined.rpca <- IntegrateData(anchorset = anchors)
DefaultAssay(combined.rpca) <- "integrated"
combined.rpca <- ScaleData(combined.rpca, verbose = FALSE)
combined.rpca <- RunPCA(combined.rpca, npcs = 30, verbose = FALSE)
combined.rpca <- RunUMAP(combined.rpca, reduction = "pca", dims = 1:30)
combined.rpca <- FindNeighbors(combined.rpca, reduction = "pca", dims = 1:30)
combined.rpca <- FindClusters(combined.rpca, resolution = 0.5)
combined[['umap_rpca']] <- CreateDimReducObject(embeddings = Embeddings(combined.rpca, reduction = "umap"),
                                                key="UMAPRPCA_")
p1 <- DimPlot(combined, reduction = "umap_css", group.by = "Sample.name")+NoLegend()
p2 <- DimPlot(combined, reduction = "umap_harmony", group.by = "Sample.name")+theme(legend.position = "bottom")
p3 <- DimPlot(combined, reduction = "umap_mnn", group.by = "Sample.name")+NoLegend()
p4 <- DimPlot(combined, reduction = "umap_rpca", group.by = "Sample.name")+NoLegend()
p1+p2+p3+p4+plot_layout(ncol=2)

# 3D umap
combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), n.components = 3, reduction.name = "umap_css_3d", reduction.key = "UMAPCSS3D_")

values <- combined$Individual.or.stem.cell.line
library(plotly)
plot_ly(as.data.frame(combined@reductions$umap_css_3d@cell.embeddings), x = ~UMAPCSS3D_1, y = ~UMAPCSS3D_2, z = ~UMAPCSS3D_3, size = 1,
        color = values, colors = setNames(scales::hue_pal()(length(unique(values))), unique(values)),
        sizes = c(10, 200), type = "scatter3d")

saveRDS(combined, file="Res_tCIO_and_tHIO_integrated_with_CSS_MNN_Harmony_rPCA_and_3dUMAPCSS_embedding.rds")

source("/links/groups/treutlein/USERS/zhisong_he/Tools/scripts/feature_plots.r")
cols <- setNames(colorRampPalette(prettyrainbow)(length(unique(values))), sort(unique(values)))
plotFeature3d_rotating(Embeddings(combined,"umap_css_3d"), values, colorPal=cols, pt_cex=1,
                       output_file="./plot_tCIO_tHIO_epi_umap-css_3d_stages_withlegend.gif",
                       do_legend = T, legend_cex = 2.5, legend_pt_cex = 4)

main_dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO"
tCIO_epi_list <- SplitObject(tCIO_epi, split.by = "Individual.or.stem.cell.line")
combined_list <- list()
for(x in names(tCIO_epi_list)){
  sub_dir <- paste(x,"with_tHIO",sep="_")
  wd <- file.path(main_dir, sub_dir)
  setwd(wd)
  seu_obj <- tCIO_epi_list[[x]]
  seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
  combined <- merge(x=tHIO_epi, y=seu_obj)
  hvg <- intersect(VariableFeatures(tHIO_epi), VariableFeatures(seu_obj))
  length(hvg)
  combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
  VariableFeatures(combined) <- hvg
  combined <- ScaleData(object = combined, verbose = T)
  combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
  usefulPCs <- 1:20
  combined <- FindNeighbors(object = combined, dims = usefulPCs)
  combined <- FindClusters(object = combined, resolution = 1)
  combined <- RunUMAP(object = combined, dims = usefulPCs)
  ## CSS integration
  combined <- cluster_sim_spectrum(combined, label_tag = "Sample.name", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
  combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
  combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
    FindClusters(resolution = 1)
  combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]
  DimPlot(combined, reduction = "umap_css", group.by = "Individual.or.stem.cell.line")
  combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), n.components = 3, reduction.name = "umap_css_3d", reduction.key = "UMAPCSS3D_")
  values <- combined$Individual.or.stem.cell.line
  cols <- setNames(colorRampPalette(prettyrainbow)(length(unique(values))), sort(unique(values)))
  plotFeature3d_rotating(Embeddings(combined,"umap_css_3d"), values, colorPal=cols, pt_cex=1,
                         output_file="./plot_tCIO_tHIO_epi_umap-css_3d_stages_withlegend.gif",
                         do_legend = T, legend_cex = 2.5, legend_pt_cex = 4)
  
  combined_list[[x]] <- combined
  saveRDS(combined, file=paste0("Res_",x,"_tCIO_and_tHIO_integrated_with_CSS.rds"))
}

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO")
combined <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO/Res_tCIO_and_tHIO_integrated_with_CSS_MNN_Harmony_rPCA_and_3dUMAPCSS_embedding.rds")
combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 0.1)
combined[[paste0("RNA_CSS_snn_res.", 0.1)]] <- combined[[paste0("RNA_snn_res.", 0.1)]]

p1 <- DimPlot(combined, reduction = "umap_css", group.by = "Sample.name")+theme(legend.position = "bottom")
p2 <- DimPlot(combined, reduction = "umap_css", group.by = "RNA_CSS_snn_res.0.6", label = T)+NoLegend()
p1+p2
saveRDS(combined, file="Res_tCIO_and_tHIO_integrated_with_CSS_MNN_Harmony_rPCA_and_3dUMAPCSS_embedding.rds")

tCIO_epi$Cl_id_with_tHIO <- combined@meta.data[colnames(tCIO_epi), "RNA_CSS_snn_res.0.6"]
p3 <- DimPlot(tCIO_epi, reduction="umap_css", group.by = "Cl_id_with_tHIO", label=T)+NoLegend()
p4 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", label=T)+NoLegend()
p1+p2+p3+p4+plot_layout(ncol=2)
n <- sapply(sort(unique(tCIO_epi$Cl_id_with_tHIO)), function(x){
  sapply(sort(unique(tCIO_epi$Sample.name)), function(y){
    sum(tCIO_epi$Sample.name==y & tCIO_epi$Cl_id_with_tHIO==x)
  })
})
colnames(n) <- paste0("C", sort(unique(tCIO_epi$Cl_id_with_tHIO)))
p <- t(n/rowSums(n))

# annotate tCIO epithelial cell types based on integrated data with tHIO
## option 1: label transfer
query <- tCIO_epi
anchors <- FindTransferAnchors(reference = tHIO_epi, query = query, 
                               dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = tHIO_epi$Cell_type, 
                            dims = 1:30)
tCIO_epi$labelTransfer_tHIO_cell_type <- predictions[,"predicted.id"]
p1 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "labelTransfer_tHIO_cell_type")

## option 2: kNN search on CSS space
ref_mat <- combined@reductions$css@cell.embeddings[combined$Species=="Human",]
que_mat <- combined@reductions$css@cell.embeddings[combined$Species=="Chimp",]
knn <- RANN::nn2(ref_mat, que_mat, k = 20)$nn.idx
dim(knn)
ref_idx <- combined@meta.data[rownames(ref_mat), "Cell_type"]
nn_idx <- matrix(ref_idx[as.vector(knn)], nrow=nrow(knn))
pred_id <- apply(nn_idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
tCIO_epi@meta.data[rownames(que_mat), "kNN_tHIO_cell_type"] <- pred_id
p2 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "kNN_tHIO_cell_type")
p1+p2
p3 <- FeaturePlot(tCIO_epi, reduction = "umap_css", features = c("GSTA1","GSTA2","AADAC", "REEP6"), order = T)
p4 <- DimPlot(tCIO_epi, reduction="umap_css", group.by = "Cl_id_with_tHIO", label=T)+NoLegend()

human_file <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi/Res_fetal_8-19_PCW_hg38_epi_with_CSS_integration.rds"
fetal_human <- readRDS(human_file)
FeaturePlot(fetal_human, reduction = "umap_css", features = c("ISL1", 'CHGA'), order = T)
query <- tCIO_epi
anchors <- FindTransferAnchors(reference = fetal_human, query = query, 
                               dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = fetal_human$Cell_type, 
                            dims = 1:30)
tCIO_epi$labelTransfer_fetal_human_cell_type <- predictions[,"predicted.id"]
p1 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "labelTransfer_tHIO_cell_type", label=T)
p2 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "labelTransfer_fetal_human_cell_type", label=T)
p1+p2

# option 3: co-clustering
tCIO_epi$Cl_id_with_tHIO <- combined@meta.data[colnames(tCIO_epi), "RNA_CSS_snn_res.1"]
n <- sapply(sort(unique(combined$RNA_CSS_snn_res.1)), function(i){
  sapply(sort(unique(combined$Cell_type)), function(x){
    length(which(combined$RNA_CSS_snn_res.1==i & combined$Cell_type==x))
  })
})
colnames(n) <- paste0("C", sort(unique(combined$RNA_CSS_snn_res.1)))
ct_anno <- setNames(rownames(n)[apply(n, 2, which.max)], colnames(n))
tCIO_epi$coCluster_tHIO_cell_type <- ct_anno[paste0("C", tCIO_epi$Cl_id_with_tHIO)]
tCIO_epi$pred_tHIO_cell_type <- tCIO_epi$coCluster_tHIO_cell_type
tCIO_epi$pred_tHIO_cell_type <- tCIO_epi$kNN_tHIO_cell_type
tCIO_epi$pred_tHIO_cell_type[which(tCIO_epi$labelTransfer_tHIO_cell_type=="Paneth_cell")] <- "Paneth_cell"
combined@meta.data[colnames(tCIO_epi), "Cell_type"] <- tCIO_epi$pred_tHIO_cell_type
saveRDS(tCIO_epi, file="Res_tCIO_epi_integrated_with_CSS_with_pred_tHIO_cell_type.rds")
saveRDS(combined, file="Res_tCIO_and_tHIO_epi_integrated_with_CSS.rds")

p1 <- DimPlot(combined, reduction = "umap_css", group.by = "Sample.name")+theme(legend.position = "bottom")
p2 <- DimPlot(combined, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", label=T)+NoLegend()
p3 <- DimPlot(combined, reduction = "umap_css", group.by = "Cell_type")+theme(legend.position = "bottom")
p4 <- DimPlot(tCIO_epi, reduction = "umap_css", group.by = "pred_tHIO_cell_type")
p1+p2+p3+p4+plot_layout(ncol=2)
n <- sapply(c("Stem_cell","Early_enterocyte","Enterocyte"), function(x){
  sapply(sort(unique(combined$Sample.name)), function(y){
    sum(combined$Cell_type==x & combined$Sample.name==y)
  })
})
p <- t(n/rowSums(n))
gCols <- setNames(colorRampPalette(prettyrainbow)(nrow(p)), rownames(p))
df <- data.frame("Sample"=rep(colnames(p), each=nrow(p)),
                 "Cell_type"=rep(rownames(p), ncol(p)),
                 "Proportion"=as.vector(p),
                 stringsAsFactors = F)
ggplot(df, aes(x=Sample, y=Proportion, fill=Cell_type))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=gCols)+
  theme(axis.text.x = element_text(angle =45, hjust=1))
# extract cells on stem-cell-to-enterocyte trajectory
cells <- colnames(combined)[which(combined$Cell_type %in% c("Stem_cell","Early_enterocyte","Enterocyte"))]
s2e <- subset(combined, cells=cells)
DimPlot(s2e, reduction = "umap_css", group.by = "Cell_type", shuffle = T)
# identify stem-cell-to-enterocyte trajectory Pt-dependent genes in each sample separately
seu_obj_list <- SplitObject(s2e, split.by="Sample.name")
s2e$"Pt_by_sample" <- NA
pt_genes <- list()
pt_expr_list <- list()
for(x in names(seu_obj_list)){
  print(paste(x,"start"))
  print("Reconstruct pseudotime")
  seu_obj <- seu_obj_list[[x]]
  input <- seu_obj@reductions$css@cell.embeddings
  dm <- DiffusionMap(input, k=20)
  vec <- rank(dm$DC1)/ncol(seu_obj)
  aa <- median(vec[which(seu_obj$Cell_type=="Stem_cell")]) > median(vec[which(seu_obj$Cell_type=="Enterocyte")])
  if(aa){
    seu_obj$Pt_by_sample <- rank(-dm$DC1)/ncol(seu_obj)
    
  }else{
    seu_obj$Pt_by_sample <- rank(dm$DC1)/ncol(seu_obj)
  }
  s2e@meta.data[colnames(seu_obj),"Pt_by_sample"] <- seu_obj$Pt_by_sample
  seu_obj_list[[x]] <- seu_obj
  
  print("Get pseudotime-dependent genes")
  expr_mat <- as.matrix(seu_obj@assays$RNA@data)
  pt_vec <- seu_obj$Pt_by_sample
  registerDoParallel(50)
  res <- foreach(k=seq(nrow(expr_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.vector(expr_mat[k,])
    m0 <- lm(e ~ 1)
    m1 <- lm(e ~ ns(pt_vec, df=6))
    a0 <- anova(m0)
    a1 <- anova(m1)
    p_anova <- anova(m1,m0)$Pr[2]
    p_resi <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], df1=a0["Residuals", "Df"], df2=a1["Residuals","Df"], lower.tail = F)
    return(c(p_anova, p_resi))
  }
  rownames(res) <- rownames(expr_mat)
  colnames(res) <- c("p_ANOVA", "p_Resi")
  df <- data.frame(res)
  df$p_ANOVA_BH <- p.adjust(df$p_ANOVA, method="BH")
  df$p_Resi_BH <- p.adjust(df$p_Resi, method="BH")
  saveRDS(df, file=paste0("Res_",x,"_stem_cell_to_enterocyte_Pt_test_pvals.rds"))
  g1 <- rownames(df)[which(df$p_ANOVA_BH<0.05)]
  g2 <- rownames(df)[which(df$p_Resi_BH<0.05)]
  print(length(g1))
  print(length(g2))
  pt_genes[[x]] <- list("ANOVA"=g1,
                        "Resi"=g2)
  
  print("Get Pt bin average expr")
  num_breaks <- 20
  pt_expr <- sapply(1:num_breaks, function(i){
    idx <- which(ceiling(pt_vec * num_breaks) == i)
    if (i == 1)
      idx <- c(idx, which(pt_vec == 0))
    if (length(idx) > 1){
      return(rowMeans(expr_mat[,idx]))
    } else if (length(idx) == 1){
      return(expr_mat[,idx])
    } else{
      return(rnorm(nrow(expr_mat)))
    }
  })
  saveRDS(pt_expr, file=paste0("Dat_",x,"_se_pt_expr.rds"))
  pt_expr_list[[x]] <- pt_expr
}
saveRDS(pt_genes, file="Res_Pt_genes_in_each_sample_of_human_and_chimp.rds")
union_pt_genes <- unique(unlist(lapply(names(pt_genes), function(x){
  pt_genes[[x]][["Resi"]]
})))

combined_expr_mat <- do.call('cbind', pt_expr_list)
rownames(combined_expr_mat) <- rownames(s2e)
saveRDS(pt_expr_list, file="Res_pt_by_sample_bin_average_expr.rds")
#saveRDS(seu_obj_list, file="Res_seurat_object_by_sample.rds")

ptg_expr_mat <- combined_expr_mat[union_pt_genes,]
hc <- hclust(as.dist(1-cor(t(ptg_expr_mat))))
plot(hc, hang=1, cex=0.1)
ptg_module_id <- cutree(hc, h=1.75)
g.list <- lapply(sort(unique(ptg_module_id)),
                 function(i){names(ptg_module_id[which(ptg_module_id==i)])})
g.size <- sapply(g.list, length)
g.size

group.cols <- setNames(c("#fcc5c0","#fa9fb5","#c51b8a","#7a0177",
                         "#7fcdbb","#253494"),
                       names(pt_expr_list))
group.vec <- rep(names(pt_expr_list), each=num_breaks)
time.vec <- rep(seq(num_breaks), length(pt_expr_list))
scale.expr = t(scale(t(ptg_expr_mat)))

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
library(tidyverse)
df <- tibble(Relative_expr = as.vector(mean.combined.mat), 
             Pt = rep(rep(seq(num_breaks), length(seu_obj_list)), ncol(mean.combined.mat)),
             Cluster = rep(seq(ncol(mean.combined.mat)), each = nrow(mean.combined.mat)),
             Sample = rep(rep(names(seu_obj_list), each=num_breaks), ncol(mean.combined.mat))) %>% 
  mutate(Cluster = factor(Cluster, levels = seq(ncol(mean.combined.mat))))

plot.name <- "Plot_ggplot_tHIO_tCIO_Pt_dependent_gene_cluster_average_expr.pdf"
pdf(plot.name, height=10, width=20)
ggplot(df, aes(x=Pt, y=Relative_expr, color = Sample)) +
  geom_line() +
  scale_color_manual(values=group.cols)+
  facet_wrap(~Cluster)+
  theme_minimal()
dev.off()


# load in fetal human multi-region intestine data
human_multi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_human_multi_intestine_region_hg19/with_colon/Res_fetal_human_multi_intestinal_region_with_CSS_integration_and_UMAP_model.rds")
p1 <- DimPlot(human_multi, reduction = "umap_css", group.by = "Tissue")
p2 <- DimPlot(human_multi, reduction = "umap_css", group.by = "Cell_type")
p1+p2
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO/intestine_regional_identity")
## load fetal atlas highly variable genes
fetal.hvg <- VariableFeatures(human_multi)
## load fetal atlas meta.data after decompression, which contains CSS-based UMAP embedding of the fetal atlas and organ identity
ref.idx <- human_multi$Tissue
## load fetal Cluster Similarity Spectrum (CSS) model
css.model <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_human_multi_intestine_region_hg19/with_colon/Res_fetal_human_multi_intestinal_region_CSS_model.rds")

## use the fetal highly variable genes that are also detected in tIO data to calculate similarity between tIO cells and fetal clusters per sample
## because here the fetal cluster average profile is reference, and the tIO cells are query, the obtained similarity spectrum is called reference similarity spectrum (RSS)
shared.genes <- intersect(rownames(combined), fetal.hvg)
length(shared.genes)
que.data <- as.matrix(combined@assays$RNA@data[shared.genes,])
rss.list <- lapply(seq(length(css.model$model$profiles)), function(i){
  ref <- css.model$model$profiles[[i]][shared.genes,]
  cor.mat <- cor(que.data, ref, method = "spearman")
  cor.z <- t(scale(t(cor.mat)))
  return(cor.z)
})
rss.mat <- do.call('cbind', rss.list)
combined[["rss"]] <- CreateDimReducObject(embeddings = rss.mat, key="RSS_", assay=DefaultAssay(combined))

## calculate the Euclidean distance between tIO cells and fetal cells on the space of similarities to fetal clusters per sample, instead of on the space of expression. correlation is potentially more robust to batch variations than expression 
## after the distance calculation, get the 20-nearest fetal cells for each tIO cell
## assign the majority of tissue identity of the 20-nearest cells to be the inferred tissue identity of tHIO cells 
css <- css.model$sim2profiles 
knn <- RANN::nn2(css, rss.mat, k = 20)$nn.idx
## map the tissue identity of fetal cells to tHIO cells
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
trans.id <- apply(nn.idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
combined@meta.data$Mapped_fetal_intestine_region <- trans.id
p3 <- DimPlot(combined, reduction = "umap_css", group.by = "Cell_type")
p4 <- DimPlot(combined, reduction = "umap_css", group.by = "Mapped_fetal_intestine_region")
p1+p2+p3+p4+plot_layout(ncol=2)
combined <- RunUMAP(combined, 
                    reduction = "rss", 
                    dims = 1:ncol(combined@reductions$rss@cell.embeddings), 
                    reduction.name = "umap_rss", 
                    reduction.key = "UMAPRSS_")

p4 <- DimPlot(combined, reduction = "umap_rss", group.by = "Cell_type")
p5 <- DimPlot(combined, reduction = "umap_rss", group.by = "Sample.name")
saveRDS(combined, file="Res_tIO_with_CSS_integration_and_fetal_human_intestine_region_identity_projection.rds")





