setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum")
library(Seurat)
library(simspec)
source("~/Work/commonScript/Script_functions.R")

# step 1. read scRNA-seq data of primary samples
library(readxl)
readout <- "scRNA-seq_respective_genome"
species <- "Mouse"
tissue <- "Fetal primary"

sample_info <- read_excel("~/Work/Endoderm/intestine_evolution/sample/HIO_CIO_metadata.xlsx", sheet=readout)
selected_sample_info <- sample_info[which(sample_info$Species==species & sample_info$Tissue %in% tissue),]
dim(selected_sample_info)

j <- 0
# read data
seu_obj_list <- list()
for(i in seq(nrow(selected_sample_info))){
  aa <- selected_sample_info$`Sample name`[i]
  print(paste("Reading", aa))
  vec <- as.vector(as.matrix(selected_sample_info[i,3:11]))
  info.mat <- cbind(names(selected_sample_info)[3:11], vec)
  path <- paste0(selected_sample_info[i,12])
  seu_obj <- prepareSeuratObject2(rawData = Read10X_h5(path), 
                                  namePrefix = paste0("MR-",i+j),
                                  species="mouse", mito.gene.list="/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/Ensembl/Mouse/List_mouse_Mt_genes.txt",
                                  additional.sample.info.mat=info.mat, seu.version=3, 
                                  mito.cutoff=Inf, min.cell.cutoff=1, 
                                  min.gene.cutoff=500, gene.num.low=500, gene.num.high=Inf)
  seu_obj_list[[i]] <- seu_obj
}

combined <- merge(x=seu_obj_list[[1]], y=seu_obj_list[-1])
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size=0)
dim(combined)
combined_sub <- subset(combined, subset = percent.mito < 0.2 & nFeature_RNA > 1000 & nCount_RNA < 75000)
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
human_mouse_gene_link <- readRDS("~/Work/Annotation/Ensembl/Human/Dat_human_mouse_one2one_ortholog.rds")
mouse_gene_to_remove <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% genes.to.remove, "Mouse_symbol"]
counts <- combined@assays$RNA@counts
counts <- counts[setdiff(rownames(counts), mouse_gene_to_remove),]
meta <- combined@meta.data[,4:13]
combined <- CreateSeuratObject(counts=counts, meta=meta)
combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 1e4)
saveRDS(combined, file="Dat_mouse_fetal_duodenum_13-18d_normalized_non-integrated_data.rds")

combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)
saveRDS(combined, file="Res_mouse_fetal_duodenum_13-18d_no_batch_correction.rds")
DimPlot(combined, reduction = "umap", group.by = "Age")+DimPlot(combined, reduction = "umap_css", group.by = "Age")

# CSS integration
combined <- cluster_sim_spectrum(combined, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(combined, file="Res_mouse_fetal_duodenum_13-18d_with_CSS_integration.rds")

DimPlot(combined, reduction = "umap_css", group.by = "Age")+DimPlot(combined, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", label=T)
# plot cell class marker expression
g1 <- c("EPCAM","COL1A2","CDH5","PTPRC","ASCL1")
mouse_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% g1, "Mouse_symbol"]
plotFeature.batch(seu.obj=combined, dr="umap_css", genes.to.plot=mouse_g1, 
                  col.num = 5, plot.name = "Plot_UMAP_CSS_major_cell_type_marker_expr.png",
                  nCols = blue.cols, cex=1)

# classify clusters into cell classes
known.markers <- read.table("~/Work/Annotation/cellTypeMarker/Intestine/Table_major_cell_type_markers_from_literature_search.txt", head=T, sep="\t",stringsAsFactors=F)  
cell.type.order <- c("Neuronal", "Epithelial", "Mesenchymal", "Endothelial", "Immune", "Erythroid", "Hepatocyte")
pan.cm <- unique(unlist(lapply(cell.type.order, function(ct){
  g1 <- known.markers$Gene[which(known.markers$Used_pan_cell_type_markers==ct)]
  mouse_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% g1, "Mouse_symbol"]
  intersect(mouse_g1, rownames(combined))
}))) 

# cell class annotation
cell_class <- list(
  "Neural" = c(6, 13, 15),
  "Epithelial" = c(2,4,9,10,18, 27,28, 32),
  "Mesenchymal" = c(0, 1, 3, 5, 7, 11, 16, 17, 19, 21, 24, 30),
  "Endothelial" = 22,
  "Immune" = c(12, 23),
  "Doublet"= c(8, 14, 20, 25, 26, 29, 31)
)

cl_vec <- combined$RNA_CSS_snn_res.1
class_vec <- rep(NA, length(cl_vec))
for(x in names(cell_class)){
  cl_x <- cell_class[[x]]
  class_vec[cl_vec %in% cl_x] <- x
}
combined$Cell_class <- class_vec
saveRDS(combined, file="Res_mouse_fetal_duodenum_13-18d_with_CSS_integration.rds")

p1 <- DimPlot(combined, reduction = "umap_css", label = T, group.by = "RNA_CSS_snn_res.1") + NoLegend()
p2 <- DimPlot(combined, reduction = "umap_css", group.by = "Age")
p3 <- DimPlot(combined, reduction = "umap_css", group.by = "Cell_class", label = T) + NoLegend()
p1+p2+p3

cl_order <- unlist(cell_class)
combined.expr <- getAveExpr(seu.obj=combined, feature.to.calc = "RNA_CSS_snn_res.1", specified.order = cl_order, genes=pan.cm)
combined.prop <- getExpressedProp(seu.obj = combined, feature.to.calc = "RNA_CSS_snn_res.1", specified.order = cl_order, genes=pan.cm)
getDotPlot(cl.expr.mat=combined.expr, cl.prop.mat=combined.prop, gene.reorder=FALSE, specified.order=NULL, cl.reorder.by.hc=FALSE, 
           genes=pan.cm, colors=c("#d9d9d9", "#252525"), point.size.factor=4.5, 
           plot.name="Plot_dot_plot_cell_class_marker_RNA.pdf", 
           plot.height=20, plot.width=30, plot.margin=c(8,12,5,8), 
           plot.cex=2, max.diag=FALSE, col.space.factor=0.5, row.space.factor=0.7)

fetal <- subset(combined, cells=colnames(combined)[which(combined$Cell_class!="Doublet")])
DimPlot(fetal, reduction = "umap_css", group.by = "Cell_class")
cell_class <- c("Endothelial",  "Epithelial", "Immune", "Mesenchymal", "Neural")
# subclustering for each cell class
fetal_cell_class_list <- list()
for(x in cell_class){
  cat(paste(x, "start\n"))
  
  # preprocess query data
  que_obj <- subset(fetal, cells=colnames(fetal)[which(fetal$Cell_class==x)])
  que_obj <- FindVariableFeatures(object = que_obj, selection.method = "vst", nfeatures = 3000)
  que_obj <- ScaleData(object = que_obj, verbose = T)
  que_obj <- RunPCA(object = que_obj, features = VariableFeatures(que_obj), verbose = F, npcs = 50)
  usefulPCs <- 1:20
  que_obj <- FindNeighbors(object = que_obj, dims = usefulPCs)
  que_obj <- FindClusters(object = que_obj, resolution = 1)
  que_obj <- RunUMAP(object = que_obj, dims = usefulPCs)
  
  # CSS integration
  que_obj <- cluster_sim_spectrum(que_obj, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
  mat <- que_obj@reductions$css@cell.embeddings
  idx <- which(colSums(is.na(mat))==0)
  mat2 <- mat[,idx]
  que_obj[["css"]] <- CreateDimReducObject(embeddings = mat2, key = "CSS_", assay = DefaultAssay(que_obj))
  que_obj <- RunUMAP(que_obj, reduction = "css", dims = 1:ncol(que_obj@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
  que_obj <- FindNeighbors(object = que_obj, reduction = "css", dims = 1:ncol(que_obj@reductions$css@cell.embeddings), force.recalc = T) %>%
    FindClusters(resolution = 1)
  que_obj[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- que_obj[[paste0("RNA_snn_res.", 0.2*5)]]
  
  fetal_cell_class_list[[x]] <- que_obj
}
saveRDS(fetal_cell_class_list, file="Res_mouse_fetal_duodenum_cell_class_subcluster_res.rds")

g1 <- c("EPCAM","COL1A2","CDH5","PTPRC","ASCL1")
mouse_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% g1, "Mouse_symbol"]
col.num=length(g1)
row.num=length(fetal_cell_class_list)
png("Plot_UMAP_CSS_fetal_mouse_duodenum_cell_class_marker_of_each_cell_class.png", height=2000*row.num, width=2000*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,10,5))
for(x in names(fetal_cell_class_list)){
  seu_obj <- fetal_cell_class_list[[x]]
  coor = Embeddings(seu_obj, reduction = "umap_css")
  for(g in mouse_g1){
    plotFeature2(coor, values=seu_obj@assays$RNA@data[g,], point.order = "sorted", main=paste(g, x, sep="@"),
                 cex=4, cex.main=8, lwd=0.3, nCols = blue.cols)
  }
}
dev.off()

row.num=1
col.num=length(fetal_cell_class_list)
png("Plot_UMAP_CSS_fetal_human_duodenum_by_cell_class.png", height=2000*row.num, width=2000*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,10,5))
for(x in names(fetal_cell_class_list)){
  seu_obj <- fetal_cell_class_list[[x]]
  coor = Embeddings(seu_obj, reduction = "umap_css")
  plotFeature2(coor, values=paste0("C", seu_obj$RNA_CSS_snn_res.1), point.order = "random", main=x,
               cex=4, add.label = T, label.cex = 8, cex.main=12, lwd=0.6)
}
dev.off()

seu_obj <- fetal_cell_class_list[["Epithelial"]]
g1 <- c("CDX2","PDX1","ONECUT2", "SATB2","SOX2", "FOXJ1", "TP63", "CLDN18", "LGALS4", "SHH","LGR5", "ASCL2", "OLFM4", "MUC2", "SPINK4", "SPDEF", "MKI67","CDK1",
        "CHGA", "ISL1", "BEST4", "SPIB", "GP2", "CA7", "MUC1","DPP4", "APOA4", "FABP2", "SI","DEFA5", "LYZ", "RGS13")
mouse_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% g1, "Mouse_symbol"]
plotFeature.batch(seu.obj = seu_obj, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_duodenum_epithelium_epi_cell_type_marker_expression.png", col.num = 8)

# C17 of the epithelial data is a doublet cluster of epithelium and mesenchyme
seu_obj_sub <- subset(seu_obj, cells=colnames(seu_obj)[which(seu_obj$RNA_CSS_snn_res.1!=17)])
FeaturePlot(seu_obj, reduction = "umap_css", features = "Col1a2", order = T)+FeaturePlot(seu_obj_sub, reduction = "umap_css", features = "Col1a2", order = T)
fetal_cell_class_list[["Epithelial"]] <- seu_obj_sub
saveRDS(fetal_cell_class_list, file="Res_mouse_fetal_duodenum_cell_class_subcluster_res.rds")

# epithelium
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium")
mouse_epi <- fetal_cell_class_list[["Epithelial"]]
DimPlot(mouse_epi, reduction = "umap_css", group.by = "orig.ident")
que_obj <- mouse_epi
que_obj <- FindVariableFeatures(object = que_obj, selection.method = "vst", nfeatures = 3000)
que_obj <- ScaleData(object = que_obj, verbose = T)
que_obj <- RunPCA(object = que_obj, features = VariableFeatures(que_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20

# CSS integration
que_obj <- cluster_sim_spectrum(que_obj, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
mat <- que_obj@reductions$css@cell.embeddings
idx <- which(colSums(is.na(mat))==0)
mat2 <- mat[,idx]
que_obj[["css"]] <- CreateDimReducObject(embeddings = mat2, key = "CSS_", assay = DefaultAssay(que_obj))
que_obj <- RunUMAP(que_obj, reduction = "css", dims = 1:ncol(que_obj@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
que_obj <- FindNeighbors(object = que_obj, reduction = "css", dims = 1:ncol(que_obj@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
que_obj[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- que_obj[[paste0("RNA_snn_res.", 0.2*5)]]
DimPlot(que_obj, reduction = "umap_css")
fetal_cell_class_list[["Epithelial"]] <- que_obj
mouse_epi <- fetal_cell_class_list[["Epithelial"]]
saveRDS(mouse_epi, file="Res_fetal_mouse_duodenum_epithelium.rds")
# GATA4 and SOX2 to mark the stomach;
# PTF1A, PROX1, NKX6-1 and ONECUT1 to mark pancreas; 
# SOX2 and NKX2-1 to mark lung; 
# SOX2 and depletion of GATA4 to mark esophagus

mouse_g1 <- c("Gata4", "Gata6", "Sox2", "Ptf1a", "Prox1", "Nkx6-1", "Onecut1","Nkx2-1", "Onecut2", "Pdx1", "Cdx2", "Satb2")
length(mouse_g1)
plotFeature.batch(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3, col.num = 6,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_duodenum_epithelium_regional_marker_expression.png")

row.num=1
col.num=5
coor = Embeddings(mouse_epi, reduction = "umap_css")
png("Plot_UMAP_CSS_fetal_human_duodenum_epithelium.png", height=2000*row.num, width=2000*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,10,5))
plotFeature2(coor, values=paste0("C", mouse_epi$RNA_CSS_snn_res.1), point.order = "random", main="Cluster",
            cex=4, add.label = T, label.cex = 8, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$Age, point.order = "random", main="Age", legend.pos = "topleft",
             cex=4, add.legend = T, legend.cex = 8, cex.main=12, lwd=0.6)
plotFeature2(coor, values=log10(mouse_epi$nCount_RNA), point.order = "random", main="nCount",
             cex=4, cex.main=12, lwd=0.6, nCols = blue.cols)
plotFeature2(coor, values=log10(mouse_epi$nFeature_RNA), point.order = "random", main="nFeature",
             cex=4, cex.main=12, lwd=0.6, nCols = blue.cols)
plotFeature2(coor, values=mouse_epi$percent.mito, point.order = "random", main="percent.mito",
             cex=4, cex.main=12, lwd=0.6, nCols = blue.cols)
dev.off()

de_res <- wilcoxauc(X=mouse_epi, group_by = "RNA_CSS_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[which(de_res$auc>0.6 & de_res$padj<0.01 & de_res$pct_diff>20 & de_res$pct_in>20 & de_res$logFC>0.1),]
saveRDS(sig_res, file="Res_fetal_mouse_duodenum_epithelial_cluster_markers.rds")
saveRDS(mouse_epi, file="Res_fetal_mouse_duodenum_epithelium.rds")

# canonical markers identified in human data
g1 <- c("CDX2","PDX1","ONECUT2", "SATB2","SOX2", "FOXJ1", "TP63", "CLDN18", "LGALS4", "SHH","LGR5", "ASCL2", "OLFM4", "MUC2", "SPINK4", "SPDEF", "MKI67","CDK1",
        "CHGA", "ISL1", "BEST4", "SPIB", "GP2", "CA7", "MUC1","DPP4", "APOA4", "FABP2", "SI","DEFA5", "LYZ", "RGS13")
mouse_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% g1, "Mouse_symbol"]
# Aviv_2017_Nature_adult_mouse_epi_cell_type_marker
mouse_g1 <- c("Mep1a","Fgf15","Clec2h","Lct","Cbr1","Ephx2","Tff3","Manf","Ccl9","Lrmp", "Dclk1","Cd24a","Lyz1","Mptx2","Defa-rs1","Chgb","Cpe","Gfra3")
# Aviv_2017_Nature_adult_mouse_scToEnt_regional_marker
mouse_g1 <- c("Gkn3", "Lgr5", "Bex1", "Sox4", "Foxm1", "Mxd3", "Batf2","Creb3l3","Gata4","Nr1i3","Osr2","Jund","Nr1h4")
# Aviv_2017_Nature_adult_mouse_scToEnt_experimentally_dissected_regional_marker
plot_name <- "Plot_UMAP_CSS_fetal_Aviv_2017_Nature_adult_mouse_scToEnt_experimentally_dissected_regional_marker_expression_in_fetal_mouse_duodenum_epithelium.png"
mouse_g1 <- c("Fabp1", "Lct", "Gstm3","Adh6a","Creb3l3","Fabp6","Mep1a","Muc3","Clec2h","Fgf15")
# cell class
mouse_g1 <- c("Epcam","Col1a2","Cdh5","Ptprc","Ascl1")
plot_name <- "Plot_UMAP_CSS_cell_class_markers_in_fetal_mouse_duodenum_epithelium.png"
plotFeature.batch(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3,
                  plot.name = plot_name, col.num = 5)

human_epi <- readRDS("~/Work/Endoderm/used_seurat_objects/Res_fetal_duodenum_and_early_intestine_epi.rds")
# C19 is a potential doublet cluster of enterocytes and goblet cells
# enterocyte markers: Lct, Apoa4, Fabp1
# goblet cell markers: Muc2, Spink4, Ccl9, Tff3 
# C21 is potential pancreatic cluster which coexpresses Ptf1a, Prox1, Nkx6-1, Onecut1
# small subset of C17 cells are also potential pancreatic cells
### pancreatic cells
ct_vec <- mouse_epi$RNA_CSS_snn_res.1
idx <- which(ct_vec==17)
scaled_expr <- scale(t(as.matrix(mouse_epi@assays$RNA@data[c("Nkx6-1", "Onecut1", "Cdx2","Lgals4"), idx])))
scaled_expr[,4] <- -scaled_expr[,4]
scaled_expr[,3] <- -scaled_expr[,3]
score_vec <- rowSums(scaled_expr)
pan_score_z <- scale(score_vec <- rowSums(scaled_expr))[,1]
plot(seq(length(score_vec)), sort(score_vec, decreasing = T))
cutoff <- qnorm(0.90)
pan_cells <- names(pan_score_z)[pan_score_z>cutoff]
abline(v=length(pan_cells))
plot(Embeddings(mouse_epi, reduction = "umap_css"), pch=16, col="#d0d0d0")
points(Embeddings(mouse_epi, reduction = "umap_css")[pan_cells,], pch=16)
length(pan_cells)

# clean up those cells
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells")
cells_to_remove <- c(pan_cells, colnames(mouse_epi)[which(mouse_epi$RNA_CSS_snn_res.1 %in% c(19, 21))])
cells <- setdiff(colnames(mouse_epi), cells_to_remove)
que_obj <- subset(mouse_epi, cells=cells)
DimPlot(que_obj, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", label=T)
mouse_g1 <- c("Gata4", "Gata6", "Sox2", "Ptf1a", "Prox1", "Nkx6-1", "Onecut1","Onecut2", "Pdx1", "Lgals4","Cdx2", "Satb2")
length(mouse_g1)
plotFeature.batch(seu.obj = que_obj, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=2, col.num = 6,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_duodenum_epithelium_regional_marker_expression.png")

que_obj <- FindVariableFeatures(object = que_obj, selection.method = "vst", nfeatures = 3000)
que_obj <- ScaleData(object = que_obj, verbose = T)
que_obj <- RunPCA(object = que_obj, features = VariableFeatures(que_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20

# CSS integration
que_obj <- cluster_sim_spectrum(que_obj, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
que_obj <- RunUMAP(que_obj, reduction = "css", dims = 1:ncol(que_obj@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
que_obj <- FindNeighbors(object = que_obj, reduction = "css", dims = 1:ncol(que_obj@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
que_obj[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- que_obj[[paste0("RNA_snn_res.", 0.2*5)]]
DimPlot(que_obj, reduction = "umap_css", label=T)
fetal_cell_class_list[["Epithelial"]] <- que_obj
mouse_epi <- fetal_cell_class_list[["Epithelial"]]
saveRDS(mouse_epi, file="Res_fetal_mouse_duodenum_epithelium.rds")

# regional_markers
mouse_g1 <- c("Gata4", "Gata6", "Sox2", "Ptf1a", "Prox1", "Nkx6-1", "Onecut1","Nkx2-1", "Onecut2", "Pdx1", "Cdx2", "Satb2")
length(mouse_g1)
plotFeature.batch(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=2, col.num = 6,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_duodenum_epithelium_regional_marker_expression.png")
mouse_g1 <- c("Gata4","Onecut2", "Pdx1", "Osr2","Satb2","Cdx2")
length(mouse_g1)
plotFeature.batch(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3, col.num = 6,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_duodenum_epithelium_subset_regional_marker_expression.png")

# cell type markers
mouse_g1 <- c("Lgr5", "Ascl2", "Olfm4",# stem cell
              "Mki67", "Cdk1", # cell cycle
              "Alpi", "Apoa1", "Apoa4", "Fabp1", "Fabp2", "Si", "Dpp4",# enterocyte
              "Muc2", "Clca3", "Tff3", "Agr2", "Spink4",# goblet
              "Lyz1", "Defa17", "Defa22", "Defa24", "Mptx2", "Defa5", # paneth cell
              "Chga", "Chgb", "Tac1", "Tph1", "Neurog3", # EEC
              "Dclk1", "Trpm5", "Gfi1b", "Il25", "Rgs13" # Tuft
              )
length(mouse_g1)
plotFeature.batch(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3, col.num = 7,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_SI_epithelium_cell_type_marker_expression.png")
mouse_g1 <- c("Shh", "Lgr5", "Ascl2", "Olfm4",# stem cell
              "Mki67", # cell cycle
              "Alpi", "Apoa4", "Fabp2", # enterocyte
              "Muc2", "Spink4",# goblet
              "Lyz1", "Defa22", "Mptx2", # paneth cell
              "Chga", "Chgb", "Neurog3", # EEC
              "Dclk1", "Rgs13" # Tuft
)
length(mouse_g1)
plotFeature.batch(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3, col.num = 6,
                  plot.name = "Plot_UMAP_CSS_fetal_mouse_SI_epithelium_cell_type_marker_expression.png")

row.num=1
col.num=5
coor = Embeddings(mouse_epi, reduction = "umap_css")
png("Plot_UMAP_CSS_fetal_mouse_SI_epithelium.png", height=2000*row.num, width=2000*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,10,5))
plotFeature2(coor, values=paste0("C", mouse_epi$RNA_CSS_snn_res.1), point.order = "random", main="Cluster",
             cex=4, add.label = T, label.cex = 8, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$Age, point.order = "random", main="Age", legend.pos = "topleft",
             cex=4, add.legend = T, legend.cex = 8, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$nCount_RNA, point.order = "random", main="nCount",
             cex=4, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$nFeature_RNA, point.order = "random", main="nFeature",
             cex=4, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$percent.mito, point.order = "random", main="percent.mito",
             cex=4, cex.main=12, lwd=0.6)
dev.off()

cl_vec <- mouse_epi$RNA_CSS_snn_res.1
ct_vec <- rep("G2M-phase ISC", length(cl_vec))
ct_anno <- list(
  "EEC" = 18,
  "Goblet" = 13,
  "G2M-phase goblet" = 17,
  "Lgr5+ / Olfm4+ ISC" = c(7, 16),
  "Enterocyte" = c(4, 20, 0, 1),
  "Enterocyte precursor" = c(15, 3),
  "Shh+ ISC" = c(6, 8, 12, 10),
  "Lgr5+/Shh+/Olfm4+ ISC" = c(9, 11, 5))
for(x in names(ct_anno)){
  cl_x <- ct_anno[[x]]
  ct_vec[which(cl_vec %in% cl_x)] <- x
}
mouse_epi$Cell_type <- ct_vec

row.num=3
col.num=1
png("Plot_UMAP_CSS_fetal_mouse_SI_epithelial_cell_type.png", height=2000*row.num, width=2000*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,10,20), xpd=TRUE)
plotFeature2(coor, values=paste0("C", mouse_epi$RNA_CSS_snn_res.1), point.order = "random", main="Cluster",
             cex=4, add.label = T, label.cex = 8, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$Cell_type, point.order = "random", main="Cell type",
             cex=4, add.label = T, label.cex = 8, cex.main=12, lwd=0.6)
plotFeature2(coor, values=mouse_epi$Age, point.order = "random", main="Age", legend.pos = "topleft",
             cex=4, add.legend = T, legend.cex = 8, cex.main=12, lwd=0.6)

dev.off()

# infer regional identity in fetal mouse epithelium data
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/infer_intestine_region_identity")
cells <- colnames(mouse_epi)[which(!mouse_epi$Cell_type %in% c("EEC", "G2M-phase goblet", "Goblet"))]
mouse_epi_sub <- subset(mouse_epi, cells=cells)

que_mat <- t(as.matrix(mouse_epi@assays$RNA@data[VariableFeatures(mouse_epi), cells]))
idx <- which(colSums(que_mat)>0)
que_mat <- que_mat[,idx]

## proximal SI score
prox_ref_pattern <- mouse_epi@assays$RNA@data["Onecut2", cells]
prox_cor_vec <- cor(que_mat, prox_ref_pattern)
prox_cor_vec <- prox_cor_vec[,1]
idx_gene <- order(prox_cor_vec, decreasing = T)[1:20]
prox_feature_genes <- names(prox_cor_vec[idx_gene])
plotFeature(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = prox_feature_genes, 
            plot.name = "Plot_UMAP_CSS_Onecut2_correlating_genes.png",
            nCols = blue.cols)

prox_score <- colMeans(mouse_epi@assays$RNA@data[prox_feature_genes,cells])


## distal SI score
dist_ref_pattern <- mouse_epi@assays$RNA@data["Osr2", cells]
dist_cor_vec <- cor(que_mat, dist_ref_pattern)
dist_cor_vec <- dist_cor_vec[,1]
idx_gene <- order(dist_cor_vec, decreasing = T)[1:20]
dist_feature_genes <- names(dist_cor_vec[idx_gene])
plotFeature(seu.obj = mouse_epi, dr="umap_css", genes.to.plot = dist_feature_genes, 
            plot.name = "Plot_UMAP_CSS_Osr2_correlating_genes.png",
            nCols = blue.cols)
distal_score <- colMeans(mouse_epi@assays$RNA@data[dist_feature_genes,cells])

mouse_epi_sub <- AddModuleScore(
  object = mouse_epi_sub,
  features = list(
    "Proximal_SI_features" = prox_feature_genes,
    "Distal_SI_features" = dist_feature_genes
    ),
  name = c("Proximal_SI_features","Distal_SI_features")
)

mouse_epi_sub$SI_region_score_diff <- mouse_epi_sub$Proximal_SI_features1 - mouse_epi_sub$Distal_SI_features2


plotFeature2(Embeddings(mouse_epi_sub, reduction = "umap_css"),
             values = mouse_epi_sub$SI_region_score_diff, point.order = "sorted")
prox_cells <- names(which(mouse_epi_sub$SI_region_score_diff>0))
dist_cells <- names(which(mouse_epi_sub$SI_region_score_diff<0))

mouse_epi_sub$Inferred_SI_region <- "Proximal"
mouse_epi_sub@meta.data[dist_cells, "Inferred_SI_region"] <- "Distal"
DimPlot(mouse_epi_sub, reduction = "umap_css", group.by = "Inferred_SI_region")

ent_cells <- colnames(mouse_epi_sub)[which(mouse_epi_sub$Cell_type == "Enterocyte")]
prog_cells <- c(colnames(mouse_epi_sub)[grep("ISC", mouse_epi_sub$Cell_type)], 
                colnames(mouse_epi_sub)[which(mouse_epi_sub$Cell_type == "Enterocyte precursor")])

ent_score <- scale(mouse_epi_sub@meta.data[ent_cells, "SI_region_score_diff"])
prox_ent_cells <- ent_cells[which(ent_score>0)]
dist_ent_cells <- ent_cells[which(ent_score<0)]

prog_score <- scale(mouse_epi_sub@meta.data[prog_cells, "SI_region_score_diff"])
prox_prog_cells <- prog_cells[which(prog_score>0)]
dist_prog_cells <- prog_cells[which(prog_score<0)]

plot(Embeddings(mouse_epi_sub, reduction = "umap_css"), pch=16, col="#d0d0d0")
points(Embeddings(mouse_epi_sub, reduction = "umap_css")[prox_ent_cells,], pch=16, col="light blue")
points(Embeddings(mouse_epi_sub, reduction = "umap_css")[dist_ent_cells,], pch=16, col="dark blue")
points(Embeddings(mouse_epi_sub, reduction = "umap_css")[prox_prog_cells,], pch=16, col="light green")
points(Embeddings(mouse_epi_sub, reduction = "umap_css")[dist_prog_cells,], pch=16, col="dark green")

# analyze each individual sample separately
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/individual_sample")
mouse_epi_by_age <- SplitObject(mouse_epi, split.by = "Age")
names(mouse_epi_by_age)
for(x in names(mouse_epi_by_age)){
  cat(paste(x,"start\n"))
  seu_obj <- mouse_epi_by_age[[x]]
  seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
  seu_obj <- ScaleData(object = seu_obj, verbose = T)
  seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
  usefulPCs <- 1:20
  seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
  seu_obj <- FindClusters(object = seu_obj, resolution = 1)
  seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
  mouse_epi_by_age[[x]] <- seu_obj
}
saveRDS(mouse_epi_by_age, file="Res_mouse_fetal_duodenum_by_age.rds")

# plot regional markers in each time point separately
for(x in names(mouse_epi_by_age)){
  seu_obj <- mouse_epi_by_age[[x]]
  plot_name <- paste0("Plot_UMAP_mouse_epi_regional_markers_at_",x,".png")
  plotFeature(seu.obj=seu_obj, dr="umap", genes.to.plot = c("Pdx1", "Osr2", "Gata4", "Onecut2", "Satb2","Cdx2", "Hoxa10"),
              nCols = blue.cols, plot.name = plot_name, col.num = 7)
  
}


# day 13.5
age <- "13.5d"
seu_obj <- mouse_epi_by_age[[age]]
DimPlot(seu_obj)
seu_obj <- AddModuleScore(
  object = seu_obj,
  features = list(
    prox_feature_genes,
    dist_feature_genes
  ),
  name = c("Proximal_SI_features","Distal_SI_features")
)
seu_obj$Prox_dist_SI_score_diff <- seu_obj$Proximal_SI_features1 - seu_obj$Distal_SI_features2
seu_obj$Inferred_SI_region_identity <- "Prox_SI"
seu_obj$Inferred_SI_region_identity[which(seu_obj$Prox_dist_SI_score_diff<0)] <- "Dist_SI"
DimPlot(seu_obj, group.by = "Inferred_SI_region_identity")
mouse_epi_by_age[[age]] <- seu_obj
saveRDS(seu_obj, file=paste0("Res_fetal_mouse_",age,"_with_inferred_SI_region_identity.rds"))

# day 14.5
age <- "14.5d"
seu_obj <- mouse_epi_by_age[[age]]
DimPlot(seu_obj)

cell_type_group <- list(
  "ISC"=setdiff(unique(seu_obj$RNA_snn_res.1),11),
  "Secretory"=11
)
cl_vec <- seu_obj$RNA_snn_res.1
group_vec <- rep(NA, length(cl_vec))
for(x in names(cell_type_group)){
  cl_x <- cell_type_group[[x]]
  group_vec[which(cl_vec%in%cl_x)] <- x
}
seu_obj$Cell_type_group <- group_vec
DimPlot(seu_obj, group.by = "Cell_type_group", label=T)

score_vec <- rep(NA, ncol(seu_obj))
names(score_vec) <- colnames(seu_obj)

for(x in names(cell_type_group)){
  seu_obj_sub <- subset(seu_obj, subset = RNA_snn_res.1%in%cell_type_group[[x]])
  seu_obj_sub <- AddModuleScore(
    object = seu_obj_sub,
    features = list(
      prox_feature_genes,
      dist_feature_genes
    ),
    name = c("Proximal_SI_features","Distal_SI_features")
  )
  score_vec[colnames(seu_obj_sub)] <- seu_obj_sub$Proximal_SI_features1 - seu_obj_sub$Distal_SI_features2
}
seu_obj$Prox_dist_SI_score_diff <- score_vec

region_vec <- rep(NA, ncol(seu_obj))
names(region_vec) <- colnames(seu_obj)
for(x in names(cell_type_group)){
  prox_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff>0 & seu_obj$Cell_type_group==x)]
  dist_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff<0 & seu_obj$Cell_type_group==x)]
  region_vec[prox_cells] <- "Prox_SI"
  region_vec[dist_cells] <- "Dist_SI"
}
seu_obj$Inferred_SI_region_identity <- region_vec
DimPlot(seu_obj, group.by = "Inferred_SI_region_identity")
mouse_epi_by_age[[age]] <- seu_obj
saveRDS(seu_obj, file=paste0("Res_fetal_mouse_",age,"_with_inferred_SI_region_identity.rds"))

# day 15.5
age <- "15.5d"
seu_obj <- mouse_epi_by_age[[age]]
DimPlot(seu_obj, label=T)

cell_type_group <- list(
  "ISC"=setdiff(unique(seu_obj$RNA_snn_res.1),c(6,10)),
  "Secretory"=c(6,10)
)
cl_vec <- seu_obj$RNA_snn_res.1
group_vec <- rep(NA, length(cl_vec))
for(x in names(cell_type_group)){
  cl_x <- cell_type_group[[x]]
  group_vec[which(cl_vec%in%cl_x)] <- x
}
seu_obj$Cell_type_group <- group_vec
DimPlot(seu_obj, group.by = "Cell_type_group", label=T)

score_vec <- rep(NA, ncol(seu_obj))
names(score_vec) <- colnames(seu_obj)

for(x in names(cell_type_group)){
  seu_obj_sub <- subset(seu_obj, subset = RNA_snn_res.1%in%cell_type_group[[x]])
  seu_obj_sub <- AddModuleScore(
    object = seu_obj_sub,
    features = list(
      prox_feature_genes,
      dist_feature_genes
    ),
    name = c("Proximal_SI_features","Distal_SI_features")
  )
  score_vec[colnames(seu_obj_sub)] <- seu_obj_sub$Proximal_SI_features1 - seu_obj_sub$Distal_SI_features2
}
seu_obj$Prox_dist_SI_score_diff <- score_vec

region_vec <- rep(NA, ncol(seu_obj))
names(region_vec) <- colnames(seu_obj)
for(x in names(cell_type_group)){
  prox_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff>0 & seu_obj$Cell_type_group==x)]
  dist_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff<0 & seu_obj$Cell_type_group==x)]
  region_vec[prox_cells] <- "Prox_SI"
  region_vec[dist_cells] <- "Dist_SI"
}
seu_obj$Inferred_SI_region_identity <- region_vec
DimPlot(seu_obj, group.by = "Inferred_SI_region_identity")
mouse_epi_by_age[[age]] <- seu_obj
saveRDS(seu_obj, file=paste0("Res_fetal_mouse_",age,"_with_inferred_SI_region_identity.rds"))

# day 16
age <- "16d"
seu_obj <- mouse_epi_by_age[[age]]
DimPlot(seu_obj, label=T)

cell_type_group <- list(
  "ISC_EP"=setdiff(unique(seu_obj$RNA_snn_res.1),c(4,7)),
  "Secretory"=c(4,7)
)
cl_vec <- seu_obj$RNA_snn_res.1
group_vec <- rep(NA, length(cl_vec))
for(x in names(cell_type_group)){
  cl_x <- cell_type_group[[x]]
  group_vec[which(cl_vec%in%cl_x)] <- x
}
seu_obj$Cell_type_group <- group_vec
DimPlot(seu_obj, group.by = "Cell_type_group", label=T)

score_vec <- rep(NA, ncol(seu_obj))
names(score_vec) <- colnames(seu_obj)

for(x in names(cell_type_group)){
  seu_obj_sub <- subset(seu_obj, subset = RNA_snn_res.1%in%cell_type_group[[x]])
  seu_obj_sub <- AddModuleScore(
    object = seu_obj_sub,
    features = list(
      prox_feature_genes,
      dist_feature_genes
    ),
    name = c("Proximal_SI_features","Distal_SI_features")
  )
  score_vec[colnames(seu_obj_sub)] <- seu_obj_sub$Proximal_SI_features1 - seu_obj_sub$Distal_SI_features2
}
seu_obj$Prox_dist_SI_score_diff <- score_vec

region_vec <- rep(NA, ncol(seu_obj))
names(region_vec) <- colnames(seu_obj)
for(x in names(cell_type_group)){
  prox_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff>0 & seu_obj$Cell_type_group==x)]
  dist_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff<0 & seu_obj$Cell_type_group==x)]
  region_vec[prox_cells] <- "Prox_SI"
  region_vec[dist_cells] <- "Dist_SI"
}
seu_obj$Inferred_SI_region_identity <- region_vec
DimPlot(seu_obj, group.by = "Inferred_SI_region_identity")
mouse_epi_by_age[[age]] <- seu_obj
saveRDS(seu_obj, file=paste0("Res_fetal_mouse_",age,"_with_inferred_SI_region_identity.rds"))

# day 17.5
seu_obj <- mouse_epi_by_age[["17.5d"]]
cell_type_group <- list(
  "ISC"=c(2,7,8,11),
  "E_and_EP"=c(0,1,3,4,5,6,9,14),
  "Secretory"=c(10,12,13)
)
cl_vec <- seu_obj$RNA_snn_res.1
group_vec <- rep(NA, length(cl_vec))
for(x in names(cell_type_group)){
  cl_x <- cell_type_group[[x]]
  group_vec[which(cl_vec%in%cl_x)] <- x
}
seu_obj$Cell_type_group <- group_vec
DimPlot(seu_obj, group.by = "Cell_type_group", label=T)

score_vec <- rep(NA, ncol(seu_obj))
names(score_vec) <- colnames(seu_obj)

for(x in names(cell_type_group)){
  seu_obj_sub <- subset(seu_obj, subset = RNA_snn_res.1%in%cell_type_group[[x]])
  seu_obj_sub <- AddModuleScore(
    object = seu_obj_sub,
    features = list(
      prox_feature_genes,
      dist_feature_genes
    ),
    name = c("Proximal_SI_features","Distal_SI_features")
  )
  score_vec[colnames(seu_obj_sub)] <- seu_obj_sub$Proximal_SI_features1 - seu_obj_sub$Distal_SI_features2
}
seu_obj$Prox_dist_SI_score_diff <- score_vec

region_vec <- rep(NA, ncol(seu_obj))
names(region_vec) <- colnames(seu_obj)
for(x in names(cell_type_group)){
  prox_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff>0 & seu_obj$Cell_type_group==x)]
  dist_cells <- colnames(seu_obj)[which(seu_obj$Prox_dist_SI_score_diff<0 & seu_obj$Cell_type_group==x)]
  region_vec[prox_cells] <- "Prox_SI"
  region_vec[dist_cells] <- "Dist_SI"
}
seu_obj$Inferred_SI_region_identity <- region_vec
DimPlot(seu_obj, group.by = "Inferred_SI_region_identity")
mouse_epi_by_age[["17.5d"]] <- seu_obj
saveRDS(seu_obj, file="Res_fetal_mouse_17.5d_with_inferred_SI_region_identity.rds")

seu_obj$Cell_type_per_region <- paste(seu_obj$Cell_type_group, seu_obj$Inferred_SI_region_identity, sep="@")
SCpubr::do_ViolinPlot(seu_obj, features = "Prox_dist_SI_score_diff", group.by = "Cell_type_per_region")



region_vec <- rep(NA, ncol(mouse_epi))
DimPlot(mouse_epi, reduction = "umap_css", group.by = "Cell_type")
names(region_vec) <- colnames(mouse_epi)
for(x in names(mouse_epi_by_age)){
  seu_obj <- mouse_epi_by_age[[x]]
  region_vec[colnames(seu_obj)] <- seu_obj$Inferred_SI_region_identity
}
mouse_epi$Inferred_SI_region_identity <- region_vec
mouse_epi$Inferred_SI_region_identity[which(mouse_epi$Cell_type %in% c("EEC", "G2M-phase goblet", "Goblet"))] <- "Undetermined"
DimPlot(mouse_epi, group.by = "Inferred_SI_region_identity", reduction = "umap_css")
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells")
saveRDS(mouse_epi, file="Res_fetal_mouse_SI_epithelium_with_inferred_SI_region_identity.rds")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/individual_sample")
# plot metadata distribution in each time point separately
row_num=3
col_num=length(mouse_epi_by_age)
pos_vec <- c("topright", "topleft", "bottomleft", "topright", "bottom")
names(pos_vec) <- sort(names(mouse_epi_by_age))
png("Plot_UMAP_mouse_epi_cell_type_by_age.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,15,5))
for(x in sort(names(mouse_epi_by_age))){
  seu_obj <- mouse_epi_by_age[[x]]
  plotFeature2(Embeddings(seu_obj, reduction = "umap"),
               values = seu_obj$Cell_type,
               point.order = "random",
               main=x,
               add.legend = T,
               cex=6,
               cex.main=12, 
               legend.cex = 8,
               legend.pos = pos_vec[x],
               lwd=1)
}
for(x in sort(names(mouse_epi_by_age))){
  seu_obj <- mouse_epi_by_age[[x]]
  plotFeature2(Embeddings(seu_obj, reduction = "umap"),
               values = paste0("C", seu_obj$RNA_snn_res.1),
               point.order = "random",
               main=x,
               add.label = T,
               cex=6,
               cex.main=12, 
               label.cex = 8,
               lwd=1)
}
for(x in sort(names(mouse_epi_by_age))){
  seu_obj <- mouse_epi_by_age[[x]]
  plotFeature2(Embeddings(seu_obj, reduction = "umap"),
               values = seu_obj$Inferred_SI_region_identity,
               point.order = "random",
               main=x,
               add.legend = T,
               cex=6,
               cex.main=12, 
               legend.cex = 8,
               lwd=1)
}
dev.off()

row_num=1
col_num=1
png("Plot_UMAP_CSS_mouse_epi_inferred_SI_region.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,15,5))
plotFeature2(Embeddings(mouse_epi, reduction = "umap_css"),
             values = mouse_epi$Inferred_SI_region_identity,
             point.order = "random",
             add.legend = T,
             cex=3,
             cex.main=12, 
             legend.cex = 8,
             legend.pos = "topleft",
             lwd=1)
dev.off()

# calculate the distance between tHIO intestinal epithelial cells and fetal intestine
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_human_multi_intestine_region_hg19/with_colon/compare_with_tHIO")
## load fetal atlas highly variable genes
fetal.hvg <- VariableFeatures(fetal_int)
## load fetal atlas meta.data after decompression, which contains CSS-based UMAP embedding of the fetal atlas and organ identity
fetal.meta <- fetal_int@meta.data
ref.idx <- fetal.meta$Corrected_tissue
## load fetal Cluster Similarity Spectrum (CSS) model
css.model <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_human_multi_intestine_region_hg19/with_colon/Res_fetal_human_multi_intestinal_region_CSS_model.rds")

# apply it to mouse data
gene_pairs <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% fetal.hvg & human_mouse_gene_link[,"Mouse_symbol"] %in% rownames(mouse_epi), c("Human_symbol", "Mouse_symbol")]

## use the fetal highly variable genes that are also detected in tIO data to calculate similarity between tIO cells and fetal clusters per sample
## because here the fetal cluster average profile is reference, and the tIO cells are query, the obtained similarity spectrum is called reference similarity spectrum (RSS)
shared.genes <- gene_pairs[,"Human_symbol"]
length(shared.genes)
que.data <- as.matrix(mouse_epi@assays$RNA@data[gene_pairs[,"Mouse_symbol"],])
rownames(que.data) <- shared.genes
rss.list <- lapply(seq(length(css.model$model$profiles)), function(i){
  ref <- css.model$model$profiles[[i]][shared.genes,]
  cor.mat <- cor(que.data, ref, method = "spearman")
  cor.z <- t(scale(t(cor.mat)))
  return(cor.z)
})
rss.mat <- do.call('cbind', rss.list)
mouse_epi[["rss"]] <- CreateDimReducObject(embeddings = rss.mat, key="RSS_", assay=DefaultAssay(mouse_epi))

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
mouse_epi@meta.data$Mapped_fetal_intestine_region_after_correction <- trans.id
DimPlot(mouse_epi, reduction = "umap_css", group.by = "Mapped_fetal_intestine_region_after_correction")
saveRDS(mouse_epi, file="Res_mouse_epi_with_CSS_integration_and_fetal_human_intestine_region_identity_projection.rds")

tIO_rna <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/C7-tCIO_H9-tHIO/Res_C7-tCIO_H9-tHIO_no_integration.rds")
plot_name <- "Plot_UMAP_tIO_non-integrated_intestine_epi_regional_markers.png"

tHIO_hg19 <- readRDS("~/Work/Endoderm/used_seurat_objects/Res_H9_ENR_tHIO.rds")
plot_name <- "Plot_UMAP_tHIO_hg19_intestine_epi_regional_markers.png"

tHIO_fetal_hg38 <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/ENR-tHIO_and_fetal_duo/epithelium/remove_mesenchymal_doublet/Res_tHIO_and_fetal_epi_with_Paneth_cells.rds")
plot_name <- "Plot_UMAP_tHIO_and_fetal_hg38_intestine_epi_regional_markers.png"

plotFeature(seu.obj=tHIO_fetal_hg38, dr="umap", genes.to.plot = c("PDX1", "OSR2", "GATA4", "ONECUT2", "SATB2","CDX2", "HOXA10"),
            nCols = blue.cols, plot.name = plot_name, col.num = 4)

plot_name <- "Plot_UMAP_mouse_epi_regional_markers.png"
plotFeature(seu.obj=mouse_epi, dr="umap", genes.to.plot = c("Pdx1", "Osr2", "Gata4", "Onecut2", "Satb2","Cdx2", "Hoxa10"),
            nCols = blue.cols, plot.name = plot_name, col.num = 4)

# check Shh, Lgr5, Olfm4 expression in each cluster across time points
# check SI regional marker expression in each cluster across time points

# focus on proximal SI stem cells-to-enterocyte differentiation
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt")
prox <- subset(mouse_epi, subset = Inferred_SI_region_identity=="Prox_SI")
DimPlot(prox, reduction = "umap_css")
prox <- FindVariableFeatures(object = prox, selection.method = "vst", nfeatures = 3000)
prox <- ScaleData(object = prox, verbose = T)
prox <- RunPCA(object = prox, features = VariableFeatures(prox), verbose = F, npcs = 50)
usefulPCs <- 1:20
prox <- FindNeighbors(object = prox, dims = usefulPCs)
prox <- FindClusters(object = prox, resolution = 1)
prox <- RunUMAP(object = prox, dims = usefulPCs)
saveRDS(prox, file="Res_fetal_mouse_proximal_SI_stemToEnterocyte_13-18d_no_batch_correction.rds")

# CSS integration
prox <- cluster_sim_spectrum(prox, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
prox <- RunUMAP(prox, reduction = "css", dims = 1:ncol(prox@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
prox <- FindNeighbors(object = prox, reduction = "css", dims = 1:ncol(prox@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
prox[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- prox[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(prox, file="Res_fetal_mouse_proximal_SI_stemToEnterocyte_13-18d_with_CSS_integration.rds")
DimPlot(prox, reduction = "umap", group.by = "Age")+DimPlot(prox, reduction = "umap_css", group.by = "Age")
plotFeature(seu.obj = prox, dr="umap", genes.to.plot = c("Pdx1", "Onecut2", "Osr2", "Satb2", "Cdx2"), col.num = 5,
            nCols = blue.cols, plot.name = "Plot_UMAP_fetal_mouse_proximal_SI_intestine_regional_markers.png")
plotFeature(seu.obj = human_epi, dr="umap_css", genes.to.plot = toupper(c("Pdx1", "Onecut2", "Osr2", "Satb2", "Cdx2")), col.num = 5,
            nCols = blue.cols, plot.name = "Plot_UMAP_fetal_human_duodenum_intestine_regional_markers.png")

DimPlot(human_epi, reduction = "umap_css")

# Harmony integration
library(harmony)
s2 <- prox
s2 <- RunHarmony(s2, "orig.ident")
s2 <- RunUMAP(s2, reduction = "harmony", dims = usefulPCs, reduction.key = "UMAPHARMONY_", reduction.name = "umap_harmony")
s2 <- FindNeighbors(object = s2, reduction = "harmony", dims = 1:20, force.recalc = T) %>% FindClusters(resolution = 0.8)
DimPlot(s2, reduction = "umap_harmony", group.by = "Cell_type")

# diffusion map for each individual sample
library(destiny)
prox_by_age <- SplitObject(prox, split.by = "Age")
for(x in names(prox_by_age)){
  seu_obj <- prox_by_age[[x]]
  seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
  seu_obj <- ScaleData(object = seu_obj, verbose = T)
  seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
  input.mat <- Embeddings(seu_obj, reduction = "pca")[,1:20]
  dm <- DiffusionMap(data=input.mat, distance="euclidean", k=20)
  dcs <- dm@eigenvectors
  dc1 <- dcs[,1]
  dc2 <- dcs[,2]
  dc3 <- dcs[,3]
  seu_obj$DC1 <- dc1
  seu_obj$DC2 <- dc2
  seu_obj$DC3 <- dc3
  prox_by_age[[x]] <- seu_obj
}

pos_vec["15.5d"] <- "bottomright"
gCols <- setNames(colorRampPalette(prettyrainbow)(length(unique(prox$Cell_type))), sort(unique(prox$Cell_type)))
row_num=3
col_num=length(prox_by_age)
png("Plot_diffusionMap_fetal_mouse_proxSI_scToEnt_by_age.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,15,5))
for(x in sort(names(prox_by_age))){
  seu_obj <- prox_by_age[[x]]
  coor <- cbind(seu_obj$DC1, seu_obj$DC2)
  plotFeature2(coor,
               values = seu_obj$Cell_type,
               point.order = "random",
               main=x,
               add.legend = T,
               cex=6,
               cex.main=12, 
               legend.cex = 8,
               legend.pos = pos_vec[x],
               lwd=1,
               gCols = gCols)
}
for(x in sort(names(prox_by_age))){
  seu_obj <- prox_by_age[[x]]
  coor <- cbind(seu_obj$DC1, seu_obj$DC3)
  plotFeature2(coor,
               values = seu_obj$Cell_type,
               point.order = "random",
               main=x,
               add.legend = T,
               cex=6,
               cex.main=12, 
               legend.cex = 8,
               legend.pos = pos_vec[x],
               lwd=1,
               gCols = gCols)
}
for(x in sort(names(prox_by_age))){
  seu_obj <- prox_by_age[[x]]
  coor <- cbind(seu_obj$DC2, seu_obj$DC3)
  plotFeature2(coor,
               values = seu_obj$Cell_type,
               point.order = "random",
               main=x,
               add.legend = T,
               cex=6,
               cex.main=12, 
               legend.cex = 8,
               legend.pos = pos_vec[x],
               lwd=1,
               gCols = gCols)
}

dev.off()

# get cell type composition across time point in proximal SI and distal SI
values <- setdiff(unique(mouse_epi$Cell_type), c("EEC", "G2M-phase goblet", "Goblet"))
gCols <- setNames(colorRampPalette(prettyrainbow)(length(values)), values)
time_point <- sort(unique(mouse_epi$Age))
n_prox <- sapply(names(gCols), function(x){
  sapply(time_point, function(t){
    sum(prox$Cell_type==x & prox$Age==t)
  })
})
p_prox <- t(n_prox/rowSums(n_prox))
n_dist <- sapply(names(gCols), function(x){
  sapply(time_point, function(t){
    sum(dist$Cell_type==x & dist$Age==t)
  })
})
p_dist <- t(n_dist/rowSums(n_dist))

mat <- cbind(p_prox, p_dist)
colnames(mat) <- paste(colnames(mat), rep(c("Prox", "Dist"), each=ncol(p_prox)), sep="_")
pdf("Plot_barplot_cell_type_composition_across_time_point.pdf", height=5, width=10)
par(mfrow=c(1,2), mar=c(8,5,3,3))
barplot(mat, col = gCols, las=2, ylab="Cell type proportion")
plot(1:10,1:10,type="n", bty="n", xaxt="n",yaxt="n",xlab="",ylab="")
legend("topright", legend = names(gCols), fill = gCols)
dev.off()

coor <- cbind(dc1, dc2)
par(mfrow=c(1,3))
x <- "Age"
plotFeature2(cbind(dc1, dc2), values = prox[[x]][,1], point.order = "random", add.legend = T, legend.pos = "bottomleft")
plotFeature2(cbind(dc1, dc3), values = prox[[x]][,1], point.order = "random", add.legend = T, legend.pos = "bottomleft")
plotFeature2(cbind(dc2, dc3), values = prox[[x]][,1], point.order = "random", add.legend = T, legend.pos = "bottomleft")


# cell type
plotFeature(seu.obj = prox, dr="umap_css", genes.to.plot = c("Shh","Lgr5", "Ascl2", "Olfm4", "Mki67", "Apoa4"), col.num = 6,
            nCols = blue.cols, plot.name = "Plot_UMAP_fetal_mouse_proximal_SI_ISC_and_enterocyte_markers.png")
plotFeature(seu.obj = dist, dr="umap_css", genes.to.plot = c("Shh","Lgr5", "Ascl2", "Olfm4", "Mki67", "Apoa4"), col.num = 6,
            nCols = blue.cols, plot.name = "Plot_UMAP_fetal_mouse_distal_SI_ISC_and_enterocyte_markers.png")



# for the distal regions
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/distal_SI_scToEnt")
dist <- subset(mouse_epi, subset = Inferred_SI_region_identity=="Dist_SI")
DimPlot(dist, reduction = "umap_css")
dist <- FindVariableFeatures(object = dist, selection.method = "vst", nfeatures = 3000)
dist <- ScaleData(object = dist, verbose = T)
dist <- RunPCA(object = dist, features = VariableFeatures(dist), verbose = F, npcs = 50)
usefulPCs <- 1:20
dist <- FindNeighbors(object = dist, dims = usefulPCs)
dist <- FindClusters(object = dist, resolution = 1)
dist <- RunUMAP(object = dist, dims = usefulPCs)
saveRDS(dist, file="Res_fetal_mouse_distal_SI_stemToEnterocyte_13-18d_no_batch_correction.rds")

# CSS integration
dist <- cluster_sim_spectrum(dist, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
dist <- RunUMAP(dist, reduction = "css", dims = 1:ncol(dist@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
dist <- FindNeighbors(object = dist, reduction = "css", dims = 1:ncol(dist@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
dist[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- dist[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(dist, file="Res_fetal_mouse_distal_SI_stemToEnterocyte_13-18d_with_CSS_integration.rds")
DimPlot(dist, reduction = "umap", group.by = "Age")+DimPlot(dist, reduction = "umap_css", group.by = "Age")
plotFeature(seu.obj = dist, dr="umap", genes.to.plot = c("Pdx1", "Onecut2", "Osr2", "Satb2", "Cdx2"), col.num = 5,
            nCols = blue.cols, plot.name = "Plot_UMAP_fetal_mouse_distal_SI_intestine_regional_markers.png")

genes <- c("Lgr5", "Ascl2", "Olfm4", "Shh", "Apoa4", "Mki67", "Pdx1", "Onecut2", "Osr2")
expr_list <- lapply(genes, function(g){
  
  res <- list()
  for(ct in names(gCols)){
    for(t in time_point){
      idx <- which(mouse_epi$Cell_type==ct & mouse_epi$Age==t)
      if(length(idx)>20){
        res[[paste(ct,t,sep="@")]] <- mouse_epi@assays$RNA@data[g, idx]
      }
    }
  }
  return(res)
  
})
names(expr_list) <- genes

row_num=3
col_num=3
pdf("Plot_boxplot_marker_gene_expr_across_cell_type_and_time_point.pdf", height=5*row_num, width=7*col_num)
par(mfrow=c(row_num, col_num), mar=c(15,5,3,3))
for(g in names(expr_list)){
  boxplot(expr_list[[g]], main=g, las=2, frame=F, ylab="Normalized expression")
}
dev.off()

# merge proximal and distal stem cells to enterocytes
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/SI_scToEnt")
combined <- subset(mouse_epi, cells=colnames(mouse_epi)[which(!mouse_epi$Cell_type %in% c("EEC", "G2M-phase goblet", "Goblet"))])
DimPlot(combined, reduction = "umap_css")
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)

# CSS integration
combined <- cluster_sim_spectrum(combined, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(combined, file="Res_fetal_mouse_combined_SI_stemToEnterocyte_13-18d_with_CSS_integration.rds")
DimPlot(combined, reduction = "umap", group.by = "Age")+DimPlot(combined, reduction = "umap_css", group.by = "Age")
plotFeature(seu.obj = combined, dr="umap", genes.to.plot = c("Pdx1", "Onecut2", "Osr2", "Satb2", "Cdx2"), col.num = 5,
            nCols = blue.cols, plot.name = "Plot_UMAP_fetal_mouse_combined_SI_intestine_regional_markers.png")

# get cell type composition across time point
values <- setdiff(unique(mouse_epi$Cell_type), c("EEC", "G2M-phase goblet", "Goblet"))
gCols <- setNames(colorRampPalette(prettyrainbow)(length(values)), values)
time_point <- sort(unique(mouse_epi$Age))
n1 <- sapply(names(gCols), function(x){
  sapply(time_point, function(t){
    sum(combined$Cell_type==x & combined$Age==t)
  })
})
p1 <- t(n1/rowSums(n1))

pdf("Plot_barplot_cell_type_composition_across_time_point.pdf", height=5, width=10)
par(mfrow=c(1,2), mar=c(8,5,3,3))
barplot(p1, col = gCols, las=2, ylab="Cell type proportion")
plot(1:10,1:10,type="n", bty="n", xaxt="n",yaxt="n",xlab="",ylab="")
legend("topright", legend = names(gCols), fill = gCols)
dev.off()

# plot cell type A-P distribution bias
n1 <- sapply(values, function(ct){
  sapply(c("Prox_SI", "Dist_SI"), function(y){
    sum(combined$Cell_type==ct & combined$Inferred_SI_region_identity==y)
  })
})
p1 <- t(t(n1)/colSums(n1))
region_cols <- setNames(c("#d9d9d9", "#505050"), c("Prox_SI", "Dist_SI"))
pdf("Plot_barplot_cell_type_AP_distribution.pdf")
par(mar=c(12,5,5,5))
barplot(p1, col = region_cols[rownames(p1)], las=2, ylab="Proportion")
dev.off()

# plot time point A-P distribution bias
n1 <- sapply(time_point, function(t){
  sapply(c("Prox_SI", "Dist_SI"), function(y){
    sum(combined$Age==t & combined$Inferred_SI_region_identity==y)
  })
})
p1 <- t(t(n1)/colSums(n1))
pdf("Plot_barplot_time_point_AP_distribution.pdf")
par(mar=c(12,5,5,5))
barplot(p1, col = region_cols[rownames(p1)], las=2, ylab="Proportion")
dev.off()

# for human-mouse comparison
## build duodenum human model
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/human_duo_scToEnt_model")
human_epi <- readRDS("~/Work/Endoderm/used_seurat_objects/Res_fetal_duodenum_and_early_intestine_epi.rds")
cells <- colnames(human_epi)[which(human_epi$Cell_type %in% c("Stem_cell", "Enterocyte", "Enterocyte_precursor"))]
human_se <- subset(human_epi, cells=cells)
DimPlot(human_se, reduction = "umap_css", group.by = "Cell_type")
human_se <- FindVariableFeatures(object = human_se, selection.method = "vst", nfeatures = 3000)
human_se <- ScaleData(object = human_se, verbose = T)
human_se <- RunPCA(object = human_se, features = VariableFeatures(human_se), verbose = F, npcs = 50)
usefulPCs <- 1:20
human_se <- FindNeighbors(object = human_se, dims = usefulPCs)
human_se <- FindClusters(object = human_se, resolution = 1)
human_se <- RunUMAP(object = human_se, dims = usefulPCs)
# CSS integration
css.model <- cluster_sim_spectrum(human_se, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, 
                                  return_seuratObj = FALSE)
saveRDS(css.model, file="Res_human_fetal_duo_se_CSS_model.rds")
css <- css.model$sim2profiles
colnames(css) <- paste("CSS", colnames(css), sep="_")
human_se[["css"]] <- CreateDimReducObject(embeddings = css, key="CSS_", assay=DefaultAssay(human_se))
human_se <- RunUMAP(human_se, reduction = "css", dims = 1:ncol(human_se@reductions$css@cell.embeddings), 
                 reduction.name = "umap_css", reduction.key = "UMAPCSS_",
                 return.model = TRUE)
DimPlot(human_se, reduction = "umap_css", group.by = "Cell_type")
saveRDS(human_se, file="Res_fetal_human_duo_se_with_models.rds")


# calculate the distance between mouse proximal s2e and human duodenum s2e
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/compare_to_human_duo/CSS_model_based")
## load fetal atlas highly variable genes
fetal.hvg <- VariableFeatures(human_se)
## load fetal atlas meta.data after decompression, which contains CSS-based UMAP embedding of the fetal atlas and organ identity
fetal.meta <- human_se@meta.data

# apply it to mouse data
gene_pairs <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"] %in% fetal.hvg & human_mouse_gene_link[,"Mouse_symbol"] %in% rownames(mouse_epi), c("Human_symbol", "Mouse_symbol")]

## use the fetal highly variable genes that are also detected in tIO data to calculate similarity between tIO cells and fetal clusters per sample
## because here the fetal cluster average profile is reference, and the tIO cells are query, the obtained similarity spectrum is called reference similarity spectrum (RSS)
shared.genes <- gene_pairs[,"Human_symbol"]
length(shared.genes)
que.data <- as.matrix(prox@assays$RNA@data[gene_pairs[,"Mouse_symbol"],])
rownames(que.data) <- shared.genes
rss.list <- lapply(seq(length(css.model$model$profiles)), function(i){
  ref <- css.model$model$profiles[[i]][shared.genes,]
  cor.mat <- cor(que.data, ref, method = "spearman")
  cor.z <- t(scale(t(cor.mat)))
  return(cor.z)
})
rss.mat <- do.call('cbind', rss.list)
prox[["rss"]] <- CreateDimReducObject(embeddings = rss.mat, key="RSS_", assay=DefaultAssay(prox))

## calculate the Euclidean distance between tIO cells and fetal cells on the space of similarities to fetal clusters per sample, instead of on the space of expression. correlation is potentially more robust to batch variations than expression 
## after the distance calculation, get the 20-nearest fetal cells for each tIO cell
## assign the majority of tissue identity of the 20-nearest cells to be the inferred tissue identity of tHIO cells 
css <- css.model$sim2profiles 
knn <- RANN::nn2(css, rss.mat, k = 20)$nn.idx
saveRDS(knn, file="Res_20nn_fetal_human_duo_stem_cell_to_enterocyte_of_mouse_proximal_s2e.rds")
## map the cell type identity of fetal cells to tHIO cells
ref.idx <- fetal.meta$Cell_type
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
trans.id <- apply(nn.idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
prox@meta.data$Mapped_human_fetal_duo_cell_type <- trans.id
## map the age of fetal cells to tHIO cells
ref.idx <- fetal.meta$Age
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
trans.id <- apply(nn.idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
prox@meta.data$Mapped_human_fetal_duo_age <- trans.id

DimPlot(prox, reduction = "umap", group.by = "Mapped_human_fetal_duo_cell_type")+DimPlot(prox, reduction = "umap", group.by = "Mapped_human_fetal_duo_age")
saveRDS(prox, file="Res_prox_with_CSS_integration_and_fetal_human_duodenum_projection.rds")

table(prox$Mapped_human_fetal_duo_age[which(prox$Age=="17.5d" & prox$Cell_type=="Lgr5+ / Olfm4+ ISC")])
# projection-based age inference seems not reliable

# extract the stem cells of different ages from human and mouse, and compare them based on transcriptome similarity
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/compare_to_human_duo/similarity_based")
load("/nas/groups/treutlein/USERS/Qianhui_Yu/Tools/scoreHIO/data/feature_gene_set.RData")
mouse_prox_sc <- subset(prox, subset = Cell_type %in% c("G2M-phase ISC", "Lgr5+ / Olfm4+ ISC", "Lgr5+/Shh+/Olfm4+ ISC", "Shh+ ISC"))
mouse_prox_sc_expr <- getAveExpr(seu.obj = mouse_prox_sc, feature.to.calc = "Age", colname.prefix = "mouse_prox")
human_duo_sc <- subset(human_epi, subset = Cell_type == "Stem_cell")
human_duo_sc_expr <- getAveExpr(seu.obj = human_duo_sc, feature.to.calc = "Age", colname.prefix = "human_duo",
                                specified.order = paste0("d", c(47,59,72,80,85,101,122,127,132)))
## use deconvolution feature genes
g1 <- feature_gene_set$deconvolution_feature_genes
## use top phase markers
g1 <- feature_gene_set$top_phase_markers
## identify DEGs between stem cells of different ages in human duodenum - this does not work well
#de_res <- wilcoxauc(human_duo_sc, group_by = "Age")
#de_res$pct_diff <- de_res$pct_in - de_res$pct_out
#sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
#library(dplyr)
#top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
#deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
#saveRDS(deg_res, file="Res_DEGs_between_stem_cells_of_different_ages_in_human_duodenum.rds")
#g1 <- unique(top_res$feature)

gene_pairs <- human_mouse_gene_link[which(human_mouse_gene_link[,"Human_symbol"] %in% g1 & human_mouse_gene_link[,"Mouse_symbol"] %in% rownames(mouse_epi)),]
ref_expr <- human_duo_sc_expr[gene_pairs[,"Human_symbol"],]
que_expr <- mouse_prox_sc_expr[gene_pairs[,"Mouse_symbol"],]
rownames(que_expr) <- gene_pairs[,"Human_symbol"]

scc_mat <- cor(que_expr, ref_expr, method="spearman")
mapped_human_age <- colnames(scc_mat)[apply(scc_mat, 1, which.max)]
names(mapped_human_age) <- rownames(scc_mat)
mapped_human_age

library(gplots)
#pdf("Plot_heatmap_SCC_human_duo_vs_mouse_prox_ISC_ages_based_on_stem_cell_maturity_phase_markers.pdf")
pdf("Plot_heatmap_SCC_human_duo_vs_mouse_prox_ISC_ages_based_on_deconvolution_feature_genes.pdf")
heatmap.2(scc_mat, col = darkBlue2Red.heatmap, trace = "none", density.info = "none", mar=c(12,12), 
          dendrogram = "none", scale="none", Rowv=FALSE, Colv = FALSE, cexRow = 1.5, cexCol = 1.5)
dev.off()

# find a reasonable way to resolve the heterogeneity of the proximal small intestine stem cell pool
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity")
## UMAP-CSS
mouse_prox_sc <- FindVariableFeatures(object = mouse_prox_sc, selection.method = "vst", nfeatures = 3000)
mouse_prox_sc <- ScaleData(object = mouse_prox_sc, verbose = T)
mouse_prox_sc <- RunPCA(object = mouse_prox_sc, features = VariableFeatures(mouse_prox_sc), verbose = F, npcs = 50)
usefulPCs <- 1:20
mouse_prox_sc <- RunUMAP(object = mouse_prox_sc, dims = usefulPCs)
DimPlot(mouse_prox_sc, reduction = "umap", group.by = "Age")+DimPlot(mouse_prox_sc, reduction = "umap_css", group.by = "Age")

# CSS integration
mouse_prox_sc <- cluster_sim_spectrum(mouse_prox_sc, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
mouse_prox_sc <- RunUMAP(mouse_prox_sc, reduction = "css", dims = 1:ncol(mouse_prox_sc@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
mouse_prox_sc <- FindNeighbors(object = mouse_prox_sc, reduction = "pca", dims = 1:ncol(mouse_prox_sc@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
mouse_prox_sc[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- mouse_prox_sc[[paste0("RNA_snn_res.", 0.2*5)]]
mouse_prox_sc <- FindNeighbors(object = mouse_prox_sc, reduction = "pca", dims = usefulPCs, force.recalc = T) %>%
  FindClusters(resolution = 1)
saveRDS(mouse_prox_sc, file="Res_mouse_fetal_prox_SC_13-18d_with_CSS_integration.rds")
DimPlot(mouse_prox_sc, reduction = "umap", group.by = "RNA_snn_res.1", label=T)+DimPlot(mouse_prox_sc, reduction = "umap", group.by = "Age")
plotFeature(seu.obj=mouse_prox_sc, 
            dr="umap",
            genes.to.plot = c("Shh", "Lgr5", "Ascl2", "Olfm4", "Mki67", "Fabp2", "Apoa4"),
            plot.name = "Plot_UMAP_proximal_SI_stem_cell_ISC_marker_expr.png",
            col.num = 4)

# identify DEGs between C7 and C9
idx1 <- which(mouse_prox_sc$RNA_snn_res.1==9)
idx2 <- which(mouse_prox_sc$RNA_snn_res.1==7)
X <- mouse_prox_sc@assays$RNA@data[,c(idx1, idx2)]
y <- rep(c("Early", "Late"), c(length(idx1), length(idx2)))
de_res <- wilcoxauc(X=X, y=y)
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[which(de_res$auc>0.6 & de_res$padj<0.01 & de_res$pct_diff>20 & de_res$pct_in>20 & de_res$logFC>0.1),]
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
mouse_decov_features <- unique(top_res$feature)
mouse_decov_ref_pattern <- cbind(rowMeans(mouse_prox_sc@assays$RNA@data[mouse_decov_features,idx1]),
                                 rowMeans(mouse_prox_sc@assays$RNA@data[mouse_decov_features,idx2]))

feature_gene <- mouse_decov_features
X <- mouse_decov_ref_pattern

que_expr <- as.matrix(mouse_prox_sc@assays$RNA@data[feature_gene,])
identity_que <- matrix(NA, nrow=ncol(que_expr), ncol=4)
rownames(identity_que) <- colnames(que_expr)
colnames(identity_que) <- c("frxn_early", "frxn_late", "Lagrangian", "Error")
for (j in seq(ncol(que_expr))){
  Y <- as.matrix(que_expr[,j])
  Rinv <- solve(chol(t(X) %*% X));
  C <- cbind(rep(1,2), diag(2))
  b <- c(1,rep(0,2))
  d <- t(Y) %*% X
  QP <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  Error <-sum(abs(Y-X%*%QP$solution))
  identity_que[j,] <- c(QP$solution[1],QP$solution[2],QP$Lagrangian[1],Error)
}
isc_score <- identity_que[,2]

plotFeature2(Embeddings(mouse_prox_sc, reduction = "umap"),
             values = isc_score,
             point.order = "sorted",
             zeroAsGray = FALSE)
mouse_prox_sc$Mouse_decov_score <- isc_score

## deconvolution based on features identified in human
load("/nas/groups/treutlein/USERS/Qianhui_Yu/Tools/scoreHIO/data/decov_ref.RData")
gene_pairs <- human_mouse_gene_link[which(human_mouse_gene_link[,"Human_symbol"]%in%rownames(decov_ref) & human_mouse_gene_link[,"Mouse_symbol"]%in%rownames(mouse_prox_sc)),]
feature_gene <- gene_pairs[,"Human_symbol"]
X <- decov_ref[feature_gene,]

sum(gene_pairs[,"Mouse_symbol"] %in% mouse_decov_features)

que_expr <- as.matrix(mouse_prox_sc@assays$RNA@data[gene_pairs[,"Mouse_symbol"],])
rownames(que_expr) <- gene_pairs[,"Human_symbol"]

identity_que <- matrix(NA, nrow=ncol(que_expr), ncol=4)
rownames(identity_que) <- colnames(que_expr)
colnames(identity_que) <- c("frxn_early", "frxn_late", "Lagrangian", "Error")
for (j in seq(ncol(que_expr))){
  Y <- as.matrix(que_expr[,j])
  Rinv <- solve(chol(t(X) %*% X));
  C <- cbind(rep(1,2), diag(2))
  b <- c(1,rep(0,2))
  d <- t(Y) %*% X
  QP <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  Error <-sum(abs(Y-X%*%QP$solution))
  identity_que[j,] <- c(QP$solution[1],QP$solution[2],QP$Lagrangian[1],Error)
}
isc_score <- identity_que[,2]
summary(isc_score)

mouse_prox_sc$Human_decov_score <- isc_score

plot(mouse_prox_sc$Human_decov_score, mouse_prox_sc$Mouse_decov_score, pch=16)
saveRDS(mouse_prox_sc, file="Res_mouse_fetal_prox_SC_13-18d_with_CSS_integration.rds")

human_score <- lapply(time_point, function(t){
  mouse_prox_sc$Human_decov_score[which(mouse_prox_sc$Age==t)]
})
names(human_score) <- time_point
beanplot(human_score, what = c(1,1,1,0))

mouse_score <- lapply(time_point, function(t){
  mouse_prox_sc$Mouse_decov_score[which(mouse_prox_sc$Age==t)]
})
names(mouse_score) <- time_point
beanplot(mouse_score, what = c(1,1,1,0))

png("Plot_UMAP_ISC_score_in_mouse_proxSI_stem_cells.png", height=2000, width=2000*2)
par(mfrow=c(1,2), mar=c(5,5,10,5))
plotFeature2(Embeddings(mouse_prox_sc, reduction = "umap"),
             values = mouse_prox_sc$Mouse_decov_score,
             point.order = "sorted",
             zeroAsGray = FALSE,
             main="Mouse-based ISC score",
             cex=4,
             cex.main=10,
             lwd=0.7)

plotFeature2(Embeddings(mouse_prox_sc, reduction = "umap"),
             values = mouse_prox_sc$Human_decov_score,
             point.order = "sorted",
             zeroAsGray = FALSE,
             main="Human-based ISC score",
             cex=4,
             cex.main=10,
             lwd=0.7)
dev.off()

pdf("Plot_ISC_score_in_mouse_proxSI_stem_cells.pdf")
plot(mouse_prox_sc$Human_decov_score, 
     mouse_prox_sc$Mouse_decov_score, 
     pch=16, 
     xlab="Human-based ISC score",
     ylab="Mouse-based ISC score",
     bty="n",
     cex.lab=1.2)
dev.off()


# add coarse cell type annotation
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells")
ct_vec <- mouse_epi$Cell_type
ct_group_vec <- ct_vec
ct_group_vec[grep("ISC", ct_group_vec)] <- "Stem_cell"
ct_group_vec[ct_group_vec=="Enterocyte precursor"] <- "Enterocyte_precursor"
ct_group_vec[ct_group_vec=="G2M-phase goblet"] <- "Goblet"
mouse_epi$Coarse_cell_type <- ct_group_vec
DimPlot(mouse_epi, group.by = "Coarse_cell_type", reduction = "umap_css")
saveRDS(mouse_epi, file="Res_mouse_SI_epi_with_coarse_cell_type_annotation.rds")

# identify cell type markers in human duodenum and mouse proximal SI
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI")
prox_epi <- subset(mouse_epi, subset = Inferred_SI_region_identity!="Dist_SI")
de_res <- wilcoxauc(prox_epi, group_by = "Coarse_cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=auc)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_fetal_mouse_proximal_SI_coarse_cell_type_markers.rds")
saveRDS(prox_epi, file="Res_fetal_mouse_proximal_SI.rds")

human_epi_sub <- subset(human_epi, cells=colnames(human_epi)[!human_epi$Cell_type %in% c("M_cell", "Tuft")]) 
human_epi_sub$Coarse_cell_type <- human_epi_sub$Cell_type
human_epi_sub$Coarse_cell_type[which(human_epi_sub$Coarse_cell_type=="Enteroendocrine")] <- "EEC"
saveRDS(human_epi_sub, file="Res_fetal_human_duodenum.rds")
de_res <- wilcoxauc(human_epi_sub, group_by = "Coarse_cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=auc)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_fetal_human_duodenum_coarse_cell_type_markers.rds")

# identify DEGs between corresponding cell types in human and mouse
ct <- c("EEC", "Enterocyte", "Enterocyte_precursor", "Goblet", "Stem_cell")
human_expr <- as.matrix(human_epi_sub@assays$RNA@data)
human_cell_type <- human_epi_sub$Coarse_cell_type
mouse_expr <- as.matrix(prox_epi@assays$RNA@data)
mouse_cell_type <- prox_epi$Coarse_cell_type

gene_pairs <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"]%in%rownames(human_expr) & human_mouse_gene_link[,"Mouse_symbol"]%in%rownames(mouse_expr),]
human_expr <- human_expr[gene_pairs[,"Human_symbol"],]
mouse_expr <- mouse_expr[gene_pairs[,"Mouse_symbol"],]
rownames(mouse_expr) <- gene_pairs[,"Human_symbol"]
de_res_by_cell_type <- list()
for(x in ct){
  print(x)
  human_idx <- which(human_cell_type==x)
  mouse_idx <- which(mouse_cell_type==x)
  X <- cbind(human_expr[,human_idx], mouse_expr[,mouse_idx])
  y <- rep(c("Human", "Mouse"), c(length(human_idx), length(mouse_idx)))
  de_res <- wilcoxauc(X=X, y=y)
  de_res$pct_diff <- de_res$pct_in - de_res$pct_out
  sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
  top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=auc)
  de_res_by_cell_type[[x]] <- list("sig_res"=sig_res, "top_res"=top_res)
}
saveRDS(de_res_by_cell_type, file="Res_fetal_human_duodenum_and_mouse_proximal_SI_DEGS_by_coarse_cell_types.rds")

lapply(names(de_res_by_cell_type), function(x){
  table(de_res_by_cell_type[[x]]$sig_res$group)
})

# identify cell type markers in day 17.5 mouse proximal SI
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5")
prox_17.5 <- subset(prox_epi, subset = Age=="17.5d")
prox_17.5 <- FindVariableFeatures(object = prox_17.5, selection.method = "vst", nfeatures = 3000)
prox_17.5 <- ScaleData(object = prox_17.5, verbose = T)
prox_17.5 <- RunPCA(object = prox_17.5, features = VariableFeatures(prox_17.5), verbose = F, npcs = 50)
usefulPCs <- 1:20
prox_17.5 <- FindNeighbors(object = prox_17.5, dims = usefulPCs)
prox_17.5 <- FindClusters(object = prox_17.5, resolution = 1)
prox_17.5 <- RunUMAP(object = prox_17.5, dims = usefulPCs)
saveRDS(prox_17.5, file="Res_mouse_fetal_prox_SI_d17.5.rds")
DimPlot(prox_17.5, reduction = "umap", group.by = "Phase", label=T)+

mouse_g1 <- c("Shh", "Lgr5", "Ascl2", "Olfm4",# stem cell
              "Mki67", "Pcna", "Top2a", "Cdk1", # cell cycle
              "Alpi", "Apoa4", "Fabp2", # enterocyte
              "Muc2", "Spink4",# goblet
              "Lyz1", "Defa22", "Mptx2", # paneth cell
              "Chga", "Chgb", "Neurog3", # EEC
              "Dclk1", "Rgs13", # Tuft
              "Gp2", "Best4", "Spib", "Ca7" # M cell
)
length(mouse_g1)
plotFeature.batch(seu.obj = prox_17.5, dr="umap", genes.to.plot = mouse_g1, nCols = blue.cols, cex=3, col.num = 6,
                  plot.name = "Plot_UMAP_fetal_mouse_d17.5_proxSI_epithelial_cell_type_marker_expression.png")

DimPlot(prox_17.5, reduction = "umap", group.by = "RNA_snn_res.1", label=T)/DimPlot(prox_17.5, reduction = "umap", group.by = "Coarse_cell_type", label=T)/DimPlot(prox_17.5, reduction = "umap", group.by = "Phase", label=T)&NoLegend()
coor=Embeddings(prox_17.5, reduction = "umap")
row_num=4
col_num=1
png("Plot_UMAP_d17.5_fetal_mouse_prox_SI.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,12,5))
plotFeature2(coor, 
             values = paste0("C", prox_17.5$RNA_snn_res.1),
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="RNA_snn_res.1",
             add.label = T,
             label.cex = 10)
plotFeature2(coor, 
             values = prox_17.5$d17.5_coarse_cell_type,
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="Refined_cell_type",
             add.label = T,
             label.cex = 10)

plotFeature2(coor, 
             values = prox_17.5$Coarse_cell_type,
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="Cell_type",
             add.label = T,
             label.cex = 10)
plotFeature2(coor, 
             values = prox_17.5$Phase,
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="Phase",
             add.legend = T,
             legend.cex = 10,
             legend.pos="bottomleft")
dev.off()

ct_anno <- list(
  "Stem_cell" = c(4, 8, 9, 10),
  "Enterocyte_precursor" = 3,
  "Enterocyte" = c(0, 1, 2, 5, 7),
  "Goblet" = c(12,6),
  "EEC" = 11
)
cl_vec <- prox_17.5$RNA_snn_res.1
ct_vec <- rep(NA, length(cl_vec))
for(x in names(ct_anno)){
  cl_x <- ct_anno[[x]]
  ct_vec[which(cl_vec%in%cl_x)] <- x
}
prox_17.5$d17.5_coarse_cell_type <- ct_vec
saveRDS(prox_17.5, file="Res_mouse_fetal_prox_SI_d17.5.rds")
de_res <- wilcoxauc(prox_17.5, group_by = "Coarse_cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=auc)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_fetal_mouse_d17.5_proximal_SI_coarse_cell_type_markers.rds")


# focusing on human duodenum 12 PCW and above
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/hg19_12w_and_above")
human_epi_old <- subset(human_epi_sub, cells=colnames(human_epi_sub)[human_epi_sub$Age_week > 11]) 
human_epi_old <- FindVariableFeatures(object = human_epi_old, selection.method = "vst", nfeatures = 3000)
human_epi_old <- ScaleData(object = human_epi_old, verbose = T)
human_epi_old <- RunPCA(object = human_epi_old, features = VariableFeatures(human_epi_old), verbose = F, npcs = 50)
usefulPCs <- 1:20
human_epi_old <- FindNeighbors(object = human_epi_old, dims = usefulPCs)
human_epi_old <- FindClusters(object = human_epi_old, resolution = 1)
human_epi_old <- RunUMAP(object = human_epi_old, dims = usefulPCs)
# CSS integration
human_epi_old <- cluster_sim_spectrum(human_epi_old, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
human_epi_old <- RunUMAP(human_epi_old, reduction = "css", dims = 1:ncol(human_epi_old@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
human_epi_old <- FindNeighbors(object = human_epi_old, reduction = "css", dims = 1:ncol(human_epi_old@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
human_epi_old[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- human_epi_old[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(human_epi_old, file="Res_fetal_human_duodenum_12PCW_and_older.rds")

de_res <- wilcoxauc(human_epi_old, group_by = "Coarse_cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=auc)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_fetal_human_duodenum_coarse_cell_type_markers.rds")

coor=Embeddings(human_epi_old, reduction = "umap_css")
row_num=3
col_num=1
png("Plot_UMAP_CSS_fetal_human_PCW12_and_above.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,12,5))
plotFeature2(coor, 
             values = paste0("C", human_epi_old$RNA_CSS_snn_res.1),
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="RNA_snn_res.1",
             add.label = T,
             label.cex = 10)
plotFeature2(coor, 
             values = human_epi_old$Coarse_cell_type,
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="Cell_type",
             add.label = T,
             label.cex = 10)
plotFeature2(coor, 
             values = human_epi_old$Age,
             point.order = "random",
             cex=5,
             lwd=0.8,
             cex.main=12,
             main="Age",
             add.legend = T,
             legend.cex = 10,
             legend.pos="bottomleft")
dev.off()


setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5")
# identify DEGs between corresponding cell types in human and mouse
ct <- c("EEC", "Enterocyte", "Enterocyte_precursor", "Goblet", "Stem_cell")
human_expr <- as.matrix(human_epi_old@assays$RNA@data)
human_cell_type <- human_epi_old$Coarse_cell_type
mouse_expr <- as.matrix(prox_17.5@assays$RNA@data)
mouse_cell_type <- prox_17.5$Coarse_cell_type

gene_pairs <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"]%in%rownames(human_expr) & human_mouse_gene_link[,"Mouse_symbol"]%in%rownames(mouse_expr),]
human_expr <- human_expr[gene_pairs[,"Human_symbol"],]
mouse_expr <- mouse_expr[gene_pairs[,"Mouse_symbol"],]
rownames(mouse_expr) <- gene_pairs[,"Human_symbol"]
de_res_by_cell_type <- list()
for(x in ct){
  print(x)
  human_idx <- which(human_cell_type==x)
  mouse_idx <- which(mouse_cell_type==x)
  X <- cbind(human_expr[,human_idx], mouse_expr[,mouse_idx])
  y <- rep(c("Human", "Mouse"), c(length(human_idx), length(mouse_idx)))
  de_res <- wilcoxauc(X=X, y=y)
  de_res$pct_diff <- de_res$pct_in - de_res$pct_out
  sig_res <- de_res[de_res$padj<0.01 & de_res$logFC>0.1 & de_res$auc>0.6 & de_res$pct_diff>10 & de_res$pct_in>20,]
  top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=auc)
  de_res_by_cell_type[[x]] <- list("sig_res"=sig_res, "top_res"=top_res)
}
saveRDS(de_res_by_cell_type, file="Res_fetal_human_duodenum_and_mouse_proximal_SI_DEGS_by_coarse_cell_types.rds")

sapply(names(de_res_by_cell_type), function(x){
  table(de_res_by_cell_type[[x]]$sig_res$group)
})

# get some examples
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/prox_SI/d17.5/examples")
mat <- de_res_by_cell_type[["EEC"]]$sig_res
#g1 <- intersect(mat$feature[mat$group=="Human"], c("ONECUT3", "ARX", "PCBP4", "CD200", "CALD1", "DPYSL2", "FABP5", "SCG5"))
g1 <- intersect(mat$feature[mat$group=="Mouse"], c("CITED4", "LMX1A", "IDS", "MBNL2", "ECM1", "SMAD5"))
length(g1)
mouse_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"]%in%g1,"Mouse_symbol"]
length(mouse_g1)
g1
row_num=length(g1)
col_num=4
png(paste0("Plot_UMAP_human-mouse_EECs_mouse_high_", paste(g1,collapse = "_"), ".png"), height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,12,5))
for(g in g1){
  mouse_g <- human_mouse_gene_link[human_mouse_gene_link[,"Human_symbol"]==g,"Mouse_symbol"]
  plotFeature2(Embeddings(prox_17.5, reduction = "umap"),
               values = prox_17.5@assays$RNA@data[mouse_g,],
               point.order = "sorted",
               nCols = blue.cols,
               cex=7,
               lwd=0.8,
               cex.main=12,
               main=paste(mouse_g,"mouse",sep="@"))
  plotFeature2(Embeddings(human_epi_old, reduction = "umap_css"),
               values = human_epi_old@assays$RNA@data[g,],
               point.order = "sorted",
               nCols = blue.cols,
               cex=7,
               lwd=0.8,
               cex.main=12,
               main=paste(g,"human",sep="@"))
  plotFeature2(Embeddings(mouse_epi, reduction = "umap_css"),
               values = mouse_epi@assays$RNA@data[mouse_g,],
               point.order = "sorted",
               nCols = blue.cols,
               cex=7,
               lwd=0.8,
               cex.main=12,
               main=paste(mouse_g,"mouse",sep="@"))
  plotFeature2(Embeddings(human_epi_sub, reduction = "umap_css"),
               values = human_epi_sub@assays$RNA@data[g,],
               point.order = "sorted",
               nCols = blue.cols,
               cex=7,
               lwd=0.8,
               cex.main=12,
               main=paste(g,"human",sep="@"))
  
}
dev.off()

# get cluster markers of stem cell pools
de_res <- wilcoxauc(X=mouse_prox_sc, group_by = "RNA_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[which(de_res$auc>0.6 & de_res$padj<0.05 & de_res$pct_diff>20 & de_res$pct_in>20 & de_res$logFC>0.1),]
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
deg_res <- list("sig_res"=sig_res,
                "top_res"=top_res)
saveRDS(deg_res, file="Res_fetal_mouse_proximal_ISC_cluster_markers_without_integration.rds")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity/de_novo_cluster_markers")
example_mat <- sig_res %>% group_by(group) %>% top_n(5, wt=logFC)
for(i in unique(sort(example_mat$group))){
  g1 <- example_mat$feature[which(example_mat$group==i)]
  plot_name <- paste0("Plot_UMAP_C",i,"_markers.png")
  plotFeature(seu.obj=mouse_prox_sc, 
              dr="umap",
              genes.to.plot = g1,
              nCols = blue.cols,
              plot.name = plot_name,
              col.num = length(g1))
}

score_mat <- sapply(sort(unique(sig_res$group)), function(i){
  g1 <- sig_res$feature[which(sig_res$group==i)]
  if(length(g1)>1){
    vec <- colSums(mouse_epi@assays$RNA@data[g1,])
  }else{
    vec <- mouse_epi@assays$RNA@data[g1,]
  }
  return(vec)
})
colnames(score_mat) <- paste0("C", sort(unique(sig_res$group)))

row_num=3
col_num=4
png("Plot_UMAP_CSS_moues_proxSI_cluster_marker_expr_in_whole_mouse_epi.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,12,5))
for(j in colnames(score_mat)){
  plotFeature2(Embeddings(mouse_epi, reduction = "umap_css"),
               values = score_mat[,j],
               point.order = "sorted",
               cex=4,
               lwd=0.8,
               cex.main=12,
               main=j)
}
dev.off()
# cluster 10 may be hemoglobin cluster, which should be excluded

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity/human_duo_sc")
human_duo_sc <- FindVariableFeatures(object = human_duo_sc, selection.method = "vst", nfeatures = 3000)
human_duo_sc <- ScaleData(object = human_duo_sc, verbose = T)
human_duo_sc <- RunPCA(object = human_duo_sc, features = VariableFeatures(human_duo_sc), verbose = F, npcs = 50)
usefulPCs <- 1:20
human_duo_sc <- FindNeighbors(object = human_duo_sc, dims = usefulPCs)
human_duo_sc <- FindClusters(object = human_duo_sc, resolution = 1)
human_duo_sc <- RunUMAP(object = human_duo_sc, dims = usefulPCs)
saveRDS(human_duo_sc, file="Res_fetal_human_duodenum_stem_cell_no_batch_correction.rds")

for(i in c(2,3,4,5,7,9)){
  g1 <- example_mat$feature[which(example_mat$group==i)]
  human_g1 <- human_mouse_gene_link[human_mouse_gene_link[,"Mouse_symbol"]%in%g1,"Human_symbol"]
  plot_name <- paste0("Plot_UMAP_mouse_C",i,"_markers_in_human_stem_cell.png")
  plotFeature(seu.obj=human_duo_sc, 
              dr="umap",
              genes.to.plot = human_g1,
              nCols = blue.cols,
              plot.name = plot_name,
              col.num = length(g1))
}



# investigate the stem cell maturation program in mouse and compare it to human
## check how stem cell DEGs change expression along stem cell maturation
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity")
human_qp_res_list <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/include_new_SI/human_endoderm_map/remove_sex_MT_ribo_genes/intestine_epi/duo_epi/stem_cell_maturation/Res_stem_cell_maturity_score_list_of_organoids_and_fetal.rds")
selected_time_points <- c("d47", "d59", "d72", "d80", "d85", "d101", "d122", "d127", "d132")
human_qp_vec <- unlist(human_qp_res_list[selected_time_points])
new_name <- sapply(names(human_qp_vec), function(vec){
  aa <- strsplit(x=vec, split = ".", fixed = T)[[1]][2]
  sub("_fetal", "", aa)
})
names(human_qp_vec) <- new_name
human_duo_sc$Human_decov_score <- human_qp_vec[colnames(human_duo_sc)]
saveRDS(human_duo_sc, file="Res_fetal_human_duodenum_stem_cell_no_batch_correction.rds")
#sc.dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/include_new_SI/human_endoderm_map/remove_sex_MT_ribo_genes/intestine_epi/duo_epi/stem_cell_maturation/HIO_to_adult/"
#human_expr_by_bin <- readRDS(paste0(sc.dir, "add_adult/Res_HIO_fetal_and_adult_stem_cell_expr_by_bin.rds"))
cell.num.per.bin=50
human_expr_by_bin <- getExprByPt(pt.vec=human_duo_sc$Human_decov_score, 
                                 expr.mat=as.matrix(human_duo_sc@assays$RNA@data),
                                 cell.num.per.bin = cell.num.per.bin)
saveRDS(human_expr_by_bin, file="Res_fetal_human_duodenum_stem_cell_expr_by_bin.rds")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity")
mouse_expr_by_bin <- getExprByPt(pt.vec=mouse_prox_sc$Mouse_decov_score, 
                                 expr.mat=as.matrix(mouse_prox_sc@assays$RNA@data),
                                 cell.num.per.bin = cell.num.per.bin)
saveRDS(mouse_expr_by_bin, file="Res_fetal_mouse_proxSI_stem_cell_expr_by_bin.rds")

g1 <- c("Shh","Ascl2","Lgr5","Olfm4","Mki67","Apoa4")
expr_mat <- mouse_expr_by_bin[g1,]
input <- t(apply(expr_mat, 1, function(vec){
  (vec-min(vec))/(max(vec)-min(vec))
}))
pdf("Plot_heatmap_ISC_marker_expr_along_stem_cell_maturation_in_mouse.pdf")
heatmap.2(input, col = darkBlue2Red.heatmap, trace = "none", density.info = "none", mar=c(12,12), 
          dendrogram = "none", scale="none", Rowv=FALSE, Colv = FALSE, cexRow = 1.5, cexCol = 1.5)
dev.off()

g1 <- c("Shh","Ascl2","Lgr5","Olfm4","Mki67","Apoa4")
pdf("Plot_barplot_ISC_marker_expr_along_stem_cell_maturation_in_mouse_and_human-2.pdf", height=3*length(g1), width=5*2)
par(mfrow=c(length(g1),2))
for(g in g1){
  barplot(mouse_expr_by_bin[g,], main=paste(g,"mouse",sep="@"))
  barplot(human_expr_by_bin[toupper(g),], main=paste(toupper(g),"human",sep="@"))
}
dev.off()

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity")
g1 <- c("Cldn6", "Hsd17b2", "Tuba1a", "Cdc20", "Cenpa", "Cdca3", "Ccnb2", "Birc5", "Ccl25", "Dmbt1", "Olfm4")
col_num=6
row_num=ceiling(2*length(g1)/col_num)
pdf("Plot_barplot_ISC_cluster_marker_expr_along_stem_cell_maturation_in_mouse_and_human-2.pdf", height=3*row_num, width=4*col_num)
par(mfrow=c(row_num, col_num))
for(g in g1){
  barplot(mouse_expr_by_bin[g,], main=paste(g,"mouse",sep="@"))
  barplot(human_expr_by_bin[toupper(g),], main=paste(toupper(g),"human",sep="@"))
}
dev.off()

# identify genes with variable expression levels along the trajectory
age.test.p <- splineBasedAgeTest(pseudotime.vec=seq(ncol(mouse_expr_by_bin)), expr.mat=mouse_expr_by_bin, df=5, mode="stringent")
padj <- p.adjust(age.test.p, method="BH")
logfc <- apply(mouse_expr_by_bin, 1, function(vec){
  max(vec)-min(vec)
})
age.test.res <- data.frame("Residual_comparison_P"=age.test.p, "BH_corrected_P"=padj, "Log(max-min)"=logfc, stringsAsFactors = F)
saveRDS(age.test.res, file="Res_age_test_along_fetal_mouse_proximal_ISC_maturation.rds")  

sc.genes <- rownames(age.test.res)[which(age.test.res$BH_corrected_P<0.01 & age.test.res$Log.max.min.>0.5)]
saveRDS(sc.genes, file="Res_fetal_mouse_proximal_ISC_maturation_dependent_genes.rds")

length(sc.genes)
sc.expr <- mouse_expr_by_bin[sc.genes,]
# hierarchical clustering for genes
hc.sc <- hclust(as.dist(1-cor(t(sc.expr), method="spearman")), method="ward.D2")
plot(hc.sc, hang=-1, cex=0.1)
abline(h=4)
g=cutree(hc.sc, h=4)
pt.gene.modules <- plotClusterExprProfile(expr=sc.expr, time.vec=seq(ncol(sc.expr)), group.vec=rep("fetal", ncol(sc.expr)), 
                                          cluster.vec=g, group.cols="#31a354", return.value=T, to.plot=T, 
                                          plot.name="Plot_fetal_mouse_proximal_ISC_maturation_related_gene_module_expr_profile_smoothed.pdf", 
                                          add.legend=F, legend.pos="topleft", cex.legend=2, col.num = 3, mean.do.smooth = T)
saveRDS(pt.gene.modules, file="Res_fetal_mouse_proximal_ISC_maturation_related_gene_module_expr_profile_in_HIO_and_fetal.rds")

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/fetal_mouse_vs_human/mouse_fetal_duodenum/epithelium/remove_doublet_and_pancreatic_cells/proximal_SI_scToEnt/resolve_sc_heterogeneity/exclude_mouse_c10")
mouse_prox_sc_sub <- subset(mouse_prox_sc, cells=colnames(mouse_prox_sc)[which(mouse_prox_sc$RNA_snn_res.1 != 10)])
mouse_pt_bin_res <- getExprByPt(pt.vec=mouse_prox_sc_sub$Mouse_decov_score, 
                                expr.mat=as.matrix(mouse_prox_sc_sub@assays$RNA@data),
                                mode="fix.bin.num",
                                bin.num=25,
                                return.idx = TRUE)
mouse_expr_by_bin <- mouse_pt_bin_res$expr.mat
barplot(mouse_expr_by_bin["Olfm4",])
mouse_prox_sc_sub$Mouse_decov_score_bin <- mouse_pt_bin_res$cell.idx.vec
# get the cell cycle phase distribution in each pseudotime bin
n1 <- sapply(sort(unique(mouse_prox_sc_sub$Mouse_decov_score_bin)), function(bin_idx){
  sapply(sort(unique(mouse_prox_sc_sub$Phase)), function(cycle_phase){
    sum(mouse_prox_sc_sub$Mouse_decov_score_bin==bin_idx & mouse_prox_sc_sub$Phase==cycle_phase)
  })
})
p1 <- t(t(n1)/colSums(n1))
phase_cols <- setNames(colorRampPalette(prettyrainbow)(3),c("G1", "G2M", "S"))
barplot(p1, col = phase_cols[rownames(p1)])

par(mfrow=c(1,2))
plotFeature2(Embeddings(mouse_prox_sc_sub, reduction = "umap"),
             values = mouse_prox_sc_sub$Phase,
             point.order = "random",
             gCols = phase_cols,
             add.legend = T,
             legend.pos = "topright",
             lwd=0.4,
             main="Cell cycle phase",
             cex.main=2)
plotFeature2(Embeddings(mouse_prox_sc_sub, reduction = "umap"),
             values = mouse_prox_sc_sub$Mapped_human_fetal_duo_cell_type,
             point.order = "random",
             add.legend = T,
             legend.pos = "topright",
             lwd=0.4,
             main="Mapped human cell type",
             cex.main=2)

human_expr_by_bin <- getExprByPt(pt.vec=human_duo_sc$Human_decov_score, 
                                 expr.mat=as.matrix(human_duo_sc@assays$RNA@data),
                                 mode="fix.bin.num",
                                 bin.num=25)
barplot(human_expr_by_bin["OLFM4",])
expr_by_bin_list <- list("human"=human_expr_by_bin,
                         "mouse"=mouse_expr_by_bin)
saveRDS(expr_by_bin_list, file="Res_fetal_mouse_proxSI_and_human_duodenum_stem_cell_expr_by_bin.rds")

g1 <- c("Shh","Ascl2","Lgr5","Olfm4","Mki67","Apoa4")
g1 <- c("Apoa1", "Dpp4", "Fabp1", "Aldob", "Sis", "Fabp2", "Apoa4")
plotFeature(seu.obj=mouse_prox_sc,
            genes.to.plot = g1,
            plot.name = "Plot_UMAP_enterocyte_marker_expr_in_moues_proxSI.png",
            col.num = 4)
pdf("Plot_barplot_enterocyte_marker_expr_along_stem_cell_maturation_in_mouse_and_human-2.pdf", height=3*length(g1), width=5*2)
par(mfrow=c(length(g1),2))
for(g in g1){
  barplot(mouse_expr_by_bin[g,], main=paste(g,"mouse",sep="@"))
  if(g=="Sis"){
    barplot(human_expr_by_bin["SI",], main="SI @ human")
  }else{
    barplot(human_expr_by_bin[toupper(g),], main=paste(toupper(g),"human",sep="@"))
  }
}
dev.off()

col_num=3
row_num=ceiling(length(g1)/col_num)
pdf("Plot_barplot_stem_cell_marker_expr_along_stem_cell_maturation_in_mouse.pdf", height=3*row_num, width=5*col_num)
par(mfrow=c(row_num,col_num))
for(g in g1){
  barplot(mouse_expr_by_bin[g,], main=g)
}
dev.off()


g1 <- feature_gene_set$deconvolution_feature_genes
g1 <- rownames(human_duo_sc)
gene_pairs <- human_mouse_gene_link[which(human_mouse_gene_link[,"Human_symbol"]%in%g1 & human_mouse_gene_link[,"Mouse_symbol"]%in%rownames(mouse_prox_sc)),]
e1 <- human_expr_by_bin[gene_pairs[,"Human_symbol"],]
e2 <- mouse_expr_by_bin[gene_pairs[,"Mouse_symbol"],]
rownames(e2) <- rownames(e1)

idx <- which(rowSums(e1)>0 & rowSums(e2)>0)
e1 <- e1[idx,]
e2 <- e2[idx,]

all_cor_vec <- sapply(rownames(e2), function(g){
  cor(e1[g,], e2[g,], method="spearman")
})

plot(density(all_cor_vec))

g1 <- c("Cldn6", "Hsd17b2", "Tuba1a", "Cdc20", "Cenpa", "Cdca3", "Ccnb2", "Birc5", "Ccl25", "Dmbt1", "Olfm4")
all_cor_vec[toupper(g1)]


# tHIO - remove weird duodenum cells that show depletion of PDX1, ONECUT2 and enrichment of SATB2, HOXA10
DimPlot(tIO_rna, reduction = "umap")

# in vitro HIO - remove WTC line
# read human fetal individuals covering different intestine regions (probably without colon) and run CSS integration on them to see whether different SI regions could be separated like in mouse

