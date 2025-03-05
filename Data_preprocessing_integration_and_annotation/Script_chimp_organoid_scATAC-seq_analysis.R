setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/signac_unified_peak_list_based")
source("~/Work/commonScript/Script_functions.R")
source("~/Work/commonScript/Script_scATAC-seq_function.R")
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(simspec)
library(dplyr)
library(presto)
library(Signac)
library(AnnotationHub)
library(GenomeInfoDb)
library(ChIPseeker)
library(reticulate)
library(igraph)
library(gplots)
library(GenomicRanges)

# session 1 - peak annotation
# read unified chimp peaks (.bed) and covert to GRange object
chimp_atac_aggr <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/chimp/Res_seurat_atac_aggr_with_peaks_called_in_individual_sample.rds")


# session 2 - analysis with chimp organoid unified peaks and re-quantified count matrix
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/signac_unified_peak_list_based")
# read unified peak data
chimp_atac <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/chimp/Res_updated_seurat_atac_aggr_with_peaks_called_in_individual_sample.rds")
seurat.atac <- chimp_atac
mat <- cbind(c("1020-UK-1-ATAC",
               "993-UK-1-ATAC",
               "UK-A1-ATAC",
               "UK-R7-JoC-D31-CIO"),
             c("C7-CIO-D89",
               "C2-CIO-D47",
               "Sandra-CIO-D32",
               "JoC-D31-CIO"))
for(i in seq(nrow(mat))){
  seurat.atac$Dataset[which(seurat.atac$Dataset==mat[i,1])] <- mat[i,2]
}
Idents(seurat.atac) <- seurat.atac$Dataset

#filtering and normalization
plot.name <- "Plot_ATAC_sample_quality_before_filtering_of.pdf"
pdf(plot.name, height=7, width=12)
par(mfrow=c(1,2))
VlnPlot(
  object = seurat.atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
)
TSSPlot(seurat.atac, group.by = 'high.tss') + NoLegend()
dev.off()

seurat.atac.subset <- subset(seurat.atac,
                             subset = peak_region_fragments > 5000 &
                               peak_region_fragments < 100000 &
                               pct_reads_in_peaks > 30 &
                               blacklist_ratio < 0.001 &
                               nucleosome_signal < 2 &
                               TSS.enrichment > 2
)

plot.name <- "Plot_ATAC_sample_quality_after_filtering.pdf"
pdf(plot.name, height=7, width=15)
par(mfrow=c(1,2))
VlnPlot(
  object = seurat.atac.subset,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
)
TSSPlot(seurat.atac.subset, group.by = 'high.tss') + NoLegend()
dev.off()

bk <- seurat.atac
seurat.atac <- seurat.atac.subset
seurat.atac <- RunTFIDF(seurat.atac)
seurat.atac <- FindTopFeatures(seurat.atac, min.cutoff = ncol(seurat.atac) * 0.01)
det_rates <- rowMeans(seurat.atac@assays$peaks_stage@counts > 0) 
VariableFeatures(seurat.atac) <- names(which(det_rates>0.01))
seurat.atac <- RunSVD(seurat.atac)
seurat.atac <- RunUMAP(seurat.atac, reduction="lsi", dims = 2:20)
seurat.atac <- FindNeighbors(seurat.atac, reduction="lsi", dims=2:20) %>% FindClusters(resolution = 1)
DimPlot(seurat.atac, group.by = "Dataset")
seurat.atac[['RNA']] <- CreateAssayObject(GeneActivity(seurat.atac))
seurat.atac <- NormalizeData(seurat.atac, assay = "RNA")
saveRDS(seurat.atac, file="Res_chimp_organoid_unified_ATAC_preprocessed.rds")

# C2-tCIO
setwd("C2_tCIO")
## extract C2-tCIO cells, which include all cell class
C2_tCIO <- subset(seurat.atac, cells=colnames(seurat.atac)[seurat.atac$Dataset=="C2-tCIO-W12"])
det_rates <- rowMeans(C2_tCIO@assays$peaks@counts > 0)
VariableFeatures(C2_tCIO) <- names(which(det_rates>0.05))
length(VariableFeatures(C2_tCIO))
C2_tCIO <- RunSVD(C2_tCIO)
C2_tCIO <- RunUMAP(C2_tCIO, reduction="lsi", dims = 2:20)
C2_tCIO <- FindNeighbors(C2_tCIO, reduction="lsi", dims=2:20) %>% FindClusters(resolution = 1)
DimPlot(C2_tCIO)
saveRDS(C2_tCIO, file="Res_chimp_C2_tCIO_organoid_ATAC_preprocessed.rds")

# C2/7-tCIO
setwd("C2_and_C7_tCIO")
## extract C2 and C7-tCIO cells, which include all cell class
tCIO <- subset(seurat.atac, cells=colnames(seurat.atac)[seurat.atac$Dataset%in%c("C2-tCIO-W12","C7-tCIO-W12")])
det_rates <- rowMeans(tCIO@assays$peaks@counts > 0)
VariableFeatures(tCIO) <- names(which(det_rates>0.05))
length(VariableFeatures(tCIO))
tCIO <- RunSVD(tCIO)
tCIO <- RunUMAP(tCIO, reduction="lsi", dims = 2:20)
tCIO <- FindNeighbors(tCIO, reduction="lsi", dims=2:20) %>% FindClusters(resolution = 1)
DimPlot(tCIO, group.by = "Dataset")
saveRDS(tCIO, file="Res_chimp_tCIO_organoid_ATAC_preprocessed.rds")


# find integration anchors
tCIO_list <- SplitObject(tCIO, split.by = "Dataset")
integration.anchors <- FindIntegrationAnchors(
  object.list = tCIO_list,
  anchor.features = rownames(tCIO),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = tCIO[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
saveRDS(integrated, file="Res_Signac_integrated_tCIO_ATAC_preprocessed.rds")
par(mfrow=c(2,3))
for(g in c("CDX2", "SOX2", "LGALS4", "ONECUT2", "PDX1")){
  plotFeature2(coor=Embeddings(integrated, reduction = "umap"),
               values = integrated@assays$RNA@data[g,],
               point.order = "sorted",
               main=g)
}


p2 <- DimPlot(integrated, group.by = "Dataset")

# CSS (pearson) integration
seurat_atac_aggr_css_pcc <- cluster_sim_spectrum(tCIO,
                                                 label_tag = "Dataset",
                                                 cluster_resolution = 0.6,
                                                 merge_spectrums = F,
                                                 verbose = T,
                                                 return_seuratObj = TRUE,
                                                 use_dr="lsi",
                                                 dims_use=2:20,
                                                 corr_method = "pearson")
seurat_atac_aggr_css_pcc <- RunUMAP(seurat_atac_aggr_css_pcc,
                                    reduction = "css",
                                    dims = 1:ncol(seurat_atac_aggr_css_pcc@reductions$css@cell.embeddings),
                                    reduction.name = "umap_css",
                                    reduction.key = "UMAPCSS_")


# CSS (spearman) integration
seurat_atac_aggr_css_scc <- cluster_sim_spectrum(tCIO,
                                                 label_tag = "Dataset",
                                                 cluster_resolution = 0.6,
                                                 merge_spectrums = F,
                                                 verbose = T,
                                                 return_seuratObj = TRUE,
                                                 use_dr="lsi",
                                                 dims_use=2:20,
                                                 corr_method = "spearman")

seurat_atac_aggr_css_scc <- RunUMAP(seurat_atac_aggr_css_scc,
                                    reduction = "css",
                                    dims = 1:ncol(seurat_atac_aggr_css_scc@reductions$css@cell.embeddings),
                                    reduction.name = "umap_css",
                                    reduction.key = "UMAPCSS_")
saveRDS(seurat_atac_aggr_css_pcc, file="Res_chimp_tCIO_seurat_atac_aggr_css_pcc.rds")
saveRDS(seurat_atac_aggr_css_scc, file="Res_chimp_tCIO_seurat_atac_aggr_css_scc.rds")

library(harmony)
seurat_atac_aggr_harmony <- RunHarmony(tCIO,
                                       group.by.vars = "Dataset",
                                       reduction="lsi",
                                       dims.use = 2:20,
                                       assay.use = "peaks_stage",
                                       project.dim = F)
seurat_atac_aggr_harmony <- RunUMAP(seurat_atac_aggr_harmony,
                                    reduction = "harmony",
                                    dims = 1:20,
                                    reduction.key = "UMAPHARMONY_",
                                    reduction.name="umap_harmony")
saveRDS(seurat_atac_aggr_harmony, file="Res_chimp_tCIO_seurat_atac_aggr_harmony.rds")

# LIGER
library(rliger)
seurat_atac_aggr_liger <- tCIO
seurat_atac_aggr_liger <- ScaleData(seurat_atac_aggr_liger, split.by = "Dataset", do.center = FALSE)
seurat_atac_aggr_liger <- RunOptimizeALS(seurat_atac_aggr_liger, k = 20, lambda = 5, split.by = "Dataset")
seurat_atac_aggr_liger <- RunQuantileNorm(seurat_atac_aggr_liger, split.by = "Dataset")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
seurat_atac_aggr_liger <- FindNeighbors(seurat_atac_aggr_liger, reduction = "iNMF", dims = 1:20)
seurat_atac_aggr_liger <- FindClusters(seurat_atac_aggr_liger, resolution = 0.3)
# Dimensional reduction and plotting
seurat_atac_aggr_liger <- RunUMAP(seurat_atac_aggr_liger, 
                                  dims = 1:ncol(seurat_atac_aggr_liger[["iNMF"]]),
                                  reduction = "iNMF",
                                  reduction.key = "UMAPLIGER_",
                                  reduction.name="umap_liger")
DimPlot(seurat_atac_aggr_liger, group.by = "Dataset", reduction = "umap_liger")
saveRDS(seurat_atac_aggr_liger, file="Res_chimp_tCIO_seurat_atac_aggr_liger.rds")

p1 <- DimPlot(tCIO, group.by = "Dataset")+ggtitle("No integration")
p2 <- DimPlot(integrated, group.by = "Dataset")+ggtitle("Seurat")
p3 <- DimPlot(seurat_atac_aggr_css_pcc, group.by = "Dataset", reduction = "umap_css")+ggtitle("CSS-PCC")
p4 <- DimPlot(seurat_atac_aggr_css_scc, group.by = "Dataset", reduction = "umap_css")+ggtitle("CSS-SCC")
p5 <- DimPlot(seurat_atac_aggr_harmony, group.by = "Dataset", reduction = "umap_harmony")+ggtitle("Harmony")
p6 <- DimPlot(seurat_atac_aggr_liger, group.by = "Dataset", reduction = "umap_liger")+ggtitle("LIGER")
p1+p2+p3+p4+p5+p6+plot_layout(ncol=3)

seurat_object_list <- list(
  "Non-integrated"=tCIO,
  "Seurat"=integrated,
  "CSS_PCC"=seurat_atac_aggr_css_pcc,
  "CSS_SCC"=seurat_atac_aggr_css_scc,
  "Harmony"=seurat_atac_aggr_harmony,
  "LIGER"=seurat_atac_aggr_liger
)

dr_vec <- setNames(c("umap", "umap", "umap_css", "umap_css", "umap_harmony", "umap_liger"), 
                   names(seurat_object_list))
for(x in names(seurat_object_list)){
  file <- paste0("Plot_selected_markers_in_",x,"_chimp_in_vitro_organoids.png")
  col_num=4
  row_num=2
  seu <- seurat_object_list[[x]]
  coor=Embeddings(seu, reduction=dr_vec[x])
  png(file,
      height=2000*row_num,
      width=2000*col_num)
  par(mfrow=c(row_num, col_num),
      mar=c(5,5,10,5))
  
  for(g in c("CDH1","CDX2", "SOX2", "LGALS4", "ONECUT2", "PDX1", "COL1A2", "ACTA2")){
    plotFeature2(coor = coor, 
                 values = seu@assays$RNA@data[g,],
                 main=g,
                 point.order = "sorted",
                 cex=4,
                 lwd=0.3,
                 cex.main=10)
    
  }
  dev.off()
}


# load gene-peak annotation
anno_df <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/panTro6/TxDB/Data_frame_peakAnno_by_ChIPseeker_2k_all.rds")
C7_tCIO_epi <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/Res_C7_tCIO_epi_with_imputed_RNA_mapped_cell_type.rds")

DimPlot(tCIO, label=T, group.by = "Cell_class")+NoLegend()
tCIO$Cell_class <- "Epithelial"
tCIO$Cell_class[which(tCIO$peaks_stage_snn_res.1%in%c(2,4,6,9,12,13))] <- "Mesenchymal"

seurat_atac_aggr_harmony <- FindNeighbors(seurat_atac_aggr_harmony, reduction = "harmony", dims = 1:20)
seurat_atac_aggr_harmony <- FindClusters(seurat_atac_aggr_harmony, resolution = 1)
DimPlot(seurat_atac_aggr_harmony, label=T, group.by = "Cell_class", reduction = "umap_harmony")+NoLegend()
seurat_atac_aggr_harmony$Cell_class <- "Epithelial"
seurat_atac_aggr_harmony$Cell_class[which(seurat_atac_aggr_harmony$peaks_stage_snn_res.1%in%c(1,6,9,11,12))] <- "Mesenchymal"
saveRDS(seurat_atac_aggr_harmony, file="Res_chimp_tCIO_seurat_atac_aggr_harmony.rds")

tCIO_epi <- subset(seurat_atac_aggr_harmony, Cell_class=="Epithelial")
det_rates <- rowMeans(tCIO_epi@assays$peaks@counts > 0)
VariableFeatures(tCIO_epi) <- names(which(det_rates>0.05))
length(VariableFeatures(tCIO_epi))
tCIO_epi <- RunSVD(tCIO_epi)
tCIO_epi <- RunUMAP(tCIO_epi, reduction="lsi", dims = 2:20)
tCIO_epi <- RunHarmony(tCIO_epi,
                       group.by.vars = "Dataset",
                       reduction="lsi",
                       dims.use = 2:20,
                       assay.use = "peaks_stage",
                       project.dim = F)
tCIO_epi <- RunUMAP(tCIO_epi,
                    reduction = "harmony",
                    dims = 1:20,
                    reduction.key = "UMAPHARMONY_",
                    reduction.name="umap_harmony")
p1 <- DimPlot(tCIO_epi, reduction ="umap_harmony", group.by = "Dataset")+ggtitle("Harmony")
p2 <- DimPlot(tCIO_epi, reduction ="umap", group.by = "Dataset")+ggtitle("Non-integrated")
p1+p2
tCIO_epi <- FindNeighbors(tCIO_epi, reduction="harmony", dims=1:20) %>% FindClusters(resolution = 1)

tCIO_epi$Mapped_RNA_cell_type <- NA
tCIO_epi@meta.data[intersect(colnames(tCIO_epi), colnames(C7_tCIO_epi)), "Mapped_RNA_cell_type"] <- C7_tCIO_epi@meta.data[intersect(colnames(tCIO_epi), colnames(C7_tCIO_epi)), "Mapped_RNA_cell_type"]
nn_mat <- tCIO_epi@graphs$peaks_stage_nn
na_cells <- colnames(tCIO_epi)[which(is.na(tCIO_epi$Mapped_RNA_cell_type))]
mapped_ct <- sort(unique(tCIO_epi$Mapped_RNA_cell_type))
nn_ident_freq <- t(sapply(na_cells, function(i){
  vec <- as.vector(as.matrix(nn_mat[i,]))
  nn_cells <- colnames(nn_mat)[which(vec==1)]
  freq <- sapply(mapped_ct, function(x){
    sum(tCIO_epi@meta.data[nn_cells, "Mapped_RNA_cell_type"]==x, na.rm = T)
  })
  return(freq) 
}))
cell_max_ct <- setNames(colnames(nn_ident_freq)[apply(nn_ident_freq, 1, which.max)],
                        rownames(nn_ident_freq))
tCIO_epi@meta.data[na_cells, "Mapped_RNA_cell_type"] <- cell_max_ct
p0 <- DimPlot(tCIO_epi, reduction ="umap", group.by = "Dataset")+ggtitle("Non-integrated")
p1 <- DimPlot(tCIO_epi, reduction ="umap_harmony", group.by = "Dataset")+ggtitle("Sample")
p2 <- DimPlot(tCIO_epi, reduction ="umap_harmony", label=T)+ggtitle("Cluster")+NoLegend()
p3 <- DimPlot(tCIO_epi, reduction ="umap_harmony", label=T, group.by = "Mapped_RNA_cell_type")+ggtitle("Mapped_RNA_cell_type")
p0+p1+p2+p3+plot_layout(ncol=2)


# get the example peak coverage plot for cell type markers
g1 <- c("LGR5", "OLFM4", "ASCL2", "FABP2", "DPP4", "APOA4", "SI","MUC2", "SPINK4","CHGA", "CHGB", "CDK1", "MKI67", "TOP2A", "PCNA")
ct <- c(rep("Stem cell", 3), rep("Enterocyte", 4), rep("Goblet cells", 2), rep("Enteroendocrine cells", 2), rep("Cycling stem cells", 4))

## get peaks annotated to selected cell type markers
p1 <- anno_df[which(anno_df$gene_symbol%in%g1),]
DefaultAssay(C7_tCIO_epi) <- "peaks_stage"
# get the peak detection rate in each mapped RNA cell type
r1 <- sapply(mapped_ct, function(x){
  idx <- which(C7_tCIO_epi$Mapped_RNA_cell_type==x)
  rowMeans(C7_tCIO_epi@assays$peaks@counts[,idx] > 0)
})
saveRDS(r1, file="Res_peak_detection_rate_per_mapped_RNA_cell_type.rds")
diff_r1 <- sapply(seq(ncol(r1)), function(j){
  r1[,j]-rowMeans(r1[,-j])
})
cm_idx <- sapply(seq(ncol(r1)), function(j){
  r1[,j]>0.2 & diff_r1[,j]>0.2
})
colnames(cm_idx) <- colnames(r1)
colnames(diff_r1) <- colnames(r1)
saveRDS(diff_r1, file="Res_peak_detection_rate_diff_per_mapped_RNA_cell_type.rds")
saveRDS(cm_idx, file="Res_cell_type_enriched_peaks.rds")

normed_expr <- getAveExpr(seu.obj=C7_tCIO_epi, assay.type = "peaks_stage", feature.to.calc = "Mapped_RNA_cell_type", colname.prefix = NULL)

order_cell_type <- c("Goblet cells", "Cycling stem cells", "Stem_cell_0", "Stem_cell_2", "Stem_cell_1", "Enterocyte", "Enteroendocrine cells")
cols_cell_type <- setNames(colorRampPalette(prettyrainbow)(length(order_cell_type)), sort(order_cell_type))
C7_tCIO_epi$Mapped_RNA_cell_type <- factor(C7_tCIO_epi$Mapped_RNA_cell_type,
                                           levels = order_cell_type)
library(ggplot2)
# cycling stem cell-enriched peaks
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/Cycling_stem_cell")
max_ct <- colnames(diff_r1)[apply(diff_r1, 1, which.max)]
mat <- diff_r1[which(max_ct=="Cycling stem cells"),]
idx <- order(mat[,"Cycling stem cells"], decreasing = T)[1:10]
p1 <- rownames(mat)[idx]
#p1 <- rownames(cm_idx)[which(cm_idx[,"Cycling stem cells"])]
anno_mat <- anno_df[which(anno_df$name%in%p1), c("name", "gene_symbol", "annotation_category")]
anno_mat <- anno_mat[which(anno_mat$gene_symbol%in%rownames(C7_tCIO_epi@assays$RNA@counts) & anno_mat$annotation_category=="Promoter"),]
for(i in seq(nrow(anno_mat))){
  p <- anno_mat$name[i]
  g <- anno_mat$gene_symbol[i]
  t <- anno_mat$annotation_category[i]
  plot_name <- paste0("Plot_C7-tCIO_",p,"_",g,"_",t,".png")
  png(plot_name)
  cp <- CoveragePlot(
    object = C7_tCIO_epi,
    region = p,
    extend.upstream = 1000,
    extend.downstream = 1000,
    features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE
  )
  print(cp & scale_fill_manual(values = cols_cell_type))
  dev.off()
}


# stem cell-enriched peaks
p1 <- rownames(cm_idx)[which(cm_idx[,"Stem_cell_0"])]
anno_mat <- anno_df[which(anno_df$name%in%p1), c("name", "gene_symbol", "annotation_category")]
anno_mat <- anno_mat[which(anno_mat$gene_symbol%in%rownames(C7_tCIO_epi@assays$RNA@counts)),]
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/Stem_cell")
for(i in seq(nrow(anno_mat))){
  p <- anno_mat$name[i]
  g <- anno_mat$gene_symbol[i]
  t <- anno_mat$annotation_category[i]
  plot_name <- paste0("Plot_C7-tCIO_",p,"_",g,"_",t,".pdf")
  pdf(plot_name)
  cp <- CoveragePlot(
    object = C7_tCIO_epi,
    region = p,
    extend.upstream = 1000,
    extend.downstream = 1000,
    features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE
  )
  cp & scale_fill_manual(values = cols_cell_type)
  dev.off()
}

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO")
cm_p_list <- lapply(c("Enterocyte", "Enteroendocrine cells", "Goblet cells"), function(x){
  intersect(anno_df$name[which(anno_df$gene_symbol%in%g1[which(ct==x)])],
            rownames(cm_idx)[which(cm_idx[,x])])
})
names(cm_p_list) <- c("Enterocyte", "Enteroendocrine cells", "Goblet cells")

anno_mat <- anno_df[which(anno_df$name%in%peaks), c("name", "gene_symbol", "annotation_category")]


for(i in seq(nrow(anno_mat))){
  p <- anno_mat$name[i]
  g <- anno_mat$gene_symbol[i]
  t <- anno_mat$annotation_category[i]
  plot_name <- paste0("Plot_C7-tCIO_",p,"_",g,"_",t,".pdf")
  pdf(plot_name)
  cp <- CoveragePlot(
    object = C7_tCIO_epi,
    group.by="Mapped_RNA_cell_type",
    region = p,
    extend.upstream = 1000,
    extend.downstream = 1000,
    features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE
  )
  print(cp & scale_fill_manual(values = cols_cell_type))
  dev.off()
}
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/selected_peaks")
selected_peaks <- c("chr11-112595349-112595978",
                    "chr20-6256539-6257472",
                    "chr14-73871176-73871463",
                    "chr11-1000115-1001145",
                    "chr1-106340566-106341765",
                    "chr4-54722597-54723637")
anno_mat <- anno_df[which(anno_df$name%in%selected_peaks), c("name", "gene_symbol", "annotation_category")]
for(i in seq(nrow(anno_mat))){
  p <- anno_mat$name[i]
  g <- anno_mat$gene_symbol[i]
  t <- anno_mat$annotation_category[i]
  plot_name <- paste0("Plot_C7-tCIO_",p,"_",g,"_",t,".pdf")
  pdf(plot_name)
  cp <- CoveragePlot(
    object = C7_tCIO_epi,
    group.by="Mapped_RNA_cell_type",
    region = p,
    extend.upstream = 1000,
    extend.downstream = 1000,
    features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE
  )
  print(cp & scale_fill_manual(values = cols_cell_type))
  dev.off()
}

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/RNA_ATAC_same_embeddings")
tCIO_panTro6_rna_matched <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/MCMF_based_integration/smoothed_SCC_based_ATAC_cell_type_annotation/only_averaging_ATAC_cells/Res_C7-tCIO_epi_panTro6_RNA_only_with_ATAC_matched_cells_seurat_object.rds")
meta_cell_info <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/MCMF_based_integration/smoothed_SCC_based_ATAC_cell_type_annotation/only_averaging_ATAC_cells/Res_meta_cell_info.rds")

# get the normalized peak intensity of constructed pseudocells
freq <- table(meta_cell_info$Pseudocell_idx)
atac_cells <- paste("CA-4",meta_cell_info$ATAC[match(names(freq)[which(freq==1)],meta_cell_info$Pseudocell_idx)],sep="_")
averaged_atac_data <- matrix(NA, nrow=nrow(C7_tCIO_epi@assays$peaks_stage@data), ncol=length(unique(meta_cell_info$Pseudocell_idx)))
rownames(averaged_atac_data) <- rownames(C7_tCIO_epi@assays$peaks_stage@data)
colnames(averaged_atac_data) <- paste("pseudocell", seq(length(unique(meta_cell_info$Pseudocell_idx))), sep="_")
averaged_atac_data[,names(freq)[which(freq==1)]] <- as.matrix(C7_tCIO_epi@assays$peaks_stage@data[,atac_cells])
for(j in names(freq)[which(freq>1)]){
  atac_cells <- paste("CA-4",meta_cell_info$ATAC[which(meta_cell_info$Pseudocell_idx==j)],sep="_")
  averaged_atac_data[,j] <- rowMeans(C7_tCIO_epi@assays$peaks_stage@data[,atac_cells])
}
saveRDS(averaged_atac_data, file="Dat_pseudocell_averaged_ATAC_normed_data.rds")


averaged_atac_gene_activity <- matrix(NA, nrow=nrow(C7_tCIO_epi@assays$RNA@data), ncol=length(unique(meta_cell_info$Pseudocell_idx)))
rownames(averaged_atac_gene_activity) <- rownames(C7_tCIO_epi@assays$RNA@data)
colnames(averaged_atac_gene_activity) <- paste("pseudocell", seq(length(unique(meta_cell_info$Pseudocell_idx))), sep="_")
averaged_atac_gene_activity[,names(freq)[which(freq==1)]] <- as.matrix(C7_tCIO_epi@assays$RNA@data[,atac_cells])
for(j in names(freq)[which(freq>1)]){
  atac_cells <- paste("CA-4",meta_cell_info$ATAC[which(meta_cell_info$Pseudocell_idx==j)],sep="_")
  averaged_atac_gene_activity[,j] <- rowMeans(C7_tCIO_epi@assays$RNA@data[,atac_cells])
}
saveRDS(averaged_atac_gene_activity, file="Dat_pseudocell_averaged_ATAC_gene_activity_normed_data.rds")

mat <- unique(meta_cell_info[, c("Pseudocell_idx", "RNA_UMAP.UMAP_1", "RNA_UMAP.UMAP_2")])
rownames(mat) <- mat$Pseudocell_idx

coor <- mat[colnames(averaged_atac_data),-1]
# finalize cell type nomenclature
update_cell_type <- list(
  "Stem_cell_0"  = "Stem cell",
  "Stem_cell_1" = "Early enterocyte",
  "Stem_cell_2" = "Progenitor",
  "Cycling stem cells" = "Transit amplifying"
)
new_ct_vec <- tCIO_panTro6_rna_matched$Cell_type
for(x in names(update_cell_type)){
  new_ct_vec[which(tCIO_panTro6_rna_matched$Cell_type==x)] <- update_cell_type[[x]]
}
tCIO_panTro6_rna_matched$Refined_cell_type <- new_ct_vec
DimPlot(tCIO_panTro6_rna_matched, group.by = "Refined_cell_type")
saveRDS(tCIO_panTro6_rna_matched, file="Res_tCIO_panTro6_rna_matched_with_updated_cell_type_name.rds")
new_ct_vec <- C7_tCIO_epi
for(x in names(update_cell_type)){
  new_ct_vec[which(tCIO_panTro6_rna_matched$Cell_type==x)] <- update_cell_type[[x]]
}
tCIO_panTro6_rna_matched$Refined_cell_type <- new_ct_vec
DimPlot(tCIO_panTro6_rna_matched, group.by = "Refined_cell_type")
saveRDS(tCIO_panTro6_rna_matched, file="Res_tCIO_panTro6_rna_matched_with_updated_cell_type_name.rds")

for(x in names(update_cell_type)){
  meta_cell_info$RNA_cell_ident[which(meta_cell_info$RNA_cell_ident==x)] <- update_cell_type[[x]]
}

cell_ident <- meta_cell_info$RNA_cell_ident[match(rownames(coor), meta_cell_info$Pseudocell_idx)]

dat <- list("RNA_coor"=coor,
            "Gene_activity"=averaged_atac_gene_activity,
            "Peak_accessibility"=averaged_atac_data,
            "RNA_cell_type"=cell_ident)
saveRDS(dat, file="Res_pseudocell_activity_and_embeddings.rds")

C7_tCIO_epi$Mapped_RNA_cell_type <- as.character(C7_tCIO_epi$Mapped_RNA_cell_type)
for(x in names(update_cell_type)){
  C7_tCIO_epi$Mapped_RNA_cell_type[which(C7_tCIO_epi$Mapped_RNA_cell_type==x)] <- update_cell_type[[x]]
}
saveRDS(C7_tCIO_epi, file="/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/Res_C7_tCIO_epi_with_imputed_RNA_mapped_cell_type.rds")

order_cell_type <- c("Goblet cells", "Transit amplifying", "Stem cell", "Progenitor", "Early enterocyte", "Enterocyte", "Enteroendocrine cells")
cols_cell_type <- setNames(colorRampPalette(prettyrainbow)(length(order_cell_type)), sort(order_cell_type))
C7_tCIO_epi$Mapped_RNA_cell_type <- factor(C7_tCIO_epi$Mapped_RNA_cell_type,
                                           levels = order_cell_type)


# plot cell type marker gene expression
genes <- c("MKI67", "OLFM4", "LGR5", "MUC2",  "APOA4", "CHGA")
row_num=1
col_num=6
png("Plot_UMAP_C7_tCIO_epi_RNA.png",
    height=2000*row_num,
    width=2000*col_num)
par(mfrow=c(row_num,col_num),
    mar=c(5,5,10,5))
for(g in genes){
  plotFeature2(Embeddings(tCIO_panTro6_rna_matched, reduction = "umap"),
               values = tCIO_panTro6_rna_matched@assays$RNA@data[g,],
               point.order = "sorted",
               nCols = beach.col,
               main=g,
               lwd=0.5,
               cex=6,
               cex.main=8)
}
dev.off()

png("Plot_UMAP_C7_tCIO_epi_ATAC_gene_activity.png",
    height=2000*row_num,
    width=2000*col_num)
par(mfrow=c(row_num,col_num),
    mar=c(5,5,10,5))
for(g in genes){
  plotFeature2(coor,
               values = averaged_atac_gene_activity[g,],
               point.order = "sorted",
               nCols = blue.cols,
               main=g,
               lwd=0.5,
               cex=6,
               cex.main=8)
}
dev.off()


png("Plot_UMAP_C7_tCIO_epi_ATAC_gene_activity_on_ATAC_embedding.png",
    height=2000*row_num,
    width=2000*col_num)
par(mfrow=c(row_num,col_num),
    mar=c(5,5,10,5))
for(g in genes){
  plotFeature2(Embeddings(C7_tCIO_epi, reduction = "umap"),
               values = C7_tCIO_epi@assays$RNA@data[g,],
               point.order = "sorted",
               nCols = blue.cols,
               main=g,
               lwd=0.5,
               cex=6,
               cex.main=8)
}
dev.off()

row_num=1
col_num=6
png("Plot_UMAP_C7_tCIO_epi_ATAC_peak_accessibility.png",
    height=2000*row_num,
    width=2000*col_num)
par(mfrow=c(row_num,col_num),
    mar=c(5,5,10,5))
for(g in anno_mat$gene_symbol){
  p <- anno_mat$name[which(anno_mat$gene_symbol==g)]
  plotFeature2(coor,
               values = averaged_atac_data[p,],
               point.order = "sorted",
               nCols = c(
                 "#feebe2",
                 "#fbb4b9",
                 "#f768a1",
                 "#c51b8a",
                 "#7a0177"),
               main=g,
               lwd=1,
               cex=8,
               cex.main=8)
}
dev.off()


png("Plot_UMAP_C7_tCIO_epi_RNA_and_ATAC_cell_type_no_label.png",
    height=2000*2,
    width=2000*1)
par(mfrow=c(2,1),xpd=TRUE,mar=c(15,15,15,15))
plotFeature2(Embeddings(tCIO_panTro6_rna_matched, reduction = "umap"),
             values = as.character(tCIO_panTro6_rna_matched$Refined_cell_type),
             point.order = "random",
             gCols = cols_cell_type,
             lwd=0.5,
             cex=3,
             cex.main=8,
             add.label = F,
             label.cex = 5)

plotFeature2(Embeddings(C7_tCIO_epi, reduction = "umap"),
             values = as.character(C7_tCIO_epi$Mapped_RNA_cell_type),
             point.order = "random",
             gCols = cols_cell_type,
             lwd=0.5,
             cex=3,
             cex.main=8,
             add.label = F,
             label.cex = 5)
dev.off()

pdf("Plot_UMAP_C7_tCIO_epi_RNA_and_ATAC_cell_type_label_only.pdf")
par(mfrow=c(2,1),xpd=TRUE)
plotFeature2(Embeddings(tCIO_panTro6_rna_matched, reduction = "umap"),
             values = as.character(tCIO_panTro6_rna_matched$Refined_cell_type),
             point.order = "random",
             gCols = cols_cell_type,
             lwd=0.5,
             cex=3,
             cex.main=8,
             label.only = T)

plotFeature2(Embeddings(C7_tCIO_epi, reduction = "umap"),
             values = as.character(C7_tCIO_epi$Mapped_RNA_cell_type),
             point.order = "random",
             gCols = cols_cell_type,
             lwd=0.5,
             cex=3,
             cex.main=8,
             label.only = T)
dev.off()



setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/selected_peaks")
selected_peaks <- c("chr11-112595349-112595978",
                    "chr20-6256539-6257472",
                    "chr14-73871176-73871463",
                    "chr11-1000115-1001145",
                    "chr1-106340566-106341765",
                    "chr4-54722597-54723637")
anno_mat <- anno_df[which(anno_df$name%in%selected_peaks), c("name", "gene_symbol", "annotation_category")]
for(i in seq(nrow(anno_mat))){
  p <- anno_mat$name[i]
  g <- anno_mat$gene_symbol[i]
  t <- anno_mat$annotation_category[i]
  plot_name <- paste0("Plot_C7-tCIO_",p,"_",g,"_",t,".pdf")
  pdf(plot_name)
  cp <- CoveragePlot(
    object = C7_tCIO_epi,
    group.by="Mapped_RNA_cell_type",
    region = p,
    extend.upstream = 1000,
    extend.downstream = 1000,
    features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE
  )
  print(cp & scale_fill_manual(values = cols_cell_type))
  dev.off()
}


tCIO_panTro6_rna_matched <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/RNA_ATAC_same_embeddings/Res_tCIO_panTro6_rna_matched_with_updated_cell_type_name.rds")
tCIO_panTro6_rna_matched <- RunTSNE(tCIO_panTro6_rna_matched, reduction="pca", dims = 1:20)
DimPlot(tCIO_panTro6_rna_matched, reduction = "umap", group.by = "Refined_cell_type", label = T)+DimPlot(tCIO_panTro6_rna_matched, reduction = "tsne", group.by = "Refined_cell_type", label = T)&NoLegend()
saveRDS(tCIO_panTro6_rna_matched, 
        file="/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6_scATAC/epi/analysis_with_RNA/based_on_chimp_organoid_unified_peak_list/C7-tCIO/RNA_ATAC_same_embeddings/Res_tCIO_panTro6_rna_matched_with_updated_cell_type_name.rds")

