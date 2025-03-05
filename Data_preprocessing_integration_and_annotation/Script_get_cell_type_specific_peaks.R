# session 1 - peak annotation
seurat_atac <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/HIO_and_tHIO/scATAC-seq/tHIO/Res_tHIO_ATAC_preprocessed_relax_filtering_0.01hvg.rds")
#seurat_rna <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/suspension_HIO_derived_tHIO/RNA_human_chimp_consensus/Res_tHIO_with_cell_class_and_epi_cell_type_annotation.rds")

# annotate peaks with ChIPseeker
gr_atac <- seurat_atac@assays$peaks@ranges

library(ChIPseeker) # need to install the most updated ChIPseeker from github, otherwise there would be errors
#library(Signac)
#library(Seurat)
library(AnnotationHub)
library(GenomeInfoDb)
ah <- AnnotationHub()
ensdb.list <- query(ah, c("EnsDb.Hsapiens"))
mat <- cbind(ensdb.list$ah_id, ensdb.list$title)
ens.version <- 93 # ensembl v93, same as the RNA ref (cellranger-ref v3.0)
id <- mat[which(mat[,2]==paste("Ensembl", ens.version, "EnsDb for Homo Sapiens")),1]
ensdb <- ensdb.list[[id]] 
seqlevelsStyle(ensdb) <- 'UCSC'
# ensdb object could not be saved

peakAnno <- annotatePeak(gr_atac,
                         tssRegion=c(-2000, 2000),
                         TxDb = ensdb,
                         annoDb=NULL,
                         overlap="all")
peakAnno@anno$annotation_category <- gsub(" \\(.+\\)$", "", peakAnno@anno$annotation)
anno_df <- as.data.frame(peakAnno@anno)
saveRDS(peakAnno, file="Res_peakAnno_by_ChIPseeker.rds")
anno_df$name <- paste(anno_df$seqnames, anno_df$start, anno_df$end, sep="-")
saveRDS(anno_df, file="Data_frame_peakAnno_by_ChIPseeker.rds")

peak_gene <- anno_df[,c("annotation", "geneId", "name", "annotation_category")]
idx <- which(peak_gene$annotation_category %in% c("Exon", "Intron"))
for(i in which(peak_gene$annotation_category =="Exon")){
  peak_gene$geneId[i] <- substr(peak_gene$annotation[i],23,37)
}
for(i in which(peak_gene$annotation_category =="Intron")){
  peak_gene$geneId[i] <- substr(peak_gene$annotation[i],25,39)
}

load("/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/Ensembl/Human/v93/ensembl.v93.hg38.RData")
mat <- unique(ensembl.v93.hg38[,c("gene_id", "gene_name")])
peak_gene <- peak_gene[which(peak_gene$geneId %in% mat$gene_id),]
gene_name_vec <- setNames(mat$gene_name, mat$gene_id)
peak_gene$gene_symbol <- gene_name_vec[peak_gene$geneId]
saveRDS(peak_gene, file="Data_frame_peak_gene_pairs.rds")

p1 <- DimPlot(seurat_atac, label=T)+NoLegend()+ggtitle(label = "ATAC cluster")
p2 <- FeaturePlot(seurat_atac, features = "EPCAM", order = T)
p3 <- FeaturePlot(seurat_atac, features = "COL1A2", order=T)
p1+p2+p3

# get the example peak coverage plot for cell type markers
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/HIO_and_tHIO/scATAC-seq/tHIO/hvg_0.05/coverage_plot")
g1 <- c("LGR5", "OLFM4", "ASCL2", "FABP2", "DPP4", "APOA4", "SI","MUC2", "SPINK4","CHGA", "CHGB", "CDK1", "MKI67", "TOP2A", "PCNA")
ct <- c(rep("Stem cell", 3), rep("Enterocyte", 4), rep("Goblet cells", 2), rep("Enteroendocrine cells", 2), rep("Cycling stem cells", 4))
p1 <- peak_gene[which(peak_gene$gene_symbol%in%g1),]
input <- p1
for(i in seq(nrow(input))){
  p <- input$name[i]
  g <- input$gene_symbol[i]
  t <- input$annotation_category[i]
  plot_name <- paste0("Plot_H9_tHIO_",p,"_",g,"_",t,".pdf")
  pdf(plot_name)
  cp <- CoveragePlot(
    object = seurat_atac,
    region = p,
    extend.upstream = 1000,
    extend.downstream = 1000,
    features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE
  )
  #print(cp & scale_fill_manual(values = gCols))
  print(cp)
  dev.off()
}


DefaultAssay(seurat_atac) <- "peaks"
# get the peak detection rate in each epithelial cluster
r1 <- sapply(c(1,2,3,7,9,12), function(x){
  idx <- which(seurat_atac$peaks_snn_res.1==x)
  rowMeans(seurat_atac@assays$peaks@counts[,idx] > 0)
})
colnames(r1) <- paste0("C", c(1,2,3,7,9,12))
saveRDS(r1, file="Res_peak_detection_rate_per_epithelial_cluster.rds")
diff_r1 <- sapply(seq(ncol(r1)), function(j){
  r1[,j]-rowMeans(r1[,-j])
})
cm_idx <- sapply(seq(ncol(r1)), function(j){
  r1[,j]>0.2 & diff_r1[,j]>0.2
})
colnames(cm_idx) <- colnames(r1)
colnames(diff_r1) <- colnames(r1)
saveRDS(diff_r1, file="Res_peak_detection_rate_diff_per_epithelial_cluster.rds")
saveRDS(cm_idx, file="Res_cell_type_enriched_peaks.rds")

cluster_enriched_peaks <- rownames(cm_idx)[which(rowSums(cm_idx)>0)]
annotated_cluster_enriched_peaks <- p1[p1$name%in%cluster_enriched_peaks,]
gCols <- setNames(colorRampPalette(prettyrainbow)(length(unique(seurat_rna$RNA_snn_res.1))), sort(unique(seurat_rna$RNA_snn_res.1)))



#normed_expr <- getAveExpr(seu.obj=C7_tCIO_epi, assay.type = "peaks_stage", feature.to.calc = "Mapped_RNA_cell_type", colname.prefix = NULL)

#order_cell_type <- c("Goblet cells", "Cycling stem cells", "Stem_cell_0", "Stem_cell_2", "Stem_cell_1", "Enterocyte", "Enteroendocrine cells")
#cols_cell_type <- setNames(colorRampPalette(prettyrainbow)(length(order_cell_type)), sort(order_cell_type))
#C7_tCIO_epi$Mapped_RNA_cell_type <- factor(C7_tCIO_epi$Mapped_RNA_cell_type,
#                                           levels = order_cell_type)
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
