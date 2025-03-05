setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/spheroids/Analysis")

tpm.mat <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/spheroids/expr_matrix/Dat_gene_type_fitlered_TPM_matrix_with_gene_symbol.rds")
#tpm.mat <- readRDS("~/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/spheroids/expr_matrix/Dat_gene_type_fitlered_TPM_matrix_with_gene_symbol.rds")

expressed.sample.num.vec <- rowSums(tpm.mat>0)
gene.num.vec <- colSums(tpm.mat>0)
length(expressed.sample.num.vec)
hist(gene.num.vec)

# only keep genes detected in at least 2 samples, samples with more 5k genes
library(Seurat)
source("~/Work/commonScript/Script_functions.R")
spheroid <- CreateSeuratObject(counts=tpm.mat, min.cells = 2, min.features = 5000) 
mito.genes <- grep(pattern = "^MT-", x = rownames(spheroid), value = TRUE)
# Calculate the percentage of genes in each cell that are mitochondrial and add number as metadata
percent.mito <- Matrix::colSums(spheroid@assays$RNA@counts[mito.genes, ])/ Matrix::colSums(spheroid@assays$RNA@counts)
spheroid$percent.mito <- percent.mito
mat <- do.call('rbind', strsplit(colnames(spheroid), "_"))
spheroid$Cell_line <- mat[,1]
spheroid$orig.ident <- mat[,3]
spheroid$Collection_date <- mat[,2]

VlnPlot(spheroid, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size=0 )

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
counts <- spheroid@assays$RNA@counts
counts <- counts[setdiff(rownames(counts), genes.to.remove),]
meta <- spheroid@meta.data[,-c(2,3)]
spheroid <- CreateSeuratObject(counts=counts, meta=meta)
spheroid <- NormalizeData(object = spheroid, normalization.method = "LogNormalize", scale.factor = 1e6)
# global variable genes
spheroid <- FindVariableFeatures(object = spheroid, selection.method = "vst", nfeatures = 3000)
spheroid <- ScaleData(object = spheroid, verbose = T)
spheroid <- RunPCA(object = spheroid, features = VariableFeatures(spheroid), verbose = F, npcs = 50)
saveRDS(spheroid, file="Res_spheroid_with_PCA.rds")
coor <- Embeddings(spheroid, reduction = "pca")[,1:2]
usefulPCs <- 1:5
spheroid <- RunUMAP(object = spheroid, dims = usefulPCs)
DimPlot(spheroid, reduction = "umap", group.by = "Cell_line")+DimPlot(spheroid, reduction = "umap", group.by = "Collection_date")
spheroid$Species <- "Human"
spheroid$Species[which(spheroid$Cell_line%in%c("JoC", "Sandra"))] <- "Chimp"
DimPlot(spheroid, reduction = "umap", group.by = "Cell_line")+DimPlot(spheroid, reduction = "umap", group.by = "Species")
saveRDS(spheroid, file="Res_spheroid_with_UMAP.rds")
plotFeature(seu.obj=spheroid, dr="umap", genes.to.plot = c("SOX2", "CLDN18","LGALS4","CDX2", "PDX1","ONECUT2", "SATB2", "GATA4", "GATA6"),
            col.num = 4, plot.name = "Plot_UMAP_human_and_chimp_spheroid_epi_markers.png")
# show the distribution of selected genes in each collection date of each species
expr.list <- list()
selected.genes <- c("SOX2", "CLDN18","LGALS4","CDX2", "TM4SF4","PDX1","ONECUT2", "GATA4", "SATB2","GUCA2A","OSR2","FZD10")
region.ident <- rep(c("Foregut", "Intestine", "Proximal intestine", "Non-distal intestine", "Distal intestine"),c(2,2,3,1,4))
names(region.ident) <- selected.genes
#for(line in c("B7-WT", "H9-WT", "JoC", "Sandra")){
for(x in c("Human", "Chimp")){
  for(g in selected.genes){
    expr.list[[x]][[g]] <- lapply(paste0("D", 3:6), function(date){
      #sample.idx <- which(spheroid$Cell_line==line & spheroid$Collection_date==date)
      sample.idx <- which(spheroid$Species==x & spheroid$Collection_date==date)
      spheroid@assays$RNA@data[g, sample.idx]
    })
    names(expr.list[[x]][[g]]) <- paste0("D", 3:6)
  }
}
row.num=length(expr.list)
col.num=length(selected.genes)
pdf("Plot_boxplot_epithelial_regional_marker_expr_group_by_species.pdf", height=5*row.num, width=4*col.num)
par(mfrow=c(row.num, col.num))
for(x in names(expr.list)){
  for(g in names(expr.list[[x]])){
    boxplot(expr.list[[x]][[g]], main=paste(paste(g, x, sep="@"), region.ident[g], sep="\n"), cex.main=2, cex.axis=2)
  }
}
dev.off()

# identify variable genes in human and chimp separately
combined <- spheroid
hs <- subset(spheroid, cells=colnames(spheroid)[which(spheroid$Species=="Human")])
cs <- subset(spheroid, cells=colnames(spheroid)[which(spheroid$Species=="Chimp")])
hs <- FindVariableFeatures(object = hs, selection.method = "vst", nfeatures = 3000)
cs <- FindVariableFeatures(object = cs, selection.method = "vst", nfeatures = 3000)
selected.hvg <- intersect(VariableFeatures(hs), VariableFeatures(cs))
length(selected.hvg)
VariableFeatures(combined) <- selected.hvg
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:10
combined <- RunUMAP(object = combined, dims = usefulPCs)
DimPlot(combined, reduction = "umap", group.by = "Cell_line")+DimPlot(combined, reduction = "umap", group.by = "Collection_date")

coor <- Embeddings(spheroid, reduction = "pca")[,c(1,2)]
par(mfrow=c(1,2))
plotFeature2(coor=coor, values=combined$Species, add.legend = T)
plotFeature2(coor=coor, values=combined$Collection_date, add.legend = T)

# analysis to individual species, and identify markers of each collection date in human and chimp separately
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/spheroids/Analysis/date_marker_by_species_based")
library(presto)
date.marker.list <- list()
for(x in c("hs", "cs")){
  seu.obj <- get(x)
  expr.mat <- as.matrix(seu.obj@assays$RNA@data)
  de.res <- wilcoxauc(X=seu.obj, group_by = "Collection_date")
  sig.res <- de.res[de.res$padj<0.1 & de.res$logFC>0.05 & de.res$pct_in>20,]
  sig.genes <- unique(sig.res$feature)
  date.marker.list[[x]] <- sig.genes
}
saveRDS(date.marker.list, file="Res_date_marker_by_species.rds")

union.sig.genes <- unique(unlist(date.marker.list))
sig.genes <- union.sig.genes
sig.genes <- intersect(date.marker.list[["hs"]], date.marker.list[["cs"]])
length(sig.genes)
seu.obj <- spheroid
expr.mat <- as.matrix(seu.obj@assays$RNA@data)
gene.idx <- which(rownames(expr.mat) %in% sig.genes)
d.mat <- 1-cor(expr.mat[gene.idx,])
dim(d.mat)
mds.coor <- cmdscale(d=d.mat, k=2)
pdf("Plot_union_date_marker_based_MDS_marker_interset_based.pdf", height=5, width=10)
par(mfrow=c(1,2))
plotFeature2(mds.coor, values = seu.obj$Collection_date, add.legend=T)
plotFeature2(mds.coor, values = seu.obj$Cell_line, add.legend = T)
dev.off()

for(x in c("hs", "cs")){
  sig.genes <- date.marker.list[[x]]
  seu.obj <- get(x)
  expr.mat <- as.matrix(seu.obj@assays$RNA@data)
  gene.idx <- which(rownames(expr.mat) %in% sig.genes)
  d.mat <- 1-cor(expr.mat[gene.idx,])
  mds.coor <- cmdscale(d=d.mat, k=2)
  if(x == "hs"){
    legend.pos="bottomright"
  }else{
    legend.pos="topright"
  }
  pdf(paste0("Plot_union_date_marker_based_MDS_for_", x, ".pdf"), height=5, width=10)
  par(mfrow=c(1,2))
  plotFeature2(mds.coor, values = seu.obj$Collection_date, add.legend=T, legend.pos = legend.pos)
  plotFeature2(mds.coor, values = seu.obj$Cell_line, add.legend = T, legend.pos = legend.pos)
  dev.off()
  
}


bg.list <- list()
for(x in c("hs", "cs")){
  seu.obj <- get(x)
  expr.mat <- getAveExpr(seu.obj=seu.obj, feature.to.calc = "Collection_date", colname.prefix = NULL)
  bg.list[[x]] <- rownames(expr.mat)[which(apply(expr.mat, 1, max)>0.05)]
}
bg.genes <- intersect(bg.list[["hs"]], bg.list[["cs"]])

n1 <- length(intersect(intersect(date.marker.list[["hs"]], date.marker.list[["cs"]]), bg.genes))
n2 <- length(intersect(setdiff(date.marker.list[["hs"]], date.marker.list[["cs"]]), bg.genes))
n3 <- length(intersect(setdiff(date.marker.list[["cs"]], date.marker.list[["hs"]]), bg.genes))
n4 <- length(intersect(setdiff(rownames(spheroid), union(date.marker.list[["cs"]], date.marker.list[["hs"]])), bg.genes))
fisher.test(matrix(c(n1,n2,n3,n4), c(2,2)), alternative = "g")
# get average expression of each collection of each species
spheroid$Date_by_species <- paste(spheroid$Collection_date, spheroid$Species, sep="_")
table(spheroid$Date_by_species)
ave.expr.mat <- getAveExpr(seu.obj=spheroid, feature.to.calc = "Date_by_species", colname.prefix = NULL)
saveRDS(ave.expr.mat, file="Dat_collection_by_species_average_expr.rds")
saveRDS(spheroid, file="../Res_spheroid_with_UMAP.rds")
head(ave.expr.mat)
human.ave.expr <- ave.expr.mat[sig.genes,grep("Human", colnames(ave.expr.mat))]
dim(human.ave.expr)
chimp.ave.expr <- ave.expr.mat[sig.genes,grep("Chimp", colnames(ave.expr.mat))]
cor.mat <- cor(human.ave.expr, chimp.ave.expr, method="spearman")
cor.mat
rownames(cor.mat)[apply(cor.mat, 2, which.max)]
library(gplots)
pdf("Plot_heatmap_SCC_spheroid_date_by_species_marker_interset_based.pdf")
heatmap.2(cor.mat, trace="none",main="", density.info="none",
          dendrogram="none", Rowv=FALSE, Colv=FALSE, scale="none", key=TRUE, 
          col=darkBlue2Red, cexRow=0.8, cexCol = 0.8)
dev.off()
