setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38")
tHIO_n_fetal <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/Res_tHIO_and_fetal_8-19_PCW_hg38_no_batch_correction.rds")
fetal_human <- subset(tHIO_n_fetal, Tissue=="Fetal primary")
fetal_human <- FindVariableFeatures(object = fetal_human, selection.method = "vst", nfeatures = 3000)
fetal_human <- ScaleData(object = fetal_human, verbose = T)
fetal_human <- RunPCA(object = fetal_human, features = VariableFeatures(fetal_human), verbose = F, npcs = 50)
usefulPCs <- 1:20
fetal_human <- FindNeighbors(object = fetal_human, dims = usefulPCs)
fetal_human <- FindClusters(object = fetal_human, resolution = 1)
fetal_human <- RunUMAP(object = fetal_human, dims = usefulPCs)
# CSS integration
fetal_human <- cluster_sim_spectrum(fetal_human, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
fetal_human <- RunUMAP(fetal_human, reduction = "css", dims = 1:ncol(fetal_human@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
fetal_human <- FindNeighbors(object = fetal_human, reduction = "css", dims = 1:ncol(fetal_human@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
saveRDS(fetal_human, file="Res_fetal_human_duodenum_integrated_with_CSS.rds")

plotFeature(seu.obj=fetal_human,
            dr="umap_css",
            genes.to.plot = c("EPCAM","COL1A2","PTPRC","PECAM1","PRPH","PHOX2B","HBA1","ALB","AFP"),
            plot.name = "Plot_UMAP_CSS_cell_class_markers.png")
DimPlot(fetal_human, group.by = "RNA_snn_res.1", label=T, reduction = "umap_css")

cell_class_anno <- list(
  "Epithelial"=c(5,7,19,24,26),
  "Mesenchymal"=c(0,1,2,3,4,11,12,13,14,17,21,30,31),
  "Immune"=c(6,8,16,18,20,22,23,28,29,33),
  "Endothelial"=c(10,25),
  "Neural"=c(9,15,27,32)
)
cl_vec <- fetal_human$RNA_snn_res.1
id_vec <- rep(NA, length(cl_vec))
for(x in names(cell_class_anno)){
  id_vec[which(cl_vec%in%cell_class_anno[[x]])] <- x
}
fetal_human$Major_cell_type <- id_vec
DimPlot(fetal_human, group.by = "Major_cell_type", reduction = "umap_css", label=T)+NoLegend()
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fetal_human <- CellCycleScoring(fetal_human, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
DimPlot(fetal_human, group.by = "Phase", reduction = "umap_css")
saveRDS(fetal_human, file="Res_fetal_human_duodenum_integrated_with_CSS.rds")

# extract epithelial cells
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/epi")
fetal_human_epi <- subset(fetal_human, Major_cell_type=="Epithelial")
DimPlot(fetal_human_epi, group.by = "RNA_snn_res.1", reduction = "umap_css", label=T)+NoLegend()+DimPlot(fetal_human_epi, group.by = "Phase", reduction = "umap_css")
fetal_human_epi <- FindVariableFeatures(object = fetal_human_epi, selection.method = "vst", nfeatures = 3000)
fetal_human_epi <- ScaleData(object = fetal_human_epi, verbose = T)
fetal_human_epi <- RunPCA(object = fetal_human_epi, features = VariableFeatures(fetal_human_epi), verbose = F, npcs = 50)
usefulPCs <- 1:20
fetal_human_epi <- FindNeighbors(object = fetal_human_epi, dims = usefulPCs)
fetal_human_epi <- FindClusters(object = fetal_human_epi, resolution = 1)
fetal_human_epi <- RunUMAP(object = fetal_human_epi, dims = usefulPCs)
# CSS integration
fetal_human_epi <- cluster_sim_spectrum(fetal_human_epi, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
fetal_human_epi <- RunUMAP(fetal_human_epi, reduction = "css", dims = 1:ncol(fetal_human_epi@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
fetal_human_epi <- FindNeighbors(object = fetal_human_epi, reduction = "css", dims = 1:ncol(fetal_human_epi@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
DimPlot(fetal_human_epi, reduction = "umap_css", label=T)+NoLegend()

plotFeature(seu.obj=fetal_human_epi,
            dr="umap_css",
            genes.to.plot = c("EPCAM","COL1A2","PTPRC","PECAM1","PRPH","PHOX2B","HBA1","ALB","AFP"),
            plot.name = "Plot_UMAP_CSS_cell_class_markers.png")

# remove potential doublet clusters with mesenchyme
# get genes correlated with known mesenchymal cell type markers
anchor_genes <- c("COL1A2","MFAP4","DCN","COL6A2","COL3A1","PDGFRA")
score <- colSums(fetal_human_epi@assays$RNA@data[anchor_genes,])
plot(seq(length(score)), sort(score))
abline(h=quantile(score, 0.9))
cells <- colnames(fetal_human_epi)[which(score>quantile(score, 0.9))]
cl_vec <- fetal_human_epi$RNA_snn_res.1
doublet_freq <- table(fetal_human_epi@meta.data[cells, "RNA_snn_res.1"])
cl_size <- table(fetal_human_epi$RNA_snn_res.1)
doublet_prop <- doublet_freq/cl_size
doublet_cl <- names(doublet_prop)[which(doublet_prop>0.9)]
doublet_to_remove <- union(cells, colnames(fetal_human_epi)[cl_vec%in%doublet_cl])


par(mfrow=c(1,2))
plotFeature2(coor=Embeddings(fetal_human_epi, reduction = "umap_css"),
             values = score,
             point.order="sorted")
plotFeature2(coor=Embeddings(fetal_human_epi, reduction = "umap_css"),
             values = score,
             point.order="sorted",
             emphasize = which(colnames(fetal_human_epi)%in%doublet_to_remove))


cells <- setdiff(colnames(fetal_human_epi), doublet_to_remove)
epi_sub <- subset(fetal_human_epi, cells=cells)
DimPlot(epi_sub, reduction = "umap_css")
epi <- epi_sub
epi <- FindVariableFeatures(object = epi, selection.method = "vst", nfeatures = 3000)
epi <- ScaleData(object = epi, verbose = T)
epi <- RunPCA(object = epi, features = VariableFeatures(epi), verbose = F, npcs = 50)
usefulPCs <- 1:20
epi <- FindNeighbors(object = epi, dims = usefulPCs)
epi <- FindClusters(object = epi, resolution = 1)
epi <- RunUMAP(object = epi, dims = usefulPCs)
# CSS integration
epi <- cluster_sim_spectrum(epi, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
epi <- RunUMAP(epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
epi <- FindNeighbors(object = epi, reduction = "css", dims = 1:ncol(epi@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
DimPlot(epi, reduction = "umap_css", label=T)+NoLegend()
g1 <- c("CDX2", "CDH1", "EPCAM", "COL1A2", "DCN", "MFAP4", "COL6A2", "COL3A1","PTPRC", "PHOX2B", "CDH5", "PECAM1")
plotFeature(seu.obj=epi,
            dr="umap_css",
            genes.to.plot = g1,
            plot.name = "Plot_UMAP_CSS_cell_class_markers_after_removing_mes_doublet_cls.png")
saveRDS(epi, file="Res_fetal_human_duo_epi_after_removing_potential_mesenchymal_doublets.rds")


g1 <- c("SOX2", "FOXJ1", "TP63", "CLDN18", "LGALS4","CDX2","LGR5", "ASCL2", "OLFM4", "MUC2", "SPINK4", "HES6", "MKI67",
        "NEUROG3","CHGA", "ISL1", "BEST4", "SPIB", "GP2", "CA7", "ICAM2","MUC1","DPP4", "APOA4", "FABP2", "DEFA5", "LYZ", "RGS13")
plot_name <- "Plot_UMAP_CSS_epithelial_cell_type_marker_expression.png"
plotFeature.batch(seu.obj = epi, dr="umap_css", genes.to.plot = g1, nCols = blue.cols, 
                  plot.name =plot_name , col.num = 7, cex = 2)

col_num=3
row_num=1
png("Plot_UMAP_CSS_8-11_PCW_fetal_epi_basic_info.png", height=2000*row_num, width=2000*col_num)
par(mfrow=c(row_num, col_num), mar=c(5,5,15,5))
plotFeature2(Embeddings(epi, reduction = "umap_css"), values = epi$Sample.name, main="Sample name @ CSS", point.order = "random",
             add.legend = T, legend.pos = "bottomleft", legend.cex = 4, cex.main=10, cex=3, lwd = 0.5)
plotFeature2(Embeddings(epi, reduction = "umap_css"), values = epi$RNA_snn_res.1, main="Cluster @ CSS", point.order = "random",
             add.label = T, label.cex = 6, cex.main=10, cex=3, lwd = 0.5)
plotFeature2(Embeddings(epi, reduction = "umap_css"), values = epi$Phase, main="Cell cycle phase @ CSS", point.order = "random",
             add.legend = T, legend.cex = 6, cex.main=10, cex=3, lwd = 0.5)

dev.off()

# extract G1-phase stem cells-to-enterocytes
g1_s2e <- subset(epi, cells=colnames(epi)[which(epi$RNA_snn_res.1%in%c(0,13,1,6,11,8,4,7,5,2,3,9) & epi$Phase=="G1")])
library(destiny)
input <- g1_s2e@reductions$css@cell.embeddings
dm <- DiffusionMap(input, k=20)
coor <- cbind(-dm$DC1, dm$DC2)
colnames(coor) <- c("DC1","DC2")
plotFeature2(coor=Embeddings(g1_s2e, reduction = "umap_css"),
             values = rank(dm$DC1))

plotFeature(seu.obj=g1_s2e, 
            dr="umap_css",
            genes.to.plot = c("LGR5","OLFM4","ASCL2","MKI67","CDK1","TOP2A","FABP1","SI","FABP2","DPP4","APOA1","APOA4"),
            plot.name = "Plot_UMAP_CSS_stem_cell_and_enterocyte_marker_expr_along_g1phase_s2E.png")


age_test_pval <- splineBasedAgeTest(pseudotime.vec=rank(dm$DC1), 
                                    expr.mat=as.matrix(g1_s2e@assays$RNA@data), 
                                    df=6, 
                                    mode="relax", 
                                    core.num=20)
padj <- p.adjust(age_test_pval, method="bonferroni")
age_genes <- names(padj)[which(padj<0.05)]
writeLines(age_genes, con="List_variable_genes_along_stem_cell_to_enterocyte_differentiation.csv")

pt_expr <- getExprByPt(pt.vec=rank(dm$DC1),
                       expr.mat = g1_s2e@assays$RNA@data,
                       mode = "fix.cell.num")
saveRDS(pt_expr, file="Res_G1phase_stemCell2Enterocyte_expr_pseudotime_bin_expr.rds")
pt_expr_fc <- rowMax(pt_expr)-rowMin(pt_expr)
names(pt_expr_fc) <- rownames(g1_s2e)

age_genes_pt_expr <- getExprByPt(pt.vec=rank(dm$DC1),
                                 expr.mat = g1_s2e@assays$RNA@data[age_genes,],
                                 mode = "fix.cell.num")
aa <- getExprByPt(pt.vec=rank(dm$DC1),
                  expr.mat = g1_s2e@assays$RNA@data[age_genes,],
                  mode = "fix.cell.num",
                  return.idx = T)
saveRDS(aa, file="Res_G1phase_stemCell2Enterocyte_pt_dependent_gene_pseudotime_bin_expr.rds")
age_genes_pt_expr <- aa$expr.mat
bin_pt <- sapply(sort(unique(aa$cell.idx.vec)), function(i){
  median(rank(dm$DC1)[which(aa$cell.idx.vec==i)])
})

plot(bin_pt, pt_expr["DPP4",])
hc <- hclust(as.dist(1-cor(t(age_genes_pt_expr))),
             method="ward.D2")
saveRDS(hc, file="Res_G1phase_stemCell2Enterocyte_pt_dependent_gene_module_hc.rds")

age_genes_module_id <- cutree(hc, h=8)
age_genes_module <- lapply(sort(unique(age_genes_module_id)),
                           function(i){names(age_genes_module_id[which(age_genes_module_id==i)])})
age_genes_module_expr <- plotClusterExprProfile(expr=age_genes_pt_expr, 
                                                time.vec=seq(ncol(age_genes_pt_expr)), 
                                                group.vec=rep("Human", ncol(age_genes_pt_expr)), 
                                                cluster.vec=age_genes_module_id, 
                                                group.cols="#31a354", 
                                                return.value=T, 
                                                to.plot=T, 
                                                plot.name="Plot_stemCell2Enterocyte_pt_dependent_gene_module_average_expr_profile.pdf", 
                                                add.legend=F, 
                                                legend.pos="topleft", 
                                                cex.legend=2, 
                                                col.num=4, 
                                                border.do.smooth=T, 
                                                mean.do.smooth=T, 
                                                df=8)
saveRDS(age_genes_module_expr, file="Res_G1phase_stemCell2Enterocyte_pt_dependent_gene_module_expr.rds")
age_genes_module_id[intersect(c("LGR5","OLFM4","ASCL2","MKI67","CDK1","TOP2A","FABP1","SI","FABP2","DPP4","APOA1","APOA4"),
                              names(age_genes_module_id))]
m6_genes <- names(age_genes_module_id)[age_genes_module_id==6]
genes <- m6_genes[order(pt_expr_fc[m6_genes], decreasing = T)[1:20]]
plotFeature(seu.obj=g1_s2e, 
            dr="umap_css",
            genes.to.plot = genes,
            plot.name = "Plot_UMAP_CSS_s2E_pt_gene_C6_expr_along_g1phase_s2E.png")


bin_hc <- hclust(as.dist(1-cor(age_genes_pt_expr)),
                 method="ward.D2")

cor_vec <- cor(t(as.matrix(g1_s2e@assays$RNA@data[which(pt_expr_fc>0),])), g1_s2e@assays$RNA@data["GSTA1",])[,1]
genes <- names(cor_vec)[order(cor_vec, decreasing = T)[1:20]]
plotFeature(seu.obj=g1_s2e, 
            dr="umap_css",
            genes.to.plot = genes,
            plot.name = "Plot_UMAP_CSS_s2E_GSTA1_correlated_gene_expr_along_g1phase_s2E.png")

# GSTA1, GSTA2, AADAC, REEP6,FABP1,RBP2,ANXA13,GSTO1,HADH,FTL,HEBP1
genes <- c("LGR5","OLFM4","ASCL2","MKI67","CDK1","TOP2A","FABP1","SI","FABP2","DPP4","APOA1","APOA4",
           "GSTA1", "GSTA2", "AADAC", "REEP6","FABP1","RBP2","ANXA13","GSTO1","HADH","FTL","HEBP1")
plotFeature(seu.obj=g1_s2e, 
            dr="umap_css",
            genes.to.plot = genes,
            plot.name = "Plot_UMAP_CSS_s2E_selected_gene_expr_along_g1phase_s2E.png")

genes <- c("LGR5","OLFM4","ASCL2",
           "GSTA1", "GSTA2", "AADAC", "REEP6",
           "FABP2","APOA4")
plotFeature(seu.obj=g1_s2e, 
            dr="umap_css",
            genes.to.plot = genes,
            plot.name = "Plot_UMAP_CSS_s2E_final_set_selected_gene_expr_along_g1phase_s2E.png",
            col.num = 3)

gene_cols <- setNames(rep(c("#d73027","#fee08b","#1a9850"),c(3,4,2)),
                      genes)
scaled_expr <- scale(t(pt_expr[genes,]))
plot(seq(nrow(scaled_expr)), scaled_expr[,1], type="n", ylab="Scaled expr", xlab="Pt bin", bty="n",
     ylim=c(min(scaled_expr), max(scaled_expr)))
for(g in colnames(scaled_expr)){
  lines(smooth.spline(seq(nrow(scaled_expr)), scaled_expr[,g], df=8), col=gene_cols[g])
}

s2e_features <- list(
  "Stem_cell"=c("LGR5","OLFM4","ASCL2"),
  "Early_enterocyte"=c("GSTA1", "GSTA2", "AADAC", "REEP6"),
  "Enterocyte"=c("FABP2","APOA4")
)
g1_s2e <- AddModuleScore(
  object = g1_s2e,
  features = s2e_features,
  ctrl = 5,
  name = 's2E_Features'
)
score_mat <- sapply(names(s2e_features), function(x){
  scale(colSums(g1_s2e@assays$RNA@data[s2e_features[[x]],]))
})
rownames(score_mat) <- colnames(g1_s2e)
id <- colnames(score_mat)[apply(score_mat, 1, which.max)]
g1_s2e$s2E_ident_scaled_expr_based <- id
p1 <- DimPlot(g1_s2e, reduction = "umap_css", group.by = "s2E_ident", shuffle = T, label = T)+NoLegend()
p2 <- DimPlot(g1_s2e, reduction = "umap_css", group.by = "s2E_ident_scaled_expr_based", shuffle = T, label = T)+NoLegend()
p1+p2

s2e <- subset(epi, cells=colnames(epi)[which(epi$RNA_snn_res.1%in%c(0,13,1,6,11,8,4,7,5,2,3,9))])
score_mat <- sapply(names(s2e_features), function(x){
  scale(colSums(s2e@assays$RNA@data[s2e_features[[x]],]))
})
rownames(score_mat) <- colnames(s2e)
id <- colnames(score_mat)[apply(score_mat, 1, which.max)]
s2e$s2E_ident_scaled_expr_based <- id

# the cell state classification result based on scaled expression levels fits better with individual gene feature plots
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


epi$Cell_type <- rep(NA, ncol(epi))
epi@meta.data[colnames(s2e), "Cell_type"] <- s2e$s2E_smoothed_ident_scaled_expr_based
cells <- colnames(epi)[which(epi$Phase=="G2M" & epi$Cell_type=="Stem_cell")]
epi@meta.data[cells, "Cell_type"] <- "Transit_amplifying_stem_cells"

par(mfrow=c(2,2))
plotFeature2(coor=Embeddings(s2e, reduction = "umap_css"),
             values = s2e$Cell_type, 
             point.order = "random")
for(g in c("OLFM4","GSTA1","APOA4")){
  plotFeature2(coor=Embeddings(s2e, reduction = "umap_css"),
               values=s2e@assays$RNA@data[g,],
               main=g,
               point.order = "sorted")
}

# cell type annotation
cell_anno <- list(
  "EEC"=10,
  "Goblet_and_Paneth"=12,
  "Tuft"=15,
  "BEST4+ enterocyte type 1"=14,
  "BEST4+ enterocyte type 2"=16,
  "Potential_doublet_EECs_and_Enterocyte"=17
)
cl_vec <- epi$RNA_snn_res.1
ct_vec <- epi$Cell_type
for(x in names(cell_anno)){
  cl_x <- cell_anno[[x]]
  ct_vec[which(cl_vec %in% cl_x)] <- x
}
epi$Cell_type <- ct_vec
DimPlot(epi, reduction = "umap_css", group.by = "Cell_type", label=T)
ref_expr <- colSums(epi@assays$RNA@data[c("DEFA5","DEFA6"), which(epi$Cell_type=="Goblet_and_Paneth")])
que_expr <- as.matrix(epi@assays$RNA@data[,which(epi$Cell_type=="Goblet_and_Paneth")])
cor_vec <- cor(t(que_expr), ref_expr)[,1]
genes <- names(cor_vec)[order(cor_vec, decreasing = T)[1:12]]
score <- colSums(epi@assays$RNA@data[genes,which(epi$Cell_type=="Goblet_and_Paneth")])
plot(seq(length(score)), sort(score))
cutoff=quantile(score, 0.8)
paneth_cells <- names(score)[which(score>cutoff)]
goblet_cells <- setdiff(colnames(epi)[which(epi$Cell_type=="Goblet_and_Paneth")],
                        paneth_cells)
epi@meta.data[paneth_cells, "Cell_type"] <- "Paneth_cell"
epi@meta.data[goblet_cells, "Cell_type"] <- "Goblet_cell"
DimPlot(epi, reduction = "umap_css", group.by = "Cell_type", label=T)+NoLegend()
saveRDS(epi, file="Res_fetal_8-19_PCW_hg38_epi_with_CSS_integration.rds")


# get the proportion of Lgr5+ cells and Olfm4+ cells in the stem cell to enterocyte population at different cell cycle phases
p1 <- sapply(c("LGR5", "OLFM4"), function(gene){
  positive_cells <- colnames(epi)[which(epi@assays$RNA@data[gene, ]>0)]
  sapply(c("G1","S","G2M"), function(phase){
    sum(epi@meta.data[positive_cells,"Phase"]==phase)/length(positive_cells)
  })
})

df1 <- data.frame("Phase"=rep(rownames(p1), ncol(p1)),
                  "Gene"=rep(colnames(p1), each=nrow(p1)),
                  "Proportion"=as.vector(p1))
pdf("Plot_fetal_human_barplot_LGR5_OLFM4_positive_cell_by_phase.pdf")
plot1 <- ggplot(data=df1, aes(x=Gene, y=Proportion, fill=Phase)) +
  geom_bar(stat="identity", width=0.5, position=position_dodge()) +
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
print(plot1)
dev.off()

