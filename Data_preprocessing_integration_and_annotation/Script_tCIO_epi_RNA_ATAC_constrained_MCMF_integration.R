library(reticulate)
source_python("matching.py")

tCIO_epi_RNA <- readRDS("Res_C7_tCIO_panTro6_epi_3_removing_mesenchyme_immune_doublet_cells.rds")
tCIO_epi_ATAC <- readRDS("Res_C7-tCIO_epi_ATAC_preprocessed.rds")
DefaultAssay(tCIO_epi_ATAC) <- "RNA"
DimPlot(tCIO_epi_ATAC, label = T, group.by = "peaks_snn_res.1")+DimPlot(tCIO_epi_RNA, label = T, group.by = "Cell_type")

# step 1 - identify the ATAC cell populations corresponding to RNA clusters
## three approaches have been tried. For all approaches, scATAC-seq data takes the gene activity assay
### approach 1 - match RNA and ATAC clusters according to cluster average gene activity similarity
### approach 2 - match ATAC cells to RNA clusters according to the similarity between RNA cluster and ATAC cell plus their nearest neighboring cells
### approach 3 - CCA based Seurat label transfer
#### specifically for this tCIO epithelium subset, approach 2 works the best, in the sense that the mapped RNA cell cluster composition is the closest to the RNA dataset
#### however, when I run the analysis on the whole dataset which include epithelial cells and mesenchymal cells, Seurat label transfer seems to work better
g1 <- intersect(VariableFeatures(tCIO_epi_RNA),
                rownames(tCIO_epi_ATAC))

rna_cl_expr <- getAveExpr(seu.obj=tCIO_epi_RNA,
                          feature.to.calc = "Cell_type",
                          colname.prefix = "RNA")
DefaultAssay(tCIO_epi_ATAC) <- "peaks"
tCIO_epi_ATAC <- FindNeighbors(tCIO_epi_ATAC, reduction="lsi", dims=2:20) %>% FindClusters(resolution = 2)
DimPlot(tCIO_epi_ATAC, label=T)
DefaultAssay(tCIO_epi_ATAC) <- "RNA"
atac_cl_expr <- getAveExpr(seu.obj=tCIO_epi_ATAC,
                           feature.to.calc = "peaks_snn_res.2",
                           colname.prefix = "ATAC")
rna_expr <- rna_cl_expr[g1,]
atac_expr <- atac_cl_expr[g1,]
cor_mat <- cor(rna_expr, atac_expr, method="spearman")
atac_cell_expr <- as.matrix(tCIO_epi_ATAC@assays$RNA@data[g1,])
cell2cl_cor_mat <- cor(atac_cell_expr, rna_expr, method="spearman")
nn_mat <- tCIO_epi_ATAC@graphs$peaks_nn
cell2cl_incl_nn_cor_mat <- t(sapply(seq(nrow(nn_mat)), function(i){
  vec <- as.vector(as.matrix(nn_mat[i,]))
  nn_cells <- colnames(nn_mat)[which(vec==1)]
  nn_cor <- colMeans(rbind(cell2cl_cor_mat[nn_cells,], cell2cl_cor_mat[i,]))
}))
rownames(cell2cl_incl_nn_cor_mat) <- rownames(nn_mat)
cell_max_ct <- setNames(colnames(cell2cl_incl_nn_cor_mat)[apply(cell2cl_incl_nn_cor_mat, 1, which.max)],
                        rownames(cell2cl_incl_nn_cor_mat))
tCIO_epi_ATAC$ATAC_cell_mapped_RNA_cl <- cell_max_ct

# step 2 - run Seurat CCA
seurat.rna <- tCIO_epi_RNA
seurat.atac <- tCIO_epi_ATAC
seurat.rna$modality <- "RNA"
seurat.atac$modality <- "ATAC"
selected.hvg <- intersect(VariableFeatures(seurat.rna), rownames(seurat.atac@assays$RNA@data))
DefaultAssay(seurat.atac) <- "RNA"
VariableFeatures(seurat.atac) <- selected.hvg
seurat.atac <- ScaleData(seurat.atac)
anchor_features <- SelectIntegrationFeatures(
  object.list = list('RNA'=seurat.rna, 'ATAC'=seurat.atac),
  nfeatures = length(selected.hvg),
  assay = rep('RNA', 2)
)
seurat.cca <- RunCCA(object1=seurat.rna, object2=seurat.atac, assay1="RNA", assay2="RNA", features = anchor_features, num.cc = 20)
seurat.cca <- RunUMAP(seurat.cca, reduction="cca", dims = 1:20)

rna_cca <- Embeddings(seurat.cca, reduction = "cca")[which(seurat.cca$modality=="RNA"),]
atac_cca <- Embeddings(seurat.cca, reduction = "cca")[which(seurat.cca$modality=="ATAC"),]
# run MCMF on the CCA space on each of the RNA and ATAC matched cell population
for(x in sort(unique(cell_max_ct))){
  rna_cells <- colnames(seurat.rna)[which(seurat.rna$Cell_type==sub("RNA_", "", x))]
  atac_cells <- colnames(seurat.atac)[which(seurat.atac$ATAC_cell_mapped_RNA_cl==x)]
  cost_graph <- get_cost_knn_graph(
    source = rna_cca[rna_cells,],
    target = atac_cca[atac_cells,],
    knn_k = 5,
    knn_n_jobs = 10,
    null_cost_percentile = 99,
    capacity_method = "uniform"
  )
  grp_idx <- setNames(data.frame(do.call(cbind, mcmf(cost_graph))+1), c("rna","atac"))
  grp_cell <- data.frame("RNA"=rna_cells[grp_idx[,1]], 
                         "ATAC"=atac_cells[grp_idx[,2]])
  idx <- which(!is.na(grp_cell[,2]))
  matched_cells <- rbind(matched_cells, grp_cell[idx,])
}
matched_cells <- unique(matched_cells)

matched_cells_atac_ref <- c()
for(x in sort(unique(cell_max_ct))){
  rna_cells <- colnames(seurat.rna)[which(seurat.rna$Cell_type==sub("RNA_", "", x))]
  atac_cells <- colnames(seurat.atac)[which(seurat.atac$ATAC_cell_mapped_RNA_cl==x)]
  cost_graph <- get_cost_knn_graph(
    source = atac_cca[atac_cells,],
    target = rna_cca[rna_cells,],
    knn_k = 5,
    knn_n_jobs = 10,
    null_cost_percentile = 99,
    capacity_method = "uniform"
  )
  grp_idx <- setNames(data.frame(do.call(cbind, mcmf(cost_graph))+1), c("atac","rna"))
  grp_cell <- data.frame("RNA"=rna_cells[grp_idx[,"rna"]], 
                         "ATAC"=atac_cells[grp_idx[,"atac"]])
  idx <- which(!is.na(grp_cell[,"RNA"]))
  matched_cells_atac_ref <- rbind(matched_cells_atac_ref, grp_cell[idx,])
}
matched_cells_atac_ref <- unique(matched_cells_atac_ref)
matched_cells_combined <- as.matrix(unique(rbind(matched_cells, matched_cells_atac_ref)))

# since there are more RNA samples than ATAC samples, so here I don't averaging the RNA cells, but only averaging the ATAC cells that are linked to the same RNA cell
matched_RNA_cells <- sort(unique(matched_cells_combined[,"RNA"]))
pseudo_cell_idx <- rep(NA, nrow(matched_cells_combined))
for(i in seq(length(matched_RNA_cells))){
  idx <- which(matched_cells_combined[,"RNA"]==matched_RNA_cells[i])
  pseudo_cell_idx[idx] <- paste0("pseudocell_",i)
}
matched_cells_atac_combined <- data.frame("RNA"=matched_cells_combined[,"RNA"],
                                          "ATAC"=matched_cells_combined[,"ATAC"],
                                          "Pseudocell_idx"=pseudo_cell_idx,
                                          stringsAsFactors = F)
saveRDS(matched_cells_atac_combined, file="Res_matched_cells_only_averaging_ATAC_cells.rds")
atac_group_embedding <- t(sapply(sort(unique(matched_cells_atac_combined$Pseudocell_idx)), function(i){
  cells <- matched_cells_atac_combined$ATAC[which(matched_cells_atac_combined$Pseudocell_idx==i)]
  if(length(cells)==1){
    return(atac_embeddings[cells,])
  }else{
    colMeans(atac_embeddings[cells,])
  }
}))
rownames(atac_group_embedding) <- sort(unique(matched_cells_atac_combined$Pseudocell_idx))
saveRDS(atac_group_embedding, file="Res_atac_group_embedding.rds")
mat <- unique(matched_cells_atac_combined[,c(1,3)])
mat <- mat[order(mat$Pseudocell_idx),]

meta_cell_info <- data.frame("Pseudocel_idx"=mat$Pseudocell_idx,
                             "RNA_cell"=mat$RNA,
                             "RNA_UMAP"=Embeddings(seurat.rna.matched, reduction = "umap")[mat$RNA,],
                             "RNA_cell_ident"=seurat.rna.matched@meta.data[mat$RNA, "Cell_type"],
                             "ATAC_UMAP"=atac_group_embedding,
                             stringsAsFactors = F)
saveRDS(meta_cell_info, file="Res_meta_cell_info.rds")
meta_cell_info <- readRDS("Res_meta_cell_info.rds")
matched_cells <- readRDS("Res_matched_cells_only_averaging_ATAC_cells.rds")
rownames(meta_cell_info) <- meta_cell_info$Pseudocel_idx
df <- data.frame(matched_cells, meta_cell_info[matched_cells$Pseudocell_idx, -c(1,2)], stringsAsFactors = F)
meta_cell_info <- df
saveRDS(meta_cell_info, 
        file="Res_meta_cell_info.rds")
