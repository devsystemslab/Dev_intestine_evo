# this is script is used to integrate the scRNA-seq data of all cells of tHIO, fetal stem cell derived enteroids and fetal tissues
setwd("~/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class")
library(Seurat)
library(Signac)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# prefix of each sample 
# HA-2 -> H9
# HA-5 -> iPSC
# fHA-1 -> fetal_Jason
# fHA-2 -> fetal_Umut
# fHA-3 -> fetal_enteroid

# get the tHIO all cell class data
HIO <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/HIO_and_tHIO/scRNA-seq/hg38/remove_low_quality_cluster/remove_doublets/Res_in_vitro_and_transplanted_human_organoids_with_CSS_and_MNN_and_Harmony_integration.rds")
tHIO_rna <- subset(HIO, Tissue=="transplanted")
DefaultAssay(tHIO_rna) <- "RNA"
tHIO_rna <- DietSeurat(
  object = tHIO_rna,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)
rm(HIO)


# get the public fetal tissue all cell class data
fetal_duo <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/Res_fetal_human_duodenum_integrated_with_CSS.rds")
DefaultAssay(fetal_duo) <- "RNA"
fetal_duo <- DietSeurat(
  object = fetal_duo,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)


# get the Jason fetal multiome all cell class data
fetal_crypt_Jason <- readRDS("~/ihb-intestine-evo/lukas_area/fetal_crypt_multiome/data.set.joined.umap.v2.rds")
DefaultAssay(fetal_crypt_Jason) <- "RNA"
fetal_crypt_Jason <- DietSeurat(
  object = fetal_crypt_Jason,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)
fetal_crypt_Jason <- RenameCells(fetal_crypt_Jason, add.cell.id = "fHA-1")

# get fetal enteroid data
fetal_enteroid <- readRDS("/projects/site/pred/ihb-intestine-evo/lukas_area/fetal_enteroid/processed/fetal_enteroid_multiome.processed.umap.annotated.rds")
fetal_enteroid <- RenameCells(fetal_enteroid, add.cell.id = "fHA-3")
DefaultAssay(fetal_enteroid) <- "RNA"
fetal_enteroid <- DietSeurat(
  object = fetal_enteroid,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)
## get fetal tissue data from Umut
fetal_umut <- readRDS("/projects/site/pred/ihb-intestine-evo/lukas_area/fetal_tissue/processed/fSI_17PCW_F_230515.processed.annotated.rds")
fetal_umut <- RenameCells(fetal_umut, add.cell.id = "fHA-2")
DefaultAssay(fetal_umut) <- "RNA"
fetal_umut <- DietSeurat(
  object = fetal_umut,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)


combined_rna <- merge(x=tHIO_rna,
                      y=list("fetal_duo"=fetal_duo,
                             "fetal_jason"=fetal_crypt_Jason,
                             "fetal_umut"=fetal_umut,
                             "fetal_enteroid"=fetal_enteroid))
combined_rna@meta.data[colnames(fetal_umut), "Sample.name"] <- "UmutSI-W17"
combined_rna@meta.data[colnames(fetal_enteroid), "Sample.name"] <- "fEnteroid-1"

# remove confounding genes before normalization
dir <- "/projects/site/pred/ihb-intestine-evo/Annotation/confound_genes/" 
all.files <- list.files(dir, pattern=".txt")
confound.genes <- lapply(seq(length(all.files)), function(i){
  file <- paste0(dir, all.files[i])
  g <- readLines(file)
  return(g)
})
names(confound.genes) <- c("Cell_cycle", "Experiment_induced", "HB", "MT", "Ribosome", "Sex_chr") 
all.confound.genes <- unique(unlist(confound.genes))

genes.to.remove <- unique(c(confound.genes[["MT"]], confound.genes[["Ribosome"]], confound.genes[["Sex_chr"]]))
counts <- combined_rna@assays$RNA@counts
counts <- counts[setdiff(rownames(counts), genes.to.remove),]
meta <- combined_rna@meta.data
fetal_rna <- CreateSeuratObject(counts=counts, meta=meta)
fetal_rna <- NormalizeData(object = fetal_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
fetal_rna$Sample.name[which(is.na(fetal_rna$Sample.name))] <- "HT−multi1−hDuo−80d"
saveRDS(fetal_rna, file="Dat_merged_fetal_tissue_enteroid_and_tHIO_scRNA-seq_data.rds")

# CCA integration
seu_obj_list <- SplitObject(fetal_rna, split.by = "Sample.name")

seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seu_obj_list)
saveRDS(features, file="Res_features_for_fetal_all_cell_class_integration.rds")
anchors <- FindIntegrationAnchors(object.list = seu_obj_list, anchor.features = features)
# this command creates an 'integrated' data assay
fetal_rna.cca <- IntegrateData(anchorset = anchors)
DefaultAssay(fetal_rna.cca) <- "integrated"
fetal_rna.cca <- ScaleData(fetal_rna.cca, verbose = FALSE)
fetal_rna.cca <- RunPCA(fetal_rna.cca, npcs = 30, verbose = FALSE)
fetal_rna.cca <- RunUMAP(fetal_rna.cca, reduction = "pca", dims = 1:30)
fetal_rna.cca <- FindNeighbors(fetal_rna.cca, reduction = "pca", dims = 1:30)
fetal_rna.cca <- FindClusters(fetal_rna.cca, resolution = 0.5)
saveRDS(fetal_rna.cca, file="Res_fetal_tissue_enteroid_and_tHIO_all_cell_class_scRNA-seq_with_CCA_integration_corrected.rds")



p1 <- SCpubr::do_DimPlot(fetal_rna.cca, 
                   reduction = "umap",
                   group.by = "Sample.name", 
                   #colors.use = epi.ct.cols, 
                   pt.size = 6, legend.icon.size = 15, font.size = 50)
png(paste0("Plot_UMAP_CCA_human_intestine_tissue_and_organoid_all_cell_class_RNA_colored_by_sample_including_off_target_cells.png"), height = 2000, width=2000)
print(p1)
dev.off()  

# extract the CCA-corrected-expression-based PCA embeddings
pca_coor <- fetal_rna.cca@reductions$pca@cell.embeddings
saveRDS(pca_coor, file="Res_CCA_corrected_expression_based_PCA_embedding.rds")
# extract the UMAP embeddings
umap_coor <- fetal_rna.cca@reductions$umap@cell.embeddings
saveRDS(umap_coor, file="Res_CCA_integrated_UMAP_embedding.rds")
fetal_rna[["CCA_PCA"]] <- CreateDimReducObject(embeddings = pca_coor, key="CCAPCA_",assay = "RNA")
fetal_rna[["cca_umap"]] <- CreateDimReducObject(embeddings = umap_coor, key="CCAUMAP_",assay = "RNA")

# sort the fetal tissues according to their ages
fetal_rna$Source <- fetal_rna$Tissue
fetal_rna$Tissue[which(fetal_rna$Sample.name%in%c("HT−multi1−hDuo−80d", "UmutSI-W17"))] <- "Fetal primary"
fetal_rna$Source[which(fetal_rna$Sample.name%in%c("HT−multi1−hDuo−80d", "UmutSI-W17"))] <- "Fetal primary"
fetal_rna$Tissue[which(fetal_rna$Sample.name=="fEnteroid-1" | fetal_rna$Tissue=="transplanted")] <- "Fetal organoid"
fetal_rna$Source[which(fetal_rna$Sample.name=="fEnteroid-1")] <- "Fetal enteroid"

fetal_rna$Age.week[which(fetal_rna$Sample.name%in%c("H9-tHIO-EGF-mesentery-W8", "H9-tHIO-EGF-kidney-W8"))] <- 11
fetal_rna$Age.week[which(fetal_rna$Sample.name%in%c("H9-tHIO-ENR-kidney-W8"))] <- 12
fetal_rna$Age.week[which(fetal_rna$Sample.name%in%c("iPSC72.3-tHIO-W10.5"))] <- 13.5
fetal_rna$Age.week[which(fetal_rna$Sample.name=="UmutSI-W17")] <- 17
fetal_rna$Age.week[which(fetal_rna$Sample.name=="HT−multi1−hDuo−80d")] <- round(80/7)
saveRDS(fetal_rna, file="Dat_merged_fetal_tissue_enteroid_and_tHIO_scRNA-seq_data_with_CCA_integration.rds")


mat <- unique(fetal_rna@meta.data[,c("Sample.name","Age.week", "Tissue", "Source")])
mat
mat <- mat[which(mat[,"Tissue"]=="Fetal primary"),]
values <- mat[order(as.numeric(mat[,"Age.week"])),"Sample.name"]

si_tissue_cols <- setNames(colorRampPalette(c("#feebe2","#fcc5c0","#fa9fb5","#f768a1","#c51b8a","#7a0177"))(length(unique(values))), values)
organoid_cols <- setNames(c("#c7e9b4","#7fcdbb","#2c7fb8","#253494","#F7DC6F"),
                          c(sort(unique(fetal_rna$Sample.name[which(fetal_rna$Source=="transplanted")])), "fEnteroid-1"))
group.cols <- c(si_tissue_cols, organoid_cols)
p1 <- SCpubr::do_DimPlot(fetal_rna, 
                         reduction = "cca_umap",
                         group.by = "Sample.name", 
                         colors.use = group.cols, 
                         pt.size = 4, 
                         legend.icon.size = 10, 
                         font.size = 30)

png("Plot_CCA_UMAP_fetal_tissue_and_organoid_all_cell_class_including_off_target_cells.png", height=2000, width=2000)
p1
dev.off()

# add tHIO cell type annotation
HIO_annotated <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_tHIO/RNA_ATAC_integration/Res_HIO_rna_with_cell_type_annotation.rds")
tHIO_cells_annotated <- colnames(HIO_annotated)[which(HIO_annotated$Tissue=="transplanted")] 
length(tHIO_cells_annotated)
sum(fetal_rna$Source=="transplanted")
fetal_rna$Cell_type <- NA
fetal_rna@meta.data[tHIO_cells_annotated, "Cell_type"] <- HIO_annotated@meta.data[tHIO_cells_annotated, "Cell_type"]
fetal_rna@meta.data[tHIO_cells_annotated, "Major_cell_type"] <- HIO_annotated@meta.data[tHIO_cells_annotated, "Cell_class"]
# add fetal_Jason cell type annotation
idx1 <- which(!is.na(fetal_rna$cell.type))
fetal_rna@meta.data[idx1,"Cell_type"] <- fetal_rna$cell.type[idx1]
# add fetal_Umut and fetal_enteroid cell type annotation
idx1 <- which(!is.na(fetal_rna$cell_type_initial))
fetal_rna@meta.data[idx1,"Cell_type"] <- fetal_rna$cell_type_initial[idx1]
# load hg19 fetal duodenum data with cell type annotation
fetal_duo_hg19 <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/Merscope_gene_panel/fetal_duodenum/Res_fetal_human_duodenum_with_CSS_integration_and_cell_type_annotation.rds")
fetal_duo <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tHIO_and_fetal_primary_RNA/fetal_duo_hg38/Res_fetal_human_duodenum_integrated_with_CSS.rds")
DefaultAssay(fetal_duo) <- "RNA"
fetal_duo <- DietSeurat(
   object = fetal_duo,
   counts = TRUE,
   data = TRUE,
   scale.data = FALSE,
   features = NULL,
   assays = "RNA",
   dimreducs = NULL,
   graphs = NULL,
   misc = FALSE
 )
# calculate the transcriptome similarity between hg38 fetal cells and hg19-based fetal cell types, and assign the most similar cell type
features <- intersect(VariableFeatures(fetal_duo_hg19), rownames(fetal_duo))
duo_hg19_cell_type_expr <- getAveExpr(seu.obj = fetal_duo_hg19, feature.to.calc = "Cell_type", genes = features, colname.prefix = NULL)
que_expr <- as.matrix(fetal_duo@assays$RNA@data[features,])
cor_mat <- cor(que_expr, duo_hg19_cell_type_expr)
cell_type_vec <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)], rownames(cor_mat))
fetal_rna@meta.data[names(cell_type_vec), "Cell_type"] <- cell_type_vec
cells <- colnames(fetal_rna)[which(!is.na(fetal_rna$Cell_type))]
fetal_rna <- subset(fetal_rna, cells = cells)

fetal_rna$Major_cell_type[which(fetal_rna$Sample.name=="fEnteroid-1")] <- "Epithelial"
ct <- sort(unique(fetal_rna$Cell_type[which(is.na(fetal_rna$Major_cell_type))]))
cell_class_list <- list(
  "Epithelial"=c("APOA1/FABP2+ mature Enterocytes",
                 "Enterocytes",
                 "GDA+ 'Enterocytes'",
                 "SI+ precursor Enterocytes",
                 "Colonocytes",
                 "Basal-like cell",
                 "BEST4+ cells",
                 "BEST4+_epithelium",
                 "CHGA/RIMBP2+ Enteroendocrine cells",
                 "ISL1+ mature Enteroendocrine cells",
                 "EECs",
                 "Goblet cells",
                 "MUC2/FCGBP+ Goblet cells",
                 "LGR5+ Stem cells",
                 "Stem cells",
                 "Secr. progenitors",
                 "TA cells",
                 "Tuft cells",
                 "Enterocyte",
                 "Enterocyte_precursor",
                 "Enteroendocrine",
                 "Goblet",
                 "M_cell",
                 "Stem_cell",
                 "Tuft"
  ),
  "Mesenchymal"=c("CD81_high",
                  "DCN+ Mesoderm",
                  "DIAPH3/TOP2A+ Cycling Mesoderm",
                  "FREM1+ Stromal 3 cells",
                  "HPSE2+ Mesoderm",
                  "Mesotheial cells",
                  "PDGFRA+ Mesoderm",
                  "Pericytes",
                  "PRKG1/ADRA1A+ Pericytes",
                  "SGCD+ Mesoderm",
                  "Smooth-muscle like cells",
                  "Stromal 1/3 cells",
                  "Stromal 2 cells",
                  "EBF2+ cells",
                  "CCL19+",
                  "CD81_high",
                  "CRABP1+",
                  "GDF10+",
                  "GPX3+_villus-core",
                  "ICC",
                  "KCNN3+",
                  "LXN+_intermediate",
                  "Mesothelial-cell",
                  "NRG1+_subepithelial",
                  "Pericyte",
                  "PITX1+_SM",
                  "Proliferative",
                  "SFRP1+_proliferative",
                  "VSM"
  ),
  "Technical_noise"=c("Doublet",
                      "Inflammed/Stressed cells"),
  "Immune"=c("AOAH+ Immune cells",
             "Monocytes/Macrophages",
             "SKAP1+ Immune cells",
             "T cells",
             "B cell",
             "Macrophage/monocyte 2",
             "T cell/NK cell 2"
  ),
  "Neural"=c("CNGB1+ Enteroendocrine cells",
             "Enteric glia cells",
             "Enteric neurons",
             "CDH19+ Enteric neurons",
             "APOE+ glia",
             "CRYAB+ glia",
             "DBH+ neuron",
             "FBP1+ PEMN",
             "GFRA3+/CRYAB+ glia",
             "GRP+ PIN",
             "ID4+/BNC2+ neuron",
             "INSM1+ glia",
             "NOS1+ PIMN",
             "NPY+ PSVN",
             "Proliferative glia",
             "SST+ PSN"
  ),
  "Endothelial"=c("CALCRL+ Venous EC",
                  "ECs",
                  "Lymphatic ECs",
                  "PODXL+ Arterial EC",
                  "RELN+ Lymphatic EC",
                  "Endothelial"
  )
)

for(x in names(cell_class_list)){
  fetal_rna$Major_cell_type[which(fetal_rna$Cell_type%in%cell_class_list[[x]] & fetal_rna$Batch=="Batch_1")] <- x
}
saveRDS(fetal_rna, file="Res_fetal_rna_with_updated_cell_class_annotation.rds")



Major_cell_type_vec <- fetal_rna$Major_cell_type
ct_vec <- fetal_rna$Cell_type
for(x in names(cell_class_list)){
  Major_cell_type_vec[which(ct_vec%in%cell_class_list[[x]])] <- x
}
fetal_rna$Major_cell_type <- Major_cell_type_vec
SCpubr::do_DimPlot(fetal_rna, reduction = "cca_umap", group.by = "Major_cell_type")
saveRDS(fetal_rna, file="Dat_merged_fetal_tissue_enteroid_and_tHIO_scRNA-seq_data_with_CCA_integration_cell_type_anno_per_sample.rds")

# check where does the potential offtarget cells (both epithelial and mesenchymal) cells locate
setwd("/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only")
# CCA on only fetal primary tissue data, exclude potential doublets
fetal_tissue_rna <- subset(fetal_rna, Tissue=="Fetal primary")
saveRDS(fetal_tissue_rna, file="Dat_fetal_tissue_rna_before_integration.rds")
rm(fetal_rna)
# CCA integration
seu_obj_list <- SplitObject(fetal_tissue_rna, split.by = "Sample.name")

seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seu_obj_list)
saveRDS(features, file="Res_features_for_fetal_tissue_all_cell_class_integration.rds")
anchors <- FindIntegrationAnchors(object.list = seu_obj_list, anchor.features = features)
# this command creates an 'integrated' data assay
fetal_tissue_rna.cca <- IntegrateData(anchorset = anchors)
DefaultAssay(fetal_tissue_rna.cca) <- "integrated"
fetal_tissue_rna.cca <- ScaleData(fetal_tissue_rna.cca, verbose = FALSE)
fetal_tissue_rna.cca <- RunPCA(fetal_tissue_rna.cca, npcs = 30, verbose = FALSE)
fetal_tissue_rna.cca <- RunUMAP(fetal_tissue_rna.cca, reduction = "pca", dims = 1:30)
fetal_tissue_rna.cca <- FindNeighbors(fetal_tissue_rna.cca, reduction = "pca", dims = 1:30)
fetal_tissue_rna.cca <- FindClusters(fetal_tissue_rna.cca, resolution = 0.5)
saveRDS(fetal_tissue_rna.cca, file="Res_fetal_tissue_all_cell_class_scRNA-seq_with_CCA_integration_corrected.rds")

# update the fetal duo cells annotation
# instead of based on correlation, try label transfer to transfer the cell type annotation from the hg19-based duo annotation to the current dataset
# load hg19 fetal duodenum data with cell type annotation
fetal_duo_hg19 <- readRDS("/Volumes/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/Merscope_gene_panel/fetal_duodenum/Res_fetal_human_duodenum_with_CSS_integration_and_cell_type_annotation.rds")

public_duo_samples <- c("HT-260-hDuo-59d",
                        "HT-274-hDuo-72d",
                        "HT-234-hDuo-80d",
                        "HT-184-hDuo-85d",
                        "HT-228-hDuo-101d",
                        "HT-158-hDuo-122d",
                        "HT-172-hDuo-127d",
                        "HT-236-hDuo-132d" )
cells <- colnames(fetal_rna)[which(fetal_rna$Sample.name%in%public_duo_samples)]
fetal_duo <- subset(fetal_rna, cells=cells)

anchors <- FindTransferAnchors(reference = fetal_duo_hg19, query = fetal_duo, 
                               dims = 1:30)
saveRDS(anchors, file="tmp_anchors.rds")
predictions <- TransferData(anchorset = anchors, refdata = fetal_duo_hg19$Cell_type, 
                            dims = 1:30)
saveRDS(predictions, file="Res_fetal_duo_hg19_to_hg38_predictions.rds")
fetal_duo$label_transfer_predicted_cell_type <- predictions$predicted.id
saveRDS(fetal_duo, file="tmp_fetal_duo_with_predicted_id.rds")
# when I updated the cell type annotation with label transfer on the public scRNA-seq dataset, I have already removed the potential doublets based on the transcriptome similarity


fetal_rna <- readRDS("/Volumes/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/Res_fetal_tissue_all_cell_class_scRNA-seq_with_CCA_integration_and_complete_data.rds")
SCpubr::do_DimPlot(fetal_rna, group.by = "Cell_type")
cells <- intersect(colnames(fetal_duo), colnames(fetal_rna))
fetal_rna@meta.data[cells, "Cell_type"] <- fetal_duo@meta.data[cells, "label_transfer_predicted_cell_type"]

# update cell class annotation
cell_class_list <- list(
  "Epithelial"=c("APOA1/FABP2+ mature Enterocytes",
                 "Enterocytes",
                 "GDA+ 'Enterocytes'",
                 "SI+ precursor Enterocytes",
                 "Colonocytes",
                 "Basal-like cell",
                 "BEST4+ cells",
                 "BEST4+_epithelium",
                 "CHGA/RIMBP2+ Enteroendocrine cells",
                 "ISL1+ mature Enteroendocrine cells",
                 "EECs",
                 "Goblet cells",
                 "MUC2/FCGBP+ Goblet cells",
                 "LGR5+ Stem cells",
                 "Stem cells",
                 "Secr. progenitors",
                 "TA cells",
                 "Tuft cells",
                 "Enterocyte",
                 "Enterocyte_precursor",
                 "Enteroendocrine",
                 "Goblet",
                 "M_cell",
                 "Stem_cell",
                 "Tuft"
  ),
  "Mesenchymal"=c("CD81_high",
                  "DCN+ Mesoderm",
                  "DIAPH3/TOP2A+ Cycling Mesoderm",
                  "FREM1+ Stromal 3 cells",
                  "HPSE2+ Mesoderm",
                  "Mesotheial cells",
                  "PDGFRA+ Mesoderm",
                  "Pericytes",
                  "PRKG1/ADRA1A+ Pericytes",
                  "Smooth-muscle like cells",
                  "Stromal 1/3 cells",
                  "Stromal 2 cells",
                  "EBF2+ cells",
                  "CCL19+",
                  "CD81_high",
                  "CRABP1+",
                  "GDF10+",
                  "GPX3+_villus-core",
                  "ICC",
                  "KCNN3+",
                  "LXN+_intermediate",
                  "Mesothelial-cell",
                  "NRG1+_subepithelial",
                  "Pericyte",
                  "PITX1+_SM",
                  "Proliferative",
                  "SFRP1+_proliferative",
                  "VSM"
  ),
  "Technical_noise"=c("Doublet",
                      "Inflammed/Stressed cells",
                      "Potential doublet of neuron and mesenchyme",
                      "Potential doublet of glia and mesenchyme",
                      "SGCD+ Mesoderm"
                      ),
  "Immune"=c("AOAH+ Immune cells",
             "Monocytes/Macrophages",
             "SKAP1+ Immune cells",
             "T cells",
             "B cell",
             "Macrophage/monocyte 2",
             "T cell/NK cell 2"
  ),
  "Neural"=c("CNGB1+ Enteroendocrine cells",
             "Enteric glia cells",
             "Enteric neurons",
             "CDH19+ Enteric neurons",
             "APOE+ glia",
             "CRYAB+ glia",
             "DBH+ neuron",
             "FBP1+ PEMN",
             "GFRA3+/CRYAB+ glia",
             "GRP+ PIN",
             "ID4+/BNC2+ neuron",
             "INSM1+ glia",
             "NOS1+ PIMN",
             "NPY+ PSVN",
             "Proliferative glia",
             "SST+ PSN"
  ),
  "Endothelial"=c("CALCRL+ Venous EC",
                  "ECs",
                  "Lymphatic ECs",
                  "PODXL+ Arterial EC",
                  "RELN+ Lymphatic EC",
                  "Endothelial"
  )
)

for(x in names(cell_class_list)){
  fetal_rna$Major_cell_type[which(fetal_rna$Cell_type%in%cell_class_list[[x]])] <- x
}
mat <- unique(fetal_rna@meta.data[,c("Cell_type", "Major_cell_type")])
mat <- mat[order(mat[,2]),]
saveRDS(fetal_rna, file="Res_fetal_rna_with_updated_cell_class_annotation.rds")

# exclude potential doublets
# i.e. those marked as potential doublets, as well as "SGCD+ Mesoderm"
setwd("/Volumes/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet")
cells <- colnames(fetal_rna)[which(fetal_rna$Major_cell_type!="Technical_noise")]
fetal_rna <- subset(fetal_rna, cells=cells)
saveRDS(fetal_rna, file="Res_fetal_tissue_RNA_exclude_potential_doublet.rds")


fetal_rna$Batch <- "Batch_1"
fetal_rna$Batch[which(fetal_rna$Sample.name=="HT−multi1−hDuo−80d")] <- "Batch_2"
fetal_rna$Batch[which(fetal_rna$Sample.name=="UmutSI-W17")] <- "Batch_3"
#saveRDS(fetal_rna, file="Res_fetal_tissue_all_cell_class_scRNA-seq_with_CCA_integration_and_complete_data.rds")
saveRDS(fetal_rna, file="Res_fetal_tissue_RNA_exclude_potential_doublet.rds")


p1 <- SCpubr::do_DimPlot(fetal_rna, 
                         reduction = "cca_umap",
                         group.by = "Major_cell_type", 
                         label = F,
                         pt.size = 6, legend.icon.size = 12, font.size = 30, label.size = 10)+NoLegend()+NoAxes()
png("Plot_UMAP_CCA_human_intestine_tissue_cell_class.png", height = 2000, width=2000)
print(p1)
dev.off()  

fetal_rna$Cell_type_per_batch <- paste(fetal_rna$Cell_type, fetal_rna$Batch, sep="@")
saveRDS(fetal_rna, file="Res_fetal_tissue_RNA_exclude_potential_doublet.rds")

# to link clusters by cluster transcriptome similarity sometimes give confusing results, where some immune cell cluster is most similar to enterocyte, but with CCA-correctedd PCA Euclcidean distance, the linking results make more sense
#source("/Volumes/pred/ihb-intestine-evo/common_script/Script_functions.R")
#cell_type_per_batch_expr <- getAveExpr(seu.obj=fetal_rna, feature.to.calc = "Cell_type_per_batch", colname.prefix = NULL, size.cutoff = 1)
#saveRDS(cell_type_per_batch_expr, file="Dat_cell_type_per_batch_expr.rds")
#
## load the feature gene set for integration
#features <- readRDS("Res_features_for_fetal_tissue_all_cell_class_integration.rds")
#ref_expr <- cell_type_per_batch_expr[,grep("Batch_1", colnames(cell_type_per_batch_expr), value = T)]
#que_expr <- cell_type_per_batch_expr[,grep("Batch_2|Batch_3", colnames(cell_type_per_batch_expr), value = T)]
#cor_mat <- cor(que_expr[features, ], ref_expr[features,])
#colnames(cor_mat) <- sub("@Batch_1","",colnames(cor_mat))
#max_ct <- cbind(colnames(cor_mat)[apply(cor_mat,1,which.max)], rownames(cor_mat))

pca_coor <- fetal_rna@reductions$CCA_PCA@cell.embeddings
ave_pca <- getAveExpr(input.type = "matrix",
                      X=t(pca_coor),
                      y=setNames(fetal_rna$Cell_type_per_batch, colnames(fetal_rna)),
                      colname.prefix = NULL,
                      size.cutoff = 1)
ref_ave_pca <- ave_pca[,grep("Batch_1",colnames(ave_pca))]
que_ave_pca <- ave_pca[,grep("Batch_2|Batch_3",colnames(ave_pca))]


cl2cl <- sapply(seq(ncol(ref_ave_pca)), function(ref_idx){
  sapply(seq(ncol(que_ave_pca)), function(que_idx){
    a <- ref_ave_pca[,ref_idx]
    b <- que_ave_pca[,que_idx]
    sqrt(sum((a - b)^2))
  })
})
colnames(cl2cl) <- colnames(ref_ave_pca)
rownames(cl2cl) <- colnames(que_ave_pca)

matched_cl_que2ref <- cbind(colnames(cl2cl)[apply(cl2cl, 1, which.min)], rownames(cl2cl))
matched_cl_ref2que <- cbind(colnames(cl2cl), rownames(cl2cl)[apply(cl2cl, 2, which.min)])
matched_cl <- unique(rbind(matched_cl_que2ref, matched_cl_ref2que))

graph <- igraph::graph_from_edgelist(matched_cl, directed=FALSE)
grouped.cells <- igraph::components(graph)
group.vec <- grouped.cells$membership
group.list <- lapply(sort(unique(group.vec)), function(i){
  names(group.vec)[which(group.vec==i)]
})
names(group.list) <- c("Macrophage/monocyte",
                       "Enterocyte",
                       "BEST4+ cell",
                       "Endothelial cell",
                       "Enteric glia",
                       "EEC",
                       "Enteric neuron subgroup 1",
                       "Mesenchyme subgroup 1",
                       "Proliferative mesenchyme",
                       "Enteric neuron subgroup 2",
                       "Mesenchyme subgroup 2",
                       "Goblet cell",
                       "Mesenchyme subgroup 3",
                       "Stem cell",
                       "Mesenchyme subgroup 4",
                       "Pericyte",
                       "Early enterocyte",
                       "T cell/NK cell",
                       "Tuft cell")
saveRDS(group.list, file="Res_matched_cell_type_list.rds")
unified_ct_vec <- rep(NA, ncol(fetal_rna))
for(x in names(group.list)){
  unified_ct_vec[fetal_rna$Cell_type_per_batch %in% group.list[[x]]] <- x
}
fetal_rna$Unified_cell_type <- unified_ct_vec
# the B cells are grouped together with Tuft cells, to correct it, manuallly take B cells out as an independent cell type
fetal_rna$Unified_cell_type[which(fetal_rna$Cell_type=="B cell")] <- "B cell"
fetal_rna$Cell_type_per_sample <- paste(fetal_rna$Cell_type, fetal_rna$Sample.name, sep="@")
saveRDS(fetal_rna, file="Res_fetal_rna_with_updated_cell_class_annotation.rds")

mat <- unique(fetal_rna@meta.data[,c("Major_cell_type", "Unified_cell_type")])
mat <- mat[order(mat[,1]),]
epi_cols <- c("#fee5d9","#fcae91","#fb6a4a","#de2d26","#a50f15")
mes_cols <- c("#edf8e9","#bae4b3","#74c476","#31a354","#006d2c")
imm_cl_cols <- setNames(c("#deebf7", "#9ecae1","#3182bd"), mat[which(mat[,1]=="Immune"),2])
neu_cl_cols <- setNames(c("#fde0dd","#fa9fb5","#c51b8a"), mat[which(mat[,1]=="Neural"),2])
epi_cl_cols <- setNames(colorRampPalette(epi_cols)(sum(mat[,1]=="Epithelial")), mat[which(mat[,1]=="Epithelial"),2])
mes_cl_cols <- setNames(colorRampPalette(mes_cols)(sum(mat[,1]=="Mesenchymal")), mat[which(mat[,1]=="Mesenchymal"),2])
endo_cols <- setNames("#8E44AD", "Endothelial cell")
cl_cols <- c(imm_cl_cols, neu_cl_cols, epi_cl_cols, mes_cl_cols, endo_cols)
p1 <- SCpubr::do_DimPlot(fetal_rna, 
                         reduction = "cca_umap",
                         group.by = "Unified_cell_type", 
                         colors.use = cl_cols,
                         label = F,
                         pt.size = 4)+NoLegend()+NoAxes()
png("Plot_UMAP_CCA_human_intestine_tissue_colored_by_unified_cell_type.png", height = 2000, width=2000)
print(p1)
dev.off()  
pdf("Plot_UMAP_CCA_human_intestine_tissue_colored_by_unified_cell_type_label_only.pdf", height=5, width=10)
plotFeature2(coor=Embeddings(fetal_rna, reduction = "cca_umap"),
             values = fetal_rna$Unified_cell_type,
             gCols = cl_cols,
             label.only = T)
dev.off()  


# get dN/dS score per cell based on all genes
dnds_df <- readRDS("/Volumes/pred/ihb-intestine-evo/evo_signature/DN_DS_ratio/data/Dat_human_to_10sp_one2one_orthologs_dnds_ratio.rds")
# get dn/ds score for all genes
setwd("/Volumes/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/exclude_potential_doublet")
source("/Volumes/pred/ihb-intestine-evo/common_script/Script_functions.R")
ave_expr <- getAveExpr(seu.obj=fetal_rna, feature.to.calc = "Cell_type_per_sample", colname.prefix = NULL, size.cutoff = 1)
saveRDS(ave_expr, file="Dat_cell_type_per_sample_average_expr.rds")

shared_genes <- intersect(rownames(dnds_df), rownames(fetal_rna))
mat <- data.frame(unique(fetal_rna@meta.data[,c("Cell_type_per_sample","Major_cell_type")]),stringsAsFactors = F)
score_list <- list()
for(sp in c("primate","mammal")){
  ave_dnds <- setNames(dnds_df[[sp]], rownames(dnds_df))
  weighted_expr <- ave_expr[shared_genes,]*ave_dnds[shared_genes]
  sum_weighted_expr <- colSums(weighted_expr, na.rm=T)
  sum_log_expr <- colSums(ave_expr[shared_genes,])
  expr_weighted_ave_dnds_score <- sum_weighted_expr/sum_log_expr
  res <- lapply(sort(unique(mat$Major_cell_type)), function(x){
    y=intersect(mat$Cell_type_per_sample[which(mat$Major_cell_type==x)],names(expr_weighted_ave_dnds_score))
    expr_weighted_ave_dnds_score[y]
  })
  names(res) <- sort(unique(mat$Major_cell_type))
  score_list[[sp]] <- res
  
}
saveRDS(score_list, file="Res_average_dnds_for_cell_class_on_cell_type_per_sample.rds")

cell_type_order <- sapply(names(score_list), function(sp){
  res <- score_list[[sp]]
  median_vec <- sapply(seq(length(res)), function(i){mean(res[[i]])})
  names(res)[order(median_vec, decreasing = T)]
})

pdf("Plot_dNdS_ratio.pdf", height=5, width=10)
par(mfrow=c(1,2), mar=c(8,6,4,6))
for(sp in names(score_list)){
  beanplot::beanplot(score_list[[sp]][cell_type_order[,sp]], what=c(1,1,1,0), frame.plot =F, main=sp, las=2, ylab="Expr. weighted average dN/dS")
  sp}
dev.off()


# generate dotplot for major cell type markers
known_markers <- read.table("/Volumes/pred/ihb-intestine-evo/Annotation/cellTypeMarker/Intestine/Table_major_cell_type_markers_from_literature_search.txt", sep="\t", stringsAsFactors = F, head=T)
cell_class <- c("Endothelial",
                "Epithelial", 
                "Immune",
                "Mesenchymal",
                "Neuronal" )
cell_class_markers <- lapply(cell_class, function(x){
  unique(known_markers$Gene[which(known_markers$Used_pan_cell_type_markers==x)])
})
names(cell_class_markers) <- cell_class
p1 <- SCpubr::do_DotPlot(fetal_rna, assay = "RNA", features = cell_class_markers, cluster.idents = TRUE, group.by = "Unified_cell_type")
pdf("Plot_dotplot_cell_class_marker.pdf", width=12)
p1
dev.off()



mat <- unique(fetal_rna@meta.data[,c("Sample.name","Age.week", "Tissue", "Source")])
values <- mat[order(as.numeric(mat[,"Age.week"])),"Sample.name"]
si_tissue_cols <- setNames(colorRampPalette(c("#feebe2","#fcc5c0","#fa9fb5","#f768a1","#c51b8a","#7a0177"))(length(unique(values))), values)
p1 <- SCpubr::do_DimPlot(fetal_rna, 
                         reduction = "cca_umap",
                         group.by = "Sample.name", 
                         colors.use = si_tissue_cols,
                         label = F,
                         pt.size = 6, legend.icon.size = 12, font.size = 30, label.size = 10)+NoLegend()+NoAxes()
png("Plot_UMAP_CCA_human_intestine_tissue_colored_by_sample.png", height = 2000, width=2000)
print(p1)
dev.off()  




# subset to fetal tissue cells with scMultiome data
setwd("/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only")
cells <- colnames(fetal_tissue_rna)[which(fetal_tissue_rna$Sample.name %in% c("HT−multi1−hDuo−80d", "UmutSI-W17"))]
fetal_tissue_multiome <- subset(fetal_tissue_rna, cells=cells)
saveRDS(fetal_tissue_multiome, file="Dat_fetal_tissue_scMultiome_snRNA_seq_portion.rds")




fetal_rna.cca <- readRDS("/Volumes/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_cell_class/fetal_primary_tissue_only/Res_fetal_tissue_all_cell_class_scRNA-seq_with_CCA_integration_corrected.rds")
# extract the CCA-corrected-expression-based PCA embeddings
pca_coor <- fetal_rna.cca@reductions$pca@cell.embeddings
saveRDS(pca_coor, file="Res_CCA_corrected_expression_based_PCA_embedding.rds")
# extract the UMAP embeddings
umap_coor <- fetal_rna.cca@reductions$umap@cell.embeddings
saveRDS(umap_coor, file="Res_CCA_integrated_UMAP_embedding.rds")
fetal_rna[["CCA_PCA"]] <- CreateDimReducObject(embeddings = pca_coor, key="CCAPCA_",assay = "RNA")
fetal_rna[["cca_umap"]] <- CreateDimReducObject(embeddings = umap_coor, key="CCAUMAP_",assay = "RNA")
SCpubr::do_DimPlot(fetal_rna, reduction = "cca_umap", group.by = "Major_cell_type") -> p1
saveRDS(fetal_rna, file="Res_fetal_tissue_all_cell_class_scRNA-seq_with_CCA_integration_and_complete_data.rds")

source("/Volumes/pred/ihb-intestine-evo/common_script/Script_functions.R")
fetal_rna$Cell_type_per_sample <- paste(fetal_rna$Cell_type, fetal_rna$Sample.name, sep="@")
ave_expr <- getAveExpr(seu.obj=fetal_rna, feature.to.calc = "Cell_type", colname.prefix = NULL, size.cutoff = 1)
saveRDS(ave_expr, file="Dat_cell_type_average_expr.rds")

fetal_duo_hg19 <- readRDS("/Volumes/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/Merscope_gene_panel/fetal_duodenum/Res_fetal_human_duodenum_with_CSS_integration_and_cell_type_annotation.rds")
features <- intersect(VariableFeatures(fetal_duo_hg19), rownames(fetal_rna))
duo_hg19_cell_type_expr <- getAveExpr(seu.obj = fetal_duo_hg19, feature.to.calc = "Cell_type", colname.prefix = NULL)
saveRDS(duo_hg19_cell_type_expr, file="Dat_fetal_human_duo_hg19_cell_type_expr.rds")
ref_expr <- duo_hg19_cell_type_expr[features,]
que_expr <- ave_expr[features,]
saveRDS(features, file="Dat_features_to_calculate_cell_type_similarity_to_duo_hg19.rds")

cor_mat <- cor(que_expr, ref_expr)
matched_cell_type <- setNames(colnames(cor_mat)[apply(cor_mat, 1, which.max)], 
                                                rownames(cor_mat))

cells <- colnames(fetal_rna)[which(fetal_rna$Cell_type=="AOAH+ Immune cells")]
SCpubr::do_DimPlot(fetal_rna, reduction = "cca_umap", cells.highlight = cells)

# gene age distribution
gene_age <- readRDS("/projects/site/pred/ihb-intestine-evo/Annotation/GenTree/hg38_Ens95_gene_age_with_gene_name.rds")
# get the gene age of cell type marker combining all samples
marker_mat <- tissue_markers
gene_age_by_cell_type <- sapply(selected_cell_type, function(ct){
  sapply(sort(unique(gene_age$branch)), function(age){
    idx <- which(gene_age$gene_name[which(gene_age$branch==age)]%in%marker_mat$feature[which(marker_mat$group1==ct)])
    length(idx)
  })
  
})

total_n <- sum(gene_age_by_cell_type)
enrichment <- sapply(seq(ncol(gene_age_by_cell_type)), function(j){
  sapply(seq(nrow(gene_age_by_cell_type)), function(i){
    a <- gene_age_by_cell_type[i,j]
    b <- sum(gene_age_by_cell_type[,j])
    c <- sum(gene_age_by_cell_type[i,])
    d <- total_n
    fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")$estimate
  })
})
colnames(enrichment) <- selected_cell_type

# get the gene age of cell type marker of each sample
## include both primary tissues and organoids
## load the cell type marker identification result per sample
per_sample_marker_res <- readRDS("/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/RNA_more_sample/all_epi_cell_type_markers/Res_combined_per_sample_epi_cell_type_marker_res.rds")

# only consider the primary tissue data
samples <- sort(unique(fetal_rna$Sample.name[which(fetal_rna$Tissue=="Fetal_primary")]))
enrichment_mat <- sapply(samples, function(x){
  marker_mat <- per_sample_marker_res[which(per_sample_marker_res$sample==x),]
  gene_age_by_cell_type <- sapply(selected_cell_type, function(ct){
    sapply(sort(unique(gene_age$branch)), function(age){
      idx <- which(gene_age$gene_name[which(gene_age$branch==age)]%in%marker_mat$feature[which(marker_mat$group1==ct)])
      length(idx)
    })
    
  })
  
  total_n <- sum(gene_age_by_cell_type)
  enrichment <- sapply(seq(ncol(gene_age_by_cell_type)), function(j){
    sapply(seq(nrow(gene_age_by_cell_type)), function(i){
      a <- gene_age_by_cell_type[i,j]
      b <- sum(gene_age_by_cell_type[,j])
      c <- sum(gene_age_by_cell_type[i,])
      d <- total_n
      fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")$estimate
    })
  })
  as.vector(enrichment)
})

rownames(enrichment_mat) <- paste(rep(selected_cell_type, each=length(ages)), 
                                  rep(seq(length(ages)), length(selected_cell_type)),
                                  sep="_")

colnames(enrichment_mat) <- samples



dev.off()
pdf("Plot_gene_age_per_cell_type_markers.pdf", height=5*2, width=5*4)
par(mfrow=c(2,4))
for(ct in selected_cell_type){
  input <- enrichment_mat[paste(ct, seq(length(ages)), sep="_"),]
  plot(x=rep(seq(length(ages)), length(samples)),
       y=as.vector(input), 
       pch=16, 
       col=paste0(epi.ct.cols[sub(pattern = "_", replacement = " ", ct)], "30"),
       bty="n",
       xlab="Gene Age",
       ylab="Enrichment",
       main=sub(pattern = "_", replacement = " ", ct))
  
  for(sample in samples){
    vec <- input[, sample]
    if (sum(vec)==0) next
    lines(smooth.spline(x=seq(length(ages)), y=vec, df=8), 
          col=paste0(epi.ct.cols[sub(pattern = "_", replacement = " ", ct)],"30"), 
          lwd=2)
  }
  lines(smooth.spline(x=seq(length(ages)), y=enrichment[,ct], df=8), 
        col=epi.ct.cols[sub(pattern = "_", replacement = " ", ct)], 
        lwd=2)
  
  abline(h=1, lty=2, col="grey", lwd=2)
  
}
dev.off()



samples <- sort(unique(fetal_rna$Sample.name[which(fetal_rna$Tissue=="Fetal_primary")]))
enrichment_mat <- sapply(samples, function(x){
  marker_mat <- per_sample_marker_res[which(per_sample_marker_res$sample==x),]
  gene_age_by_cell_type <- sapply(selected_cell_type, function(ct){
    sapply(sort(unique(gene_age$branch)), function(age){
      idx <- which(gene_age$gene_name[which(gene_age$branch==age)]%in%marker_mat$feature[which(marker_mat$group1==ct)])
      length(idx)
    })
    
  })
  
  total_n <- sum(gene_age_by_cell_type)
  enrichment <- sapply(seq(ncol(gene_age_by_cell_type)), function(j){
    sapply(seq(nrow(gene_age_by_cell_type)), function(i){
      a <- gene_age_by_cell_type[i,j]
      b <- sum(gene_age_by_cell_type[,j])
      c <- sum(gene_age_by_cell_type[i,])
      d <- total_n
      fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")$estimate
    })
  })
  as.vector(enrichment)
})

rownames(enrichment_mat) <- paste(rep(selected_cell_type, each=length(ages)), 
                                  rep(seq(length(ages)), length(selected_cell_type)),
                                  sep="_")

colnames(enrichment_mat) <- samples



dev.off()
par(mfrow=c(2,4))
for(ct in selected_cell_type){
  input <- enrichment_mat[paste(ct, seq(length(ages)), sep="_"),]
  plot(x=rep(seq(length(ages)), length(samples)),
       y=as.vector(input), 
       pch=16, 
       col=paste0(epi.ct.cols[sub(pattern = "_", replacement = " ", ct)], "30"),
       bty="n",
       xlab="Gene Age",
       ylab="Enrichment",
       main=sub(pattern = "_", replacement = " ", ct))
  
  for(sample in samples){
    vec <- input[, sample]
    if (sum(vec)==0) next
    lines(smooth.spline(x=seq(length(ages)), y=vec, df=8), 
          col=paste0(epi.ct.cols[sub(pattern = "_", replacement = " ", ct)],"30"), 
          lwd=2)
  }
  lines(smooth.spline(x=seq(length(ages)), y=enrichment[,ct], df=8), 
        col=epi.ct.cols[sub(pattern = "_", replacement = " ", ct)], 
        lwd=2)
  
  abline(h=1, lty=2, col="grey", lwd=2)
  
}




# harmonize the cell type annotation




# identify the pan cell type marker using presto and DElegate, and plot the average expression patterns 
setwd("/Volumes/pred/ihb-intestine-evo/fetal_human_duo_crypt/GRN/")
fetal_tissue_multiome <- readRDS("/Volumes/pred/ihb-intestine-evo/fetal_human_duo_crypt/GRN/Res_fetal_tissue_all_cell_class_Pando_glm_model_grn.rds")





