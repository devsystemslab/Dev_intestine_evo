setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/consensus_genome_RNA")
library(Seurat)
library(ggplot2)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# get the paths to the data
dir_list <- list(
  "iPSC72.3-tHIO-W10.5"="/projects/site/pred/ihb-intestine-evo/Data_from_Spence_lab/suspension_HIO_derived_tHIO/RNA/processed/6011-AW-1/outs",
  "H9-tHIO-ENR-kidney-W8"="/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/consensus_genome_RNA/tHIO-matrigel-8-weeks",
  "H9-tHIO-EGF-kidney-W8"="/projects/site/pred/ihb-intestine-evo/Data_from_Spence_lab/Sample_517-EH_tHIO_EGF/processed_HnC_consensus/517-EH-1/outs",
  "H9-tHIO-EGF-mesentery-W8"="/projects/site/pred/ihb-intestine-evo/Data_from_Spence_lab/Sample_517-EH_tHIO_EGF/processed_HnC_consensus/517-EH-3/outs",
  "C2-tCIO-W12"="/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/consensus_genome_RNA/C2_tCIO",
  "C7-tCIO-W12"="/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/consensus_genome_RNA/C7_tCIO"
)

# read tIO epithelial data with cell type annotation
tIO_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
df <- data.frame(unique(tIO_rna@meta.data[,c("orig.ident","Sample.name")]),
                  stringsAsFactors = F)


# read the data and create seurat object
count_list <- list()
for(sample in names(dir_list)){
  path <- dir_list[[sample]]
  file <- paste0(path, "/filtered_feature_bc_matrix.h5")
  counts <- Read10X_h5(file)
  id <- df$orig.ident[which(df$Sample.name==sample)] 
  colnames(counts) <- paste(id, colnames(counts), sep = "_")
  shared_cells <- intersect(colnames(tIO_rna), colnames(counts))
  used_counts <- counts[, shared_cells] 
  count_list[[sample]] <- used_counts
}

combined_counts <- do.call('cbind', count_list)
shared_cells <- intersect(colnames(combined_counts), colnames(tIO_rna))
tIO_rna <- subset(tIO_rna, cells=shared_cells)

# generate integrated UMAP cell embedding
# remove confounding genes
dir <- "/projects/site/pred/ihb-intestine-evo/Annotation/confound_genes/" 
all.files <- list.files(dir, pattern=".txt")
confound.genes <- lapply(seq(length(all.files)), function(i){
  file <- paste0(dir, all.files[i])
  g <- readLines(file)
  return(g)
})
names(confound.genes) <- c("Cell_cycle", "Experiment_induced", "HB", "MT", "Ribosome", "Sex_chr") 
genes.to.remove <- unique(c(confound.genes[["MT"]], confound.genes[["Ribosome"]], confound.genes[["Sex_chr"]]))

counts <- combined_counts[setdiff(rownames(counts), genes.to.remove), shared_cells]
tIO_rna[["HnC_consensus_genome_based_RNA"]] <- CreateAssayObject(counts=counts)
DefaultAssay(tIO_rna) <- "HnC_consensus_genome_based_RNA"
tIO_rna <- NormalizeData(object = tIO_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
seu_obj_list <- SplitObject(tIO_rna, split.by = "Sample.name")
# normalize and identify variable features for each dataset independently
seu_obj_list <- lapply(X = seu_obj_list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seu_obj_list)
tIO_rna <- FindVariableFeatures(object = tIO_rna, selection.method = "vst", nfeatures = 3000)
VariableFeatures(tIO_rna) <- features
tIO_rna <- ScaleData(object = tIO_rna, verbose = T)
tIO_rna <- RunPCA(object = tIO_rna, features = VariableFeatures(tIO_rna), verbose = F, npcs = 50)
usefulPCs <- 1:20
tIO_rna <- FindNeighbors(object = tIO_rna, dims = usefulPCs)
tIO_rna <- FindClusters(object = tIO_rna, resolution = 1)
tIO_rna <- RunUMAP(object = tIO_rna, dims = usefulPCs)

# CSS integration
tIO_rna <- simspec::cluster_sim_spectrum(tIO_rna, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE, reduction.name="consensus_CSS")
tIO_rna <- RunUMAP(tIO_rna, reduction = "consensus_CSS", dims = 1:ncol(tIO_rna@reductions$consensus_CSS@cell.embeddings), reduction.name = "consensus_umap_css", reduction.key = "CONUMAPCSS_")
DimPlot(tIO_rna, reduction = "consensus_umap_css", group.by = "Sample.name") + DimPlot(tIO_rna, reduction = "consensus_umap_css", group.by = "subtype_group")
saveRDS(tIO_rna, file="Res_tIO_HnC_consensus_genome_with_CSS_integration.rds")

# re-run DEG analyses
species_ct_expr <- getAveExpr(seu.obj = tIO_rna, feature.to.calc = "subtype_group_per_species", colname.prefix = NULL, assay.type = "HnC_consensus_genome_based_RNA")
saveRDS(species_ct_expr, file="Dat_tIO_epi_subtype_group_per_species_expr.rds")

subtypes <- sort(unique(tIO_rna$subtype_group))
size <- sapply(subtypes, function(x){
  sapply(c("Human","Chimp"), function(y){
    sum(tIO_rna$subtype_group==x & tIO_rna$Species==y)
  })
})
subtypes <- names(which(colSums(size>5)>1))
deg_num <- matrix(NA, ncol=2, nrow=length(subtypes))
rownames(deg_num) <- subtypes
colnames(deg_num) <- c("Human_high", "Chimp_high")
for(ct in subtypes){
  
  cells <- colnames(tIO_rna)[which(tIO_rna$subtype_group==ct)]
  seurat <- subset(tIO_rna, cells=cells)
  idx <- grep(paste0(ct,"@"), colnames(species_ct_expr), fixed = T)
  expressed_genes <- rownames(species_ct_expr)[which(rowSums(species_ct_expr[,idx])>0)]
  
  ## introduce the mapped region probability as numeric covariate for cross species comparison
  species_vec <- as.factor(seurat$Species)
  i=2 
  mapped_region_prob_mat <- seurat@meta.data[,paste(paste0("Intestine_region_L",i), 
                                                    paste0("prob.",c("Colon","Prox_SI")), 
                                                    sep=":")]
  
  registerDoParallel(20)
  species_test_res <- foreach(k=expressed_genes, .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.numeric(as.vector(seurat@assays$HnC_consensus_genome_based_RNA@data[k,]))
    m0 <- lm(e ~ mapped_region_prob_mat[,1]+mapped_region_prob_mat[,2])
    m1 <- lm(e ~ mapped_region_prob_mat[,1]+mapped_region_prob_mat[,2]+species_vec)
    a0 <- anova(m0)
    a1 <- anova(m1)
    p_anova <- anova(m1,m0)$Pr[2]
    p_resi <- pf(a0["Residuals","Mean Sq"] / a1["Residuals","Mean Sq"], 
                 df1 = a0["Residuals","Df"], df2 = a1["Residuals","Df"], lower.tail = F)
    coef <- coef(m1)[4]
    return(c(p_anova, p_resi ,coef))
  }
  stopImplicitCluster()
  rownames(species_test_res) <- expressed_genes
  colnames(species_test_res) <- c("p_ANOVA", "p_Resi", "Coef")
  resi_p_adj <- p.adjust(species_test_res[,2], method = "BH")
  anova_p_adj <- p.adjust(species_test_res[,1], method = "BH")
  species_test_res <- data.frame(species_test_res,
                                 "BH_corrected_p_Resi"=resi_p_adj,
                                 "BH_corrected_p_ANOVA"=anova_p_adj)
  species_test_res$Significant <- species_test_res$BH_corrected_p_ANOVA<0.05 & abs(species_test_res$Coef)>0.1
  human_high_genes <- rownames(species_test_res)[which(species_test_res$Significant & species_test_res$Coef>0)]
  chimp_high_genes <- rownames(species_test_res)[which(species_test_res$Significant & species_test_res$Coef<0)]
  deg_num[ct, "Human_high"] <- length(human_high_genes)
  deg_num[ct, "Chimp_high"] <- length(chimp_high_genes)
  saveRDS(species_test_res, file=paste0("Res_",ct,"_species_test_res.rds"))
}
saveRDS(deg_num, file="Res_tIO_human_chimp_DEG_number_per_cell_type.rds")

enterocyte_res <- readRDS("Res_Enterocyte_species_test_res.rds")


p_list <- lapply(c("SLC5A12", "IL32", "GAS6", "NPY"), function(g){
  SCpubr::do_BoxPlot(tIO_rna, assay = "HnC_consensus_genome_based_RNA", feature = g, group.by = "Unified_cell_type_per_species")
})
library(gridExtra)
do.call("grid.arrange", c(p_list, ncol=4))


# compare the consensus genome based and respective genome based human-chimp logFC in each epithelial cell type of tIOs
## load the respective genome based gene expression levels
expr_own_genome <- getAveExpr(seu.obj = tIO_rna, feature.to.calc = "subtype_group_per_species", colname.prefix = NULL, assay.type = "RNA")
expr_consensus_genome <- readRDS("Dat_tIO_epi_subtype_group_per_species_expr.rds")
saveRDS(expr_own_genome, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_RNA_high_resolution_cell_type_average_expr.rds")
saveRDS(expr_consensus_genome, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_HnC-based_tIO_epi_RNA_high_resolution_cell_type_average_expr.rds")
shared_genes <- intersect(rownames(expr_own_genome), rownames(expr_consensus_genome))
shared_cell_type <- intersect(colnames(expr_own_genome), colnames(expr_consensus_genome))
shared_cell_type <- grep(paste0(subtypes,collapse = "|"), shared_cell_type, value = T)
expr_own_genome <- expr_own_genome[shared_genes, shared_cell_type]
expr_consensus_genome <- expr_consensus_genome[shared_genes, shared_cell_type]
mat <- do.call('rbind', strsplit(shared_cell_type, split = "@"))
combined_fc_mat <- matrix(NA, nrow=nrow(expr_own_genome), ncol=ncol(expr_own_genome))
rownames(combined_fc_mat) <- rownames(expr_own_genome)
colnames(combined_fc_mat) <- paste(rep(unique(mat[,1]), each=2), rep(c("Respective_genome", "Consensus_genome"), length(unique(mat[,2]))), sep="@")
for(ct in unique(mat[,1])){
  human_id <- paste0(ct, "@Human")
  chimp_id <- paste0(ct, "@Chimp")
  fc_own <- expr_own_genome[,human_id] - expr_own_genome[,chimp_id]
  fc_con <- expr_consensus_genome[,human_id] - expr_consensus_genome[,chimp_id]
  
  combined_fc_mat[names(fc_own),paste(ct,"Respective_genome",sep="@")] <- fc_own
  combined_fc_mat[names(fc_con),paste(ct,"Consensus_genome",sep="@")] <- fc_con
}
saveRDS(combined_fc_mat, file='Dat_expr_H_vs_C_logFC_respective_genome_and_consensus_genome_based_quantification.rds')
# generate scatter plot to compare the logFC
p_list <- list()
for(ct in unique(mat[,1])){
  df <- data.frame(combined_fc_mat[,paste(ct,c("Respective_genome", "Consensus_genome" ),sep="@")])
  names(df) <- c("X", "Y")
  p <- ggplot(df, aes(x=X, y=Y))+
    geom_point()+
    labs(title=ct,
         x ="Respective genome (H - C)", 
         y = "Consensus genome (H - C)")+
    theme_minimal()
  p_list[[ct]] <- p
    
}
pdf("Plot_scatter_plot_gene_expr_logFC.pdf", height=5*2, width=5*length(p_list)/2)
do.call("grid.arrange", c(p_list, ncol=length(p_list)/2))
dev.off()

# examine how many respecitve genome based DEG are also detected in the consensus genome based analysis
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/consensus_genome_RNA/")
selected_cell_type <- c("Stem_cell", "Enterocyte", "EC-cell", "non-EC-EEC", "Goblet_cell", "BEST4+_epithelium") 
# load the DEG list from the respective-genome-quantification-based analysis
tIO_deg <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_one2one_orhologs_DEGs_with_input_gene_list.rds")


# load the DEG list from the consensus-genome-quantification-based analysis
human_chimp_orthologs <- readRDS("/Volumes/pred/ihb-intestine-evo/used_object/orthologs/List_human_chimp_one2one_orthologs.rds")
deg_per_cell_type <- list()
all_deg <- c()
for(ct in selected_cell_type){
  species_test_res <- readRDS(paste0("Res_",ct,"_species_test_res.rds"))
  con_human_high_genes <- intersect(rownames(species_test_res)[which(species_test_res$Significant & species_test_res$Coef>0)], human_chimp_orthologs)
  con_chimp_high_genes <- intersect(rownames(species_test_res)[which(species_test_res$Significant & species_test_res$Coef<0)], human_chimp_orthologs)
  deg_per_cell_type[[ct]][["Human_high"]] <- con_human_high_genes
  deg_per_cell_type[[ct]][["Chimp_high"]] <- con_chimp_high_genes
  all_deg <- c(all_deg, con_human_high_genes, con_chimp_high_genes)
}
all_deg <- sort(unique(all_deg))
deg_mat <- matrix(NA, nrow=length(all_deg), ncol=length(selected_cell_type))
rownames(deg_mat) <- all_deg
colnames(deg_mat) <- selected_cell_type
for(ct in names(deg_per_cell_type)){
  con_human_high_genes <- deg_per_cell_type[[ct]][["Human_high"]]
  con_chimp_high_genes <- deg_per_cell_type[[ct]][["Chimp_high"]]
  deg_mat[con_human_high_genes, ct] <- "Human_high"
  deg_mat[con_chimp_high_genes, ct] <- "Chimp_high"
}
write.table(deg_mat, file="Table_consensus_genome_based_tIO_DEG.txt", sep="\t", quote=F)
write.table(deg_mat, file="/Volumes/pred/ihb-intestine-evo/used_object/supplementary_tables/Table_consensus_genome_based_tIO_DEG.txt", sep="\t", quote=F)



num_mat <- t(sapply(selected_cell_type, function(ct){
  species_test_res <- readRDS(paste0("Res_",ct,"_species_test_res.rds"))
  con_human_high_genes <- intersect(rownames(species_test_res)[which(species_test_res$Significant & species_test_res$Coef>0)], human_chimp_orthologs)
  con_chimp_high_genes <- intersect(rownames(species_test_res)[which(species_test_res$Significant & species_test_res$Coef<0)], human_chimp_orthologs)
  res_human_high_genes <- rownames(tIO_deg)[which(tIO_deg[,ct]=="Human_high")]
  res_chimp_high_genes <- rownames(tIO_deg)[which(tIO_deg[,ct]=="Chimp_high")]
  intersect_human_high_genes <- intersect(con_human_high_genes, res_human_high_genes)
  intersect_chimp_high_genes <- intersect(con_chimp_high_genes, res_chimp_high_genes)
  res <- c(length(con_human_high_genes), length(con_chimp_high_genes), length(res_human_high_genes), length(res_chimp_high_genes), length(intersect_human_high_genes), length(intersect_chimp_high_genes))
  return(res)
}))
colnames(num_mat) <- c("Consensus_genome_human_high", "Consensus_genome_chimp_high", "Respective_genome_human_high", "Respective_genome_chimp_high","Intersect_human_high", "Intersect_chimp_high")
rownames(num_mat) <- gsub("_", " ", rownames(num_mat))
saveRDS(num_mat, file="Res_tIO_epi_DEG_number_comparison_between_respective_and_consensus_genome_based_quantification.rds")

num_mat <- readRDS("Res_tIO_epi_DEG_number_comparison_between_respective_and_consensus_genome_based_quantification.rds")
# get the proportion of concordant DEGs in the respecitve genome identified DEGs
p_human <- num_mat[,"Intersect_human_high"]/num_mat[,"Respective_genome_human_high"]
p_chimp <- num_mat[,"Intersect_chimp_high"]/num_mat[,"Respective_genome_chimp_high"]
n1 <- sum(as.vector(num_mat[,grep("Respective",colnames(num_mat))]))
n2 <- sum(as.vector(num_mat[,grep("Intersect",colnames(num_mat))]))


input <- num_mat[,c("Consensus_genome_human_high", "Consensus_genome_chimp_high","Intersect_human_high", "Intersect_chimp_high")]
input[,grep("chimp",colnames(input))] <- -input[,grep("chimp",colnames(input))]
# source colors
source("/projects/site/pred/ihb-intestine-evo/colors/colors.R")
pdf("Plot_DEG_number_consensus_genome_and_intersect.pdf", height=5, width=5)
barplot(input[,"Consensus_genome_human_high"],beside=T, col=epi.ct.cols[rownames(input)], ylim=c(min(input), max(input)), ylab="DEG number", names="")
barplot(input[,"Intersect_human_high"],beside=T, density=20, add=T, names="", border="#303030", angle=45)
barplot(input[,"Consensus_genome_chimp_high"],beside=T, col=paste0(epi.ct.cols[rownames(input)], "A0"), add=T, names="")
barplot(input[,"Intersect_chimp_high"],beside=T, density=20, add=T, names="", border="#303030", angle=135)
dev.off()


