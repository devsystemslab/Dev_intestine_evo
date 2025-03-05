setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/DE_vs_DA")
library(Seurat)
library(Signac)
library(ggplot2)
library(ggrepel)

# load the tIO scRNA-seq gene expression data
tIO_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff/Dat_tIO_epi_species_ct_expr.rds")
# load the tIO scATAC-seq genomic region accessibility data
tIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/Dat_tIO_epi_ATAC_hg38_coor_cell_type_per_species_detection_rate.rds")
saveRDS(tIO_atac, file="/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/Dat_tIO_epi_ATAC_merged_peak_id_cell_type_per_species_detection_rate.rds")
rownames(tIO_atac) <- rownames(peak_anno)[match(rownames(tIO_atac), peak_anno$peak_id)]
saveRDS(tIO_atac, file="/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/Dat_tIO_epi_ATAC_hg38_coor_cell_type_per_species_detection_rate.rds")
# load multi layer annotation
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/Res_tIO_merged_peak_annotation-2.rds")
peak_id <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/Res_combined_peaks_in_panTro6_hg38_coor.rds")

peak_anno$peak_id <- peak_id$peak_ID[match(rownames(peak_anno), peak_id$combined_peaks_hg38_coor)]
peak_anno$panTro6_coor <- peak_id$combined_peaks_panTro6_coor[match(rownames(peak_anno), peak_id$combined_peaks_hg38_coor)]
saveRDS(peak_anno, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation.rds")
saveRDS(peak_anno, file="/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/Res_tIO_merged_peak_annotation-2.rds")

# load human scATAC-seq data
human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/ATAC/Res_tHIO_and_fetal_combined_scATAC-seq_with_fragment_files.rds")
human_atac$Cell_type_per_tissue <- paste(human_atac$Unified_cell_type, human_atac$Tissue, sep="@")
saveRDS(human_atac, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_human_and_tHIO_scATAC_with_fragment_file.rds")

# load chimp scATAC-seq data
tCIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/CIO_and_tCIO/scATAC-seq/tCIO/add_fragment_file/Dat_tCIO_epi_atac_with_fragment_file.rds")
saveRDS(tCIO_atac, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tCIO_scATAC_with_fragment_file.rds")

epi_rna_seurat <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_updated_human_chimp_mouse_epi_integrated_with_CSS_seurat_obj.rds")

papers_pos_selection <- c("Selective Sweep - Apes - Cagan 2016",
                                 "Human-specific segmental duplications - Dennis 2017",
                                 "Human-accelerated regions (combined) - Doan 2016",
                                 "Human-accelerated DNase I hypersensitive sites - Gittelman 2015",
                                 "Human-accelerated regions (zooHARs) - Keough 2023",                
                                 "Human-specific structural variants - Kronenberg 2018",
                                 "Hhuman ancestor quickly evolved regions - Mangan 2022",           
                                 "Human-specific deletions - McLean 2011",
                                 "Selective Sweep - Modern Humans - Peyregne 2017",                  
                                 "Human-accelerated regions - Pollard 2006",
                                 "Human-accelerated conserved non-coding sequences - Prabhakar 2006",
                                 "Human-specific deletions - Xue 2023")
har_papers <- c("Human-accelerated regions (combined) - Doan 2016",
                "Human-accelerated DNase I hypersensitive sites - Gittelman 2015",
                "Human-accelerated regions (zooHARs) - Keough 2023",                
                "Hhuman ancestor quickly evolved regions - Mangan 2022",           
                "Human-accelerated regions - Pollard 2006",
                "Human-accelerated conserved non-coding sequences - Prabhakar 2006")
regions_under_pos_selection <- rownames(peak_anno)[which(rowSums(peak_anno[,papers_pos_selection])>0)]
har_regions <- rownames(peak_anno)[which(rowSums(peak_anno[,har_papers])>0)]
selected_groups <- setdiff(intersect(colnames(tIO_rna), colnames(tIO_atac)), c("Early_enterocyte@Chimp",  "Early_enterocyte@Human"))
mat <- do.call(
  'rbind', strsplit(selected_groups, split="@")
)


input_list <- list()
for(ct in unique(mat[,1])){
  idx <- which(peak_anno[,paste0("Expressed:",ct,"@fetal")] & peak_anno$Closest_gene %in% rownames(tIO_rna))
  genes <- peak_anno$Closest_gene[idx]
  peaks <- peak_anno$peak_id[idx]
  rna_logFC <- tIO_rna[genes,paste(ct,"Human",sep="@")] - tIO_rna[genes,paste(ct,"Chimp",sep="@")]
  atac_logFC <- tIO_atac[peaks,paste(ct,"Human",sep="@")] - tIO_atac[peaks,paste(ct,"Chimp",sep="@")]
  
  atac_human_ct_enrichment <- (tIO_atac[peaks,paste(ct,"Human",sep="@")] - rowMeans(tIO_atac[peaks, grep("@Human", colnames(tIO_atac))]))+1
  
  dar<- c(rownames(peak_anno)[which(peak_anno[,paste0("Human_high:",ct,"@fetal")])],
          rownames(peak_anno)[which(peak_anno[,paste0("Chimp_high:",ct,"@fetal")])])
  if(ct != "EEC"){
    regions_linked_to_deg <- rownames(peak_anno)[which(peak_anno[,paste0("H_and_C:",ct,"@fetal")])]
  }else{
    regions_linked_to_deg <- rownames(peak_anno)[which(rowSums(peak_anno[,paste0("H_and_C:",c("EC-cell","non-EC-EEC","EEC_progenitor"),"@fetal")])>0)]
  }
  
  
  #if(ct != "BEST4+_epithelium"){
  #  regions_linked_to_deg <- c(regions_linked_to_deg,
  #                             rownames(peak_anno)[which(peak_anno[,paste0("H_and_M:",ct,"@fetal")])],
  #                             rownames(peak_anno)[which(peak_anno[,paste0("Human_specific:",ct,"@fetal")])])
  #}
  dar_linked_to_deg <- intersect(dar, regions_linked_to_deg)
  dar_only <- setdiff(dar, regions_linked_to_deg)
  linked_to_deg_only <- setdiff(regions_linked_to_deg, dar)
  group_vec <- rep("Detected", length(idx))
  group_vec[which(rownames(peak_anno)[idx]%in%dar_only)] <- "DAR_only"
  group_vec[which(rownames(peak_anno)[idx]%in%linked_to_deg_only)] <- "Linked_to_DEG_only"
  group_vec[which(rownames(peak_anno)[idx]%in%dar_linked_to_deg)] <- "DAR_linked_to_DEG"
  input_list[[ct]] <- data.frame("X"=rna_logFC,
                                 "Y"=atac_logFC,
                                 "size"=atac_human_ct_enrichment,
                                 "name"=rownames(peak_anno)[idx],
                                 "group"=group_vec)
}
saveRDS(input_list, file="Dat_plotting_DA_DE_input.rds")


region_group_cols <- setNames(c("#696969","#F4D03F","#3498DB","#16A085"),
                              c("Detected","DAR_only","Linked_to_DEG_only","DAR_linked_to_DEG"))

p_list <- list()
for(ct in names(input_list)){
  df <- input_list[[ct]]
  p_list[[ct]] <- ggplot(df, aes(x=X, y=Y, color=group, size=size^4*10)) +
    geom_point()+
    scale_color_manual(values=region_group_cols)+
    labs(title=ct,
         x="Expression (Human-Chimp)", 
         y = "Accessibility (Human-Chimp)",
         size="Cell type enrichment", 
         color="Groups")+
    theme_minimal()+
    theme(text = element_text(size = 80))+
    guides(color = guide_legend(override.aes = list(size=10)))
}
                                   
#library("gridExtra")                                             
row_num=1
col_num=5
png("Plot_scatter_DA_vs_DE.png", height=2000*row_num, width=2500*col_num)
do.call("grid.arrange", c(p_list, ncol = 5))
dev.off()

# scatter plot per cell type, x and y axis represents region detection rate difference and gene expression difference in human and chimp
# highlight the regions that are DA, linked to DE, cell type enriched and overlapped with positive selection signatures
selected_marker_da_and_de_idx_by_cell_type <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/group_DAR/Res_positive_selected_cell_type_marker_DA_regions_linked_to_DE_by_cell_type.rds")
input_list <- list()
p_list <- list()
selected_cell_type <- colnames(selected_marker_da_and_de_idx_by_cell_type)[which(colSums(selected_marker_da_and_de_idx_by_cell_type)>0)]
for(ct in selected_cell_type){
  idx <- which(peak_anno[,paste0("Expressed:",ct,"@fetal")] & peak_anno$Closest_gene %in% rownames(tIO_rna))
  genes <- peak_anno$Closest_gene[idx]
  peaks <- peak_anno$peak_id[idx]
  rna_logFC <- tIO_rna[genes,paste(ct,"Human",sep="@")] - tIO_rna[genes,paste(ct,"Chimp",sep="@")]
  atac_logFC <- tIO_atac[peaks,paste(ct,"Human",sep="@")] - tIO_atac[peaks,paste(ct,"Chimp",sep="@")]

  regions_to_highlight <- rownames(selected_marker_da_and_de_idx_by_cell_type)[which(selected_marker_da_and_de_idx_by_cell_type[,ct]>0)]
  df <- data.frame("X"=rna_logFC,
                   "Y"=atac_logFC,
                   "name"=rownames(peak_anno)[idx])
  input_list[[ct]] <- df
  p_list[[ct]] <- ggplot(df, aes(x=X, y=Y, label=name)) +
    geom_point(color="#d0d0d0")+
    geom_label_repel(data = subset(df, df$name%in%regions_to_highlight),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     max.overlaps = Inf)+
    labs(title=ct,
         x="Expression (Human-Chimp)", 
         y = "Accessibility (Human-Chimp)")+
    geom_point(data=subset(df, df$name%in%regions_to_highlight),
               color="#303030")+
    theme_minimal()
  
}

#library("gridExtra")                                             
row_num=1
col_num=length(p_list)
pdf("Plot_scatter_DA_vs_DE_for_ct_marker_and_positive_selected_regions.pdf", height=10*row_num, width=10*col_num)
do.call("grid.arrange", c(p_list, ncol = length(p_list)))
dev.off()

# generate coverage plots and gene expression feature plots
tIO_epi <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
source("/projects/site/pred/ihb-intestine-evo/colors/colors.R")
selected_cell_type <- c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell","Stem cell")
vec1 <- rep(selected_cell_type, each=2)
vec2 <- rep(c("Fetal primary", "transplanted"), length(selected_cell_type))
selected_groups <- paste(vec1,
                         vec2,
                         sep="@")
alpha <- setNames(c("30","FF"), c("Fetal primary", "transplanted"))
tissue.epi.ct.cols <- setNames(paste0(epi.ct.cols[vec1],alpha[vec2]),
                               selected_groups)

for(ct in colnames(selected_marker_da_and_de_idx_by_cell_type)){
  regions <- rownames(selected_marker_da_and_de_idx_by_cell_type)[which(selected_marker_da_and_de_idx_by_cell_type[,ct]>0)]
  if(length(regions)==0) next
  path <- paste0("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/DE_vs_DA/",ct)
  if(!dir.exists(path)){
    dir.create(path)
  }
  setwd(path)
  for(p_human in regions){
    
    g <- peak_anno[p_human,"Closest_gene"]
    p_chimp <- peak_anno[p_human, "panTro6_coor"]
    
    # coverage plot in chimp
    p1 <- CoveragePlot(
      object = tCIO_atac,
      assay="CIO_unified_peaks",
      group.by = "Cell_type",
      region = p_chimp,
      idents = selected_cell_type,
      extend.upstream = 100,
      extend.downstream = 100,
      #features = g,
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = TRUE,
      window = 500
    )
    p1 <- p1 & scale_fill_manual(values = epi.ct.cols)
    
    # coverage plot in human
    p2 <- CoveragePlot(
      object = human_atac,
      assay="HIO_unified_peaks",
      group.by = "Cell_type_per_tissue",
      idents = selected_groups,
      region = p_human,
      extend.upstream = 100,
      extend.downstream = 100,
      #features = g,
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = TRUE,
      window = 500
    )
    p2 <- p2 & scale_fill_manual(values = tissue.epi.ct.cols) 
    
    # gene expression feature plot in human, chimp and mouse combined data
    p3 <- SCpubr::do_FeaturePlot(tIO_epi, 
                                 reduction = "umap_css", 
                                 features = g, 
                                 order = T,
                                 pt.size = 5,
                                 font.size = 80,
                                 legend.framewidth = 2,
                                 legend.width = 5,
                                 legend.length = 50)
    
    plot_name <- paste0("Plot_coveragePlot_hg38_",p_human,"_panTro6_",p_chimp,"_",g,"_10k_smoothed.pdf")
    pdf(plot_name, height=7, width=12)
    print(p1+p2)
    dev.off()
    plot_name <- paste0("Plot_UMAP_CSS_tIO_epi_",g,"_expr.png")
    png(plot_name, height=2000, width=2000)
    print(p3)
    dev.off()
  }
}



# generate scatter plot showing the accessibility difference between human and chimp in stem cells and enterocytes
stem_cell_atac_logFC <- tIO_atac[,"Enterocyte@Human"] - tIO_atac[,"Enterocyte@Chimp"]
enterocyte_atac_logFC <- tIO_atac[,"Stem_cell@Human"] - tIO_atac[,"Stem_cell@Chimp"]
eec_atac_logFC <- tIO_atac[,"EEC@Human"] - tIO_atac[,"EEC@Chimp"]
selected_groups <- c("Human-specific segmental duplications - Dennis 2017",               
                     "Human-accelerated regions (combined) - Doan 2016",                 
                     "Human-accelerated DNase I hypersensitive sites - Gittelman 2015",   
                     "Human-accelerated regions (zooHARs) - Keough 2023",                
                     "Human-specific structural variants - Kronenberg 2018",             
                     "Hhuman ancestor quickly evolved regions - Mangan 2022",            
                     "Human-specific deletions - McLean 2011",                            
                     "Selective Sweep - Modern Humans - Peyregne 2017",                  
                     "Human-accelerated regions - Pollard 2006",                          
                     "Human-accelerated conserved non-coding sequences - Prabhakar 2006",
                     "Human-specific deletions - Xue 2023")
pos_selected_regions <- rownames(peak_anno)[which(rowSums(peak_anno[,selected_groups])>0)]
df <- data.frame("X"=enterocyte_atac_logFC,
                 "Y"=eec_atac_logFC,
                 "name"=rownames(peak_anno))
ggplot(df, aes(x=X,y=Y))+
  geom_point(color="#d0d0d0")+
  labs(title="Accessibility (Human - Chimp)",
       x="Enterocyte", 
       y = "EEC")+
  geom_point(data=subset(df, df$name%in%pos_selected_regions),
             color="#252525")+
  theme_minimal()


regions_to_highlight <- rownames(selected_marker_da_and_de_idx_by_cell_type)[which(selected_marker_da_and_de_idx_by_cell_type[,ct]>0)]
df <- data.frame("X"=rna_logFC,
                 "Y"=atac_logFC,
                 "name"=rownames(peak_anno)[idx])
input_list[[ct]] <- df
p_list[[ct]] <- ggplot(df, aes(x=X, y=Y, label=name)) +
  geom_point(color="#d0d0d0")+
  geom_label_repel(data = subset(df, df$name%in%regions_to_highlight),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = Inf)+
  labs(title=ct,
       x="Expression (Human-Chimp)", 
       y = "Accessibility (Human-Chimp)")+
  geom_point(data=subset(df, df$name%in%regions_to_highlight),
             color="#303030")+
  theme_minimal()











# scatter plot to show gene expression logFC in human-chimp tIO comparison vs human-mouse developing tissue comparisons
setwd("/projects/site/pred/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/pairwise_DEGs/DE_vs_DE")
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
epi_expr <- getAveExpr(seu.obj=epi_rna_seurat, feature.to.calc = "Cell_type_per_group", colname.prefix = NULL) 
saveRDS(epi_expr, file="Dat_human_chimp_mouse_cell_type_per_group_expr.rds")

# get human-chimp DEGs
hc_degs <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_DEGs_with_input_gene_list.rds")
# get human-mouse DEGs
hm_degs <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_fetal_human_and_mouse_ANCOVA_DEG_per_cell_type_res.rds")

gene_group_cols <- setNames(c("#696969","#F4D03F","#3498DB","#16A085", "#C0392B"),
                              c("Expressed","HC_only","HM_only","Both", "Highlight"))

selected_cell_type <- c("EC-cell", "Enterocyte", "Goblet cell", "non-EC-EEC", "Stem cell")

hc_logFC_mat <- sapply(selected_cell_type, function(ct){
  epi_expr[,paste0(ct,"@Human@transplanted")] - epi_expr[,paste0(ct,"@Chimp@transplanted")]
})
hm_logFC_mat <- sapply(selected_cell_type, function(ct){
  epi_expr[,paste0(ct,"@Human@Fetal primary")] - epi_expr[,paste0(ct,"@Mouse@Fetal primary")]
})
hm_logFC_ct_enrichment <- t(apply(abs(hm_logFC_mat), 1, function(vec){
  if(sum(vec)==0){
    return(rep(0,ncol(hm_logFC_mat)))
  }else{
    return((vec-min(vec))/(max(vec)-min(vec)))
  }
}))
colnames(hm_logFC_ct_enrichment) <- colnames(hm_logFC_mat)
hc_logFC_ct_enrichment <- t(apply(abs(hc_logFC_mat), 1, function(vec){
  if(sum(vec)==0){
    return(rep(0,ncol(hm_logFC_mat)))
  }else{
    return((vec-min(vec))/(max(vec)-min(vec)))
  }
}))
colnames(hc_logFC_ct_enrichment) <- colnames(hc_logFC_mat)
saveRDS(hc_logFC_mat, file="Dat_tIO_human_chimp_expr_logFC.rds")
saveRDS(hm_logFC_mat, file="Dat_fetal_human_and_mouse_expr_logFC.rds")
saveRDS(hm_logFC_ct_enrichment, file="Dat_fetal_human_and_mouse_logFC_cell_type_enrichment.rds")
saveRDS(hc_logFC_ct_enrichment, file="Dat_tIO_logFC_cell_type_enrichment.rds")

highlight_genes <- list(
  "Enterocyte"=c("SLC26A3","SLC5A12","SLC46A1"),
  "Stem_cell"=c("ETV4","MYC"),
  "EC-cell"="FEV"
)

# read disease-associated gene list
disease_genes <- read.table("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/Intestine_disease_associated/with_data_cm/v2/Table_manually_curated_GI_disease_related_genes.txt",
                            sep="\t",
                            stringsAsFactors = F,
                            fill=T,
                            quote="",
                            head=T)
disease_gene_symbols <- sort(unique(disease_genes$Gene))




p_list <- list()
input_list <- list()
for(ct in selected_cell_type){
  ct2 <- sub(pattern=" ", replacement = "_", ct)
  hc_genes <- rownames(hc_degs)[which(hc_degs[,ct2]%in%c("Human_high", "Chimp_high"))]
  hm_genes <- rownames(hm_degs)[which(hm_degs[,ct2]%in%c("Human_high", "Mouse_high"))]
  shared_degs <- intersect(intersect(hc_genes, hm_genes), rownames(epi_expr))
  hc_specific_degs <- intersect(setdiff(hc_genes, hm_genes), rownames(epi_expr))
  hm_specific_degs <- intersect(setdiff(hm_genes, hc_genes), rownames(epi_expr))
  group_vec <- setNames(rep("Expressed", nrow(epi_expr)), rownames(epi_expr))
  group_vec[shared_degs] <- "Both"
  group_vec[hc_specific_degs] <- "HC_only"
  group_vec[hm_specific_degs] <- "HM_only"
  
  
  df <- data.frame("Human_Chimp"=hc_logFC_mat[,ct],
                   "Human_Mouse"=hm_logFC_mat[,ct],
                   "group"=group_vec,
                   "size"=hc_logFC_ct_enrichment[,ct]*0.2,
                   "name"=rownames(hc_logFC_mat))
  input_list[[ct]] <- df
  
  p_list[[ct]] <- ggplot(df, aes(x=Human_Chimp, y=Human_Mouse, color=group, size=size, label=name)) +
    geom_point()+
    scale_color_manual(values=gene_group_cols)+
    labs(title=ct)+
    theme_minimal()+
    #theme(text = element_text(size = 8))+
    #geom_text(aes(label=ifelse(df$group=="Highlight",as.character(name),'')),
    #          hjust=0,
    #          vjust=0,
    #          color="#252525")+
    geom_label_repel(data = subset(df, df$name%in%disease_gene_symbols),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     max.overlaps = Inf)
    #guides(color = guide_legend(override.aes = list(size=10)))
  
}
saveRDS(input_list, file="Dat_hc_vs_hm_DE_plotting_input.rds")

row_num=1
col_num=length(p_list)
pdf("Plot_scatter_DE_vs_DE.pdf", height=10*row_num, width=10*col_num)
do.call("grid.arrange", c(p_list, ncol = 5))
dev.off()



# generate scatter plot for genes, using cell type expression logFC, dn/ds ratio as X and Y axis
selected_cell_type <- c("BEST4+_epithelium", "EEC", "Enterocyte", "Goblet_cell", "Stem_cell")
tIO_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff/Dat_tIO_epi_species_ct_expr.rds")
idx <- which(rowSums(tIO_rna)>0)
length(idx)
nrow(tIO_rna)
tIO_rna <- tIO_rna[idx,]


rna_logFC_per_cell_type <- sapply(selected_cell_type, function(ct){
  tIO_rna[,paste(ct,"Human",sep="@")] - tIO_rna[,paste(ct,"Chimp",sep="@")]
})
rna_logFC <- apply(rna_logFC_per_cell_type, 1, function(vec){
  vec[which.max(abs(vec))]
})
dnds <- readRDS("/projects/site/pred/ihb-intestine-evo/evo_signature/DN_DS_ratio/Human_mouse/Res_human_chimp_ensemblv99_dnds.rds")
dnds <- dnds[which(dnds<5)]
shared_genes <- intersect(names(dnds), names(rna_logFC))
length(shared_genes)

# load overlapping between gene regions (gene body +/- 10k) and selection signatures
evo_sig_overlap_mat <- readRDS("/projects/site/pred/ihb-intestine-evo/evo_signature/overlap_with_gene_10k/Res_ensembl93_human_gene_10k_overlap_with_evolution_signatures.rds")
pos_selected_region_nearby_genes <- sort(unique(rownames(evo_sig_overlap_mat)[which(rowSums(evo_sig_overlap_mat[,-1])>0)]))
length(pos_selected_region_nearby_genes)
hc_deg_info <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_DEGs_with_input_gene_list.rds") 
union_deg <- unique(unlist(lapply(seq(ncol(hc_deg_info)), function(j){
  rownames(hc_deg_info)[which(hc_deg_info[,j]%in%c("Human_high", "Chimp_high"))]
})))
pos_selected_deg <- intersect(pos_selected_region_nearby_genes, union_deg)
length(pos_selected_deg)

# load cell type markers
ct_markers <- readRDS("/projects/site/pred/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/cell_type_marker/Res_cell_type_marker_per_species_res_list.rds")
human_mat <- ct_markers[["Human"]]
human_ct_markers <- unique(human_mat$feature)

chimp_mat <- ct_markers[["Chimp"]]
chimp_ct_markers <- unique(chimp_mat$feature)
union_ct_markers <- union(human_ct_markers, chimp_ct_markers)
length(union_ct_markers)
pos_selected_marker_deg <- intersect(union_ct_markers, pos_selected_deg)
length(pos_selected_marker_deg)

df <- data.frame("X"=rna_logFC[shared_genes],
                 "Y"=dnds[shared_genes],
                 "name"=shared_genes)
genes <- sort(intersect(df$name[which(abs(df$X)>0.8)], pos_selected_marker_deg))
length(genes)

pdf("Plot_scatter_plot_expression_dnds.pdf")
ggplot(df, aes(x=X,y=Y, label=name))+
  geom_point(color="#d0d0d0")+
  labs(x="Expression logFC", 
       y = "dN/dS")+
  geom_point(data=subset(df, df$name%in%genes), color="#252525")+
  geom_label_repel(data = subset(df, df$name%in%genes),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = Inf)+
  theme_minimal()+
  theme(text = element_text(size = 20))
dev.off()

# generate dot plot to show the gene expression in tIO cell type
selected_groups <- paste(rep(colnames(hc_deg_info), 2),
                        rep(c("Human","Chimp"), each=ncol(hc_deg_info)),
                        sep="@")
cells <- colnames(tIO_epi)[which(tIO_epi$Cell_type_per_species%in%selected_groups)]
length(cells)
tIO_epi_subset <- subset(tIO_epi, cells=cells)
SCpubr::do_DotPlot(tIO_epi_subset, features = genes, group.by = "Cell_type_per_species", cluster.idents = T)


expr_mat <- tIO_rna[genes, selected_groups]





heagenes <- intersect(c("NPC1L1",
           "APOB48",
           "APOA1",
           "ABCA1",
           "ABCG5",
           "ABCG8",
           "SCARB1",
           "PLAGL2"),
           rownames(epi_rna_seurat))
col_num=4
row_num=ceiling(length(genes)/col_num)
png("Plot_UMAP_CSS_human_chimp_mouse_cholesterol_gene_expression.png",
    height=2000*row_num,
    width=2000*col_num)
SCpubr::do_FeaturePlot(epi_rna_seurat, 
                       reduction = "umap_css", 
                       features = genes, 
                       order = T, 
                       ncol = col_num, 
                       pt.size = 5,
                       font.size = 80,
                       legend.framewidth = 2,
                       legend.width = 5,
                       legend.length = 50)
dev.off()


# perform GO enrichment on the DEGs without considering the directions
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/DE_vs_DA/GO")
human_chimp_filtered_DEGs <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff/Res_human_chimp_filtered_DEGs_with_input_gene_list.rds")
idx <- (rowSums(human_chimp_filtered_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Human_high", na.rm=T)+rowSums(human_chimp_filtered_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Chimp_high", na.rm = T))>0
sig_mat <- sapply(c("BEST4+_epithelium", "Enterocyte", "Goblet_cell", "Stem_cell"), function(ct){
  human_chimp_filtered_DEGs[,ct]%in%c("Human_high", "Chimp_high")
})
rownames(sig_mat) <- rownames(human_chimp_filtered_DEGs)
sig_df <- data.frame(sig_mat, "EEC"=idx, stringsAsFactors = F)
colnames(sig_df)[1] <- "BEST4+_epithelium"
saveRDS(sig_df, file="Res_tIO_filtered_DEG_without_considering_direction.rds")


idx <- (rowSums(human_chimp_filtered_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Human_high", na.rm=T)+
  rowSums(human_chimp_filtered_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Chimp_high", na.rm = T)+
    rowSums(human_chimp_filtered_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Expressed", na.rm = T))>0
expressed_mat <- sapply(c("BEST4+_epithelium", "Enterocyte", "Goblet_cell", "Stem_cell"), function(ct){
  human_chimp_filtered_DEGs[,ct]%in%c("Human_high", "Chimp_high", "Expressed")
})
rownames(expressed_mat) <- rownames(human_chimp_filtered_DEGs)
expressed_df <- data.frame(expressed_mat, "EEC"=idx, stringsAsFactors = F)
colnames(expressed_df)[1] <- "BEST4+_epithelium"
saveRDS(expressed_df, file="Res_tIO_expressed_genes.rds")


GO_anno <- read.csv("/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v109/Ensembl_v109_GO.csv")
idx <- which(GO_anno$GO.domain!="")
GO_anno <- GO_anno[idx,]
terms <- c("reactive oxygen species metabolic process",
           "carbohydrate metabolic process",
           "sterol metabolic process",
           "NF-kappaB transcription factor activity",
           "retinoid metabolic process")
#for(group in sort(unique(GO_anno$GO.domain))){
group <- "biological_process"
print(paste(group,"start"))
selected_gs <- sort(unique(GO_anno$GO.term.name[which(GO_anno$GO.domain==group)]))

pval <- matrix(NA, nrow=length(selected_gs), ncol=ncol(sig_df))
rownames(pval) <- selected_gs
colnames(pval) <- colnames(sig_df)

or <- matrix(NA, nrow=length(selected_gs), ncol=ncol(sig_df))
rownames(or) <- selected_gs
colnames(or) <- colnames(sig_df)

prop <- matrix(NA, nrow=length(selected_gs), ncol=ncol(sig_df))
rownames(prop) <- selected_gs
colnames(prop) <- colnames(sig_df)

n1 <- matrix(NA, nrow=length(selected_gs), ncol=ncol(sig_df))
rownames(n1) <- selected_gs
colnames(n1) <- colnames(sig_df)


for(gs in selected_gs){
  go_genes <- GO_anno$HGNC.symbol[which(GO_anno$GO.term.name==gs)]
  for(x in colnames(sig_df)){
    cluster_genes <- rownames(sig_df)[which(sig_df[,x])]
    all_genes <- rownames(expressed_df)[which(expressed_df[,x])]
    genes <- intersect(go_genes, all_genes)
    a <- sum(cluster_genes%in%genes)
    b <- length(genes)
    c <- length(cluster_genes)
    d <- length(all_genes)
    #pval[gs,x] <- fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")$p.value
    #prop[gs,x] <- a/c
    n1[gs,x] <- a
    #or[gs,x] <- fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")$estimate
  }
}

padj <- apply(pval, 2, p.adjust, method="BH")


res <- list("Hypogeometric_nominal_P"=pval,
            "BH_corrected_P"=padj,
            "Prop"=prop,
            "Odds_ratio"=or,
            "Hit_number"=n1)
saveRDS(res, file="Res_tIO_filtered_DEG_no_direction_GO_BP_enrichment_res.rds")


deg_enriched_term_list <- lapply(seq(ncol(pval)), function(j){
  intersect(rownames(pval)[which(pval[,j]<0.05)], rownames(or)[which(or[,j]>1.5)])
})
names(deg_enriched_term_list) <- colnames(pval)
sapply(deg_enriched_term_list, length)
saveRDS(deg_enriched_term_list, file=paste0("Res_tIO_filtered_DEG_no_direction_GO_",group,"_deg_enriched_term_list.rds"))


deg_enriched_term_list <- lapply(seq(ncol(pval)), function(j){
  rownames(pval)[which(pval[,j]<0.05)]
})
names(deg_enriched_term_list) <- colnames(pval)
sapply(deg_enriched_term_list, length)
saveRDS(deg_enriched_term_list, file=paste0("Res_tIO_filtered_DEG_no_direction_GO_",group,"_deg_enriched_term_list_pval_filtering_only.rds"))



# run GREAT analysis on DAP per cell type without considering directions
# read species DA peak list
ddp_by_cell_type <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/Res_DDP_by_cell_type_2.rds")
orth_peaks <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/Res_combined_peaks_in_panTro6_hg38_coor.rds")
# run GREAT analysis
library(rGREAT)
great_res_list <- list()
for(x in colnames(sig_df)){
  print(x)
  mat <- ddp_by_cell_type[[x]]
  rownames(mat) <- orth_peaks$combined_peaks_hg38_coor[match(rownames(mat), orth_peaks$peak_ID)]
  peaks <- union(rownames(mat)[which(mat$Human_high)], rownames(mat)[which(mat$Chimp_high)])
  gr <- StringToGRanges(peaks)
  job = submitGreatJob(gr, species = "hg38")
  tbl = getEnrichmentTables(job, download_by = "tsv")
  great_res_list[[x]] <- tbl

}
saveRDS(great_res_list, file="Res_online_rGREAT_res_per_cell_type_DA_peaks_without_considering_direction.rds")

dar_enriched_term_list <- list()
for(ct in names(great_res_list)){
  
    res <- great_res_list[[ct]][["GO Biological Process"]]
    enriched_terms <- res$Desc[which(res$BinomFdrQ<0.05 & res$RegionFoldEnrich>1.5)]
    dar_enriched_term_list[[ct]] <- enriched_terms
  
}
saveRDS(dar_enriched_term_list, file="Res_tIO_DAR_no_direction_enriched_GO_BP_terms_by_cell_type.rds")


dar_enriched_term_list <- list()
for(ct in names(great_res_list)){
  
  res <- great_res_list[[ct]][["GO Biological Process"]]
  enriched_terms <- res$Desc[which(res$BinomFdrQ<0.05)]
  dar_enriched_term_list[[ct]] <- enriched_terms
  
}
saveRDS(dar_enriched_term_list, file="Res_tIO_DAR_no_direction_enriched_GO_BP_terms_by_cell_type_pval_filtering_only.rds")


sapply(names(dar_enriched_term_list), function(x){
  sapply(names(dar_enriched_term_list), function(y){
    length(intersect(dar_enriched_term_list[[x]], dar_enriched_term_list[[y]]))
  })
})


sapply(names(deg_enriched_term_list), function(x){
  sapply(names(deg_enriched_term_list), function(y){
    length(intersect(deg_enriched_term_list[[x]], deg_enriched_term_list[[y]]))
  })
})


overlapped_enriched_terms <- lapply(names(dar_enriched_term_list), function(ct){
  intersect(dar_enriched_term_list[[ct]], deg_enriched_term_list[[ct]])
})
names(overlapped_enriched_terms) <- names(dar_enriched_term_list)
sapply(overlapped_enriched_terms, length)
saveRDS(overlapped_enriched_terms, file="Res_tIO_RNA_and_ATAC_overlapped_enriched_GO_terms.rds")


union_enriched_terms <- lapply(names(dar_enriched_term_list), function(ct){
  intersect(union(dar_enriched_term_list[[ct]], deg_enriched_term_list[[ct]]), intersect(great_res_list[[ct]]$`GO Biological Process`$Desc, selected_gs))
})
names(union_enriched_terms) <- names(dar_enriched_term_list)
sapply(union_enriched_terms, length)
saveRDS(union_enriched_terms, file="Res_tIO_RNA_or_ATAC_enriched_GO_terms.rds")

intersect(GO_anno$GO.term.name[which(GO_anno$HGNC.symbol=="SLC5A12")], overlapped_enriched_terms[["Enterocyte"]])

# use the pval of DEGs and BinomFdrQ of DARs as the X and Y axis for the GO terms that are enriched in either modality and get involved in both modalities  
p_list <- list()
for(ct in names(union_enriched_terms)){
  terms <- union_enriched_terms[[ct]]
  gene_pval <- -log10(pval[terms,ct])
  dar_res <- great_res_list[[ct]][["GO Biological Process"]]
  region_pval <- -log10(dar_res$BinomP[match(terms, dar_res$Desc)])
  df <- data.frame("X"=gene_pval,
                   "Y"=region_pval,
                   "name"=terms)
  p_list[[ct]] <- ggplot(df, aes(x=X, y=Y, label=name)) +
    geom_point(color="#d0d0d0")+
    geom_label_repel(data = subset(df, df$name%in%overlapped_enriched_terms[[ct]]),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     max.overlaps = Inf)+
    labs(title=ct,
         x="DEG -log10(P)", 
         y = "DAR -log10(P)")+
    geom_point(data=subset(df, df$name%in%overlapped_enriched_terms[[ct]]),
               color="#303030")+
    theme_minimal()
}


#library("gridExtra")                                             
row_num=1
col_num=length(p_list)
pdf("Plot_scatter_DA_DE_union_enriched_GO_on_pval_and_or.pdf", height=7*row_num, width=7*col_num)
do.call("grid.arrange", c(p_list, ncol = 5))
dev.off()

terms <- c("steroid metabolic process",
           "defense response to virus",
           "cholesterol metabolic process")
# get the DEGs and DARs involved in the selected terms 
# generate scatter plot using cell type enrichment and p-value as x and y axis, size by odds ratio

tIO_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
ct_expr <- getAveExpr(seu.obj=tIO_rna, feature.to.calc = "Unified_cell_type_per_species", colname.prefix = NULL)
saveRDS(ct_expr, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_species_cell_type_average_expr.rds")

# get the cell type enrichment in human
human_ct_expr <- ct_expr[,grep("@Human", colnames(ct_expr))]
scaled_human_ct_expr <- t(scale(t(human_ct_expr)))
colnames(scaled_human_ct_expr) <- sub("@Human","",colnames(scaled_human_ct_expr))
saveRDS(scaled_human_ct_expr, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tHIO_epi_scaled_cell_type_average_expr.rds")


# calculate the cell type enrichment of DEGs of each cell type
cell_type_specificity <- t(sapply(selected_gs, function(gs){
  go_genes <- GO_anno$HGNC.symbol[which(GO_anno$GO.term.name==gs)]
  sapply(colnames(sig_df), function(x){
    cluster_genes <- rownames(sig_df)[which(sig_df[,x])]
    all_genes <- rownames(expressed_df)[which(expressed_df[,x])]
    genes <- intersect(cluster_genes, intersect(go_genes, all_genes))
    mean(abs(scaled_human_ct_expr[genes,x]))
  })
}))
saveRDS(cell_type_specificity, file="Res_per_cell_type_DEG_involved_GO_term_average_cell_type_specificity.rds")


p_list <- list()
for(ct in colnames(pval)){
  df <- data.frame("X"=cell_type_specificity[,ct],
                   "Y"=-log10(pval[,ct]),
                   "name"=rownames(pval))
  p_list[[ct]] <- ggplot(df, aes(x=X, y=Y, label=name)) +
    geom_point(color="#d0d0d0")+
    geom_label_repel(data = subset(df, df$name%in%overlapped_enriched_terms[[ct]]),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     max.overlaps = Inf)+
    labs(title=ct,
         y="DEG -log10(P)", 
         x = "Human cell type specificity")+
    geom_point(data=subset(df, df$name%in%overlapped_enriched_terms[[ct]]),
               color="#303030")+
    theme_minimal()
}


#library("gridExtra")                                             
row_num=1
col_num=length(p_list)
pdf("Plot_scatter_gene_pval_and_gene_cell_type_enrichment_for_DA_DE_overlapped_enriched_GO.pdf", height=7*row_num, width=7*col_num)
do.call("grid.arrange", c(p_list, ncol = 5))
dev.off()


# get SLC5A12 and SLC26A3 involved GO terms and check its enrichment on DEGs and DARs
g1 <- "SLC5A12"
gs <- GO_anno$GO.term.name[GO_anno$HGNC.symbol==g1 & GO_anno$GO.domain=="biological_process"]

# load the GO enrichment for DEGs classifying the species differential direction
deg_GO_res <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff/Res_GO_biological_process_enrichment_res.rds")
pval_v1 <- deg_GO_res$Hypogeometric_nominal_P
or_v1 <- deg_GO_res$Odds_ratio
prop_v1 <- deg_GO_res$Prop
pval_v1[gs,]
# SLC5A12 is slightly enriched in lactate transmembrane transport if classifying the change direction - human-high, nominal hypergeometric test P=0.12

## exclude goblet cell progenitor, combine EEC subtypes
## take the union sets of the chimp high and human high feature enriched GO terms
pval_v2 <- matrix(NA, nrow=nrow(pval_v1), ncol=5)
rownames(pval_v2) <- rownames(pval_v1)
colnames(pval_v2) <- c("BEST4+_epithelium", "EEC","Enterocyte","Goblet_cell", "Stem_cell")

or_v2 <- matrix(NA, nrow=nrow(or_v1), ncol=5)
rownames(or_v2) <- rownames(or_v1)
colnames(or_v2) <- c("BEST4+_epithelium", "EEC","Enterocyte","Goblet_cell", "Stem_cell")

pval_mat <- pval_v1[,grep("EC", colnames(pval_v1))]
idx <- apply(pval_mat, 1, which.min)
pval_v2[,"EEC"] <- apply(pval_mat, 1, min)
or_mat <- or_v1[,grep("EC", colnames(or_v1))]
or_v2[,"EEC"] <- sapply(seq(length(idx)), function(i){
  or_mat[i,idx[i]]
})

for(ct in c("BEST4+_epithelium", "Enterocyte", "Goblet_cell", "Stem_cell")){
  id <- paste(ct, c("Human","Chimp"), sep=":")
  pval_mat <- pval_v1[,id]
  or_mat <- or_v1[,id]
  idx <- apply(pval_mat, 1, which.min)
  pval_v2[,ct] <- apply(pval_mat, 1, min)
  or_mat <- or_v1[,grep("EC", colnames(or_v1))]
  or_v2[,ct] <- sapply(seq(length(idx)), function(i){
    or_mat[i,idx[i]]
  })
}

res_v2 <- list("Hypergeometric_nominal_P"=pval_v2,
               "Odds_ratio"=or_v2)
saveRDS(res_v2, file="Res_tIO_filtered_DEG_GO_biological_process_enrichment_res_union_on_species_change_direction.rds")
# get enriched GO terms with relaxing cutoffs
deg_enriched_term_list <- lapply(seq(ncol(pval_v2)), function(j){
  intersect(rownames(pval_v2)[which(pval_v2[,j]<0.05)], rownames(or_v2)[which(or_v2[,j]>1.5)])
})
names(deg_enriched_term_list) <- colnames(pval_v2)
sapply(deg_enriched_term_list, length)
saveRDS(deg_enriched_term_list, file="Res_tIO_filtered_DEG_with_direction_GO_BP_deg_enriched_term_list_union_on_species_change_direction.rds")

# load the GO enrichment for DARs classifying the species differential direction
great_res_list <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/DAP/control_regional_difference/rGREAT/Res_online_rGREAT_res_per_cell_type_DA_peaks.rds")

pval_dar <- matrix(NA, nrow=nrow(pval_v2), ncol=5)
rownames(pval_dar) <- rownames(pval_v2)
colnames(pval_dar) <- colnames(pval_v2)

dar_enriched_term_list <- list()
for(ct in names(great_res_list)){
  for(sp in c("Human", "Chimp")){
    res <- great_res_list[[ct]][[sp]][["GO Biological Process"]]
    enriched_terms <- res$Desc[which(res$BinomFdrQ<0.05 & res$RegionFoldEnrich>1.5)]
    enriched_terms <- intersect(enriched_terms, rownames(pval_dar))
    if(length(enriched_terms)<1) next
    if(length(dar_enriched_term_list[[ct]])<1){
      dar_enriched_term_list[[ct]] <- enriched_terms
    }else{
      dar_enriched_term_list[[ct]] <- union(enriched_terms,dar_enriched_term_list[[ct]])
    }
    
    pval_vec <- setNames(res$BinomFdrQ[match(enriched_terms, res$Desc)], enriched_terms)
    
    na_idx <- which(is.na(pval_dar[enriched_terms,ct]))
    pval_dar[enriched_terms[na_idx], ct] <- pval_vec[na_idx]
    
    non_na_idx <- which(!is.na(pval_dar[enriched_terms,ct]))
    pval_dar[enriched_terms[non_na_idx], ct] <- apply(cbind(pval_dar[enriched_terms[non_na_idx],ct], pval_vec[non_na_idx]), 1, min) 
    
  }
}
saveRDS(pval_dar, file="Res_tIO_DAR_pval_mat.rds")
saveRDS(dar_enriched_term_list, file="Res_tIO_DAR_with_direction_enriched_GO_BP_terms_by_cell_type.rds")






overlapped_enriched_terms <- lapply(names(dar_enriched_term_list), function(ct){
  intersect(dar_enriched_term_list[[ct]], deg_enriched_term_list[[ct]])
})
names(overlapped_enriched_terms) <- names(dar_enriched_term_list)
sapply(overlapped_enriched_terms, length)
saveRDS(overlapped_enriched_terms, file="Res_tIO_RNA_and_ATAC_overlapped_enriched_GO_terms_classifying_direction.rds")


union_enriched_terms <- lapply(names(dar_enriched_term_list), function(ct){
  union(dar_enriched_term_list[[ct]], deg_enriched_term_list[[ct]])
})
names(union_enriched_terms) <- names(dar_enriched_term_list)
sapply(union_enriched_terms, length)
saveRDS(union_enriched_terms, file="Res_tIO_RNA_or_ATAC_enriched_GO_terms_classifying_direction.rds")

# get the enterocyte DEGs and DARs involved in the following terms  
terms <- overlapped_enriched_terms[["Enterocyte"]]
x <- "Enterocyte"
go_genes <- GO_anno$HGNC.symbol[which(GO_anno$GO.term.name%in%terms)]
hc_degs <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_filtered_DEGs_with_input_gene_list.rds")
cluster_genes <- rownames(hc_degs)[which(hc_degs[,x]%in%c("Human_high","Chimp_high"))]
all_genes <- rownames(hc_degs)[which(hc_degs[,x]%in%c("Human_high","Chimp_high","Expressed"))]
genes <- intersect(cluster_genes, intersect(go_genes, all_genes))
length(genes)
plotFeature(seu.obj = tIO_rna, dr="umap_css", genes.to.plot = genes, plot.name = "Plot_UMAP_CSS_tIO_epi_lipid_related_DEG_expr.png")

# get the overlapped enriched GO terms classifying the species differential directions
overlapped_enriched_terms <- list()
for(ct in names(dar_enriched_term_list)){
  
}
overlapped_enriched_terms <- lapply(names(dar_enriched_term_list), function(ct){
  intersect(dar_enriched_term_list[[ct]], deg_enriched_term_list[[ct]])
})
names(overlapped_enriched_terms) <- names(dar_enriched_term_list)
sapply(overlapped_enriched_terms, length)
saveRDS(overlapped_enriched_terms, file="Res_tIO_RNA_and_ATAC_overlapped_enriched_GO_terms.rds")


setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/cell_type_vs_species_diff")
# load tIO scRNA-seq data
tIO_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
# load tIO cell type per species expression
sp_ct_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_species_cell_type_average_expr.rds")
scaled_sp_ct_expr <- t(scale(t(sp_ct_expr)))
saveRDS(scaled_sp_ct_expr, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_scaled_species_cell_type_average_expr.rds")

## highlight DEGs in stem cells and enterocytes
cell_type_diff <- scaled_sp_ct_expr[,"Enterocyte@Human"] - scaled_sp_ct_expr[,"Stem_cell@Human"]
stem_cell_high_features <- names(cell_type_diff)[which(cell_type_diff<=0)]
enterocyte_high_features <- names(cell_type_diff)[which(cell_type_diff>=0)]

stem_cell_sp_diff <- scaled_sp_ct_expr[,"Stem_cell@Human"] - scaled_sp_ct_expr[,"Stem_cell@Chimp"]
enterocyte_sp_diff <- scaled_sp_ct_expr[,"Enterocyte@Human"] - scaled_sp_ct_expr[,"Enterocyte@Chimp"]

# stem cell vs enterocyte in human
idx <- which(tIO_rna$Unified_cell_type_per_species%in%c("Stem_cell@Human","Enterocyte@Human"))
X <- tIO_rna@assays$RNA@data[,idx]
y <- tIO_rna$Unified_cell_type[idx]
stem_cell_enterocyte_deg_res <- presto::wilcoxauc(X=X, y=y)
stem_cell_enterocyte_deg_res$pct_diff <- stem_cell_enterocyte_deg_res$pct_in - stem_cell_enterocyte_deg_res$pct_out
stem_cell_enterocyte_sig_res <- stem_cell_enterocyte_deg_res[which(stem_cell_enterocyte_deg_res$auc>0.7 & stem_cell_enterocyte_deg_res$pct_diff>25 & stem_cell_enterocyte_deg_res$padj<0.01 & stem_cell_enterocyte_deg_res$logFC>0.3),]
table(stem_cell_enterocyte_sig_res$group)
genes <- stem_cell_enterocyte_sig_res$feature
output <- stem_cell_enterocyte_deg_res[which(stem_cell_enterocyte_deg_res$group=="Enterocyte"),]

# generate boxplot to show individual gene expression across cell types between species
selected_genes <- c("C1QTNF12","C3orf85","SLC2A2","RBP2", "AFP", "MTTP", "SMLR1", "SLC26A3", "SLC5A12", "C11orf86", "APOA4", "APOC3","SLC7A7", "CIDEB",  "IL32", "DNASE1", "CREB3L3", "ENPP7", "PPP1R14A","SULT2A1","NR1H4","MAX",
                    "HMGCS2", "IMPDH2","SLC12A2","NAP1L1","TKT","STMN1","MYC","ETV4")
output$Selected <- output$feature%in%selected_genes
saveRDS(output, file="Res_tHIO_stem_cell_vs_enterocyte_presto_output.rds")

groups <- sort(unique(tIO_rna$Unified_cell_type_per_species))
mat <- do.call('rbind', strsplit(groups, split = "@"))
mat[,1] <- sub("_", " ", mat[,1])
alpha <- setNames(c("30","FF"), c("Chimp", "Human"))
source("/projects/site/pred/ihb-intestine-evo/colors/colors.R")
tissue.epi.ct.cols <- setNames(paste0(epi.ct.cols[mat[,1]],alpha[mat[,2]]),
                               groups)

p_list <- list()
for(g in selected_genes){
  p_list[[g]] <- SCpubr::do_BoxPlot(tIO_rna, feature = g, group.by = "Unified_cell_type_per_species", colors.use = tissue.epi.ct.cols)
}

length(selected_genes)
row_num=5
col_num=6
pdf("Plot_boxplot_tHIO_stem_cell_and_enterocyte_markers.pdf", height=5*row_num, width=5*col_num)
do.call("grid.arrange", c(p_list, ncol = col_num))
dev.off()


# size the point by -log10(P) in enterocyte vs stem cell DA test
idx <- which(output$auc<0.5)
output$auc[idx] <- 1-output$auc[idx]
size <- setNames(output$auc^4*2, output$feature)
summary(size)

coor_enterocyte <- data.frame(
  "Enterocyte_vs_stem_cell"=cell_type_diff[enterocyte_high_features],
  "Human_vs_chimp"=enterocyte_sp_diff[enterocyte_high_features],
  "name"=enterocyte_high_features,
  "size"=size[enterocyte_high_features],
  stringsAsFactors = F
)

coor_stem_cell <- data.frame(
  "Enterocyte_vs_stem_cell"=cell_type_diff[stem_cell_high_features],
  "Human_vs_chimp"=stem_cell_sp_diff[stem_cell_high_features],
  "name"=stem_cell_high_features,
  "size"=size[stem_cell_high_features],
  stringsAsFactors = F
)

saveRDS(coor_stem_cell, file="Dat_coor_cell_type_vs_species_stem_cell_high_feature_genes.rds")
saveRDS(coor_enterocyte, file="Dat_coor_cell_type_vs_species_enterocyte_high_feature_genes.rds")


library(ggplot2)
library(ggrepel)
p_stem_cell <- ggplot(coor_stem_cell, aes(x=Enterocyte_vs_stem_cell, y=Human_vs_chimp, size=size, label=name))+
  geom_point(col="#d9d9d960")+
  geom_point(data=subset(coor_stem_cell, coor_stem_cell$name%in%selected_genes), 
             col=paste0(epi.ct.cols["Stem cell"],"60"))+
  geom_label_repel(data=subset(coor_stem_cell, coor_stem_cell$name%in%selected_genes),
                   max.overlaps = 9999)+
  theme_minimal()

p_enterocyte <- ggplot(coor_enterocyte, aes(x=Enterocyte_vs_stem_cell, y=Human_vs_chimp, size=size, label=name))+
  geom_point(col="#d9d9d960")+
  geom_point(data=subset(coor_enterocyte, coor_enterocyte$name%in%selected_genes), 
             col=paste0(epi.ct.cols["Enterocyte"],"60"))+
  geom_label_repel(data=subset(coor_enterocyte, coor_enterocyte$name%in%selected_genes),
                   max.overlaps = 9999)+
  theme_minimal()


library(gridExtra)
library(grid)
pdf("Plot_scatter_gene_expression_cell_type_vs_species_difference.pdf", height=7, width=7*2)
grid.arrange(p_stem_cell,p_enterocyte,ncol=2)
dev.off()


linked_regions <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/cor_DEG_and_DAR/Res_tIO_DAR_linked_with_DEG_1MB.rds")

sig_df <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/DE_vs_DA/GO/Res_tIO_filtered_DEG_without_considering_direction.rds")
selected_genes <- c("C1QTNF12","C3orf85","SLC2A2","RBP2", "AFP", "MTTP", "SMLR1", "SLC26A3", "SLC5A12", "C11orf86", "APOA4", "APOC3","SLC7A7", "CIDEB",  "IL32", "DNASE1", "CREB3L3", "ENPP7", "PPP1R14A","SULT2A1","NR1H4","MAX",
                    "HMGCS2", "IMPDH2","SLC12A2","NAP1L1","TKT","STMN1","MYC","ETV4")
hc_enterocyte_deg <- intersect(c("C1QTNF12","C3orf85","SLC2A2","RBP2", "AFP", "MTTP", "SMLR1", "SLC26A3", "SLC5A12", "C11orf86", "APOA4", "APOC3","SLC7A7", "CIDEB",  "IL32", "DNASE1", "CREB3L3", "ENPP7", "PPP1R14A","SULT2A1","NR1H4","MAX"), rownames(sig_df)[which(sig_df$Enterocyte)])
hc_stem_cell_deg <- intersect(c("HMGCS2", "IMPDH2","SLC12A2","NAP1L1","TKT","STMN1","MYC","ETV4"), rownames(sig_df)[which(sig_df$Stem_cell)])

# get DARs in enterocytes
human_chimp_dar <- readRDS("~/ihb-intestine-evo/used_object/differential_features/Res_tIO_DDP_by_cell_type_2.rds")
enterocyte_dar_df <- human_chimp_dar[["Enterocyte"]]
enterocyte_human_high_regions <- rownames(enterocyte_dar_df)[which(enterocyte_dar_df$Human_high)]

# get human stem cell and enterocyte marker regions
tHIO_stem_cell_enterocyte_res <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/species_vs_cellType/Res_tHIO_enterocyte_vs_stem_cell_differential_detection_res_no_wilcoxauc_filtering-2.rds")
tHIO_enterocyte_marker_regions <- rownames(tHIO_stem_cell_enterocyte_res)[which(tHIO_stem_cell_enterocyte_res$Enterocyte_high)]

human_high_enteroyte_marker_regions <- intersect(tHIO_enterocyte_marker_regions, enterocyte_human_high_regions)
length(human_high_enteroyte_marker_regions)

# get enterocyte DARs that correlate with selected enterocyte DE markers
enterocyte_link_df <- linked_regions[which(linked_regions$gene%in%hc_enterocyte_deg & linked_regions$peak_ID%in%human_high_enteroyte_marker_regions),]

df_enterocyte <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/species_vs_cellType/Dat_df_cell_type_vs_species_enterocyte_high_regions.rds")
enterocyte_regions <- intersect(df_enterocyte$name[which(df_enterocyte$X>0.25 & df_enterocyte$Y>0.25)], unique(enterocyte_link_df$peak_ID))
unique(enterocyte_link_df$gene[which(enterocyte_link_df$peak_ID%in%enterocyte_regions)])
length(enterocyte_regions)

# generate coverage plot for selected regions
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation.rds")
tCIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tCIO_scATAC_with_fragment_file.rds")
human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_human_and_tHIO_scATAC_with_fragment_file.rds")
selected_cell_type <- c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell","Stem cell")
vec1 <- rep(selected_cell_type, each=2)
vec2 <- rep(c("Fetal primary", "transplanted"), length(selected_cell_type))
selected_groups <- paste(vec1,
                         vec2,
                         sep="@")
alpha <- setNames(c("30","FF"), c("Fetal primary", "transplanted"))
tissue.epi.ct.cols <- setNames(paste0(epi.ct.cols[vec1],alpha[vec2]),
                               selected_groups)

selected_regions <- c("mergedPeak408441", "mergedPeak414224", "mergedPeak452851", "mergedPeak316300")
for(id in selected_regions){
  p_human <- rownames(peak_anno)[which(peak_anno$peak_id==id)]
  p_chimp <- peak_anno[p_human, "panTro6_coor"]
  
  # coverage plot in chimp
  p1 <- CoveragePlot(
    object = tCIO_atac,
    assay="CIO_unified_peaks",
    group.by = "Cell_type",
    region = p_chimp,
    idents = selected_cell_type,
    extend.upstream = 15000,
    extend.downstream = 20000,
    #features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE,
    window = 500
  )
  p1 <- p1 & scale_fill_manual(values = epi.ct.cols)
  
  # coverage plot in human
  p2 <- CoveragePlot(
    object = human_atac,
    assay="HIO_unified_peaks",
    group.by = "Cell_type_per_tissue",
    idents = selected_groups,
    region = p_human,
    extend.upstream = 15000,
    extend.downstream = 20000,
    #features = g,
    annotation = TRUE,
    peaks = TRUE,
    tile = FALSE,
    links = TRUE,
    window = 500
  )
  p2 <- p2 & scale_fill_manual(values = tissue.epi.ct.cols) 
  
  plot_name <- paste0("Plot_coveragePlot_",id,"_hg38_",p_human,"_panTro6_",p_chimp,"_15k_smoothed.pdf")
  pdf(plot_name, height=7, width=12)
  print(p1+p2)
  dev.off()
}

# 
slc5a12_nearby_selected_regions <- peak_anno$peak_id[which(peak_anno$Closest_gene=="SLC5A12" & rowSums(peak_anno[,41:51])>0)]
StringToGRanges(rownames(peak_anno)) -> gr_regions
target_regions <- StringToGRanges(c("chr11-26715000-26735000"))
res <- findOverlaps(query=target_regions, gr_regions)
matched_regions <- rownames(peak_anno)[res@to]
peak_anno[matched_regions,c(41:51,ncol(peak_anno))]

linked_regions[which(linked_regions$gene=="SLC5A12" & linked_regions$peak_ID%in%slc5a12_nearby_regions),]


snc_1 <- read.csv("/projects/site/pred/ihb-intestine-evo/lukas_area/evo_signatures/mike_pruefer_2013.hSNCs.hg38.csv",
                  row.names = 1)
snc_2 <- read.csv("/projects/site/pred/ihb-intestine-evo/lukas_area/evo_signatures/mike_unpublished.hSNCs.hg38.csv",
                  row.names = 1)
main_chr <- paste0("chr",c(seq(22),"X"))
idx1 <- which(snc_1$seqnames%in%main_chr)
idx2 <- which(snc_2$seqnames%in%main_chr)
snc <- union(paste(snc_1$seqnames[idx1], snc_1$start[idx1], snc_1$end[idx1], sep="-"),
             paste(snc_2$seqnames[idx2], snc_2$start[idx2], snc_2$end[idx2], sep="-"))
gr_snc_main <- StringToGRanges(snc)

tIO_peaks <- rownames(peak_anno)
tIO_peaks_main <- tIO_peaks[!grepl("chrY", tIO_peaks)]
gr_tIO_main <- StringToGRanges(tIO_peaks_main)

overlap_evo_sig <- findOverlaps(query = gr_snc_main, subject = gr_tIO_main)


# categorize the DA regions
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/group_DAR")
## summarize the information
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation.rds")
#colSums(peak_anno[,c(41:51,ncol(peak_anno))])
selected_cell_types <- c("BEST4+_epithelium","Enterocyte","Goblet_cell","Stem_cell","EEC")
# get DA regions
da_idx <- sapply(selected_cell_types, function(ct){
  human_vec <- as.numeric(peak_anno[,paste0("Human_high:", ct, "@fetal")])
  chimp_vec <- as.numeric(peak_anno[,paste0("Chimp_high:", ct, "@fetal")])
  idx <- ifelse(human_vec>chimp_vec, 1, ifelse(human_vec<chimp_vec, -1, 0))
  return(idx)
})
rownames(da_idx) <- rownames(peak_anno)
saveRDS(da_idx, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_DA_regions_per_cell_type.rds")

# get regions that are linked to DEGs
### load human-chimp DEGs
human_chimp_DEGs <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff/Res_human_chimp_DEGs_with_input_gene_list.rds")
idx <- (rowSums(human_chimp_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Human_high", na.rm=T)+rowSums(human_chimp_DEGs[,c("EC-cell", "non-EC-EEC", "EEC_progenitor")]=="Chimp_high", na.rm = T))>0
sig_mat <- sapply(c("BEST4+_epithelium", "Enterocyte", "Goblet_cell", "Stem_cell"), function(ct){
  human_chimp_DEGs[,ct]%in%c("Human_high", "Chimp_high")
})
rownames(sig_mat) <- rownames(human_chimp_DEGs)
sig_df <- data.frame(sig_mat, "EEC"=idx, stringsAsFactors = F)
colnames(sig_df)[1] <- "BEST4+_epithelium"
saveRDS(sig_df, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_DEG_without_considering_direction.rds")

## regions are linked to both the most nearby gene and genes located within 1MB range with significant correlation in tIO epi bimodal data
### load DA-DE linkage
linked_regions <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/cor_DEG_and_DAR/Res_tIO_DAR_linked_with_DEG_1MB.rds")
### load the ChIP-seeker annotation results
cs_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/human_peak_annotation/Data_frame_peak_gene_pairs.rds")
region_gene_pairs <- data.frame("gene"=c(linked_regions$gene, cs_anno$gene_symbol), 
                                "region"=c(linked_regions$peak, cs_anno$name))
region_gene_pairs <- region_gene_pairs[which(region_gene_pairs$gene%in%rownames(sig_df)),]
saveRDS(region_gene_pairs, file="/projects/site/pred/ihb-intestine-evo/used_object/Dat_tIO_region_gene_pairs.rds")


de_idx <- sig_df[region_gene_pairs$gene,]
region_linked_to_de_idx <- sapply(colnames(de_idx), function(ct){
  regions <- unique(region_gene_pairs$region[which(de_idx[,ct])])
  rownames(peak_anno)%in%regions
})
rownames(region_linked_to_de_idx) <- rownames(peak_anno)
saveRDS(region_linked_to_de_idx, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_regions_linked_to_DEGs_per_cell_type.rds")

dar_linked_to_deg_idx <- da_idx*region_linked_to_de_idx

## get detected peak lists
expressed_idx <- rowSums(peak_anno[,grep("Expressed:", colnames(peak_anno))])>0

## get cell type marker regions
marker_idx_by_cell_type <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/group_DAR/Res_cell_type_marker_index_by_cell_type.rds")
da_idx_by_cell_type <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/group_DAR/Res_DA_index_by_cell_type.rds")

## 
papers <- list(
  "UC"=c("Ultraconserved elements - Bejerano 2004",
         "Ultraconserved elements (zooUCEs) - Christmas 2023"),
  "HAR"=c("Human-accelerated regions (combined) - Doan 2016",
          "Human-accelerated DNase I hypersensitive sites - Gittelman 2015",
          "Human-accelerated regions (zooHARs) - Keough 2023",
          "Hhuman ancestor quickly evolved regions - Mangan 2022",
          "Human-accelerated regions - Pollard 2006",
          "Human-accelerated conserved non-coding sequences - Prabhakar 2006"),
  "HSD"=c("Human-specific deletions - McLean 2011", 
          "Human-specific deletions - Xue 2023"),
  "HSSV"=c("Human-specific segmental duplications - Dennis 2017",
           "Human-specific structural variants - Kronenberg 2018"),
  "Sweep"="Selective Sweep - Modern Humans - Peyregne 2017"
)

features <- data.frame("Expressed"=expressed_idx,
                       "DA"=rowSums(da_idx_by_cell_type)>0,
                       "Fixed"=peak_anno$Pruefer_et_al_fixed_SNC>0,
                       "HAR"=rowSums(peak_anno[,papers[["HAR"]]])>0,
                       "Sweep"=peak_anno[,papers[["Sweep"]]]>0,
                       "HSD"=rowSums(peak_anno[,papers[["HSD"]]])>0,
                       "HSSV"=rowSums(peak_anno[,papers[["HSSV"]]])>0,
                       "UC"=rowSums(peak_anno[,papers[["UC"]]])>0,
                       "Cell_type_enriched"=ifelse(rowSums(marker_idx_by_cell_type)>0,1,-1),
                       "DE"=ifelse(rowSums(region_linked_to_de_idx)>0,1,-1),
                       "DA_Human_high"=ifelse(rowSums(da_idx)>0, 1,-1))
saveRDS(features, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_region_features.rds")

# only plot DARs,otherwise too many dot and lines to plot
features <- features[which(features$DA),]
mat <- features[,c("DA_Human_high","Cell_type_enriched","DE")]
each_step <- sapply(seq(ncol(mat)), function(j){
  mat[,j]*2^(2-j)
})
score <- 4.5+sapply(seq(ncol(mat)), function(j){
  if(j==1){
    return(each_step[,j])
  }else{
    return(rowSums(each_step[,1:j]))
  }
})


y_coor <- cbind(rep(4.5, nrow(features)), score, matrix(rep(score[,ncol(score)], each=6), byrow = T, ncol=6))
colnames(y_coor) <- c("DA", "DA_Human_high", "Cell_type_enriched", "DE", "Fixed", "HAR", "Sweep", "HSD", "HSSV", "UC")
rownames(y_coor) <- rownames(features)
x_coor <- matrix(rep(seq(ncol(y_coor)), nrow(y_coor)),
                 byrow = T,
                 nrow=nrow(y_coor))
rownames(x_coor) <- rownames(y_coor)
colnames(x_coor) <- colnames(y_coor)

y_resi <- rnorm(nrow(score), 0, 0.1)
x_resi <- rnorm(nrow(score), 0, 0.05)

y_pos <- y_coor+y_resi
x_pos <- x_coor+x_resi

# highlight cell type enriched DA regions that are linked to DEGs
highlight_idx <- marker_idx_by_cell_type*da_idx_by_cell_type*region_linked_to_de_idx
saveRDS(highlight_idx, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_cell_type_enriched_DAR_linked_to_DEG.rds")

regions_to_highlight <- rownames(highlight_idx)[which(rowSums(highlight_idx)>0)]

plot_input <- list(
  "x_pos"=x_pos,
  "y_pos"=y_pos
)
saveRDS(plot_input, file="/projects/site/pred/ihb-intestine-evo/used_object/plotting_input/Res_humand_chimp_organoid_detected_region_catergory_plot_input.rds")

dar_group <- list(
  "x_coor"=x_coor,
  "y_coor"=y_coor
)
saveRDS(dar_group, file="/projects/site/pred/ihb-intestine-evo/used_object/plotting_input/Res_humand_chimp_organoid_detected_region_catergory_noJitter.rds")

png("Plot_DA_region_category.png", height=2000, width=2000)
plot(0,0,xlim=c(min(x_pos,na.rm=T),max(x_pos,na.rm=T)), ylim=c(min(y_pos,na.rm=T),max(y_pos,na.rm=T)), type="n", bty="n", yaxt="n", xaxt="n", xlab=NA, ylab=NA)
#axis(side = 1, at = seq(ncol(x_pos)), labels=F)
#mtext(side = 1, text = colnames(x_pos), at = seq(ncol(x_pos)), line = 2)
for(i in 1:nrow(x_pos)) lines(x_pos[i,1:4], y_pos[i,1:4], col=ifelse(rownames(x_pos)%in%regions_to_highlight, "#dedede30", "#dedede10"), lwd=0.5)
for(i in 1:4){
  points(x_pos[,i], y_pos[,i], cex=0.5, pch=16, col="#bdbdbd30")
  points(x_pos[rownames(x_pos)%in%regions_to_highlight,i], y_pos[rownames(x_pos)%in%regions_to_highlight,i], cex=0.6, pch=16, col="#30303070")
}
for(i in c("Fixed","HSSV","HAR","Sweep","HSD","UC")){
  points(x_pos[,i], y_pos[,i], cex=0.5, pch=16, col="#bdbdbd30")
  points(x_pos[which(features[rownames(x_pos),i]),i], 
         y_pos[which(features[rownames(x_pos),i]),i], 
         cex=0.6, 
         pch=16, 
         col=ifelse(i=="UC", "#0c2c8450", "#7a017790"))
}
dev.off()

pdf("Plot_DA_region_category_label_only.pdf")
plot(0,0,xlim=c(min(x_pos,na.rm=T),max(x_pos,na.rm=T)), ylim=c(min(y_pos,na.rm=T),max(y_pos,na.rm=T)), type="n", bty="n", yaxt="n", xaxt="n", xlab=NA, ylab=NA)
axis(side = 1, at = seq(ncol(x_pos)), labels=F)
mtext(side = 1, text = colnames(x_pos), at = seq(ncol(x_pos)), line = 2, las=2)
dev.off()

linked_regions <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/species_vs_cellType/Res_tIO_DAR_linked_with_DEG_100MB.rds")
selected_regions <- c("mergedPeak408441", "mergedPeak414224", "mergedPeak452851", "mergedPeak316300")
selected_genes <- c("PPP1R14A", "CREB3L3", "IL32", "DNASE1")
linked_regions
