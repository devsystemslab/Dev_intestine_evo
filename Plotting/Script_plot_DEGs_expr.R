setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/s2E_trajectory_and_across_cell_type/")

# box plot for gene expression in different species
setwd("/home/yuq22/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells")
fetal <- readRDS("/home/yuq22/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/Res_updated_human_chimp_mouse_epi_integrated_with_CSS_seurat_obj.rds")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(Seurat)
library(ggpubr)
# Basic box plot
idx <-  which(fetal$subtype_group=="Enterocyte")
df <- data.frame("expr"=fetal@assays$RNA@data["SLC5A12",idx],
"sp"=fetal$Species[idx],
stringsAsFactors=F)

query.features <- c("SLC5A12", "CD36")
cell_type_colors <- c("Enterocyte" = "#C31F3A")

for(query.feature in query.features){
    
    fetch_mat <- FetchData(fetal, c("subtype_group", "Species", query.features), slot="data")

    fetch_mat_df <- fetch_mat %>% 
      as.data.frame() %>% 
      rownames_to_column(var="cell") %>%
      pivot_longer(-c("cell", "subtype_group", "Species"), names_to="gene", values_to="expression")
      fetch_mat_df <- fetch_mat_df %>% 
      dplyr::filter(subtype_group=="Enterocyte") %>%
      dplyr::filter(Species!="Mouse")
    
    p1 <- ggplot(fetch_mat_df %>% 
              dplyr::filter(gene%in%query.feature), aes(x=Species, y=expression, fill=subtype_group, alpha=Species)) + 
            stat_boxplot(width=0.2, geom="errorbar") +
            geom_boxplot(outlier.shape = NA) + 
            coord_cartesian(ylim=as.vector(quantile(fetch_mat_df$expression[which(fetch_mat_df$gene==query.feature)], c(0.1, 0.9))))+
            theme_pubr() + 
            scale_alpha_manual(values=c(0.5, 1)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            scale_fill_manual(values = cell_type_colors) + 
            labs(x="", y="Normalized Expression") + 
            theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face="bold"))

    ggsave(paste0("/home/yuq22/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/", query.feature, ".pdf"), p1, width = 4, height = 4) 

}
# tHIOs vs tCIOs
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/s2E_trajectory_and_across_cell_type/tIO")
combined_expr_mat <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO/align_trajectory/aligned_to_human_combined_allow_flip/Res_combined_Pt_aligned_expr_mat.rds")
mat <- do.call('rbind', strsplit(colnames(combined_expr_mat), split = "_")) 
mat[,1]-> group_vec
mat[,2]-> time_vec
group.cols <- setNames(c("#fcc5c0","#fa9fb5","#c51b8a","#7a0177",
                         "#7fcdbb","#253494"),
                       unique(group_vec))


#g1 <- c("SI","SLC2A2", "SLC26A3", "CTSA")
#g1 <- c("MAX","CREB3L3","BHLHE40","AFF1",    "TULP4",   "CREB3L2", "TBX3",    "ZBTB4")
g1 <- c("SLC5A12", "ETV4", "MYC", 'SLC26A3')
g1 <- c("NR1H4", "CEBPA", "MAX", "MAF", "ZNF503", "ASCL2", "ONECUT2", "SOX9", "HMGA1", "HMGB2")
g1 <- c("CD36", "NR1H4", "CYP3A4")
row.num=1
col.num=length(g1)
plot.name="Plot_tIO_s2e_CD36_NR1H4_CYP3A4.pdf"
pdf(plot.name, height=5*row.num, width=5*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,5,5))
for(gene in g1){
  mean.vec <- combined_expr_mat[gene,]
  
  plot(time_vec, mean.vec, pch=16, col=group.cols[group_vec], xlab="Aligned Pt bin",ylab="Normed expr.", main=gene, bty="n", cex.main=2, cex.axis=2, cex.lab=2)
  for(group in unique(group_vec)){
    g.idx <- which(group_vec==group)
    g.time <- time_vec[g.idx]
    g.mean <- mean.vec[g.idx]
    lines(g.time, smooth.spline(g.time, g.mean, df=6)$y, col=group.cols[group], lwd=3)
  }
}
plot(time_vec, mean.vec, type="n", xlab="",ylab="", xaxt="n",yaxt="n", bty="n")
legend("topleft", legend=names(group.cols), text.col=group.cols, bty="n", cex=2)
dev.off()

# plot gene expression in human, chimp and mouse
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/s2E_trajectory_and_across_cell_type/human_chimp_and_mouse")
pt_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tIO_and_fetal_proxSI/tHIO_tCIO_and_developed_fetal_human_and_mouse/stem_cell_to_enterocyte/se_pt_alignment/Res_combined_Pt_aligned_expr_mat.rds")
mat <- do.call('rbind', strsplit(split="_", x=colnames(pt_expr)))
group.vec <- mat[,1]
time.vec <- as.numeric(mat[,2])
group.cols <- setNames(c("#c51b8a","#7fcdbb","#253494"),
                       c("Human", "Chimp", "Mouse"))
genes <- c("ETV4", "MYC", "SLC2A2", "MAX", "CREB3L3", "NR1H4", "MAF")
genes <- intersect(c("SLC27A5", "GAS6", "IL32", "SLC5A12", "ETV4", "MYC", "SLC2A2", "MAX", "CREB3L3", "NR1H4", "MAF"), rownames(pt_expr))
for(g in genes){
  expr.vec <- pt_expr[g,]
  pdf(paste0("Plot_",g,".pdf"))
  par(mar=c(5,5,5,5))
  plot(time.vec, expr.vec, type="n", xlab="Aligned Pt bin",ylab="Normed expr.", 
       main=g, bty="n", cex.main=2, cex.axis=2, cex.lab=2)
  for(group in unique(group.vec)){
    g.idx <- which(group.vec==group)
    g.time <- time.vec[g.idx]
    g.mean <- expr.vec[g.idx]
    lines(g.time, smooth.spline(g.time, g.mean, df=6)$y, col=group.cols[group], lwd=3)
  }
  legend("topright", legend=names(group.cols), text.col=group.cols, bty="n", cex=2)
  dev.off()
}


# plot gene expression in each individual sample of human, chimp and mouse scRNA-seq data
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/s2E_trajectory_and_across_cell_type/human_chimp_and_mouse")
by_sample_pt_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tIO_and_fetal_proxSI/tHIO_tCIO_and_developed_fetal_human_and_mouse/stem_cell_to_enterocyte/se_pt_alignment/Res_pt_by_sample_bin_average_aligned_expr.rds")
combined_pt_expr <- do.call('cbind', by_sample_pt_expr)
colnames(combined_pt_expr)[grep("Mm_E17-5", colnames(combined_pt_expr))] <- paste("Mm-E17.5-SI",seq(20), sep="_")
mat <- data.frame(do.call('rbind', strsplit(split="_", x=colnames(combined_pt_expr))),
                  stringsAsFactors = F)
colnames(mat) <- c("Sample", "Bin_idx")
species_vec <- rep("Dev.mouse", nrow(mat))
species_vec[grep("tHIO", mat[,1])] <- "tHIO"
species_vec[grep("tCIO", mat[,1])] <- "tCIO"
species_vec[grep("hDuo", mat[,1])] <- "Dev.human"
mat$Source <- species_vec
group.vec <- mat[,3]
time.vec <- as.numeric(mat[,2])
group.cols <- setNames(c("#cb4335","#f39c12","#16a085","#2980b9"),
                       c("Dev.human", "tHIO", "tCIO","Dev.mouse"))
cols <- group.cols[group.vec]
mat$Color <- cols
saveRDS(mat, file="Dat_human_chimp_and_mouse_per_sample_combined_aligned_pt_bin_meta_data.rds")

genes <- c("ETV4", "MYC", "SLC2A2", "MAX", "CREB3L3", "NR1H4", "MAF")
genes <- intersect(c("SLC27A5", "GAS6", "IL32", "SLC5A12", "ETV4", "MYC", "SLC2A2", "MAX", "CREB3L3", "NR1H4", "MAF", "SLC17A5", "SLC27A5"), rownames(combined_pt_expr))
for(g in genes){
  expr.vec <- combined_pt_expr[g,]
  pdf(paste0("Plot_",g,"_per_sample.pdf"))
  par(mar=c(5,5,5,5))
  plot(time.vec, expr.vec, type="n", xlab="Aligned Pt bin",ylab="Normed expr.", 
       main=g, bty="n", cex.main=2, cex.axis=2, cex.lab=2)
  for(sample in unique(mat[,1])){
    g.idx <- which(mat[,1]==sample)
    g.time <- time.vec[g.idx]
    g.mean <- expr.vec[g.idx]
    smooth.lines(g.time, smooth.spline(g.time, g.mean, df=6)$y, col=mat[g.idx,4], lwd=3)
  }
  legend("topright", legend=names(group.cols), text.col=group.cols, bty="n", cex=2)
  dev.off()
}
saveRDS(combined_pt_expr, file="Dat_human_chimp_and_mouse_per_sample_combined_aligned_pt_bin_expr.rds")


setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/s2E_trajectory/tIO")
# load GO annotation
GO_anno <- read.table("/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/GO/Ensembl_v109_hg38_GO.txt", sep="\t", fill = T, quote="", head=T)
# load human-chimp DEG list and extract those upregulated in human enterocytes or stem cells
human_chimp_deg <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_DEGs_with_input_gene_list.rds")
human_high_in_stem_and_enterocyte_genes <- rownames(human_chimp_deg)[which(human_chimp_deg[,"Enterocyte"]=="Human_high" | human_chimp_deg[,"Stem_cell"]=="Human_high")]

# RA and vitamin related
ra_vitamin_terms <- sort(unique(grep("retino|vitamin", GO_anno$GO.term.name, value = T)))
human_high_ra_vitamin_genes <- intersect(human_high_in_stem_and_enterocyte_genes, GO_anno$HGNC.symbol[which(GO_anno$GO.term.name%in%ra_vitamin_terms)])
length(human_high_ra_vitamin_genes)

# lipid related
lipid_terms <- sort(unique(grep("lipid|fatty acid|stero|ceramide|bile acid|triglyceride", GO_anno$GO.term.name, value = T)))
length(lipid_terms)
human_high_lipid_genes <- intersect(human_high_in_stem_and_enterocyte_genes, GO_anno$HGNC.symbol[which(GO_anno$GO.term.name%in%lipid_terms)])
length(human_high_lipid_genes)

# signaling genes
pattern_genes <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Annotation/Public_scRNA/Han_curated_signaling_gene_list_2020_Nat_comm/Dat_Nat_comm_Han_curated_signal_gene_list_human_symbol.rds")

# KEGG pathway
qusage::read.gmt("/projects/site/pred/ihb-intestine-evo/Annotation/MSigDB/v7.1/c2.cp.kegg.v7.1.symbols.gmt") -> kegg_anno
selected_terms <- grep("STARTCH|PPAR|STEROID|FATTY_ACID|LIPID", names(kegg_anno), value = T)
human_high_kegg_genes <- sort(intersect(human_high_in_stem_and_enterocyte_genes, unlist(kegg_anno[selected_terms])))
length(human_high_kegg_genes)


# plot their profiles along the s2E trajectory
combined_expr_mat <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/Chimp_tCIO_panTro6/C2_and_C7_tCIO/epi/with_tHIO/align_trajectory/aligned_to_human_combined_allow_flip/Res_combined_Pt_aligned_expr_mat.rds")
mat <- do.call('rbind', strsplit(colnames(combined_expr_mat), split = "_")) 
mat[,1]-> group_vec
mat[,2]-> time_vec
group.cols <- setNames(c("#fcc5c0","#fa9fb5","#c51b8a","#7a0177",
                         "#7fcdbb","#253494"),
                       unique(group_vec))
g1 <-human_high_lipid_genes 
# TFs
library(Pando)
data("motif2tf")
g1 <- intersect(sort(unique(motif2tf$tf)), human_high_in_stem_and_enterocyte_genes)
# pattern_genes
g1 <- intersect(pattern_genes$Human_gene_symbol, human_high_in_stem_and_enterocyte_genes)
#g1 <- c("TTR", "RBP2", "ADH1C", "AQP3", "DKK1", "AKR1C1", "CYP2C18", "GAS6")
g1 <- c("AGT", "ACAA1", "TPP1", "EBPL", "APOH", "ENPP7")
g1 <- c("STMN1", "HMGCS2", "NAP1L1", "IMPDH2", "TKT", "SLC12A2", "C11orf86", "C3orf85", "SLC2A2", "CREB3L3", "IL32", "PPP1R14A", "MTTP","AFP","SMLR1","SLC7A7","DNASE1","C1QTNF12","SULT2A1","CIDEB")
g1 <- human_high_kegg_genes
g1 <- c("ACAA1", "MSMO1", "SULT1E1", "SLC27A5", "UGT2B15", "STMN1", "C11orf86", "C3orf85", "IL32", "SLC2A2","DNASE1")
col.num=min(c(8,length(g1)))
row.num=ceiling(length(g1)/col.num)
plot.name="Plot_tIO_s2e_selected_human_high_stem_cell_enterocyte_selected_DEGs.pdf"
pdf(plot.name, height=5*row.num, width=5*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,5,5))
for(gene in g1){
  mean.vec <- combined_expr_mat[gene,]
  
  plot(time_vec, mean.vec, pch=16, col=group.cols[group_vec], xlab="Aligned Pt bin",ylab="Normed expr.", main=gene, bty="n", cex.main=2, cex.axis=2, cex.lab=2)
  for(group in unique(group_vec)){
    g.idx <- which(group_vec==group)
    g.time <- time_vec[g.idx]
    g.mean <- mean.vec[g.idx]
    lines(g.time, smooth.spline(g.time, g.mean, df=6)$y, col=group.cols[group], lwd=3)
  }
}
plot(time_vec, mean.vec, type="n", xlab="",ylab="", xaxt="n",yaxt="n", bty="n")
legend("topleft", legend=names(group.cols), text.col=group.cols, bty="n", cex=2)
dev.off()



pt_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/tIO_and_fetal_proxSI/tHIO_tCIO_and_developed_fetal_human_and_mouse/stem_cell_to_enterocyte/se_pt_alignment/Res_combined_Pt_aligned_expr_mat.rds")
mat <- do.call('rbind', strsplit(split="_", x=colnames(pt_expr)))
group.vec <- mat[,1]
time.vec <- as.numeric(mat[,2])
group.cols <- setNames(c("#c51b8a","#7fcdbb","#253494"),
                       c("Human", "Chimp", "Mouse"))
genes <- c("ETV4", "MYC", "SLC2A2", "MAX", "CREB3L3", "NR1H4", "MAF")

row.num=1
col.num=length(g1)
plot.name="Plot_human_chimp_mouse_s2e_selected_DEGs.pdf"
pdf(plot.name, height=5*row.num, width=5*col.num)
par(mfrow=c(row.num, col.num), mar=c(5,5,5,5))
for(gene in g1){
  expr.vec <- pt_expr[g,]
  plot(time.vec, expr.vec, type="n", xlab="Aligned Pt bin",ylab="Normed expr.", 
       main=g, bty="n", cex.main=2, cex.axis=2, cex.lab=2)
  for(group in unique(group.vec)){
    g.idx <- which(group.vec==group)
    g.time <- time.vec[g.idx]
    g.mean <- expr.vec[g.idx]
    lines(g.time, smooth.spline(g.time, g.mean, df=6)$y, col=group.cols[group], lwd=3)
  }
}
plot(time_vec, expr.vec, type="n", xlab="",ylab="", xaxt="n",yaxt="n", bty="n")
legend("topleft", legend=names(group.cols), text.col=group.cols, bty="n", cex=2)

dev.off()

# generate boxplot to show individual gene expression across cell types between species
tIO_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
selected_genes <- c("MAX", "CREB3L3", "ETV4", "MYC",
                    "AGT", "ACAA1", "TPP1", "EBPL", "APOH", "ENPP7",
                    "SLC5A12", "SLC2A2", "SLC26A3", "SLC2A2","SLC27A5",
                    "IL32","DNASE1","TTR","RBP2", "ADH1C","AQP3",
                    "DKK1","AKR1C1","CYP2C18","GAS6","C11orf86","C3orf85",
                    "ACAA1","MSMO1","SULT1E1","UGT2B15","STMN1"
                    )

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
library(gridExtra)
pdf("Plot_boxplot_tIO_selected_stem_cell_and_enterocyte_DEGs.pdf", height=5*row_num, width=5*col_num)
do.call("grid.arrange", c(p_list, ncol = col_num))
dev.off()

# get the DA regions that are close and show significant correlation with the above selected DEGs
tIO_epi_subset <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/cor_DEG_and_DAR/Res_tIO_epi_subset_DEG_with_Linked_DARs_1M.rds")
region_link <- tIO_epi_subset@assays$ATAC_DAR_only@links

# load human-chimp DARs
human_chimp_dar <- readRDS("~/ihb-intestine-evo/used_object/differential_features/Res_tIO_DDP_by_cell_type_2.rds")
# get the stem cell and enterocyte DARs
dar <- unique(unlist(lapply(c("Stem_cell", "Enterocyte"), function(j){
  mat <- human_chimp_dar[[j]]
  rownames(mat)[which(mat$Human_high)]
})))
orth_peaks <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/Res_combined_peaks_in_panTro6_hg38_coor.rds")
dar_hg38 <- orth_peaks$combined_peaks_hg38_coor[which(orth_peaks$peak_ID%in%dar)]

deg_cor_dar <- unique(region_link$peak[which(region_link$gene%in%selected_genes & region_link$peak%in%dar_hg38)])


# generate coverage plot for selected regions
setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/s2E_trajectory_and_across_cell_type/tIO/coverage_plot")
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

selected_regions <- deg_cor_dar
selected_regions <- c("chr1-25500668-25501944",
                      "chr1-230171443-230172509",
                      "chr3-170764926-170766002",
                      "chr9-33445385-33447963",
                      "chr10-5609541-5611177",
                      "chr11-66350609-66351906",
                      "chr13-113938482-113940344",
                      "chr13-113274989-113275898",
                      "chr13-113836816-113838180",
                      "chr13-113860046-113861235",
                      "chr13-113848640-113849203",
                      "chr13-113938482-113940344",
                      
                      
                      "chr16-2799175-2800554",
                      "chr16-2903931-2906161",
                      "chr16-3086151-3087637",
                      "chr19-57810453-57811319",
                      "chr10-5609541-5611177")


selected_regions <- c(
                      
                      "chr16-2799175-2800554",
                      "chr16-2903931-2906161",
                      "chr16-3086151-3087637",
                      "chr19-57810453-57811319",
                      "chr10-5609541-5611177")

IDs <- c("mergedPeak322039",
        "mergedPeak375040",
        "mergedPeak375174",
        "mergedPeak375193",
        "mergedPeak406714",
        "mergedPeak406763",
        "mergedPeak287960")

selection_signature <- c("Human-specific segmental duplications - Dennis 2017",               "Human-accelerated regions (combined) - Doan 2016",                 
                         "Human-accelerated DNase I hypersensitive sites - Gittelman 2015",   "Human-accelerated regions (zooHARs) - Keough 2023",                
                         "Human-specific structural variants - Kronenberg 2018",              "Hhuman ancestor quickly evolved regions - Mangan 2022",            
                         "Human-specific deletions - McLean 2011",                            "Selective Sweep - Modern Humans - Peyregne 2017",                  
                         "Human-accelerated regions - Pollard 2006",                          "Human-accelerated conserved non-coding sequences - Prabhakar 2006",
                         "Human-specific deletions - Xue 2023", "Pruefer_et_al_fixed_SNC")

regions <- union(region_link$peak[which(region_link$gene %in% selected_genes)], rownames(peak_anno)[which(peak_anno$Closest_gene%in%selected_genes)])
mat <- peak_anno[regions, selection_signature]
freq <- rowSums(mat)
names(freq)[which(freq>0)]



library(ggplot2)
for(id in IDs){
  p_human <- rownames(peak_anno)[which(peak_anno$peak_id==id)]
  p_chimp <- peak_anno[p_human, "panTro6_coor"]
  id <- peak_anno[p_human, "peak_id"]
  # coverage plot in chimp
  p1 <- CoveragePlot(
    object = tCIO_atac,
    assay="CIO_unified_peaks",
    group.by = "Cell_type",
    region = p_chimp,
    idents = selected_cell_type,
    extend.upstream = 15000,
    extend.downstream = 15000,
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
    extend.downstream = 15000,
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


id <- c("mergedPeak375174", # GAS6 # overlapped with fixed SNC
        "mergedPeak406763", # IL32 # overlapped with fixed SNC
        "mergedPeak316304"  # SLC5A12 # overlapped with fixed SNC
        ) 
selected_regions <- rownames(peak_anno)[which(peak_anno$peak_id%in%id)]

# the TFBS motif of which TFs that located in the above regions overlapped with fixed SNCs
idx_mat <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/motif_analysis/on_combined_peak_lists/Dat_overlapping_of_each_TFBS_of_each_human_scATAC-seq_region_with_published_fixed_SNC.rds")


# accessibility profile of regions overlapped with different categories of selection signature
# get the per sample per cell type detection rate of each of those regions
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
det_rate <- getExpressedProp(seu.obj=tIO_atac, feature.to.calc = "High_resolution_cell_type_per_species", colname.prefix = NULL, assay.type = "peaks_species")
det_rate <- det_rate[,setdiff(colnames(det_rate), "Putative_secretory_progenitors@Chimp")]
saveRDS(det_rate, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_atac_high_resolution_cell_type_detection_rate.rds")

# cluster the profile of those regions















