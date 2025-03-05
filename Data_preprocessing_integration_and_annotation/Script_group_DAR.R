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
