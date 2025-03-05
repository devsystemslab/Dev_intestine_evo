setwd("/home/yuq22/ihb-intestine-evo/examples/lipid_metabolism")
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
g1 <- c("CD36", "NR1H4")
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

# load human-chimp DEG list
human_chimp_DEG <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_filtered_DEGs_with_input_gene_list.rds")
enterocyte_deg <- rownames(human_chimp_DEG)[which(human_chimp_DEG[,"Enterocyte"] %in% c("Human_high" , "Chimp_high"))]

# load GRN inferred from human/chimp organoid data to get the NR1H4 targets
tIO_grn <- readRDS("/home/yuq22/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/Pando_on_tIO_ct_markers_and_DEGs/Res_tIO_epi_Pando_glm_model_grn.rds")
# extract the inferred module data
grn_meta <- tIO_grn@grn@networks$glm_network@modules@meta
saveRDS(grn_meta, file="/home/yuq22/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/Pando_on_tIO_ct_markers_and_DEGs/Res_tIO_epi_Pando_glm_model_grn_meta.rds")
grn_meta$target[which(grn_meta$tf=="NR1H4" & grn_meta$target%in%enterocyte_deg)]

species_ct_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/adult_primate_intestine_atlas/primate_duo/analysis/across_species/human_chimp_epi_respective_genome/all_cell_type_DE/Dat_species_ct_normalized_expr.rds")
species_ct_expr[c("SLC5A12", "CD36"),grep("Enterocyte", colnames(species_ct_expr))]

# get regulatory regions associated with CD36 and NR1H4 and generate coverage plots in developing data
library(Signac)
library(ggplot2)
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation-2.2.rds")
#g1 <- "CD36"
#regions <- rownames(peak_anno)[which(peak_anno$Closest_gene==g1 & (peak_anno$"Linked_by_cor@fetal" | peak_anno$"Linked_by_Pando@fetal") & peak_anno$"Union_cell_type_marker:Enterocyte")] 

regions <- "chr7-99792617-99793308"

tCIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tCIO_epi_scATAC_with_fragment_file.rds")
human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_human_and_tHIO_scATAC_with_fragment_file.rds")
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
      extend.upstream = 10000,
      extend.downstream = 10000,
      region.highlight = StringToGRanges(p_chimp),
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
      extend.upstream = 10000,
      extend.downstream = 10000,
      region.highlight = StringToGRanges(p_human),
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = TRUE,
      window = 500
    )
    p2 <- p2 & scale_fill_manual(values = tissue.epi.ct.cols) 
    
    plot_name <- paste0("Plot_coveragePlot_hg38_",p_human,"_panTro6_",p_chimp,".pdf")
    pdf(plot_name, height=7, width=12)
    print(p1+p2)
    dev.off()
 }

# get example human open chromatin regions emerged at evolutionary young nodes and overlap with acceleration/selection genomic signatures
setwd("/home/yuq22/ihb-intestine-evo/examples/pos_select_enterocyte_regions")
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation-2.2.rds")
table(peak_anno$node_assignment_singleCon)
selected_nodes <- c("chimp/bonobo", "gorilla", "macaque", "marmoset")
primate_regions <- rownames(peak_anno)[which(peak_anno$node_assignment_singleCon %in% selected_nodes)]

# get SNP positions under positive selection
data <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Dat_SNPs_under_positive_selection_hg38.rds")
## need to remove ancestral regions
pos_select <- data[which(data$"Genome-wide CMS" != "ancestral"),]
saveRDS(data, file="/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Table_S6_hg38_coor.rds")
saveRDS(pos_select, file="/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Dat_SNPs_under_positive_selection_hg38.rds")

gr_region <- StringToGRanges(rownames(peak_anno))
overlaps <- IRanges::findOverlaps(gr_region, pos_select)

aa <- data.frame(gr_region[overlaps@from], pos_select[overlaps@to])
aa$organoid_regions_hg38 <- rownames(peak_anno)[overlaps@from]
aa$Closest_gene <- peak_anno$Closest_gene[overlaps@from]
saveRDS(aa, file="/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Res_human_chimp_iPSC_organoid_region_under_pos_selection.rds")

# get number of positively selected SNPs in each human and chimp organoid regions
freq <- table(overlaps@from)
names(freq) <- rownames(peak_anno)[as.numeric(names(freq))]
peak_anno$CMS_pos_select <- 0
peak_anno[names(freq), "CMS_pos_select"] <- freq
saveRDS(peak_anno, file="/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation-2.2_with_pos_selection.rds")

# get regions that are enriched in enterocytes and under positive selection
regions <- rownames(peak_anno)[which(peak_anno$"Union_cell_type_marker:Enterocyte" & peak_anno$CMS_pos_select>0)]
output <- peak_anno[regions, ]
write.table(output, file="Table_human_or_chimp_enterocyte_enriched_regions_under_positive_selection.txt", sep="\t", quote=F)
output <- aa[match(regions, aa$organoid_regions_hg38),]
write.table(output, file="Table_human_or_chimp_enterocyte_enriched_regions_under_positive_selection_with_CMS_statistics.txt", sep="\t", quote=F, row.names=F)

table(peak_anno[regions, "Closest_gene"])
table(peak_anno[regions, "node_assignment_singleCon"])
primate_pos_regions <- intersect(regions, primate_regions)
length(primate_pos_regions)
table(peak_anno[primate_pos_regions, "Closest_gene"])

library(Signac)
library(ggplot2)

tCIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tCIO_epi_scATAC_with_fragment_file.rds")
human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_human_and_tHIO_scATAC_with_fragment_file.rds")
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

 for(p_human in primate_pos_regions){
    
    g <- peak_anno[p_human,"Closest_gene"]
    p_chimp <- peak_anno[p_human, "panTro6_coor"]
    
    # coverage plot in chimp
    p1 <- CoveragePlot(
      object = tCIO_atac,
      assay="CIO_unified_peaks",
      group.by = "Cell_type",
      region = p_chimp,
      idents = selected_cell_type,
      extend.upstream = 10000,
      extend.downstream = 10000,
      region.highlight = StringToGRanges(p_chimp),
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
      extend.upstream = 10000,
      extend.downstream = 10000,
      region.highlight = StringToGRanges(p_human),
      ranges=pos_select,
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = TRUE,
      window = 500
    )
    p2 <- p2 & scale_fill_manual(values = tissue.epi.ct.cols) 
    
    plot_name <- paste0("Plot_coveragePlot_hg38_",p_human,"_panTro6_",p_chimp,".pdf")
    pdf(plot_name, height=7, width=12)
    print(p1+p2)
    dev.off()
 }


# get the per cell type differential accessibiliy levels for each human/chimp organoid peaks
setwd("/home/yuq22/ihb-intestine-evo/annotate_human_chimp_organoid_peaks")
# load tIO scATAC-seq data
library(Signac)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
tIO_atac <- readRDS("/home/yuq22/ihb-intestine-evo/used_object/Res_tIO_epi_ATAC_seurat_object_exclude_putative_secretory_progenitors_and_early_enterocyte.rds")
cell_type_acc <- getAveExpr(seu.obj=tIO_atac, feature.to.calc="High_resolution_cell_type_per_species", assay.type="peaks_species", data.type="data", colname.prefix=NULL)
saveRDS(cell_type_acc, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_atac_per_species_cell_type_average_accessibility_level.rds")

acc_diff <- sapply(sort(unique(tIO_atac$High_resolution_cell_type)), function(ct){
  cell_type_acc[,paste0(ct,"@Human")] - cell_type_acc[,paste0(ct,"@Chimp")]
})
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation-2.2_with_pos_selection.rds")
rownames(acc_diff) <- rownames(peak_anno)
saveRDS(acc_diff, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_atac_Human-Chimp_per_cell_type_average_accessibility_diff.rds")

max_sp_diff <- apply(abs(acc_diff), 1, max)
peak_anno$Max_sp_diff <- max_sp_diff

# get human cell type accesibility specificity
tau_score <- calc_tau(cell_type_acc[,grep("Human",colnames(cell_type_acc))], byRow = TRUE)
peak_anno$tHIO_epi_tau <- tau_score
saveRDS(peak_anno, "/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_tIO_merged_peak_annotation-2.2_with_pos_selection.rds")

snp <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Table_S6_hg38_coor.rds")
gr_region <- StringToGRanges(rownames(peak_anno))
overlaps <- IRanges::findOverlaps(gr_region, snp)

aa <- data.frame(gr_region[overlaps@from], snp[overlaps@to])
aa$organoid_regions_hg38 <- rownames(peak_anno)[overlaps@from]
aa$Closest_gene <- peak_anno$Closest_gene[overlaps@from]
aa$Max_sp_diff <- peak_anno$Max_sp_diff[overlaps@from]
aa$tHIO_epi_tau  <- peak_anno$tHIO_epi_tau[overlaps@from]
saveRDS(aa, file="/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Res_human_chimp_iPSC_organoid_region_with_CMS_score.rds")

idx <- which(aa$Genome.wide.CMS!="ancestral")
pdf("Plot_association.pdf", height=5, width=15)
par(mfrow=c(1,3))
plot(aa$Normalized.CMS, aa$Max_sp_diff, pch=16)
points(aa$Normalized.CMS[idx], aa$Max_sp_diff[idx], col="dark red", pch=16)

plot(aa$Normalized.CMS, aa$tHIO_epi_tau, pch=16)
points(aa$Normalized.CMS[idx], aa$tHIO_epi_tau[idx], col="dark red", pch=16)

plot(aa$tHIO_epi_tau, aa$Max_sp_diff, pch=16)
points(aa$tHIO_epi_tau[idx], aa$Max_sp_diff[idx], col="dark red", pch=16)
dev.off()


# visualize Gene Ontology relationships
setwd("/home/yuq22/ihb-intestine-evo/annotate_human_chimp_organoid_peaks/GO_graph")
go_url <-
    "https://purl.obolibrary.org/obo/go/go-basic.obo"
path <- "go-basic_Feb2024.obo"
httr::GET(go_url, httr::write_disk(path, overwrite = TRUE))
obo <- OmnipathR::obo_parser(path, tables = TRUE)
df <- data.frame(obo$relations)
saveRDS(df, file="Dat_GO_term_relationships.rds")
freq <- sapply(seq(nrow(df)), function(i){
  length(unlist(df[i,3]))
})
all_go_pairs <- data.frame(
  "from"=rep(df$term, freq), 
  "to"=unlist(df$parents),
  stringsAsFactors=F
)
saveRDS(all_go_pairs, file="Dat_directed_GO_term_pairs.rds")
igraph_obj <- graph_from_data_frame(all_go_pairs)
saveRDS(igraph_obj, file="Res_directed_GO_term_igraph_obj.rds")

library(igraph)
library(ggraph)
library(tidygraph)
# biological process - GO:0008150
# cellular component - GO:0005575
# molecular function - GO:0003674
root_terms <- c("GO:0008150", "GO:0005575", "GO:0003674")

# test-run on DA enriched terms
GO_anno <- read.csv("/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v109/Ensembl_v109_GO.csv")
el <- readRDS("/home/yuq22/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/ATAC/peak_annotation/DE_vs_DA/GO/Res_DAR_enriche_GO_terms_igraph_input_edge_list.rds")
selected_terms <- intersect(df$term, GO_anno$GO.term.accession[which(GO_anno$GO.term.name%in%el[,1])])
node_list <- lapply(selected_terms, function(x){
  paths <- all_shortest_paths(graph=igraph_obj, from=x, to=root_terms[1])
  unique(names(unlist(paths$vpaths)))
})
kept_nodes <- unlist(node_list)
kept_pairs <- all_go_pairs[which(all_go_pairs$from %in% kept_nodes & all_go_pairs$to %in% kept_nodes),]

graph_obj <- graph_from_data_frame(kept_pairs)
layout <- layout.reingold.tilford(graph_obj, flip.y=T, root=which(V(graph_obj)$name %in% root_terms), mode="all")
saveRDS(layout, file="Res_DAR_enriched_BP_term_graph.rds")
plot(layout, pch=16)
rownames(layout) <- V(graph_obj)$name
node_df <- data.frame(
  "name"=rownames(layout),
  "X"=layout[,1],
  "Y"=layout[,2],
  #"size_pval"=pval_vec,
  #"size_or"=or_vec,
  #"group"=group_vec,
  stringsAsFactors = F
)

ggraph_obj <- as_tbl_graph(graph_obj)
ggraph_obj <- ggraph_obj %>% 
  tidygraph::activate(nodes) %>%
  left_join(node_df, by = c("name" = "name")) %>% 
  tidygraph::activate(edges) %>%
  mutate(from_name = (.N()$name[from])) %>%
  mutate(to_name = (.N()$name[to])) %>%
  mutate(pair_name = paste(from_name, to_name, sep=":")) %>%
  tidygraph::activate(nodes)

p1 <- ggraph(ggraph_obj, x=X, y=Y) +
  geom_edge_diagonal(color="#969696", width=0.5) +
  geom_node_point(aes(x=X, y=Y, filter=name%in%selected_terms),
                  shape = 19, color='darkgrey') +
  geom_node_label(aes(label=name, filter = name%in%root_terms),
                 size=3, repel=T, max.overlaps = 13) +
  theme_void()

png("Plot_ggraph_test.png", height=2000, width=2000)
p1
dev.off()

#p1 <- ggraph(ggraph_obj, x=X, y=Y) +
#  geom_edge_diagonal(aes(alpha=or, color=factor(direction)),
#                     width=0.5, arrow = arrow(length = unit(1,"mm"))) +
#  scale_edge_alpha_continuous(range=c(0.1,0.8), guide = "none") +
#  scale_edge_color_manual(values = c('-1'='#7FB3D5', '1'='#EC7063')) +
#  geom_node_point(aes(size = size_or, fill = group),
#                  shape = 21, color='darkgrey') +
#  scale_fill_manual(values = group_cols) +
#  scale_size_continuous(range = c(1,5), trans = "sqrt") +
#  geom_node_point(aes(size = size_or, filter = name%in%top_terms),
#                  shape = 1, color='#303030') +
#  geom_node_label(aes(label=name, filter = group=="Cell_type_label"),
#                 size=3, repel=T, max.overlaps = 13) +
#  geom_node_text(aes(label=name, filter = name %in% top_terms),
#                  min.segment.length = unit(0, 'lines'), size=3, repel=T, max.overlaps = 99999) +
#  theme_void()

# fetal human epi enriched regions
setwd("/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/ATAC/GO_enrichment")
human_enriched_regions <- readRDS("/pstore/data/ihb-intestine-evo/lukas_area/for_qianhui/fetal_tissue_pan_and_ct_enriched_peaks.rds")
selected_cell_types <- c(
  "Epithelial",
  "Stem_cell",
  "Enterocyte",           
  "Goblet cell",
  "EEC",                   
  "Early_enterocyte"
  )
# run GREAT analysis
library(rGREAT)
library(Signac)
great_res_list <- list()
gene_peak_pair_list <- list()
for(x in selected_cell_types){
  print(x)
  human_peaks <- human_enriched_regions$feature[which(human_enriched_regions$group==x)]
  gr <- StringToGRanges(human_peaks)
  job = submitGreatJob(gr, species = "hg38")
  tbl = getEnrichmentTables(job, download_by = "tsv")
  res = plotRegionGeneAssociationGraphs(job)
  great_res_list[[x]] <- tbl

}
saveRDS(great_res_list, file="Res_online_rGREAT_res_fetal_human_per_cell_type_enriched_peaks.rds")

# trial 1: get enriched GO terms all epi cell types
# trial 2: get epi and enterocyte enriched peaks
great_res_list <- readRDS("Res_online_rGREAT_res_fetal_human_per_cell_type_enriched_peaks.rds")
great_res_list <- great_res_list[c("Epithelial", "Enterocyte")]

pval_list <- list()
or_list <- list()

all_top_terms <- c()
for(ct in names(great_res_list)){
  res <- do.call('rbind', great_res_list[[ct]])
  all_top_terms <- union(all_top_terms, res$ID)
  pval_list[[ct]] <- setNames(res$BinomFdrQ, res$ID)
  or_list[[ct]] <-  setNames(res$RegionFoldEnrich, res$ID)
}

pval_mat <- matrix(NA, nrow=length(all_top_terms), ncol=length(pval_list))
rownames(pval_mat) <- all_top_terms
colnames(pval_mat) <- names(pval_list)

or_mat <- matrix(NA, nrow=length(all_top_terms), ncol=length(pval_list))
rownames(or_mat) <- all_top_terms
colnames(or_mat) <- names(pval_list)

for(ct in names(pval_list)){
  pval_mat[names(pval_list[[ct]]), ct] <- pval_list[[ct]]
  or_mat[names(pval_list[[ct]]), ct] <- or_list[[ct]]
}

# fill the NA in pval matrix with 1, and in odds ratio matrix with 0
pval_mat <- apply(pval_mat, 2, function(vec){
  vec[is.na(vec)] <- 1
  return(vec)
})
or_mat <- apply(or_mat, 2, function(vec){
  vec[is.na(vec)] <- 0
  return(vec)
})
saveRDS(pval_mat, file="Res_fetal_human_enterocyte_enriched_region_GREAT_pval_mat_filled_with_1.rds")
saveRDS(or_mat, file="Res_fetal_human_enterocyte_enriched_region_GREAT_odds_ratio_mat_filled_with_0.rds")

# get enriched GO terms
sig_idx <- (pval_mat<0.05) * (or_mat>2)
saveRDS(sig_idx, file="Res_fetal_human_enterocyte_enriched_region_GO_enrichment_sig_idx.rds")

selected_terms <- rownames(sig_idx)[which(rowSums(sig_idx)>0)]
saveRDS(selected_terms, file="Res_fetal_human_enterocyte_enriched_region_enriched_GO_terms.rds")

## visualize GO terms based on GO ontonlogy relationships
setwd("/home/yuq22/ihb-intestine-evo/Annotation/GO/Gene_ontology_relationship")
go_url <-
    "https://purl.obolibrary.org/obo/go/go-basic.obo"
path <- "go-basic_Feb2024.obo"
httr::GET(go_url, httr::write_disk(path, overwrite = TRUE))
obo <- OmnipathR::obo_parser(path, tables = TRUE)
df <- obo$names
saveRDS(df, file="Dat_GO_ID_and_names.rds")
#df <- data.frame(obo$relations, stringsAsFactors=F)
#saveRDS(df, file="Dat_GO_term_relationships.rds")
## only preserve "is_a" and "part_of" relationships, trim out "regulate" relationship
#df <- df[which(df$relation %in% c("is_a", "part_of")),]
#saveRDS(df, file="Dat_GO_term_isA_partOf_relationships.rds")
#
#freq <- sapply(seq(nrow(df)), function(i){
#  length(unlist(df[i,3]))
#})
#all_go_pairs <- data.frame(
#  "from"=c(as.character(rep(df$term, freq)), root_terms), 
#  "to"=c(as.character(unlist(df$parents)), rep("ROOT", 3)), # link all root terms to a fake root
#  stringsAsFactors=F
#)
#saveRDS(all_go_pairs, file="Dat_directed_GO_term_pairs.rds")
#saveRDS(all_go_pairs, file="/home/yuq22/ihb-intestine-evo/used_object/Dat_directed_GO_term_pairs.rds")
#igraph_obj <- graph_from_data_frame(all_go_pairs)
#saveRDS(igraph_obj, file="Res_directed_GO_term_igraph_obj.rds")
#saveRDS(igraph_obj, file="/home/yuq22/ihb-intestine-evo/used_object/Res_directed_GO_term_igraph_obj.rds")


setwd("/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/integrate_fetal_multiome_and_tHIO/ATAC/GO_enrichment")
library(igraph)
library(ggraph)
library(tidygraph)
# biological process - GO:0008150
# cellular component - GO:0005575
# molecular function - GO:0003674
all_go_pairs <- readRDS("/home/yuq22/ihb-intestine-evo/used_object/Dat_directed_GO_term_pairs.rds")
igraph_obj <- readRDS("/home/yuq22/ihb-intestine-evo/used_object/Res_directed_GO_term_igraph_obj.rds")
root_terms <- "GO:0008150" # focus on BP
selected_terms <- readRDS("Res_fetal_human_enterocyte_enriched_region_enriched_GO_terms.rds")
aa <- intersect(selected_terms, as.character(unlist(all_go_pairs)))
saveRDS(aa, file="Res_fetal_human_enterocyte_enriched_region_enriched_GO_terms_in_graph.rds")
length(aa) # 1360 out of the 1458 terms are preserved, for example, regulation were removed
length(selected_terms)
bk <- selected_terms
selected_terms <- aa
length(selected_terms)

node_list <- lapply(selected_terms, function(x){
  paths <- all_shortest_paths(graph=igraph_obj, from=x, to=root_terms) 
  unique(names(unlist(paths$vpaths))) 
})
kept_nodes <- unlist(node_list)
kept_pairs <- all_go_pairs[which(all_go_pairs$from %in% kept_nodes & all_go_pairs$to %in% kept_nodes),]

graph_obj <- graph_from_data_frame(kept_pairs)
layout <- layout.reingold.tilford(graph_obj, flip.y=T, root=which(V(graph_obj)$name %in% root_terms), mode="all")
rownames(layout) <- V(graph_obj)$name
saveRDS(layout, file="Res_fetal_human_enterocyte_enriched_BP_term_graph.rds")
plot(layout, pch=16)
df <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/GO/Gene_ontology_relationship/Dat_GO_ID_and_names.rds")
node_df <- data.frame(
  "name"=rownames(layout),
  "X"=layout[,1],
  "layer"=layout[,2],
  "Y"=layout[,2]+rnorm(nrow(layout), sd=0.1),
  "desc"=df$name[match(rownames(layout), df$term)],
  #"size_pval"=pval_vec,
  #"size_or"=or_vec,
  #"group"=group_vec,
  stringsAsFactors = F
)
saveRDS(node_df, file="Res_enterocyte_enriched_BP_term_node_df.rds")

ggraph_obj <- as_tbl_graph(graph_obj)
ggraph_obj <- ggraph_obj %>% 
  tidygraph::activate(nodes) %>%
  left_join(node_df, by = c("name" = "name")) %>% 
  tidygraph::activate(edges) %>%
  mutate(from_name = (.N()$name[from])) %>%
  mutate(to_name = (.N()$name[to])) %>%
  mutate(pair_name = paste(from_name, to_name, sep=":")) %>%
  tidygraph::activate(nodes)

p1 <- ggraph(ggraph_obj, x=X, y=Y) +
  geom_edge_diagonal(color="#969696", width=0.1) +
  geom_node_point(aes(x=X, y=Y, filter=name%in%selected_terms),
                  shape = 16, size=0.5, color='darkgrey') +
  geom_node_label(aes(label=desc, filter = layer>8 | layer<3),
                 size=2, repel=T, max.overlaps = 13) +
  theme_void()

pdf("Plot_ggraph_fetal_human_enterocyte_enriched_BP.pdf")
p1
dev.off()

#p1 <- ggraph(ggraph_obj, x=X, y=Y) +
#  geom_edge_diagonal(aes(alpha=or, color=factor(direction)),
#                     width=0.5, arrow = arrow(length = unit(1,"mm"))) +
#  scale_edge_alpha_continuous(range=c(0.1,0.8), guide = "none") +
#  scale_edge_color_manual(values = c('-1'='#7FB3D5', '1'='#EC7063')) +
#  geom_node_point(aes(size = size_or, fill = group),
#                  shape = 21, color='darkgrey') +
#  scale_fill_manual(values = group_cols) +
#  scale_size_continuous(range = c(1,5), trans = "sqrt") +
#  geom_node_point(aes(size = size_or, filter = name%in%top_terms),
#                  shape = 1, color='#303030') +
#  geom_node_label(aes(label=name, filter = group=="Cell_type_label"),
#                 size=3, repel=T, max.overlaps = 13) +
#  geom_node_text(aes(label=name, filter = name %in% top_terms),
#                  min.segment.length = unit(0, 'lines'), size=3, repel=T, max.overlaps = 99999) +
#  theme_void()

library(Signac)
setwd("/home/yuq22/ihb-intestine-evo/annotate_human_fetal_tissue_and_organoid_peaks")
# load CMS data
pos_select <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Dat_SNPs_under_positive_selection_hg38.rds")
# load evolutionary dating results
peak_age <- readRDS("/projects/site/pred/ihb-intestine-evo/lukas_area/for_qianhui/fetal_tissue_peaks_dated.rds")
# load scATAC-seq data
human_epi_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_tissue_tHIO_and_enteroid_epi_scATAC-seq_data_with_cell_type_annotation_and_fragment_files.rds")
human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_fetal_tissue_atac_with_all_cell_class_annotation.rds")
all_peaks <- rownames(human_atac@assays$ATAC@data)
gr <- human_atac@assays$ATAC@ranges
# get positively selected loci that overlap with human intestine open chromatin regions 
library(IRanges)
overlap <- findOverlaps(subject=gr, query=pos_select)
gr$name <- all_peaks

# get number of positively selected SNPs in each human open chromatin region
freq <- table(subjectHits(overlap))
names(freq) <- all_peaks[as.numeric(names(freq))]
gr$CMS_pos_select <- 0
idx <- match(names(freq),gr$name)
gr$CMS_pos_select[idx] <- freq
saveRDS(gr, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_human_fetal_tissue_and_organoid_peak_annotation.rds")
aa <- data.frame(gr[subjectHits(overlap)], pos_select[queryHits(overlap)])
saveRDS(aa, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_human_fetal_tissue_and_organoid_peak_under_pos_selection.rds")

# perform GREAT enrichment analysis on the overlapping SNPs
library(rGREAT)
library(Signac)

great_input <- pos_select[unique(queryHits(overlap))]
job = submitGreatJob(great_input, species = "hg38")
tbl = getEnrichmentTables(job, download_by = "tsv")
pdf("Plot_fetal_human_intestine_region_under_selection.pdf", height=5, width=15)
plotRegionGeneAssociationGraphs(job)
dev.off()
great_res <- tbl
great_res <- do.call('rbind', tbl)
# get enriched GO terms
great_res$sig_idx <- (great_res$BinomFdrQ<0.05) & (great_res$RegionFoldEnrich>2)
saveRDS(great_res, file="Res_online_rGREAT_res_fetal_human_peaks_overlap_with_positively_selected_loci.rds")

pos_select_terms <- great_res$ID[which(great_res$sig_idx)]
pos_select_term_names <- great_res$Desc[which(great_res$sig_idx)]

# get the path from the enriched terms back to BP
all_go_pairs <- readRDS("/home/yuq22/ihb-intestine-evo/used_object/Dat_directed_GO_term_pairs.rds")
igraph_obj <- readRDS("/home/yuq22/ihb-intestine-evo/used_object/Res_directed_GO_term_igraph_obj.rds")
root_terms <- "GO:0008150" # focus on BP
selected_terms <- pos_select_terms
aa <- intersect(selected_terms, as.character(unlist(all_go_pairs)))
saveRDS(aa, file="Res_fetal_human_intestine_CMS_selection_region_enriched_GO_terms_in_graph.rds")
length(aa) # 1360 out of the 1458 terms are preserved, for example, regulation were removed
length(selected_terms)
bk <- selected_terms
selected_terms <- aa
length(selected_terms)

node_list <- lapply(selected_terms, function(x){
  paths <- all_shortest_paths(graph=igraph_obj, from=x, to=root_terms) 
  unique(names(unlist(paths$vpaths))) 
})
kept_nodes <- unlist(node_list)
saveRDS(kept_nodes, file="Res_nodes_from_fetal_human_intestine_CMS_selection_region_enriched_GO_terms_to_BP_root_in_graph.rds")

# get the overlapping terms between enterocyte terms and CMS terms
node_df$CMS_path <- node_df$name%in%kept_nodes
length(node_df$desc[node_df$CMS_path & node_df$layer>6])

# load Mike's unpublished modern human specific fixed SNC
modern_human_snc <- read.table("/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Mike_unpublished_fixed_unique_cat_4col.bed", sep="\t", stringsAsFactors=F) # hg19
names(modern_human_snc) <- c("chr", "start", "end", "name")
gr_modern_snc <- makeGRangesFromDataFrame(modern_human_snc,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",,
                         start.field="start",
                         end.field="end", 
                         starts.in.df.are.0based=TRUE)
# convert from hg19 to hg38 coor
library(rtracklayer)
path = "/home/yuq22/ihb-intestine-evo/Annotation/liftOver_chain_file/hg19ToHg38.over.chain"
ch <- import.chain(path)
seqlevelsStyle(gr_modern_snc) <- "UCSC"  # necessary
gr_modern_snc_hg38 <- liftOver(gr_modern_snc, ch)
gr_modern_snc_hg38 <- unlist(gr_modern_snc_hg38)
saveRDS(gr_modern_snc_hg38, file="/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Mike_unpublished_modern_human_specific_fixed_SNC_hg38_grange_obj.rds")

# load Prufer's fixed SNC data
fixed_snc <- read.table("/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Prüfer_et_al_HumanDerived_SNC_bothgq30.all_combined_maxsco_ranked.tsv", sep="\t", comment.char = "", head=T)
fixed_snc$hg19_pos <- paste(fixed_snc$X.Chrom, fixed_snc$Pos, sep=":")
# modern human specific
fixed_snc$modern_human_specific_and_fixed_0.9 <- 0
fixed_snc$modern_human_specific <- 0
fixed_snc$modern_human_specific_and_fixed_0.995 <- 0

fixed_snc$modern_human_specific_and_fixed_0.9[which(fixed_snc$Archaic_States.AltaiNea.Denisova.=="A/A,A/A" & fixed_snc$Human_Maj!=fixed_snc$Anc)] <- 1
fixed_snc$modern_human_specific[which(fixed_snc$Archaic_States.AltaiNea.Denisova.=="A/A,A/A")] <- 1
fixed_snc$modern_human_specific_and_fixed_0.995[which(fixed_snc$Archaic_States.AltaiNea.Denisova.=="A/A,A/A" & fixed_snc$Human_Maj!=fixed_snc$Anc & fixed_snc$Human_Maj_1000G.Global.>0.995)] <- 1

gr_snc_hg19  <- makeGRangesFromDataFrame(fixed_snc,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="X.Chrom",,
                         start.field="Pos",
                         end.field="Pos", 
                         starts.in.df.are.0based=FALSE)
seqlevelsStyle(gr_snc_hg19) <- "UCSC"  # necessary
gr_snc_hg38 <- liftOver(gr_snc_hg19, ch)
gr_snc_hg38 <- unlist(gr_snc_hg38)
gr_snc_hg38$hg38_pos <- GRangesToString(gr_snc_hg38)
saveRDS(gr_snc_hg38, file="/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Prufer_et_all_SNC_hg38_grange_obj.rds")


#prufer_modern_human_snc_hg19 <- fixed_snc[idx2,]
#gr_prufer_modern_human_snc_hg19  <- makeGRangesFromDataFrame(prufer_modern_human_snc_hg19,
#                         keep.extra.columns=TRUE,
#                         ignore.strand=TRUE,
#                         seqinfo=NULL,
#                         seqnames.field="X.Chrom",,
#                         start.field="Pos",
#                         end.field="Pos", 
#                         starts.in.df.are.0based=FALSE)
#seqlevelsStyle(gr_prufer_modern_human_snc_hg19) <- "UCSC"  # necessary
#gr_prufer_modern_human_snc_hg38 <- liftOver(gr_prufer_modern_human_snc_hg19, ch)
#gr_prufer_modern_human_snc_hg38 <- unlist(gr_prufer_modern_human_snc_hg38)
#saveRDS(gr_prufer_modern_human_snc_hg38, file="/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Prufer_et_al_modern_human_specific_fixed_SNC_hg38_grange_obj.rds")

# get overlaps between human intestine open chromatin regions and SNCs
overlap <- findOverlaps(subject=gr, query=gr_snc_hg38)
pairs <- data.frame(
  "peaks"=gr$name[subjectHits(overlap)], 
  "SNC"=gr_snc_hg38$hg38_pos[queryHits(overlap)],
  stringsAsFactors=F
)
pairs$Prufer_modern_human_specific_and_fixed_0.9 <- pairs$SNC%in%gr_snc_hg38$hg38_pos[which(gr_snc_hg38$modern_human_specific_and_fixed_0.9>0)]
pairs$Prufer_modern_human_specific_and_fixed_0.995 <- pairs$SNC%in%gr_snc_hg38$hg38_pos[which(gr_snc_hg38$modern_human_specific_and_fixed_0.995>0)]
pairs$Prufer_modern_human_specific <- pairs$SNC%in%gr_snc_hg38$hg38_pos[which(gr_snc_hg38$modern_human_specific>0)]
saveRDS(pairs, file="Res_overlap_between_human_intestine_open_chromatin_regions_and_Prufer_SNCs.rds")

gr <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_human_fetal_tissue_and_organoid_peak_annotation.rds")
pairs <- readRDS("Res_overlap_between_human_intestine_open_chromatin_regions_and_Prufer_SNCs.rds")
# get number of SNCs in each human intestine open chromatin region
df <- data.frame(gr)
for(x in colnames(pairs)[-c(1,2)]){
  vec <- pairs$peaks[which(pairs[,x])]
  freq <- table(vec)
  df[,x] <- 0
  idx <- match(names(freq), df$name)
  df[idx,x] <- freq  
}
saveRDS(df, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_human_fetal_tissue_and_organoid_peak_annotation_df.rds")

# get overlap between 
gr_modern_snc_hg38 <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Mike_unpublished_modern_human_specific_fixed_SNC_hg38_grange_obj.rds")


saveRDS(df, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_human_fetal_tissue_and_organoid_peak_annotation.rds")


