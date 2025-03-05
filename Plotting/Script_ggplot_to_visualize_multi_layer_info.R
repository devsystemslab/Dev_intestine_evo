setwd("/home/yuq22/ihb-intestine-evo/examples/lipid_metabolism")
library(Seurat)
library(Signac)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
source("/projects/site/pred/ihb-intestine-evo/colors/colors.R")

# load data
## peak annotation files
peak_anno <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_human_fetal_tissue_and_organoid_peak_annotation_df.rds")
## scATAC-seq data
tCIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tCIO_epi_scATAC_with_fragment_file.rds")
#human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/annotate_human_fetal_tissue_and_organoid_peaks/characterize_epi_markers/coveragePlot/Dat_fetal_tissue_tHIO_epi_scATAC-seq_with_fragment_files.rds")
#saveRDS(human_atac, file="/projects/site/pred/ihb-intestine-evo/used_object/Dat_fetal_tissue_tHIO_epi_scATAC-seq_with_fragment_files.rds")
human_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Dat_fetal_tissue_tHIO_epi_scATAC-seq_with_fragment_files.rds")

## developing intestine scRNA-seq data
#fetal <- readRDS("/home/yuq22/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/Res_updated_human_chimp_mouse_epi_integrated_with_CSS_seurat_obj.rds")
#cell_type_expr <- getAveExpr(seu.obj=fetal, feature.to.calc="Cell_type_per_group", colname.prefix=NULL)
#cell_type_expr <- cbind(cell_type_expr, rep(NA, nrow(cell_type_expr)))
#colnames(cell_type_expr)[ncol(cell_type_expr)] <- "BEST4+ epithelium@Mouse@Fetal primary" 
#mat <- do.call('rbind', strsplit(colnames(cell_type_expr),split="@"))
#mat <- cbind(mat, paste(mat[,2], mat[,3], sep="@"))
#groups <- sort(unique(mat[,4]))
#eec_expr <- sapply(groups, function(group){
#  rowMeans(cell_type_expr[,paste(c("EC-cell","non-EC-EEC"), group, sep="@")])
#})
#colnames(eec_expr) <- paste("EEC", colnames(eec_expr),sep="@")
#cell_type_expr <- cbind(cell_type_expr, eec_expr)
#selected_cell_type <- sort(c("BEST4+ epithelium", "EEC","Enterocyte", "Goblet cell","Stem cell"))
#selected_id <- paste(rep(selected_cell_type, each=length(groups)), rep(groups, length(selected_cell_type)), sep="@")
#cell_type_expr_used <- cell_type_expr[,selected_id]
#saveRDS(cell_type_expr_used, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_developed_human_and_mouse_tIO_epi_subset_cell_type_data.rds")
#cell_type_expr_scaled <- t(scale(t(cell_type_expr_used)))
#saveRDS(cell_type_expr_scaled, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_developed_human_and_mouse_tIO_epi_subset_cell_type_data_scaled.rds")
cell_type_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_developed_human_and_mouse_tIO_epi_subset_cell_type_data.rds")
cell_type_expr_scaled <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_developed_human_and_mouse_tIO_epi_subset_cell_type_data_scaled.rds")
mat <- do.call('rbind', strsplit(colnames(cell_type_expr_scaled),split="@"))
mat <- cbind(mat, paste(mat[,2], mat[,3], sep="@"))
groups <- unique(mat[,4])
group_short_name <- setNames(c("C","Ht","Ho","M"), groups)
selected_cell_type <- unique(mat[,1])

# load tIO epi scRNA-seq data
tIO_cell_type_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_species_cell_type_average_expr.rds")
mat <- do.call('rbind', strsplit(colnames(tIO_cell_type_expr),split="@"))
idx <- which(!mat[,1]%in%c("Early_enterocyte", "Paneth_cell", "Tuft_cell"))
tIO_cell_type_expr_sub <- tIO_cell_type_expr[,paste(mat[idx,1], mat[idx,2],sep="@")]
saveRDS(tIO_cell_type_expr_sub, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_species_cell_type_average_expr_sub.rds")
# combine tIO data with fetal human tissue data
human_tissue_cell_type_expr <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_fetal_human_tissue_scRNA_epi_cell_type_average_data.rds")
colnames(human_tissue_cell_type_expr) <- sub("_", " ", colnames(human_tissue_cell_type_expr))
human_tissue_cell_type_expr_sub <- human_tissue_cell_type_expr[,selected_cell_type]
colnames(human_tissue_cell_type_expr_sub) <- paste0(colnames(human_tissue_cell_type_expr_sub), "@Ht")
saveRDS(human_tissue_cell_type_expr_sub, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_fetal_human_tissue_scRNA_epi_cell_type_average_data_sub.rds")

shared_genes <- intersect(rownames(human_tissue_cell_type_expr), rownames(tIO_cell_type_expr_sub))
combined_expr <- cbind(tIO_cell_type_expr_sub[shared_genes,], human_tissue_cell_type_expr[shared_genes,])

tIO_cell_type_expr_sub <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_species_cell_type_average_expr_sub.rds")
tIO_cell_type_expr_scaled <- t(scale(t(tIO_cell_type_expr_sub)))
saveRDS(tIO_cell_type_expr_scaled, file="/projects/site/pred/ihb-intestine-evo/used_object/cell_type_average/Dat_tIO_epi_cell_type_expr_scaled")

# load DEG list
tIO_deg <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_DEGs_with_input_gene_list.rds")

# get peak region in interest
## genes in interest
g <- "GRAMD1B"
data <- peak_anno[which(peak_anno$Closest_gene==g1 & peak_anno$Composite_analysis_score>1),]
gr_peak <- StringToGRanges(rownames(data))
gr_snc_hg38 <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Prufer_et_all_SNC_hg38_grange_obj.rds")
overlaps <- findOverlaps(gr_snc_hg38, gr_peak)
peaks <- rownames(data)[unique(subjectHits(overlaps))]
# peaks associated with GRAMD1B
peaks <- c("chr11-123586576-123587368", "chr11-123587643-123588761")
overlap_pairs <- findOverlapPairs(gr_snc_hg38, StringToGRanges(peaks)) # to double check

g <- "SLC5A12"
data <- peak_anno[which(peak_anno$Closest_gene==g & peak_anno$Human_Chimp_DAP),]
data[,c("Closest_gene","evo_branch", "CE", "SNC", "HAR", "GWAS","Human_Chimp_DAP", "annotation_category")]
peaks <- rownames(data)
overlap_pairs <- findOverlapPairs(gr_snc_hg38, StringToGRanges(peaks)) # to double check
# "chr11-26730301-26731333" +/- 1M

peaks <- "chr16-3086151-3087659" # close to IL32, promoter of "AC108134.3",

g <- "CYP3A4"
data <- peak_anno[which(peak_anno$Closest_gene==g & peak_anno$Human_Chimp_DAP),]
data[,c("Closest_gene","evo_branch", "CE", "SNC", "HAR", "GWAS","Human_Chimp_DAP", "annotation_category")]
peaks <- rownames(data)
overlap_pairs <- findOverlapPairs(gr_snc_hg38, StringToGRanges(peaks)) # to double check
overlap_pairs
peaks <- "chr7-99792734-99793569" # only take the one that overlaps with SNC
peak_anno[peaks,]

# generate coverage plot for selected regions
selected_cell_type <- sort(c("BEST4+ cell", "EEC", "Enterocyte", "Goblet cell","Stem cell"))
vec1 <- rep(selected_cell_type, each=2)
vec2 <- rep(c("fetal_primary", "transplanted"), length(selected_cell_type))
selected_groups <- paste(vec1,
                         vec2,
                         sep="@")
tissue.epi.ct.cols <- setNames(epi.ct.cols[vec1], selected_groups)


# load SNC, HAR, PS and GWAS region
evo_list <- readRDS("/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_evolutionary_signature_regions.rds")
## 
gr_snc_hg38 <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/fixed_SNC/Prufer_et_all_SNC_hg38_grange_obj.rds")
evo_list[["All_human_SNC"]] <- gr_snc_hg38
evo_col <- setNames(c("#154360","#154360","#303030","#7B241C"), c("SNC","All_human_SNC","GWAS","CMS"))
saveRDS(evo_list, file="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_evolutionary_signature_regions_plus_all_Pruefer_SNC.rds")

composite_plot_list <- list()
for(p_human in peaks){
    print(p_human)
    g <- peak_anno[p_human,"Closest_gene"]
    p_chimp <- peak_anno[p_human, "panTro6_coor"]
    

    # 1. genomic signature
    p_list <- list()
    
    ## plot the overlapped sites
    p_human_vec <- strsplit(p_human, split="-")[[1]]
    extend_up <- 10000
    extend_down <- 10000
    region <- GRanges(
        seqnames=p_human_vec[1],
        ranges=IRanges(
          start=as.numeric(p_human_vec[2])-extend_up,
          end=as.numeric(p_human_vec[3])+extend_down
        )
    )
    start.pos <- start(x = region)
    end.pos <- end(x = region)
    chromosome <- seqnames(x = region)

    ## get the genomic coordinate of the genomic signatures overlapping with the human peak
    overlapped_coor <- lapply(names(evo_col), function(x){
      overlaps <- IRanges::findOverlaps(query=evo_list[[x]], subject=region)
      evo_list[[x]][queryHits(overlaps)]
    })
    names(overlapped_coor) <- names(evo_col)
    
    for(x in c("SNC","All_human_SNC", "GWAS","CMS","CE")){ # HAR will be in the coverage plot, here only for the single nucleotide stuff
      if(length(overlapped_coor[[x]])>0){
        start_vec <- as.numeric(do.call('rbind', strsplit(GRangesToString(overlapped_coor[[x]]),split="-"))[,2])
        peak.intersect <- subsetByOverlaps(x = overlapped_coor[[x]], ranges = region)
        peak.df <- as.data.frame(x = peak.intersect)
        peak.df$start[peak.df$start < start.pos] <- start.pos # clipped to the same x ranges
        peak.df$end[peak.df$end > end.pos] <- end.pos

        p_sig <- ggplot(
          data = peak.df
        ) 
        for(pos in start_vec){
          p_sig <- p_sig + geom_vline(xintercept=pos, color=evo_col[x], linetype=1, size=0.5)
        }

      }else{
        p_sig <-  ggplot(
            data=data.frame("seqnames"=chromosome, "start"=start.pos, "end"=end.pos, stringsAsFactors=F)
        )
      }
      p_sig <- p_sig + 
        theme_classic()+
        ylab(label = x) +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        xlab(label = paste0(chromosome, " position (bp)")) +
        xlim(c(start.pos, end.pos))+
        theme(legend.position = "none") # remove legend
      p_list[[x]] <- p_sig
    }

    # merge the signature plots
    p_sig_merged <- plot_grid(plotlist = p_list, align = "v",ncol = 1, axis = "t")
    
    
    # 2. coverage plot in human
    cov_human <- CoveragePlot(
      object = human_atac,
      assay="ATAC",
      group.by = "Unified_cell_type_per_tissue",
      idents = selected_groups,
      region = p_human,
      extend.upstream = 10000,
      extend.downstream = 10000,
      region.highlight = StringToGRanges(p_human),
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = TRUE,
      window = 500,
      ranges=overlapped_coor[["HAR"]],
      ranges.title="HAR"
    )
    cov_human <- cov_human & scale_fill_manual(values = tissue.epi.ct.cols) 

    # 3. overage plot in chimp
    plot_chimp <- CoveragePlot(
      object = tCIO_atac,
      assay="CIO_unified_peaks",
      group.by = "Cell_type",
      region = p_chimp,
      idents = selected_cell_type,
      extend.upstream = 10000,
      extend.downstream = 10000,
      region.highlight = StringToGRanges(p_chimp),
      annotation = FALSE,
      peaks = FALSE,
      tile = FALSE,
      links = FALSE,
      window = 500
    )
    plot_chimp <- plot_chimp & scale_fill_manual(values = epi.ct.cols)

    # 4. human, chimp, mouse cell type expression levels
    if(sum(g%in%rownames(cell_type_expr_scaled))>0){
      data <- data.frame(matrix(cell_type_expr_scaled[g,], byrow=T, nrow=length(selected_cell_type)))
      colnames(data) <- group_short_name
      rownames(data) <- selected_cell_type
    }else{
      data <- data.frame(matrix(tIO_cell_type_expr_scaled[g,], byrow=T, nrow=length(selected_cell_type)))
      colnames(data) <- c("C", "H")
      rownames(data) <- selected_cell_type
    }
    expr_dev <- ggheatmap::ggheatmap(data, color=darkBlue2Red.heatmap, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", border="000000")+
                scale_y_discrete(limits=rev(rownames(data)))+
                theme(legend.position = "none") # remove legend

  plot_combined <- CombineTracks(
    plotlist = list(cov_human, p_sig_merged, plot_chimp),
    expression.plot=expr_dev,
    heights = c(8, 6, 6),
    widths = c(10, 2)
  )
  composite_plot_list[[p_human]] <- plot_combined

 }

ncol=length(composite_plot_list)
nrow=1
p_merged <- plot_grid(plotlist = composite_plot_list, align = "hv",ncol = ncol)
plot_name <- paste0("Plot_coveragePlot_for_selected_regions_",g,".pdf")
pdf(plot_name, height=15, width=10*ncol)
p_merged
dev.off()

plot <- composite_plot_list[[1]]
ncol=1
plot_name <- paste0("Plot_coveragePlot_for_",names(composite_plot_list)[4],"_",g,".pdf")
pdf(plot_name, height=15, width=10*ncol)
plot
dev.off()
