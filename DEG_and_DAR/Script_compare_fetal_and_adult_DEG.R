setwd("/projects/site/pred/ihb-intestine-evo/fetal_and_adult")
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(gridExtra)
adult_enteroid_de_res <- readRDS("/projects/site/pred/ihb-intestine-evo/lukas_area/for_qianhui/adult.enteroids.human.marmoset.mouse.delegate.res.rds")
adult_enteroid_deg_group <- adult_enteroid_de_res$test_deg_group_mat
tIO_deg <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_DEGs_with_input_gene_list.rds")
fetal_deg_group <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/differential_features/Res_human_chimp_mouse_pairwise_DEG_combined_mat.rds")

fetal_id <- setNames(c("Stem_cell","Enterocyte", "Goblet_cell", "EEC", "EEC"),
               c("Stem_cell","Enterocyte", "Goblet_cell", "EC-cell", "non-EC-EEC"))
adult_id <- setNames(c("Stem_cell","Enterocyte", "Goblet_cell", "EEC"),
                     c("Stem cells","Enterocytes", "Goblet cells", "EECs"))

fetal_deg_list <- list()
for(group in c("human_specific_high", "human_specific_low", "conserved")){
  for(ct in c("Stem_cell","Enterocyte", "Goblet_cell", "EC-cell", "non-EC-EEC")){
    id <- paste(group, fetal_id[ct], sep="@")
    if(is.null(fetal_deg_list[[id]])){
      fetal_deg_list[[id]] <- rownames(fetal_deg_group)[which(fetal_deg_group[,paste0("H_C_M:",ct)]==group)]
    }else{
      fetal_deg_list[[id]] <- union(fetal_deg_list[[id]],
                                    rownames(fetal_deg_group)[which(fetal_deg_group[,paste0("H_C_M:",ct)]==group)])
    }
  }
}

adult_deg_list <- list()
for(group in c("Human_specific_high", "Human_specific_low", "conserved")){
  for(ct in colnames(adult_enteroid_deg_group)){
    id <- paste(tolower(group), adult_id[ct], sep="@")
    if(is.null(adult_deg_list[[id]])){
      adult_deg_list[[id]] <- rownames(adult_enteroid_deg_group)[which(adult_enteroid_deg_group[,ct]==group)]
    }else{
      adult_deg_list[[id]] <- union(adult_deg_list[[id]],
                                    rownames(adult_enteroid_deg_group)[which(adult_enteroid_deg_group[,ct]==group)])
    }
  }
}

shared_deg_list <- list()
for(x in names(fetal_deg_list)){
  shared_deg_list[[x]] <- intersect(fetal_deg_list[[x]], adult_deg_list[[x]])
}


for(ct in unique(fetal_id)){
  shared_deg_list[[paste("divergent",ct,sep="@")]] <- union(shared_deg_list[[paste("human_specific_high",ct,sep="@")]], shared_deg_list[[paste("human_specific_low",ct,sep="@")]])
}
deg_list <- list(
  "fetal"=fetal_deg_list,
  "adult_enteroid"=adult_deg_list,
  "shared"=shared_deg_list
)
saveRDS(deg_list, file="Dat_overlap_of_fetal_and_adult_DEG_list.rds")

# get adult cell type per species average expression levels
adult_enteroid <- readRDS("/projects/site/pred/ihb-intestine-evo/lukas_area/for_qianhui/adult.enteroids.human.marmoset.mouse.css.integrated.rds")
mat <- do.call('rbind', strsplit(adult_enteroid$cell_type, split = " (", fixed = T))
adult_enteroid$coarse_cell_type <- mat[,1]
adult_enteroid$coarse_cell_type <- adult_id[adult_enteroid$coarse_cell_type]
adult_enteroid$coarse_cell_type_species <- paste(adult_enteroid$coarse_cell_type, adult_enteroid$species, sep="@")
saveRDS(adult_enteroid, file="Res_adult_enteroids_human_marmoset_mouse_css_integrated_with_corase_cell_type.rds")
adult_enteroid_expr <- getAveExpr(seu.obj=adult_enteroid, feature.to.calc = "coarse_cell_type_species", colname.prefix = NULL)
saveRDS(adult_enteroid_expr, file="Dat_adult_enteroid_coarse_cell_type_species_expr.rds")
saveRDS(adult_enteroid_expr, file="/projects/site/pred/ihb-intestine-evo/Dat_adult_enteroid_coarse_cell_type_species_expr.rds")

# get fetal cell type per species average expression levels
fetal_rna <- readRDS("/projects/site/pred/ihb-intestine-evo/tHIO_tCIO_and_developed_fetal_human_and_mouse/update_annotation/exclude_distal_SI_mouse_cells/Res_human_chimp_and_mouse_subseted_epi_distinct_cell_types.rds")
# for non-EC-EEC and EC cells, take the average of the two subtypes
fetal_expr <- getAveExpr(seu.obj=fetal_rna, feature.to.calc = "Cell_type_per_species", colname.prefix = NULL)
expr_1 <- fetal_expr[,grep("Enterocyte|Goblet_cell|Stem_cell", colnames(fetal_expr))]
expr_2 <- sapply(sort(unique(fetal_rna$Species)), function(sp){
  rowMeans(fetal_expr[,paste(c("EC-cell","non-EC-EEC"),sp,sep="@")])
})
colnames(expr_2) <- paste("EEC", colnames(expr_2), sep="@")
bk <- fetal_expr
fetal_expr <- cbind(expr_1, expr_2)
saveRDS(fetal_expr, file="Dat_fetal_coarse_cell_type_species_expr.rds")
saveRDS(fetal_expr, file="/projects/site/pred/ihb-intestine-evo/Dat_fetal_coarse_cell_type_species_expr.rds")

fetal_rna <- subset(fetal_rna, subtype_group != "BEST4+ epithelium")
fetal_rna$Coarse_cell_type <- fetal_rna$subtype_group
fetal_rna$Coarse_cell_type[which(fetal_rna$subtype_group%in%c("EC-cell","non-EC-EEC"))] <- "EEC"
fetal_rna$Coarse_cell_type_species <- paste(fetal_rna$Coarse_cell_type, fetal_rna$Species, sep="@")
saveRDS(fetal_rna, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_human_chimp_mouse_non_BEST4_cell_distinct_cell_type.rds")

# generate violin plots for the shared human-specific (divergent) genes

# get human-mouse logFC in fetal
mat_fetal <- do.call('rbind', strsplit(colnames(fetal_expr), split="@"))
fetal_hm_logFC <- sapply(sort(unique(mat_fetal[,1])), function(ct){
  fetal_expr[,paste(ct,"Human",sep="@")] - fetal_expr[,paste(ct,"Mouse",sep="@")]
})
saveRDS(fetal_hm_logFC, file="Dat_fetal_human_mouse_logFC.rds")

# get human-mouse logFC in adult
mat_adult_enteroid <- do.call('rbind', strsplit(colnames(adult_enteroid_expr), split="@"))
adult_enteroid_hm_logFC <- sapply(sort(unique(mat_adult_enteroid[,1])), function(ct){
  adult_enteroid_expr[,paste(ct,"Human",sep="@")] - adult_enteroid_expr[,paste(ct,"Mouse",sep="@")]
})
saveRDS(adult_enteroid_hm_logFC, file="Dat_adult_enteroid_human_mouse_logFC.rds")

fetal_na_genes <- rownames(fetal_hm_logFC)[which(rowSums(is.na(fetal_hm_logFC))>0)]
adult_enteroid_na_genes <- rownames(adult_enteroid_hm_logFC)[which(rowSums(is.na(adult_enteroid_hm_logFC))>0)]
shared_genes <- setdiff(intersect(rownames(fetal_hm_logFC), rownames(adult_enteroid_hm_logFC)), union(fetal_na_genes, adult_enteroid_na_genes))
length(shared_genes)

colnames(adult_enteroid_hm_logFC) <- paste(colnames(adult_enteroid_hm_logFC), "adult_enteroid", sep="@")
colnames(fetal_hm_logFC) <- paste(colnames(fetal_hm_logFC), "fetal", sep="@")
df <- data.frame(shared_genes, adult_enteroid_hm_logFC[shared_genes,], fetal_hm_logFC[shared_genes,], stringsAsFactors = F)
saveRDS(df,  file="Dat_fetal_and_adult_enteroid_human_mouse_logFC_for_shared_genes.rds")


mat <- do.call('rbind', strsplit(colnames(df)[-1],split = ".", fixed = T))
plot_list <- list()
cell_types <- sort(unique(mat[,1]))
# only higlight top 10 logFC of either direction
sum_logfc <- sapply(cell_types, function(ct){
  rowSums(df[,grep(ct, colnames(df))])
})

highlight_gene_list <- lapply(colnames(sum_logfc), function(j){
  human_high_genes <- deg_list[["shared"]][[paste("human_specific_high",ct,sep="@")]]
  human_low_genes <- deg_list[["shared"]][[paste("human_specific_low",ct,sep="@")]]
  vec <- sum_logfc[human_high_genes,j]
  g_high <- human_high_genes[order(vec, decreasing = T)[1:min(c(10, length(human_high_genes)))]]
  vec <- sum_logfc[human_low_genes,j]
  g_low <- human_low_genes[order(vec, decreasing = T)[1:min(c(10, length(human_low_genes)))]]
  return(union(g_high, g_low))
  
})
names(highlight_gene_list) <- colnames(sum_logfc)
saveRDS(highlight_gene_list, file="Dat_highlight_gene_list.rds")


highlight_color <- "#A93226" 
for(ct in cell_types){
  df_ct <- df[,c("shared_genes", grep(ct, colnames(df), value = T))]
  genes_to_highlight <- highlight_gene_list[[ct]]
  colnames(df_ct)[-1] <- c("x","y")
  p <- ggplot(df_ct, aes(x=x, y=y)) + 
    geom_point(color = "grey50")+
    geom_point(data = df_ct[which(df_ct$shared_genes %in% genes_to_highlight),],
               color = highlight_color)+
    geom_smooth(method=lm, color=highlight_color, se=FALSE, linetype="dashed")+
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(df_ct$shared_genes %in% genes_to_highlight, shared_genes, "")),
                    aes(label = label), 
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    ylim = c(-Inf, Inf),
                    xlim = c(-Inf, Inf),
                    point.size = NA,
                    color=highlight_color)+ 
    geom_hline(yintercept=0, linetype="dashed")+
    geom_vline(xintercept=0, linetype="dashed")+
    labs(title=ct,
         x ="Adult ( Human - Mouse )", 
         y = "Developing ( Human - Mouse )")+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
    theme_minimal()
  plot_list[[ct]] <- p
}

pdf("Plot_scatter_plot_human-mouse_logFC_highlight_shared_human_specific_genes.pdf", height=7, width=28)
do.call("grid.arrange", c(plot_list, ncol = length(cell_types)))
dev.off()







# perform GO enrichment on conserved and genes with human-specific changes
expressed_genes <- intersect(rownames(fetal_deg_group), rownames(adult_enteroid_deg_group))
GO_anno <- read.csv("/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v109/Ensembl_v109_GO.csv")
idx <- which(GO_anno$GO.domain!="" & GO_anno$HGNC.symbol%in%expressed_genes)
GO_anno <- GO_anno[idx,]
all_genes <- intersect(expressed_genes, GO_anno$HGNC.symbol)

library(doParallel)
registerDoParallel(cores=10)
#for(group in sort(unique(GO_anno$GO.domain))){
  
  group <- "biological_process"
  print(paste(group,"start"))
  selected_gs <- sort(unique(GO_anno$GO.term.name[which(GO_anno$GO.domain==group)]))
  length(selected_gs)
  
  go_res <- list()
  for(deg_group in names(shared_deg_list)){
    res <- foreach(term=selected_gs, .multicombine = T, .combine = 'rbind', .maxcombine = length(selected_gs))%dopar%{
      degs <- shared_deg_list[[deg_group]]
      go_genes <- intersect(GO_anno$HGNC.symbol[which(GO_anno$GO.term.name==term)], all_genes)
      go_deg <- intersect(degs, go_genes)
      a <- length(go_deg)
      b <- length(degs)
      c <- length(go_genes)
      d <- length(all_genes)
      res <- fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")
      return(c(res$p.value, res$estimate, a, a/b))
    }
    res <- data.frame(res, stringsAsFactors = F)
    rownames(res) <- selected_gs
    names(res) <- c("Hypogeometric_test_nominal_P",
                    "Odds_ratio",
                    "Hit_count",
                    "Prop"
                  )
    res$"BH_corrected_P"=p.adjust(res$"Hypogeometric_test_nominal_P", method="BH")
    go_res[[deg_group]] <- res
  }
  saveRDS(go_res, file=paste0("Res_GO_",group,"_enrichment_res.rds"))
  
  # only visualize the enriched terms of human-specific genes in each cell type
  selected_go_terms <- lapply(cell_types, function(ct){
    res <- go_res[[paste0("divergent@",ct)]]
    sig_res <- res[which(res$Hypogeometric_test_nominal_P<0.01),]
    top_res <- sig_res[order(sig_res$Hit_count, decreasing = T)[1:min(c(10, nrow(sig_res)))],]
    top_res$terms <- rownames(top_res)
    top_res$cell_type <- ct
    return(top_res)
  })
  input <- do.call('rbind', selected_go_terms)
  input$cell_type <- sub("_", " ", input$cell_type)
  input$logP <- -log10(input$Hypogeometric_test_nominal_P)
  
  p1 <- ggdotchart(input, x = "terms", y = "logP",
             color = "cell_type",                                # Color by groups
             palette = epi.ct.cols[unique(input$cell_type)], # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             rotate = FALSE,                                # Rotate vertically
             group = "cell_type",                                # Order by groups
             dot.size = 15,                                 # Large dot size
             label = round(input$Hit_count),                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 15, 
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_pubr(),                        # ggplot2 theme
             ylab="-log10(P)",
             xlab="Enriched GO terms"
  )
  pdf("Plot_top_enriched_GO_terms.pdf", height=10, width=20)
  par(xpd=TRUE)
  p1
  dev.off()
  
  
  enriched_term_list <- lapply(names(go_res), function(x){
    mat <- go_res[[x]]
    rownames(mat)[which(mat$Hypogeometric_test_nominal_P<0.01)]
  })
  names(enriched_term_list) <- names(go_res)
  all_enriched_terms <- unique(unlist(enriched_term_list))
  
  or_mat <- sapply(names(go_res), function(x){
    go_res[[x]]$Odds_ratio
  })
  rownames(or_mat) <- selected_gs
  sig_num <- apply(or_mat, 1, function(vec){
    high_cutoff <- mean(vec)+2*sd(vec)
    sum(vec>high_cutoff)
  })
  
  
  
  all_enriched_terms <- unique(unlist(enriched_terms))
  length(all_enriched_terms)
  
  
  getDotPlot(cl.expr.mat=t(-log10(pval)), 
             cl.prop.mat=sqrt(t(or)), 
             gene.reorder=TRUE,
             cl.reorder.by.hc=TRUE, 
             genes=all_enriched_terms,
             point.size.factor=0.4,
             plot.cex=1,
             plot.margin = c(45,10,10,10),
             row.space.factor = 0.2,
             col.space.factor = 1,
             colors = beach.col,
             plot.name = paste0("Plot_dotplot_",group,"_enrichment.pdf") )
  
#}







