setwd("/home/yuq22/ihb-intestine-evo/for_revision/spheroids/spheroid_stage_markers")
library(Seurat)
library(dplyr)
library(ggplot2)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

hio <- readRDS("/home/yuq22/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/used_seurat_objects/Res_updated_H9-HIO_before_transplantation_timecourse_SCC_based_MG_mapping.rds")

deg_res <- presto::wilcoxauc(hio, group_by="Cell_type")
deg_res$pct_diff <- deg_res$pct_in - deg_res$pct_out
sig_res <- deg_res %>% filter(padj < 0.05 & pct_diff > 10 & pct_in > 10 & logFC > 0.1)
saveRDS(sig_res, file="Res_H9_timecourse_cell_type_markers.rds")

deg_res <- presto::wilcoxauc(hio, group_by="Age")
deg_res$pct_diff <- deg_res$pct_in - deg_res$pct_out
sig_res <- deg_res %>% filter(padj < 0.05 & pct_diff > 10 & pct_in > 10 & logFC > 0.1)
saveRDS(sig_res, file="Res_H9_timecourse_time_point_markers.rds")


# load single spheroid (ss) data
ss <- readRDS("/home/yuq22/ihb-intestine-evo/for_revision/spheroids/Res_spheroid_with_UMAP_and_diffusion_map_and_age_marker_and_age_scores.rds")
SCpubr::do_DimPlot(ss, reduction="mds", group.by="Species", pt.size=3)
expr <- readRDS("/home/yuq22/ihb-intestine-evo/for_revision/spheroids/Dat_Date_by_species_average_expression.rds")
sig_res <- readRDS("Res_H9_timecourse_time_point_markers.rds")
top_res <- sig_res %>% group_by(group) %>% top_n(50, wt=logFC)
features <- intersect(unique(top_res$feature), rownames(expr))

# compare the spheroid data to the time course data
hio_expr <- getAveExpr(seu.obj=hio, feature.to.calc="Age", colname.prefix=NULL)
saveRDS(hio_expr, file="Dat_H9_timecourse_age_expression.rds")

ref_expr <- hio_expr[features,paste0("d", c(-5, 0, 3, 7, 14, 28, 30))]
que_expr <- expr[features,]
cor_mat <- cor(que_expr, ref_expr, method="spearman")
pdf("Plot_time_point_marker_SCC_in_spheroids_vs_timecourse.pdf", width=5, height=5)
gplots::heatmap.2(cor_mat, trace="none",main="", density.info="none",
          scale="none", Colv=FALSE, key=TRUE, 
          col=darkBlue2Red.heatmap, cexRow=1, cexCol = 1)
dev.off()


# compare the human spheroid data to chimp ones
features <- intersect(rownames(expr), sig_res$feature[which(sig_res$group=="d0")])
cor_mat <- cor(expr[features, grep("Human", colnames(expr))], expr[features, grep("Chimp", colnames(expr))])
pdf("Plot_stage_marker_PCC_in_spheroids.pdf", width=5, height=5)
gplots::heatmap.2(cor_mat, trace="none",main="", density.info="none",
          dendrogram="none", Rowv=FALSE, Colv=FALSE, scale="none", key=TRUE, 
          col=darkBlue2Red.heatmap, cexRow=1, cexCol = 1)
dev.off()

res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC) %>% filter(group=="d0")
genes <- intersect(res$feature, rownames(ss))


p1 <- SCpubr::do_DotPlot(ss, features=genes, group.by="Date_by_species")
p1
ggsave(p1, filename="Plot_selected_spheroid_marker_expr_in_spheroids.pdf", width=5, height=5)


p2 <- SCpubr::do_DotPlot(hio, features=genes, group.by="Age")
p2

p1 <- SCpubr::do_DimPlot(ss, reduction="mds", group.by="Species", pt.size=3)
p2 <- SCpubr::do_DimPlot(ss, reduction="mds", group.by="Collection_date", pt.size=3)
p1+p2
ggsave(p1+p2, filename="Plot_spheroid_mds_by_species_and_date.pdf", width=10, height=5)
