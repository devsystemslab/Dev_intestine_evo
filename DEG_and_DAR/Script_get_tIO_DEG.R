full_script_path <- "/home/yuq22/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/Script_get_tIO_regional_identity.R"
# return the mapped probability to different regions
## choose resolution level 2: proximal SI, ileum and colon
knn <- readRDS("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/Res_fetal_hg19_knn_for_tIO.rds")
i <- 2
ref.idx <- human_multi@meta.data[,paste0("Intestine_region_L",i)]
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
trans.res <- data.frame(t(apply(nn.idx, 1, function(vec){
  n1 <- sapply(sort(unique(ref.idx)), function(x){
    sum(vec==x)
  })
  p1 <- n1/sum(n1)
  return(c(names(which.max(n1)), p1))
})), stringsAsFactors = F)
colnames(trans.res) <- paste(paste0("Intestine_region_L",i), 
                             c("pred.id", paste0("prob.",sort(unique(ref.idx)))), 
                             sep=":")
for(j in 2:ncol(trans.res)){
  trans.res[,j] <- as.numeric(trans.res[,j])
}
tIO_epi@meta.data[,colnames(trans.res)] <- trans.res
saveRDS(tIO_epi, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")

# identify human-chimp DEGs in each cell type
## taking into account the regional variability
tIO_epi$subtype_group_per_species <- paste(tIO_epi$subtype_group, tIO_epi$Species, sep='@')
saveRDS(tIO_epi, file="/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_RNA_seurat_object.rds")
species_ct_expr <- getAveExpr(seu.obj = tIO_epi, feature.to.calc = "subtype_group_per_species", colname.prefix = NULL)
saveRDS(species_ct_expr, file="Dat_tIO_epi_subtype_group_per_species_expr.rds")
expressed_genes <- rownames(species_ct_expr)[which(rowSums(species_ct_expr)>0)]

setwd("/projects/site/pred/ihb-intestine-evo/iPSC_intestinal_organoid/HIO_and_CIO/tIO_epi/RNA/species_and_region_diff/organoid_regional_identity/species_DE_controlling_region_diff")
for(ct in sort(unique(tIO_epi$subtype_group))){

  cells <- colnames(tIO_epi)[which(tIO_epi$subtype_group==ct)]
  seurat <- subset(tIO_epi, cells=cells)
  idx <- grep(paste0(ct,"@"), colnames(species_ct_expr), fixed = T)
  expressed_genes <- rownames(species_ct_expr)[which(rowSums(species_ct_expr[,idx])>0)]
  
  ## introduce the mapped region probability as numeric covariate for cross species comparison
  species_vec <- as.factor(seurat$Species)
  mapped_region_prob_mat <- seurat@meta.data[,paste(paste0("Intestine_region_L",i), 
                                                    paste0("prob.",c("Colon","Prox_SI")), 
                                                    sep=":")]
  
  registerDoParallel(20)
  species_test_res <- foreach(k=expressed_genes, .multicombine = T, .combine = 'rbind')%dopar%{
    e <- as.numeric(as.vector(seurat@assays$RNA@data[k,]))
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
  saveRDS(species_test_res, file=paste0("Res_",ct,"_species_test_res-2.rds"))
}
