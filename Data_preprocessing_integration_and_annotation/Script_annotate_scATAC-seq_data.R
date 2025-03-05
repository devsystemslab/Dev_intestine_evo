# this script is used to annotate scATAC-seq data

## identify cell type markers in the scRNA-seq data
cell_type_size <- table(HIO_rna$Cell_type_group)
selected_cell_type <- names(cell_type_size)[which(cell_type_size>5)]
de_res <- presto::wilcoxauc(HIO_rna, group_by="Cell_type_group")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res[which(de_res$auc>0.6 & de_res$padj<0.05 & de_res$logFC>0.1 & de_res$pct_diff > 20),]
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=auc)
deg_res <- list("sig_res"=sig_res,
                "top_res"=top_res)
saveRDS(deg_res, file="Res_HIO_rna_cell_type_group_marker_genes.rds")

## load scATAC-seq data
HIO_atac <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/HIO_and_tHIO/scATAC-seq/inVitro_and_inVivo/integration/Res_Signac_integrated_HIO_atac_ATAC_preprocessed.rds")

## load scRNA- and scATAC-seq cell matching result
merged_meta_cell_info <- readRDS("../Res_HIO_and_tHIO_merged_meta_cell_info.rds")

paired_cells <- merged_meta_cell_info[,c("RNA_cell", "ATAC_cell")]
c1 <- intersect(merged_meta_cell_info$RNA_cell, colnames(HIO_rna)) # get the cells that are present in the final annotated scRNA-seq data and scRNA/scATAC matching results
paired_cells$Mapped_RNA_cell_cluster <- NA 
paired_cells$Mapped_RNA_cell_cluster[match(c1, paired_cells$RNA_cell)] <- as.character(HIO_rna@meta.data[c1, "Cell_type_group"]) # transfer the RNA final annotation to the matched cells
mat <- unique(paired_cells[,c("ATAC_cell", "Mapped_RNA_cell_cluster")])
# only transfer the RNA identity to an ATAC cell if the ATAC cell is uniquely matched to one RNA identity
# if an ATAC cell is matched to more than one RNA cells with different cell type identity, than do not transfer the RNA identity directly 
freq <- table(mat$ATAC_cell)
mixed_cells <- names(freq)[which(freq>1)] 
unique_cells <- names(freq)[which(freq==1)]
length(unique_cells)
HIO_atac$Mapped_RNA_cell_type <- NA
HIO_atac@meta.data[unique_cells, "Mapped_RNA_cell_type"] <- paired_cells[match(unique_cells, paired_cells$ATAC_cell), "Mapped_RNA_cell_cluster"]
saveRDS(HIO_atac, file="Dat_HIO_atac_before_annotation_refinement.rds")

DimPlot(HIO_atac, reduction = "umap", group.by = "Mapped_RNA_cell_type", label=T)+NoLegend()
unique_atac_mat <- paired_cells[match(unique_cells, paired_cells$ATAC_cell),] 
# annotate the scATAC-seq data on each cluster of each sample
## aggregate scaled top marker gene activity scores across ATAC clusters, check whether certain ATAC cluster show significantly higher expression of RNA cluster markers
## check whether overlap with enriched mapped RNA identity - if yes, assign to that identity
DefaultAssay(HIO_atac) <- "peaks_stage"
seu_obj_list <- SplitObject(HIO_atac, split.by = "Sample.name")
prop_cutoff <- 0.7
p_cutoff <- 0.1
for(x in names(seu_obj_list)){
  print(x)
  # de novo cluster on each scATAC-seq sample
  seu_obj <- seu_obj_list[[x]]
  det_rates <- rowMeans(seu_obj@assays$peaks_stage@counts > 0) 
  top <- names(which(det_rates>0.05 & det_rates<0.9))
  seu_obj <- FindTopFeatures(seu_obj, min.cutoff = ncol(seu_obj) * 0.05)
  VariableFeatures(seu_obj) <- top
  seu_obj <- RunSVD(seu_obj)
  seu_obj <- RunUMAP(seu_obj, reduction="lsi", dims = 2:20)
  seu_obj <- FindNeighbors(seu_obj, reduction="lsi", dims=2:20) %>% FindClusters(resolution=2)
  
  # aggregate scaled top marker gene activity scores across ATAC clusters, check whether certain ATAC cluster show significantly higher expression of RNA cluster markers
  ave_gene_activity <- getAveExpr(seu.obj=seu_obj, feature.to.calc="peaks_stage_snn_res.2",  genes = intersect(unique(top_res$feature), rownames(seu_obj@assays$RNA@data)), size.cutoff = 1)
  score_mat <- sapply(sort(unique(top_res$group)), function(ct){
    genes <- intersect(top_res$feature[which(top_res$group==ct)], rownames(seu_obj@assays$RNA@data))
    colSums(ave_gene_activity[genes,])
  })
  z_score <- scale(t(score_mat))
  cluster_enrichment <- z_score > qnorm(0.9)
  
  
  cl_to_test <- sort(unique(seu_obj$peaks_stage_snn_res.2))
  
  # check whether overlap with enriched mapped RNA identity - if yes, assign to that identity
  ## calculate the proportion of mapped RNA identity in each ATAC cluster
  prop <- sapply(cl_to_test, function(i){
    sapply(setdiff(sort(unique(seu_obj$Mapped_RNA_cell_type)), NA), function(x){
      a <- sum(seu_obj$peaks_stage_snn_res.2==i & seu_obj$Mapped_RNA_cell_type==x, na.rm=T)
      b <- sum(seu_obj$peaks_stage_snn_res.2==i & !is.na(seu_obj$Mapped_RNA_cell_type), na.rm=T)
      a/b
    })
  })
  colnames(prop) <- paste("Cluster", cl_to_test, sep="_")
  
  ## calculate the p-value of over-representation of RNA identity in each ATAC cluster
  pval <- sapply(cl_to_test, function(i){
    sapply(sort(unique(seu_obj$Mapped_RNA_cell_type)), function(x){
      a <- sum(seu_obj$peaks_stage_snn_res.2==i & seu_obj$Mapped_RNA_cell_type==x, na.rm=T)
      b <- sum(seu_obj$peaks_stage_snn_res.2==i, na.rm=T)
      c <- sum(seu_obj$Mapped_RNA_cell_type==x, na.rm=T)
      d <- ncol(seu_obj)
      fisher.test(matrix(c(a,b,c,d), c(2,2)), alternative = "g")$p.value
    })
  })
  colnames(pval) <- paste("Cluster", cl_to_test, sep="_")
  
  enriched_cell_type_num <- colSums(pval<p_cutoff)+colSums(prop>prop_cutoff)
  
  # if an ATAC cluster is not dominated by any mapped RNA identity and no significant enrichment
  # then annotate the ATAC cluster based on gene activity similarity to RNA reference
  for(cl_idx in which(enriched_cell_type_num==0)){
    cl <- cl_to_test[cl_idx]
    cluster_cells <- colnames(seu_obj)[which(seu_obj$peaks_stage_snn_res.2==cl)]
    cluster_activity_enriched_ct <- rownames(cluster_enrichment)[which(cluster_enrichment[,paste("Cluster",cl,sep="_")])]
    mapped_cells_identity <- unique(seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"])
    enriched_ct_overlap <- intersect(mapped_cells_identity, cluster_activity_enriched_ct) # check whether there is overlap between MCMF-based mapped RNA identity on certain cells and cluster marker gene body activity score on the cluster level  
    
    if(length(enriched_ct_overlap)==1){ # assign the overlapped identity to the whole cluster if there is one overlap 
      seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- enriched_ct_overlap
      seu_obj_list[[x]] <- seu_obj
      HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- enriched_ct_overlap
      
    }else if(length(enriched_ct_overlap)>1){ # assign the overlapped identities with the maximal score to each individual cell if there are more than one overlapping identities 
      per_cl_score_mat <- sapply(enriched_ct_overlap, function(ct){
        genes <- intersect(top_res$feature[which(top_res$group==ct)], rownames(seu_obj@assays$RNA@data))
        colSums(seu_obj@assays$RNA@data[genes, cluster_cells])
      })
      pred_id <- colnames(per_cl_score_mat)[apply(per_cl_score_mat, 1, which.max)]
      seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- pred_id
      HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- pred_id
      seu_obj_list[[x]] <- seu_obj
    }else if(length(enriched_ct_overlap)==0){ # if there is no overlap between MCMF-based mapped RNA identity on certain cells and cluster marker gene body activity score on the cluster level
      if(length(cluster_activity_enriched_ct)==1){ # assign the RNA cluster marker gene body activity score based identity to the whole cluster if only one cluster shows significantly higher score
        seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- cluster_activity_enriched_ct
        seu_obj_list[[x]] <- seu_obj
        HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- cluster_activity_enriched_ct
        
      }else if(length(cluster_activity_enriched_ct)>1){ # assign the mapped RNA cluster identities with the maximal score to each individual cell if there are more than one clusters showing significantly higher score
        per_cl_score_mat <- sapply(cluster_activity_enriched_ct, function(ct){
          genes <- intersect(top_res$feature[which(top_res$group==ct)], rownames(seu_obj@assays$RNA@data))
          colSums(seu_obj@assays$RNA@data[genes, cluster_cells])
        })
        pred_id <- colnames(per_cl_score_mat)[apply(per_cl_score_mat, 1, which.max)]
        seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- pred_id
        HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- pred_id
        seu_obj_list[[x]] <- seu_obj
      }
    }
  }
  
  # if an ATAC cluster is dominated by any mapped RNA identity or shows significant enrichment, diffuse the over-represented RNA identity from the matched cells to unmatched cells
  for(cl_idx in which(enriched_cell_type_num>0)){
    
    cl <- cl_to_test[cl_idx]
    cluster_cells <- colnames(seu_obj)[which(seu_obj$peaks_stage_snn_res.2==cl)]
    enriched_cell_type <- union(rownames(pval)[which(pval[,cl_idx]<p_cutoff)], 
                                rownames(prop)[which(prop[,paste0("Cluster_",cl)]>prop_cutoff)])
    if(length(enriched_cell_type)==1){ # assign the overlapped identity to the whole cluster if there is only one over-represented RNA identity 
      seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- enriched_cell_type
      seu_obj_list[[x]] <- seu_obj
      HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- enriched_cell_type
    }else{ # diffuse the overlapped identity from the matched cells to the unmatched cells based on neighboring cells on lsi space (exclude the first dimension) if there is more than one over-represented RNA identities 
      dr <- Embeddings(seu_obj, reduction = "lsi")[cluster_cells,-1]
      que_cells <- intersect(colnames(seu_obj)[which(is.na(seu_obj$Mapped_RNA_cell_type) | !seu_obj$Mapped_RNA_cell_type%in%enriched_cell_type)], rownames(dr))
      ref_cells <- intersect(colnames(seu_obj)[which(seu_obj$Mapped_RNA_cell_type %in% enriched_cell_type)], rownames(dr)) # take the ATAC cells with mapped RNA identity that is overrepresented as reference cells
      seu_obj@meta.data[que_cells, "Mapped_RNA_cell_type"] <- NA
      if(length(ref_cells)<5){ # if the reference ATAC cell number is too small, check the overlap with RNA cluster marker gene activity score enrichment, follow the rules when there is no over-represented identity
        cluster_activity_enriched_ct <- rownames(cluster_enrichment)[which(cluster_enrichment[,paste("Cluster",cl,sep="_")])]
        enriched_ct_overlap <- intersect(cluster_activity_enriched_ct, enriched_cell_type)
        if(length(enriched_ct_overlap)==1){
          seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- enriched_ct_overlap
          seu_obj_list[[x]] <- seu_obj
          HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- enriched_ct_overlap
          
        }else if(length(enriched_ct_overlap)>1){
          per_cl_score_mat <- sapply(cluster_activity_enriched_ct, function(ct){
            genes <- intersect(top_res$feature[which(top_res$group==ct)], rownames(seu_obj@assays$RNA@data))
            colSums(seu_obj@assays$RNA@data[genes, cluster_cells])
          })
          pred_id <- colnames(per_cl_score_mat)[apply(per_cl_score_mat, 1, which.max)]
          seu_obj@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- pred_id
          HIO_atac@meta.data[cluster_cells, "Mapped_RNA_cell_type"] <- pred_id
          seu_obj_list[[x]] <- seu_obj
        }else if(length(enriched_ct_overlap)==0){
          seu_obj@meta.data[que_cells, "Mapped_RNA_cell_type"] <- NA
          HIO_atac@meta.data[que_cells, "Mapped_RNA_cell_type"] <- NA
          seu_obj_list[[x]] <- seu_obj
        }
        
      }else{ # if the reference ATAC cell number is not too small, then diffuse the over-represented identity from the matched cells to the unmatched cells based on maximal votes of neighboring cells 
        que_mat <- dr[que_cells,]
        ref_mat <- dr[ref_cells,]
        knn <- RANN::nn2(ref_mat, que_mat, k = 5)$nn.idx
        ref_idx <- seu_obj@meta.data[ref_cells, "Mapped_RNA_cell_type"]
        nn_idx <- matrix(ref_idx[as.vector(knn)], nrow=nrow(knn))
        pred_id <- apply(nn_idx, 1, function(vec){
          freq <- table(vec)
          names(which.max(freq))
        })
        seu_obj@meta.data[que_cells, "Mapped_RNA_cell_type"] <- pred_id
        HIO_atac@meta.data[que_cells, "Mapped_RNA_cell_type"] <- pred_id
        seu_obj_list[[x]] <- seu_obj
      }
      
    }
  }
  
}
saveRDS(HIO_atac, file="Res_HIO_atac_with_per_sample_cluster_based_annotation_refined_v3.rds")
saveRDS(seu_obj_list, file="Res_HIO_atac_seu_obj_list_with_per_sample_cluster_based_annotation_refined_v3.rds")
