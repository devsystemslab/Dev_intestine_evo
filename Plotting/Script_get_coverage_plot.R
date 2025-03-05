seu_obj_list <- SplitObject(fetal_atac, split.by = "Sample.name")
for(x in names(seu_obj_list)){
  wd <- paste("/projects/site/pred/ihb-intestine-evo/fetal_human_duo_crypt/troubleShooting",x,sep="/")
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)
  seu_obj <- seu_obj_list[[x]]
  for(i in seq(nrow(mat))){
    p <- mat$name[i]
    g <- mat$gene_symbol[i]
    plot_name <- paste0("Plot_coveragePlot_hg38_",p,"_",g,"_10k_smoothed_in_",x,".pdf")
    p2 <- CoveragePlot(
      object = seu_obj,
      assay="ATAC",
      group.by = "Unified_cell_type",
      region = p,
      extend.upstream = 10000,
      extend.downstream = 10000,
      #features = g,
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = TRUE,
      window = 500,
      region.highlight = StringToGRanges(p)
    )
    p2 <- p2 & scale_fill_manual(values = epi.ct.cols) 
    pdf(plot_name, height=7, width=12)
    print(p2)
    dev.off()
    
  }
  
}
