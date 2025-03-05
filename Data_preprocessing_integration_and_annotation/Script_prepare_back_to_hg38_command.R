setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/Human_without_adult_enteroid/liftOver")
chain_file <- readLines("hg38_to_others.txt")
chain_folder <- "~/Work/Annotation/liftOver_chain_file"
out <- sapply(chain_file, function(x){
  sp <- strsplit(x, split=".", fixed = T)[[1]][1]
  sp2 <- strsplit(sp, split="To")[[1]][2]
  aa <- strsplit(sp2, split="", fixed=T)[[1]]
  sp2 <- paste0(tolower(aa[1]), paste(aa[-1], collapse = ""))
  id <- paste0(sp2,"ToHg38")
  back_chain_file <- file.path(chain_folder, paste0(id,".over.chain.gz"))
  input <- paste0("Table_human_organoid_unified_peaks_",sp,".bed")
  output1 <- paste0("Table_human_organoid_unified_peaks_",id,".bed")
  output2 <- paste0("unmapped_",id,".bed")
  command <- paste("~/Work/Tools/liftOver", input, back_chain_file, output1, output2, "-minMatch=0.1")
  return(command)
})
writeLines(out, con="run_liftOver_back_to_hg38.sh")