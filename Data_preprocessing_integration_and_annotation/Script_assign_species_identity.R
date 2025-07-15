wd <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/demultiplex"
wd <- "/home/yuq22/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/demultiplex"
setwd(wd)

folders <- c("1020-UK-1-ATAC", "1020-UK-1-RNA", "993-UK-1-RNA", "993-UK-1-ATAC", "UK-R1-RNA")
folders <- c("UK-A1-ATAC")

for(x in folders){
  print(x)
  setwd(x)
  data <- read.table(paste0(x,"_demultiplexed.tsv"), sep="\t", stringsAsFactors = F, head=T, row.names = 1)
  head(data)
  sum.vec <- rowSums(data)
  human.prop <- data[,1]/sum.vec
  chimp.prop <- data[,2]/sum.vec
  pdf("Plot_barplot_demultiplex.pdf", height=5, width=8)
  par(mfrow=c(1,2))
  hist(chimp.prop)
  hist(human.prop)
  dev.off()
  human.cell.idx <- names(human.prop)[which(human.prop>0.8)]
  chimp.cell.idx <- names(chimp.prop)[which(chimp.prop>0.8)]
  n1 <- length(chimp.cell.idx)
  n2 <- length(human.cell.idx)
  n3 <- length(sum.vec)-n1-n2
  print(paste0("Chimp cell number: ", n1,"; Human cell number: ",n2,"; Unassigned cell number: ",n3))
  sum(human.cell.idx%in%chimp.cell.idx)
  write.table(human.cell.idx, file=paste0("List_",x,"_human_cell_barcodes.txt"), quote=F, row.names = F, col.names = F)
  write.table(chimp.cell.idx, file=paste0("List_",x,"_chimp_cell_barcodes.txt"), quote=F, row.names = F, col.names = F)
  setwd(wd)
}

