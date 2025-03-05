# for human peaks
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/human/liftOver")
peaks <- readRDS("Res_human_unified_peaks_detection_in_each_sample.rds")
peaks <- rownames(peaks)
df <- do.call('rbind', strsplit(split="\\-|:", x=peaks))
name <- apply(df[,1:3], 1, function(vec){
  paste(vec, collapse="-")
})
input <- data.frame("chr"=df[,1], 
                    "start"=as.integer(df[,2])-1, 
                    "end"=as.integer(df[,3]),
                    "name"=name, 
                    "score"=945, 
                    "strand"=".",
                    stringsAsFactors = F)
write.table(input, file="Table_human_organoid_unified_peaks_hg38.bed", 
            sep="\t", row.names = F, col.names = F, quote=F)


# run liftOver in bash
## hg38 to panTro6
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToPanTro6.over.chain.gz Table_human_organoid_unified_peaks_hg38ToPanTro6.bed unmapped_hg38ToPanTro6.bed -minMatch=0.1
## panTro6 back to hg38
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38ToPanTro6.bed ~/Work/Annotation/liftOver_chain_file/panTro6ToHg38.over.chain.gz Table_human_organoid_unified_peaks_panTro6BackToHg38.bed unmapped_panTro6ToHg38.bed -minMatch=0.1

## hg38 to calJac4
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToCalJac4.over.chain.gz Table_human_organoid_unified_peaks_hg38ToCalJac4.bed unmapped_hg38ToCalJac4.bed -minMatch=0.1
## calJac4 back to hg38
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38ToCalJac4.bed ~/Work/Annotation/liftOver_chain_file/calJac4ToHg38.over.chain.gz Table_human_organoid_unified_peaks_calJac4BackToHg38.bed unmapped_calJac4ToHg38.bed -minMatch=0.1

## hg38 to mm10
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToMm10.over.chain.gz Table_human_organoid_unified_peaks_hg38ToMm10.bed unmapped_hg38ToMm10.bed -minMatch=0.1
## mm10 back to hg38
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38ToMm10.bed ~/Work/Annotation/liftOver_chain_file/mm10ToHg38.over.chain.gz Table_human_organoid_unified_peaks_mm10BackToHg38.bed unmapped_mm10ToHg38.bed -minMatch=0.1


# for chimp peaks
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/chimp/liftOver")
peak_mat <- readRDS("../Res_chimp_peaks_called_in_each_sample.rds")
end <- as.integer(peak_mat@ranges@start)+as.integer(peak_mat@ranges@width)-1
peaks <- paste(peak_mat@seqnames, peak_mat@ranges@start, end, sep="-")
input <- data.frame("chr"=peak_mat@seqnames, 
                    "start"=as.integer(as.integer(peak_mat@ranges@start)-1), 
                    "end"=as.integer(end),
                    "name"=peaks, 
                    "score"=945, 
                    "strand"=".",
                    stringsAsFactors = F)
write.table(input, file="Table_chimp_organoid_unified_peaks_panTro6.bed", 
            sep="\t", row.names = F, col.names = F, quote=F)


# run liftOver in bash
## panTro6 to hg38
#~/Work/Tools/liftOver  Table_chimp_organoid_unified_peaks_panTro6.bed ~/Work/Annotation/liftOver_chain_file/panTro6ToHg38.over.chain.gz Table_chimp_organoid_unified_peaks_panTro6ToHg38.bed unmapped_panTro6ToHg38.bed -minMatch=0.1
## hg38 back to panTro6
#~/Work/Tools/liftOver  Table_chimp_organoid_unified_peaks_panTro6ToHg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToPanTro6.over.chain.gz Table_chimp_organoid_unified_peaks_hg38BackToPanTro6.bed unmapped_hg38ToPanTro6.bed -minMatch=0.1

## panTro6 to mm10
#~/Work/Tools/liftOver  Table_chimp_organoid_unified_peaks_panTro6.bed ~/Work/Annotation/liftOver_chain_file/panTro6ToMm10.over.chain.gz Table_chimp_organoid_unified_peaks_panTro6ToMm10.bed unmapped_panTro6ToMm10.bed -minMatch=0.1
## mm10 back to panTro6
#~/Work/Tools/liftOver  Table_chimp_organoid_unified_peaks_panTro6ToMm10.bed ~/Work/Annotation/liftOver_chain_file/mm10ToPanTro6.over.chain.gz Table_chimp_organoid_unified_peaks_mm10BackToPanTro6.bed unmapped_mm10ToPanTro6.bed -minMatch=0.1

# for marmost peaks
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/marmoset_enteroid_calJac4_ATAC")
peak_mat <- read.table("/nas/groups/treutlein/DATA/sequencing/20210330_UMUT_intestine_organoid_data_from_IOB/marmoset-1_ATAC/processed_calJac4/322_marmoset-1_ATAC/outs/peaks.bed",
                       sep="\t", stringsAsFactors = F)
peaks <- paste(peak_mat[,1], peak_mat[,2]+1, peak_mat[,3], sep="-")
input <- data.frame("chr"=peak_mat[,1], 
                    "start"=as.integer(peak_mat[,2]), 
                    "end"=as.integer(peak_mat[,3]),
                    "name"=peaks, 
                    "score"=945, 
                    "strand"=".",
                    stringsAsFactors = F)
write.table(input, file="Table_marmoset_organoid_unified_peaks_calJac4.bed", 
            sep="\t", row.names = F, col.names = F, quote=F)

## calJac4 to hg38
#~/Work/Tools/liftOver  Table_marmoset_organoid_unified_peaks_calJac4.bed ~/Work/Annotation/liftOver_chain_file/calJac4ToHg38.over.chain.gz Table_marmoset_organoid_unified_peaks_calJac4ToHg38.bed unmapped_calJac4ToHg38.bed -minMatch=0.1
## hg38 back to calJac4
#~/Work/Tools/liftOver  Table_marmoset_organoid_unified_peaks_calJac4ToHg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToCalJac4.over.chain.gz Table_marmoset_organoid_unified_peaks_hg38BackToCalJac4.bed unmapped_hg38ToCalJac4.bed -minMatch=0.1

## calJac4 to mm10
#~/Work/Tools/liftOver  Table_marmoset_organoid_unified_peaks_calJac4.bed ~/Work/Annotation/liftOver_chain_file/calJac4ToMm10.over.chain.gz Table_marmoset_organoid_unified_peaks_calJac4ToMm10.bed unmapped_calJac4ToMm10.bed -minMatch=0.1
## mm10 back to calJac4
#~/Work/Tools/liftOver  Table_marmoset_organoid_unified_peaks_calJac4ToMm10.bed ~/Work/Annotation/liftOver_chain_file/mm10ToCalJac4.over.chain.gz Table_marmoset_organoid_unified_peaks_mm10BackToCalJac4.bed unmapped_mm10ToCalJac4.bed -minMatch=0.1

# for mouse peaks
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/mouse_enteroid_mm10_ATAC")
peak_mat <- read.table("/nas/groups/treutlein/DATA/sequencing/20190916_IOB_UMUT_Mouse_adult_intestine_organoid_scATAC_seq_NB551561_0007_AH5WGWBGXB/processed/aMIO_6d/outs/peaks.bed",
                       sep="\t", stringsAsFactors = F)
peaks <- paste(peak_mat[,1], peak_mat[,2]+1, peak_mat[,3], sep="-")
input <- data.frame("chr"=peak_mat[,1], 
                    "start"=as.integer(peak_mat[,2]), 
                    "end"=as.integer(peak_mat[,3]),
                    "name"=peaks, 
                    "score"=945, 
                    "strand"=".",
                    stringsAsFactors = F)
write.table(input, file="Table_mouse_organoid_unified_peaks_mm10.bed", 
            sep="\t", row.names = F, col.names = F, quote=F)

## mm10 to hg38
#~/Work/Tools/liftOver  Table_mouse_organoid_unified_peaks_mm10.bed ~/Work/Annotation/liftOver_chain_file/mm10ToHg38.over.chain.gz Table_mouse_organoid_unified_peaks_mm10ToHg38.bed unmapped_mm10ToHg38.bed -minMatch=0.1
## hg38 back to mm10
#~/Work/Tools/liftOver  Table_mouse_organoid_unified_peaks_mm10ToHg38.bed ~/Work/Annotation/liftOver_chain_file/hg38Tomm10.over.chain.gz Table_mouse_organoid_unified_peaks_hg38BackTomm10.bed unmapped_hg38Tomm10.bed -minMatch=0.1

