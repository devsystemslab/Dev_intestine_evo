setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs")

get_peak_orthologous_region <- function(p1_sp2_to_sp1_bed_file = NULL,
                                        p1_sp1_to_sp2_bed_file = NULL,
                                        peaks_to_subset=main_chr_peaks,
                                        name_idx=4,
                                        name_split="-",
                                        remove_peak_ID=TRUE,
                                        peak_ID_sep="_",
                                        chr_idx_1 = 1,
                                        chr_idx_2 = 7,
                                        start_idx_1 = 2,
                                        start_idx_2 = 8,
                                        end_idx_1 = 3,
                                        end_idx_2 = 9,
                                        length_idx = 10,
                                        back_map_prop_cutoff = 0.95
){
  
  ## NOTE
  # All input bed files should use 0-based start position, 1-based end position
  # But in the name field, i.e. the fourth column, both the start and end position are 1-based
  # The peak_ortholog_output will be 1-based, to fit the peak format of Signac object  
  
  ## step 1. Get p1 mapped back to original region by reciprocal liftOver 
  input_mat <- read.table(p1_sp2_to_sp1_bed_file, sep="\t", stringsAsFactors = F)
  colnames(input_mat) <- c("sp1_chr", "sp1_start", "sp1_end", "p1_names", "score", "strand")
  
  ### 1.1 split name to region
  if(!is.null(peaks_to_subset)){
    input_mat <- input_mat[input_mat[,name_idx]%in%peaks_to_subset,]
  }
  vec_to_split = input_mat[,name_idx]
  if(remove_peak_ID){
    vec_to_split <- sapply(vec_to_split, function(x){
      strsplit(x=x, split = peak_ID_sep)[[1]][2]
    }) 
  }
  dat <- do.call('rbind', strsplit(x=vec_to_split, split = name_split))
  input_mat$original_chr <- dat[,1]
  input_mat$original_start <- as.integer(dat[,2])
  input_mat$original_end <- as.integer(dat[,3])
  input_mat$original_region_length <- input_mat$original_end - input_mat$original_start+1
  
  ### 1.2 get overlap region length and proportion
  overlapped_region_length <- apply(input_mat, 1, function(vec){
    if(vec[chr_idx_1]==vec[chr_idx_2]){
      start <- max(as.integer(vec[start_idx_1])+1, as.integer(vec[start_idx_2]))
      end <- min(as.integer(vec[end_idx_1]), as.integer(vec[end_idx_2]))
      if(start<end){
        length <- end-start+1
        return(length)
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  })
  overlapped_region_prop <- overlapped_region_length/input_mat[,length_idx]
  
  ### 1.3 only keep peaks with overlapped proportion that pass the cutoff
  p1 <- input_mat[which(overlapped_region_prop > back_map_prop_cutoff), name_idx]
  
  ## step 2. Get orthologous region of p1
  p1_sp1_to_sp2_mat <- read.table(p1_sp1_to_sp2_bed_file, sep="\t", stringsAsFactors = F)
  p1_sp1_to_sp2_mat <- p1_sp1_to_sp2_mat[which(p1_sp1_to_sp2_mat[,4] %in% p1),]
  colnames(p1_sp1_to_sp2_mat) <- c("sp2_chr", "sp2_start", "sp2_end", "p1_names", "score", "strand")
  p1_sp1_to_sp2_mat$sp2_start <- p1_sp1_to_sp2_mat$sp2_start+1
  
  ## output peak orthologous regions, by default we used tsv format, which is 1-based for genomic coordinates
  p1_orth_region <- paste(p1_sp1_to_sp2_mat$sp2_chr, p1_sp1_to_sp2_mat$sp2_start, p1_sp1_to_sp2_mat$sp2_end, sep ="-")
  dat <- data.frame("p1_sp2_orth_region"=p1_orth_region, "p1_original_peak_name"=p1_sp1_to_sp2_mat$p1_names, stringsAsFactors = F)
  return(dat)
}

setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs")
# get orthologous regions of chimp peaks in hg38
input_1 <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/chimp/liftOver/Table_chimp_organoid_unified_peaks_hg38BackToPanTro6.bed"
input_2 <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/chimp/liftOver/Table_chimp_organoid_unified_peaks_panTro6ToHg38.bed"
chimp_peaks_in_human <- get_peak_orthologous_region(p1_sp2_to_sp1_bed_file = input_1,
                                                    p1_sp1_to_sp2_bed_file = input_2)
names(chimp_peaks_in_human) <- c("chimp_peaks_in_human_coor", "chimp_peaks_in_chimp_coor")
saveRDS(chimp_peaks_in_human, file="Res_chimp_unified_peaks_in_hg38_coor.rds")

# get orthologous regions of marmoset peaks in hg38
dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/marmoset_enteroid_calJac4_ATAC/"
input_1 <- paste0(dir, "Table_marmoset_organoid_unified_peaks_hg38BackToCalJac4.bed")
input_2 <- paste0(dir, "Table_marmoset_organoid_unified_peaks_calJac4ToHg38.bed")
marmoset_peaks_in_human <- get_peak_orthologous_region(p1_sp2_to_sp1_bed_file = input_1,
                                                       p1_sp1_to_sp2_bed_file = input_2)
names(marmoset_peaks_in_human) <- c("marmoset_peaks_in_human_coor", "marmoset_peaks_in_marmoset_coor")
saveRDS(marmoset_peaks_in_human, file="Res_marmoset_unified_peaks_in_hg38_coor.rds")

# get orthologous regions of mouse peaks in hg38
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs")
dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/mouse_enteroid_mm10_ATAC/"
input_1 <- paste0(dir, "Table_mouse_organoid_unified_peaks_hg38BackToMm10.bed")
input_2 <- paste0(dir, "Table_mouse_organoid_unified_peaks_mm10ToHg38.bed")
mouse_peaks_in_human <- get_peak_orthologous_region(p1_sp2_to_sp1_bed_file = input_1,
                                                    p1_sp1_to_sp2_bed_file = input_2)
names(mouse_peaks_in_human) <- c("mouse_peaks_in_human_coor", "mouse_peaks_in_mouse_coor")
saveRDS(mouse_peaks_in_human, file="Res_mouse_unified_peaks_in_hg38_coor.rds")

# merge peaks detected in human, chimp, marmoset and mouse
# all are in hg38 coordinate (both start and end coordinates are 1-based)
mouse_peaks <- mouse_peaks_in_human$mouse_peaks_in_human_coor
marmoset_peaks <- marmoset_peaks_in_human$marmoset_peaks_in_human_coor
chimp_peaks <- chimp_peaks_in_human$chimp_peaks_in_human_coor
human_peaks <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/human/Res_human_unified_peaks_detection_in_each_sample.rds")
human_peaks <- rownames(human_peaks)

all_peaks <- unique(c(human_peaks, chimp_peaks, marmoset_peaks, mouse_peaks))
dat <- do.call('rbind', strsplit(x=all_peaks, split = "-"))
dat <- as.data.frame(dat)
colnames(dat) <- c("chr", "start", "end")
# convert to genomic ranges
gr <- makeGRangesFromDataFrame(dat)

# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = gr)

# Filter out bad peaks based on length
peakwidths <- width(combined_peaks)
combined_peaks <- combined_peaks[peakwidths  < 10000 & peakwidths > 20]
saveRDS(combined_peaks, file="Res_merged_peaks_in_human.rds")

peaks <- GRangesToString(combined_peaks)
df <- do.call('rbind', strsplit(split="-", x=peaks))
name <- paste(paste0("peak", seq(length(peaks))), peaks, sep="_")
input <- data.frame("chr"=df[,1], 
                    "start"=as.integer(as.integer(df[,2])-1), 
                    "end"=as.integer(df[,3]),
                    "name"=name, 
                    "score"=945, 
                    "strand"=".",
                    stringsAsFactors = F)
write.table(input, file="Table_human-chimp-marmoset-mouse_organoid_merged_peaks_hg38.bed", 
            sep="\t", row.names = F, col.names = F, quote=F)


setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/combined_liftOver")
# run liftOver in bash
## hg38 to panTro6
#~/Work/Tools/liftOver  Table_human-chimp-marmoset-mouse_organoid_merged_peaks_hg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToPanTro6.over.chain.gz Table_human_organoid_unified_peaks_hg38ToPanTro6.bed unmapped_hg38ToPanTro6.bed -minMatch=0.1
## panTro6 back to hg38
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38ToPanTro6.bed ~/Work/Annotation/liftOver_chain_file/panTro6ToHg38.over.chain.gz Table_human_organoid_unified_peaks_panTro6BackToHg38.bed unmapped_panTro6ToHg38.bed -minMatch=0.1

## hg38 to calJac4
#~/Work/Tools/liftOver  Table_human-chimp-marmoset-mouse_organoid_merged_peaks_hg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToCalJac4.over.chain.gz Table_human_organoid_unified_peaks_hg38ToCalJac4.bed unmapped_hg38ToCalJac4.bed -minMatch=0.1
## calJac4 back to hg38
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38ToCalJac4.bed ~/Work/Annotation/liftOver_chain_file/calJac4ToHg38.over.chain.gz Table_human_organoid_unified_peaks_calJac4BackToHg38.bed unmapped_calJac4ToHg38.bed -minMatch=0.1

## hg38 to mm10
#~/Work/Tools/liftOver  Table_human-chimp-marmoset-mouse_organoid_merged_peaks_hg38.bed ~/Work/Annotation/liftOver_chain_file/hg38ToMm10.over.chain.gz Table_human_organoid_unified_peaks_hg38ToMm10.bed unmapped_hg38ToMm10.bed -minMatch=0.1
## mm10 back to hg38
#~/Work/Tools/liftOver  Table_human_organoid_unified_peaks_hg38ToMm10.bed ~/Work/Annotation/liftOver_chain_file/mm10ToHg38.over.chain.gz Table_human_organoid_unified_peaks_mm10BackToHg38.bed unmapped_mm10ToHg38.bed -minMatch=0.1

# get orthologous regions of combined peaks in other species
setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/combined_ortholog")
input <- read.table("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/Table_human-chimp-marmoset-mouse_organoid_merged_peaks_hg38.bed", 
                    sep="\t", stringsAsFactors = F)
# only consider peaks on main chromosomes
main_chr <- paste0("chr", c(seq(22), "X", "Y", "2A", "2B"))
input <- input[which(input[,1] %in% main_chr),]
main_chr_peaks <- input[,4]
saveRDS(main_chr_peaks, file="Res_combined_peaks_in_hg38_coor.rds")

## in chimpanzee (panTro6)
dir <- "/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/combined_liftOver/"
input_1 <- paste0(dir, "Table_human_organoid_unified_peaks_panTro6BackToHg38.bed")
input_2 <- paste0(dir, "Table_human_organoid_unified_peaks_hg38ToPanTro6.bed")
combined_peaks_in_chimp <- get_peak_orthologous_region(p1_sp2_to_sp1_bed_file = input_1,
                                                       p1_sp1_to_sp2_bed_file = input_2,
                                                       peaks_to_subset=main_chr_peaks)
names(combined_peaks_in_chimp) <- c("combined_peaks_in_chimp_coor", "combined_peaks_in_human_coor")
peak_ID <- sapply(combined_peaks_in_chimp[,2], function(x){strsplit(x=x, split="_")[[1]][1]})
combined_peaks_in_chimp[,1] <- paste(peak_ID, combined_peaks_in_chimp[,1], sep="_")
# only keep peaks located on main chromosomes
chimp_chr <- sapply(combined_peaks_in_chimp[,1], function(x){
  vec <- strsplit(x=x, split="-")[[1]][1]
  chr <- paste(strsplit(x=vec, split="_")[[1]][-1], collapse = "_")
})
idx <- which(chimp_chr %in% main_chr)
combined_peaks_in_chimp_2 <- combined_peaks_in_chimp[idx,]
saveRDS(combined_peaks_in_chimp_2, file="Res_combined_peaks_in_panTro6_hg38_coor.rds")

## in marmoset (calJac4)
input_1 <- paste0(dir, "Table_human_organoid_unified_peaks_calJac4BackToHg38.bed")
input_2 <- paste0(dir, "Table_human_organoid_unified_peaks_hg38ToCalJac4.bed")
combined_peaks_in_marmoset <- get_peak_orthologous_region(p1_sp2_to_sp1_bed_file = input_1,
                                                          p1_sp1_to_sp2_bed_file = input_2,
                                                          peaks_to_subset=main_chr_peaks)
names(combined_peaks_in_marmoset) <- c("combined_peaks_in_marmoset_coor", "combined_peaks_in_human_coor")
peak_ID <- sapply(combined_peaks_in_marmoset[,2], function(x){strsplit(x=x, split="_")[[1]][1]})
combined_peaks_in_marmoset[,1] <- paste(peak_ID, combined_peaks_in_marmoset[,1], sep="_")
# only keep peaks located on main chromosomes
marmoset_chr <- sapply(combined_peaks_in_marmoset[,1], function(x){
  vec <- strsplit(x=x, split="-")[[1]][1]
  chr <- paste(strsplit(x=vec, split="_")[[1]][-1], collapse = "_")
})
idx <- which(marmoset_chr %in% main_chr)
combined_peaks_in_marmoset_2 <- combined_peaks_in_marmoset[idx,]
saveRDS(combined_peaks_in_marmoset_2, file="Res_combined_peaks_in_calJac4_hg38_coor.rds")

## in mouse (mm10)
input_1 <- paste0(dir, "Table_human_organoid_unified_peaks_mm10BackToHg38.bed")
input_2 <- paste0(dir, "Table_human_organoid_unified_peaks_hg38ToMm10.bed")
combined_peaks_in_mouse <- get_peak_orthologous_region(p1_sp2_to_sp1_bed_file = input_1,
                                                       p1_sp1_to_sp2_bed_file = input_2,
                                                       peaks_to_subset=main_chr_peaks)
names(combined_peaks_in_mouse) <- c("combined_peaks_in_mouse_coor", "combined_peaks_in_human_coor")
peak_ID <- sapply(combined_peaks_in_mouse[,2], function(x){strsplit(x=x, split="_")[[1]][1]})
combined_peaks_in_mouse[,1] <- paste(peak_ID, combined_peaks_in_mouse[,1], sep="_")
# only keep peaks located on main chromosomes
mouse_chr <- sapply(combined_peaks_in_mouse[,1], function(x){
  vec <- strsplit(x=x, split="-")[[1]][1]
  chr <- paste(strsplit(x=vec, split="_")[[1]][-1], collapse = "_")
})
idx <- which(mouse_chr %in% main_chr)
combined_peaks_in_mouse_2 <- combined_peaks_in_mouse[idx,]
saveRDS(combined_peaks_in_mouse_2, file="Res_combined_peaks_in_mm10_hg38_coor.rds")

# combine the list
peak_list <- list("hg38"=main_chr_peaks, 
                  "panTro6"=combined_peaks_in_chimp_2[,1],
                  "calJac4"=combined_peaks_in_marmoset_2[,1],
                  "mm10"=combined_peaks_in_mouse_2[,1])
saveRDS(peak_list, file="Res_combined_peak_list_of_each_species.rds")

peak_coor_list <- lapply(names(peak_list), function(x){
  vec <- peak_list[[x]]
  mat <- data.frame(do.call('rbind', strsplit(x=vec, split="_")), stringsAsFactors = F)
  coor <- mat[,2]
  names(coor) <- mat[,1]
  return(coor)
})
names(peak_coor_list) <- names(peak_list)
peak_orth_mat <- matrix(NA, nrow=length(peak_list[["hg38"]]), ncol=4)
rownames(peak_orth_mat) <- names(peak_coor_list[["hg38"]])
colnames(peak_orth_mat) <- names(peak_coor_list)
for(x in names(peak_coor_list)){
  vec <- peak_coor_list[[x]]
  peak_orth_mat[names(vec),x] <- vec
}
freq <- rowSums(!is.na(peak_orth_mat))
peak_orth_mat[head(which(freq==4)),]
saveRDS(peak_orth_mat, file="Res_organoid_orthologous_peak_mat.rds")

