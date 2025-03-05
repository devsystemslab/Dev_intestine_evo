setwd("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_aggregate/signac_based/Human_without_adult_enteroid")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(dplyr)
source("~/Work/commonScript/Script_functions.R")

# step 1. create common peak list across human samples
# read paths to human scATAC-seq processed data
library(readxl)
readout <- "scATAC-seq"
species <- "Human"

sample_info <- read_excel("~/Work/Endoderm/intestine_evolution/sample/HIO_CIO_metadata.xlsx", sheet=readout)
sample_to_exclude <- c(
  "409B2-WT", "409B2-GLI3-KO", "B7-WT", "ms-2-adultIntestine-Adult", 
  "CT-1-cDuo-Adult", "CT-1-cIle-Adult", "CT-1-cCol-Adult",
  "Duo81-hEnteroids-Adult-enteroid", "Col87-hEnteroids-Adult-enteroid"
)
selected_sample_info <- sample_info[which(sample_info$Species==species & !sample_info$`Sample name` %in% sample_to_exclude),]

all_sample_name <- c()
seurat_atac_list <- list()
for(i in seq(nrow(selected_sample_info))){
  h5_folder <- paste0(selected_sample_info$`Cellranger output folder`[i], "/outs/")
  counts <- Read10X_h5(filename = paste0(h5_folder,"filtered_peak_bc_matrix.h5"))
  metadata <- read.csv(
    file = paste0(h5_folder,"singlecell.csv"),
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0(h5_folder,"fragments.tsv.gz"),
    min.cells = 1,
    min.feature=1
  )
  
  seurat_atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  if(selected_sample_info$`Demultiplexed cell barcode`[i]!="NA"){
    bc_folder <- selected_sample_info$`Demultiplexed cell barcode`[i]
    vec <- strsplit(bc_folder, "/")[[1]]
    sample_name <- vec[length(vec)]
    cat(paste("Get pooled sample", sample_name, "\n"))
    file_name <- paste0("List_", sample_name, "_", tolower(species), "_cell_barcodes.txt")
    bc <- readLines(paste(bc_folder, file_name, sep="/"))
    cat(paste("Subset", species, "cells\n"))
    seurat_atac <- subset(seurat_atac, cells = bc)
  }else{
    sample_name <- selected_sample_info$`Sample name`[i]
    cat(paste("Get individually sequenced sample", sample_name, "\n"))
  }
  
  cat("Add meta_data\n")
  sample_name <- gsub("_", "-", sample_name)
  all_sample_name <- c(all_sample_name, sample_name)
  seurat_atac$Dataset <- sample_name
  seurat_atac$Sample_idx <- paste0("HA-",i)
  seurat_atac$Species <- selected_sample_info$`Species`[i]
  seurat_atac$Individual_or_cell_line <- selected_sample_info$`Individual or stem cell line`[i]
  seurat_atac$Age <- selected_sample_info$Age[i]
  seurat_atac$Age_week <- selected_sample_info$Age_week[i]
  seurat_atac$Tissue <- selected_sample_info$Tissue[i]
  seurat_atac_list[[sample_name]] <- seurat_atac
}

names(all_sample_name) <- paste0("HA-", seq(length(all_sample_name)))
# merge all datasets, adding a cell ID to make sure cell names are unique
# NOTE: in the assay "peaks", the human-chimp pooled samples have peaks from chimp cells; 
seurat_atac_aggr <- merge(
  x = seurat_atac_list[[1]],
  y = seurat_atac_list[-1],
  add.cell.ids = names(all_sample_name)
)
saveRDS(seurat_atac_aggr, file="Res_seurat_atac_aggr.rds")

# call peaks in each individual sample
peaks <- CallPeaks(
  object = seurat_atac_aggr,
  group.by="Sample_idx",
  macs2.path = "/home/yuq/anaconda3/bin/macs2"
)
saveRDS(peaks, file="Res_human_peaks_called_in_each_sample.rds")

# quantify counts in each peaks of the union list
counts_atac_aggr <- FeatureMatrix(seurat_atac_aggr@assays$peaks@fragments,
                                  features = peaks,
                                  cells = colnames(seurat_atac_aggr))
assay_peaks_stage_atac <- CreateChromatinAssay(counts_atac_aggr, fragments = seurat_atac_aggr@assays$peaks@fragments)

# get annotation file necessary for ATAC samples
annotations <- readRDS("~/Work/Annotation/EnsDB_for_ATAC/data.ens93_annot_for_atac.rds")

Annotation(assay_peaks_stage_atac) <- annotations
seurat_atac_aggr[['peaks_stage']] <- assay_peaks_stage_atac
saveRDS(seurat_atac_aggr, file="Res_seurat_atac_aggr_with_peaks_called_in_individual_sample.rds")

# NOTE: in the assay "peaks", the human-chimp pooled samples have peaks from chimp cells; 
# in the assay "peaks_stage", all peaks were called from human samples

# step 2. preprocessing
DefaultAssay(seurat_atac_aggr) <- "peaks_stage"
seurat_atac_aggr$log_nCount_peaks <- log10(seurat_atac_aggr$nCount_peaks)
seurat_atac_aggr <- NucleosomeSignal(object = seurat_atac_aggr) # compute nucleosome signal score per cell
seurat_atac_aggr <- TSSEnrichment(object = seurat_atac_aggr, fast = FALSE) # compute TSS enrichment score per cell
seurat_atac_aggr$high.tss <- ifelse(seurat_atac_aggr$TSS.enrichment > 2, 'High', 'Low')
seurat_atac_aggr$pct_reads_in_peaks <- seurat_atac_aggr$peak_region_fragments / seurat_atac_aggr$passed_filters * 100
seurat_atac_aggr$blacklist_ratio <- seurat_atac_aggr$blacklist_region_fragments / seurat_atac_aggr$peak_region_fragments
seurat_atac_aggr[['peaks']] <- NULL
saveRDS(seurat_atac_aggr, file="Res_seurat_atac_aggr_with_peaks_called_in_individual_sample.rds")
