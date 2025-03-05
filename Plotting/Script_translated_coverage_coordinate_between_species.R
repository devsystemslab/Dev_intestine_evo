# step 1. liftOver the chimp coverage data to human coordinate using liftOver
# run the following in command line
#!/bin/bash
# BSUB -J C2_tCIO_fragment_liftOver 
# BSUB -n 10 
# BSUB -R "span[hosts=1]"
# BSUB -q long 
# BUSB -M 20GB 
#/projects/site/pred/ihb-intestine-evo/Tools/liftOver -minMatch=0.5 fragments.bed /projects/site/pred/ihb-intestine-evo/Annotation/liftOver_chain_file/panTro6ToHg38.over.chain.gz fragments_panTro6ToHg38.bed unmap.bed 

#sort -k1,1 -k2,2n fragments_panTro6ToHg38.bed > sorted_fragments_panTro6ToHg38.bed
#awk '{print $1"\t"($2+1)"\t"$3"\tCA-3_"$4"\t"$5}' sorted_fragments_panTro6ToHg38.bed > C2_tCIO_fragments.tsv

# change conda environment before job submission
# conda activate SCG
# compress the file
bgzip -@ 10 C2_tCIO_fragments.tsv

# build index
tabix -p bed C2_tCIO_fragments.tsv.gz

# add fragment files to species combined Seurat object
tIO_atac <- readRDS("/projects/site/pred/ihb-intestine-evo/used_object/Res_tIO_epi_ATAC_seurat_object.rds")

# step 2. add fragment files to human and chimp data
# add fragment files (chimp coverage has been liftOver to human)
seu_obj_list <- SplitObject(tIO_atac, split.by = "Sample.name")
path_list <- list(
  "H9-tHIO-W15.5+4"="/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/merged_region_list_on_all_human_data/merged_fragment_file/H9_tHIO_fragments.tsv.gz",
  "iPSC72.3-tHIO-W10.5"="/home/yuq22/ihb-intestine-evo/fetal_human_duo_crypt/merged_region_list_on_all_human_data/merged_fragment_file/iPSC_tHIO_fragments.tsv.gz",
  "C2-CIO-W12.5+4"="/home/yuq22/ihb-intestine-evo/iPSC_intestinal_organoid/tIO_scRNA_scATAC/fragment_files/C2_tCIO/liftOver/C2_tCIO_fragments.tsv.gz",
  "C7-CIO-W12.5+4"="/home/yuq22/ihb-intestine-evo/iPSC_intestinal_organoid/tIO_scRNA_scATAC/fragment_files/C7_tCIO/liftOver/C7_tCIO_fragments.tsv.gz"
)

annotations <- readRDS("/projects/site/pred/ihb-intestine-evo/Annotation/EnsDB_for_ATAC/data.ens93_annot_for_atac.rds")
peak_oths <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/intestine_evolution/peak_orthologs/HIO_CIO_orthologs/relax_and_merge/Res_combined_peaks_in_panTro6_hg38_coor.rds")

for(x in names(seu_obj_list)){
  seu_obj <- seu_obj_list[[x]]
  fpath <- path_list[[x]]
  counts <- seu_obj@assays$peaks_species@counts
  rownames(counts) <- peak_oths$combined_peaks_hg38_coor
  meta_data <- seu_obj@meta.data
  chrom_assay <- CreateChromatinAssay(counts=counts,
                                      fragments = fpath)
  Annotation(chrom_assay) <- annotations
  seu_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = meta_data
  )
  seu_obj_list[[x]] <- seu_obj
}

combined <- merge(x=seu_obj_list[[1]],
                  y=seu_obj_list[-1])
combined <- RunTFIDF(combined)

coor <- tIO_atac@reductions$umap@cell.embeddings[colnames(combined),]
combined[["LSI_integrated_UMAP"]] <- CreateDimReducObject(
  embeddings = coor,
  assay = "ATAC",
  key = "INUMAP_"
)
saveRDS(combined, file="Dat_merged_seu_obj_with_normalization.rds")

# step 3. prepare the input to get signal track
# run the following in R
# split the fragment file according to cell type per species, and get .bed files
SplitFragments(
  combined,
  assay = "ATAC",
  group.by = "Cell_type_per_species",
  idents = NULL,
  outdir = getwd(),
  file.suffix = "_tIO_scATAC",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)

# get per group scaling factor
cell_number <- table(combined$Cell_type_per_species)
seq_depth <- sapply(names(cell_number), function(x){
  mean(combined$nCount_ATAC[which(combined$Cell_type_per_species==x)])/100000
})
scale_factor <- cell_number*seq_depth
write.table(scale_factor, file="List_tIO_epi_scATAC_seq_scale_factor.tsv",row.names = F,col.names = F,quote=F,sep="\t")

# run the following command in command line

# load bed tools
#ml bedtools
## convert bed file to bedGraph
#for input in `ls *.bed`
#do
#  ID=`echo $input | cut -d'.' -f1`
#  output=$ID.bedGraph
#  bedtools genomecov -i $input -bg -g /projects/site/pred/ihb-intestine-evo/genomes/hg38_standard/hg38.chrom.sizes > $output
#  sort -k 1,1 -k 2,2n $output > sorted.bedGraph
#  mv sorted.bedGraph  $output
#done

# normalize bedGraph by number of single-cells
#while IFS=$'\t' read -r ID factor
#do
#       input=${ID}_tIO_scATAC.bedGraph
#       output=${ID}_normed.bedGraph
#       awk 'OFS="\t"{print $1,$2,$3,($4/'$factor')*100}' $input > $output 
#       #echo 'OFS="\t"{print $1,$2,$3,($4/'$factor')*100}'
#done < List_tIO_epi_scATAC_seq_scale_factor.tsv

# convert bedGraphToBigWig
#for file in `ls *normed.bedGraph`
#do
#ID=`echo $file | cut -d'.' -f1`
#output1=${ID}_noOverlaps.bedGraph
#output2=${ID}.bw
##echo $output
#bedops --partition $file | bedmap --echo --echo-map-id-uniq --delim '\t' - $file | awk '{ n = split($4, a, ";"); max = a[1]; for(i = 2; i <= n; i++) { if (a[i] > max) { max = a[i];} } print $1"\t"$2"\t"$3"\t"max; }' - > $output1
#
#/projects/site/pred/ihb-intestine-evo/Tools/bedGraphToBigWig $output1 /projects/site/pred/ihb-intestine-evo/genomes/hg38_standard/hg38.chrom.sizes $output2
#done
## normalize bedGraph by number of single-cells
##while IFS=$'\t' read -r ID factor
##do
##       input=${ID}_tIO_scATAC.bedGraph
##       output=${ID}_normed.bedGraph
##       awk 'OFS="\t"{print $1,$2,$3,($4/'$factor')*100}' $input > $output 
##       #echo 'OFS="\t"{print $1,$2,$3,($4/'$factor')*100}'
##done < List_tIO_epi_scATAC_seq_scale_factor.tsv
#
## convert bedGraphToBigWig
#for file in `ls *normed.bedGraph`
#do
#ID=`echo $file | cut -d'.' -f1`
#output1=${ID}_noOverlaps.bedGraph
#output2=${ID}.bw
##echo $output
#bedops --partition $file | bedmap --echo --echo-map-id-uniq --delim '\t' - $file | awk '{ n = split($4, a, ";"); max = a[1]; for(i = 2; i <= n; i++) { if (a[i] > max) { max = a[i];} } print $1"\t"$2"\t"$3"\t"max; }' - > $output1
#
#/projects/site/pred/ihb-intestine-evo/Tools/bedGraphToBigWig $output1 /projects/site/pred/ihb-intestine-evo/genomes/hg38_standard/hg38.chrom.sizes $output2
#done
#
#



