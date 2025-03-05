setwd("/Volumes/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v99/")

library(biomaRt)
## This gives the full list of all marts, including archives
archives <- listEnsemblArchives()

ensemblhsapiens = useMart(host = 'https://jan2020.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL',
                          dataset = "hsapiens_gene_ensembl")

filters = listFilters(ensemblhsapiens)
attributes = listAttributes(ensemblhsapiens)

selected_att <- c("ensembl_gene_id",
                  "external_gene_name",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  
                  # zebrafish orthologs
                  "drerio_homolog_ensembl_gene",                    
                  "drerio_homolog_associated_gene_name",                      
                  "drerio_homolog_chromosome",                      
                  "drerio_homolog_chrom_start",                    
                  "drerio_homolog_chrom_end",                    
                  "drerio_homolog_orthology_type",                              
                  "drerio_homolog_orthology_confidence",
                  
                  # frog orthologs
                  "xtropicalis_homolog_ensembl_gene",                    
                  "xtropicalis_homolog_associated_gene_name",                      
                  "xtropicalis_homolog_chromosome",                      
                  "xtropicalis_homolog_chrom_start",                    
                  "xtropicalis_homolog_chrom_end",                    
                  "xtropicalis_homolog_orthology_type",                              
                  "xtropicalis_homolog_orthology_confidence",
        
                  # chicken orthologs
                  "ggallus_homolog_ensembl_gene",                    
                  "ggallus_homolog_associated_gene_name",                      
                  "ggallus_homolog_chromosome",                      
                  "ggallus_homolog_chrom_start",                    
                  "ggallus_homolog_chrom_end",                    
                  "ggallus_homolog_orthology_type",                              
                  "ggallus_homolog_orthology_confidence",
                  
                  # rat
                  "rnorvegicus_homolog_ensembl_gene",                    
                  "rnorvegicus_homolog_associated_gene_name",                      
                  "rnorvegicus_homolog_chromosome",                      
                  "rnorvegicus_homolog_chrom_start",                    
                  "rnorvegicus_homolog_chrom_end",                    
                  "rnorvegicus_homolog_orthology_type",                              
                  "rnorvegicus_homolog_orthology_confidence",
                  
                  # mouse orthologs
                  "mmusculus_homolog_ensembl_gene",                    
                  "mmusculus_homolog_associated_gene_name",                      
                  "mmusculus_homolog_chromosome",                      
                  "mmusculus_homolog_chrom_start",                    
                  "mmusculus_homolog_chrom_end",                    
                  "mmusculus_homolog_orthology_type",                              
                  "mmusculus_homolog_orthology_confidence",
                  
                  # marmoset
                  "cjacchus_homolog_ensembl_gene",                    
                  "cjacchus_homolog_associated_gene_name",                      
                  "cjacchus_homolog_chromosome",                      
                  "cjacchus_homolog_chrom_start",                    
                  "cjacchus_homolog_chrom_end",                    
                  "cjacchus_homolog_orthology_type",                              
                  "cjacchus_homolog_orthology_confidence",
                  
                  # macaque
                  "mmulatta_homolog_ensembl_gene",                    
                  "mmulatta_homolog_associated_gene_name",                      
                  "mmulatta_homolog_chromosome",                      
                  "mmulatta_homolog_chrom_start",                    
                  "mmulatta_homolog_chrom_end",                    
                  "mmulatta_homolog_orthology_type",                              
                  "mmulatta_homolog_orthology_confidence",
                  
                  # chimp
                  "ptroglodytes_homolog_ensembl_gene",                    
                  "ptroglodytes_homolog_associated_gene_name",                      
                  "ptroglodytes_homolog_chromosome",                      
                  "ptroglodytes_homolog_chrom_start",                    
                  "ptroglodytes_homolog_chrom_end",                    
                  "ptroglodytes_homolog_orthology_type",                              
                  "ptroglodytes_homolog_orthology_confidence")
att_mat <- getBM(attributes=selected_att, mart=ensemblhsapiens,
                 filters = "external_gene_name",
                 values = "LGR5")
att_mat <- getBM(attributes=selected_att, mart=ensemblhsapiens)
saveRDS(att_mat, file="Dat_biomaRt_v99_human_multispecies_orthologs.rds")
