---
title: "coloc_eQTLcatalogue"
author: "AMCalejandro"
date: "18/2/2022"
output: html_document
---




# Load libraries
.libPaths("/mnt/rreal/RDS/acarrasco/R_libs/")

library("dplyr")
library("ggplot2")
library("rlang")
library("readr")
library("coloc")
library("data.table")
library("colochelpR")
library("GenomicRanges")
library("seqminer")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("AnnotationDbi")
library("rtracklayer")

# Load gwas
gwas = fread("/mnt/rreal/RDS/acarrasco/ANALYSES_WORKSPACE/GWAS_SURVIVAL/METAANALYSIS/METAL_AUG2022/METAL_ALL_BestModel/QC_metal_tpd.opdc.ppmi.pdstat.pdbp_BESTMODEL_MOTORBL_stderr_1.tbl")
gene_coords = fread("/mnt/rreal/RDS/acarrasco/ANALYSES_WORKSPACE/GWAS_SURVIVAL/POST_GWAS/V2/eQTLCatalogue/NCBI37.3.gene.loc.extendedMHCexcluded.txt")


# Get tabix paths for eQTL Catalogue data
tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", 
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", 
                                  sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Define some utils

## eQTL Catalogue util to import data
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  if(verbose){
    print(ftp_path)
  }
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::select(-rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics 
}


## GWAS QC utils
GWAS_filtering_coloc = function(df, variant, expand_left = 0.5, expand_right = 0.5) {
  myVariant_metadata = df %>% filter(MarkerName == .env[["variant"]])
  myCHR = myVariant_metadata$CHR
  myPOS = myVariant_metadata$BP
  
  gwas_filtered = df %>% 
    filter(CHR == myCHR) %>%
    filter(BP %in% ((myPOS-(expand_left*1e6)):(myPOS+(expand_right*1e6))) )
  gwas_filtered
}

remove_dups = function(gwasdf, snp_name) {
  gwasdf %>%
    dplyr::distinct() %>% #rsid duplicates
    dplyr::group_by({{snp_name}}) %>%
    dplyr::arrange(MAF) %>%
    dplyr::filter(row_number()==1) %>% ungroup()
}

## Util tu run coloc
run_coloc <- function(eqtl_sumstats, gwas_sumstats){
  eQTL_dataset = list(beta = eqtl_sumstats$beta,
                      varbeta = eqtl_sumstats$se^2,
                      N = (eqtl_sumstats$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                      MAF = eqtl_sumstats$maf, 
                      type = "quant", 
                      snp = eqtl_sumstats$id)
  gwas_dataset = list(beta = gwas_sumstats$Effect,
                      varbeta = gwas_sumstats$StdErr^2, 
                      type = "quant",
                      snp = gwas_sumstats$snpName_hg38,
                      MAF = gwas_sumstats$MAF, 
                      N = mean(gwas_sumstats$TotalSampleSize, na.rm=T))
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}




# 

# gwas_qc = gwas %>%
#   convert_loc_to_rs(., dbSNP = dbSNP) %>%
#   dplyr::mutate(MAF= ifelse(Freq1 > 0.5, 1-Freq1, Freq1),
#                 Effect = ifelse(Freq1>0.5, -1*Effect, Effect )) 
# 
# 
# length(unique(gwas_filt_list[[1]]$MarkerName))



# # Do some GWAS QC

# Get gwas data across +-0.5Mb from LID GWAS lead snps
# Then deal with duplicates
# Then get maf
hg19tohg38file = import.chain("/mnt/rreal/RDS/acarrasco/ANALYSES_WORKSPACE/GWAS_SURVIVAL/POST_GWAS/V2/eQTLCatalogue/hg19ToHg38.over.chain")
gwas_filt_list = lapply(c("4:32435284", "16:17044975", "1:53778300"), function(x) {
  
  # Makr minor QC
  gwasFilt = GWAS_filtering_coloc(df = gwas, variant = x,
                       expand_left = 1.5,
                       expand_right = 0) %>%
    dplyr::mutate(MAF= ifelse(Freq1 > 0.5, 1-Freq1, Freq1),
                  Effect = ifelse(Freq1>0.5, -1*Effect, Effect),
                  CHR = paste0("chr",CHR)) %>%
    remove_dups(., MarkerName)
  
  granges_tmp = 
    makeGRangesFromDataFrame(gwasFilt,
                             seqnames.field="CHR",
                             start.field="BP",
                             end.field="BP", 
                             keep.extra.columns = T)
  
  gwasLifted = as_tibble(rtracklayer::liftOver(granges_tmp, hg19tohg38file))
  
  # Do some processing in the outpout
  gwasLifted = gwasLifted %>%
    dplyr::select(-c(group,group_name, end, strand)) %>%
    dplyr::rename(CHR = seqnames, BP = start) %>%
    dplyr::mutate(CHR = as.integer(gsub("chr", "", CHR)),
                  snpName_hg38 = paste(CHR,BP, sep =":")) %>%
    dplyr::relocate(snpName_hg38, .before = CHR)
  
  gwasLifted
    
  
})
gwas_filt_list[[2]]


## After we looked in genomebrowser we found this gene called PCDH7 that is a bit more than 1Mb downsream from the lead SNP in chromosome 5.
## We are going to perform eQTL analysis against this gene

# This is the gene range in hg19 : chr4:30,721,991-31,148,422
# This is the enssembl id: ENSG00000169851


# The question is. Do the significant gwas locus with a span of +-0.5 colocalices with  PCDH7 ?




##########################################################################
############## Uniformly processed RNA-seq datasets ######################
###########################################################################
rnaseq_df = dplyr::filter(tabix_paths, quant_method == "ge") %>%
  dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list = setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
#Extract column names from first file
column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))
#Wrap the download function around purrr::safely to avoid catch erros
safe_import = purrr::safely(import_eQTLCatalogue)

#Import summmary stats

# ge data for PCDH7 from ge studies
region = "4:30720369-32141411"
summary_list_PCDH7= purrr::map(ftp_path_list, 
                               ~safe_import(., 
                                            region, 
                                            selected_gene_id = "ENSG00000169851",
                                            column_names))

region = "4:30720369-32400000"
summary_list_PCDH7= purrr::map(ftp_path_list, 
                               ~safe_import(., 
                                            region, 
                                            selected_gene_id = "ENSG00000169851",
                                            column_names))

tissue_noDATA = function(df) {
  is.null(df) || nrow(df) == 0
}
summary_list_PCDH7 = purrr::map(summary_list_PCDH7, ~.$result)
summary_list_PCDH7 = 
  summary_list_PCDH7[!unlist(purrr::map(summary_list_PCDH7, tissue_noDATA))]

# Run coloc against the df ge studies
coloc_df_rnaseq_PCDH7 = purrr::map_df(summary_list_PCDH7, 
                                     ~run_coloc(., gwas_filt_list[[1]]),
                                     .id = "qtl_id")




microarray_df = dplyr::filter(tabix_paths, quant_method == "tx") %>%
  dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list = setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id)

column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))

summary_list = purrr::map(ftp_path_list, 
                          ~import_eQTLCatalogue(., 
                                                region,
                                                selected_gene_id = "ENST00000265854", 
                                                column_names))

coloc_df_microarray = purrr::map_df(summary_list, ~run_coloc(., myGWAS_filtered_nodup), .id = "qtl_id")
