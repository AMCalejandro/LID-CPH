### Candidate polymorphism exploration ###


# LOAD LIBRARIES
.libPaths("/mnt/rreal/RDS/acarrasco/R_libs/")
library(data.table)
library(tidyverse)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(colochelpR)
library(LDlinkR)

# LOAD DATA

lidVariants <- read.delim("/mnt/rreal/RDS/acarrasco/ANALYSES_WORKSPACE/GWAS_SURVIVAL/POST_GWAS/V2/CANDIDATE_VARIANT_ANALYSIS/Panel_LID.tsv")
lidMetanalysis_basic = fread("/mnt/rreal/RDS/acarrasco/ANALYSES_WORKSPACE/GWAS_SURVIVAL/METAANALYSIS/METAL_AUG2022/METAL_ALL_BasicModel_LongDisDurExcluded/metal_tpd.opdc.ppmi.pdstat.pdbp_modelConfounders_NOLONGDUR_stderr_1.tbl")
lidMetanalysis_best = fread("/mnt/rreal/RDS/acarrasco/ANALYSES_WORKSPACE/GWAS_SURVIVAL/METAANALYSIS/METAL_AUG2022/METAL_ALL_BestModel/metal_tpd.opdc.ppmi.pdstat.pdbp_BESTMODEL_MOTORBL_stderr_1.tbl")


# QUERY LIDPD DATABASE 

# SNP id transfromation
dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37
lidVariants_snpUpdate = lidVariants %>% 
  filter(!dbSNP.ID == "-") %>%
  colochelpR::convert_rs_to_loc(., SNP = "dbSNP.ID", dbSNP=dbSNP) %>%
  dplyr::rename(chrpos = loc)

lidVariants_lidmetanalysis_best = lidVariants_snpUpdate %>% 
  inner_join(lidMetanalysis_best %>% 
               select(MarkerName,Freq1,BETA=Effect,
                      SE=StdErr,Pgwas=`P-value`,
                      Direction, Nsamples = TotalSampleSize), by=c("chrpos"="MarkerName")) %>%
  dplyr::arrange(Pgwas) %>%
  dplyr::relocate(c(1:2,4:6), .after=Nsamples)

lidVariants_lidmetanalysis_basic = lidVariants_snpUpdate %>% 
  inner_join(lidMetanalysis_basic %>% 
               select(MarkerName,Freq1,BETA=Effect,
                      SE=StdErr,Pgwas=`P-value`,
                      Direction, Nsamples = TotalSampleSize), by=c("chrpos"="MarkerName")) %>%
  dplyr::arrange(Pgwas) %>%
  
  dplyr::mutate(MAF = ifelse(Freq1 < 0.5, Freq1, 1-Freq1),
                BETA = ifelse(Freq1 < 0.5, BETA, -1*BETA)) %>%
  
  dplyr::select(-c(chrpos, Freq1)) %>%
  dplyr::rename(rsID = dbSNP.ID ) %>%
  dplyr::relocate(c(1:2,4:6), .after = Nsamples) %>% 
  dplyr::relocate(MAF,  .after = rsID) 
  
fwrite(lidVariants_lidmetanalysis_basic, "LIDPD_LIDGWAS_query.csv", quote= F, sep="," ,col.names=T, row.names=F)


# SEARCH FOR G20192 AND MAOA MAOB LID RELATED POLYMORPHISMS PROXIES

# check genes for which the risk variants was not present in our meta-analysis
genes_missing = setdiff(unique(lidVariants_snpUpdate$HGNC.symbol), unique(lidVariants_lidmetanalysis_basic$HGNC.symbol))

# Get variants in high LD, and check in our meta-analysios
# We do not have info for chromosome X. Exclude COMT
variants =  lidVariants_snpUpdate  %>% 
  filter( (HGNC.symbol %in% genes_missing) & (!dbSNP.ID == "-")
          & (!grepl("X:", chrpos))) %>%
  mutate(chrpos = paste0("chr", chrpos)) %>%
  pull(chrpos)

mytoken = "32e6b695b2aa"
proxies.lrrk2 = lapply(variants, FUN = function(x)
  LDproxy(x, pop = "CEU",
          r2d="r2",
          token = mytoken)
)

# We did not find good proxies for any LRRK2 variant in LD link using an 
#European population as the reference panel









  
  

