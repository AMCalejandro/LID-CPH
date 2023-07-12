########################################################################
###### Utils to automate the generation of summary statistics  #########
########################################################################

my_summary <- function(data, cohort_name, event, ...) {
  # get meansd
  meansd = data %>% 
    dplyr::summarise(across(c(...), get_meansd)) %>%
    as.data.frame
  # get proportions of LID
  prop_LID = get_propCases(data, enquo(event))
  
  # Assume there is a gender variable called sex
  propmales = ((data %>% count(SEX))$n[1] %>% as.numeric() / dim(data)[1])
  
  mydf = cbind(cohort_name,
               propmales,
               N = dim(data)[1],
               meansd,
               prop_LID)
}

get_meansd = function(var) {
  mymean = round(mean(var, na.rm = T), digits = 2)
  mysd = round(sd(var, na.rm = T), digits = 2)
  data.frame(mymean, mysd)
}
get_propCases = function(df, var) {
  # Defusing and injecting the event_var
  
  n_lid = df %>% dplyr::filter(!!var == 1) %>% 
    pull(ID) %>% length
  
  
  return(round(n_lid / dim(df)[1], digits = 2 ))
}





########################################################################
################ Util to get forestplots  ##############################
########################################################################

myforestplot = function(df, snpname = "SNP", n_studies) {
  
  # Three studiesn = 3
  numbers = c(1, n_studies + 3, n_studies + 4) %>% as.character()
  mylist = setNames(vector(mode = "list", length = 3), numbers)
  
  lines = lapply(mylist, function(x) {
    x = gpar(lwd = 1, columns = 1:4, col = "#444444") 
  })
  
  df %>% 
    forestplot::forestplot(
      labeltext = c(study, effect, ci), 
      is.summary = summary,
      clip = c(-1, 2), 
      xlog = FALSE,
      graph.pos = 4,
      colgap = unit(0.03, "npc"),
      boxsize = 0.25,
      xticks=c(-1.5,-1, -0.5, 0, 0.5, 1, 1.5),
      #txt_gp = fpTxtGp(cex=0.75),
      hrzl_lines = lines,
      col = fpColors(box = "royalblue",
                     line = "darkblue",
                     summary = "royalblue"),
      xlab = "Effect size and CI",
      title = paste0(snpname, " forest plot"),
      txt_gp = fpTxtGp(ticks=gpar(cex=0.8),
                       xlab = gpar(col = "red", lty = "solid", lwd = 3, fontsize = 20)))
  
}





########################################################################
## Util to wrangle assessments data and have in roch and long format ###
########################################################################

get_assesmentLong = function(df = NULL, outcome_names = NULL,
                             format = "wide", myid = "ID") {
  
  pattern_vector  = outcome_names 
  lapply(pattern_vector, function(pattern) {
    
    get_wide_measure = colnames(df)[grep(pattern, colnames(df))]
    df = df %>% 
      dplyr::select(all_of(c(myid, get_wide_measure)))
    
    if (format == "wide") {
      df = df %>%
        tidyr::gather(measure, score, contains(pattern)) %>%
        dplyr::mutate(visit_number = str_extract(measure, regex("(\\d+)")))
    }
    else {
      # Check the input long format has the required columns
      ###TODO
    }
    
    # Cahnge all column names to lower case to avoid issues
    setnames(df, tolower(names(df)))
    df = data.table::setnames(
      df %>%
        mutate(measure = tolower(measure)),
      tolower(names(df))
    )
    myid = tolower(myid)
    
    
    moca_exclude = c("face", "velvet", "church", "red", "daisy",
                     "diag_adj", "screen_adj", "cat_screen", "cat_diag", 
                     "total")
    
    # Get all columns that do not have the total values
    df <- df %>%
      dplyr::select(!ends_with("total")) %>%
      dplyr::filter(!str_detect(.data[["measure"]], paste(moca_exclude, collapse = "|"))) %>%
      distinct(across(c(all_of(myid), measure)), .keep_all = T)
    
    measure_name = paste0(tolower(pattern), "_total")
    # Get the total value of the clinical assessment
    df_final <- df %>%
      dplyr::group_by(across( c({{ myid }}, "visit_number") )) %>% 
      #dplyr::mutate(moca_total= sum(score)) %>%
      dplyr::mutate("{measure_name}" := sum(score)) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::select(-c(measure, score)) %>%
      ungroup()
    
    df_final
  })
  
}

visityears_convert <- function(var) {
  var = as.character(var)
  key <- c("4", "7", "9", "10", "11")
  value <- c(1.5, 3, 4.5, 7, 9.5)
  hash <- new.env(hash = TRUE, parent = emptyenv(), size = 100L)
  assign_hash(key, value, hash)
  get_hash(var, hash)[[var]]
}

get_ratechange = function(data, year_var = "visit", idvar = ID, outcome_vars = NULL) {
  mydf = data.frame()
  for (index in seq_len(length(outcome_vars))) {
    outcome = outcome_vars[index]
    tmp = data %>%
      dplyr::mutate(ratechange = 
                      (.data[[outcome]] - .data[[outcome]][visit_number == "1"]) / as.numeric(year_var),
                    years = year_var,
                    outcome = outcome) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::select({{ idvar }}, years, outcome, ratechange) %>%
      ungroup()
    
    mydf = rbind(mydf, tmp)
  }
  mydf
}

outcome_select = function(data, pattern) {
  data %>%
    select(contains(pattern)) %>%
    colnames
}
filter_df <- function(data, timevar, idvar) {

  tmp = data %>%
    dplyr::filter( (visit_number == {{ timevar }}) | (visit_number == 1)) %>%
    #dplyr::group_by(ID) %>%
    dplyr::group_by({{ idvar }}) %>%
    dplyr::arrange(desc(visit_number))
  tmp
}




########################################################################
###### Util to get assesment total score data in long format  ##########
########################################################################

get_assesmentLong = function(df = NULL, outcome_names = NULL,
                             format = "wide", myid = "ID") {
  
  pattern_vector  = outcome_names 
  lapply(pattern_vector, function(pattern) {
    
    get_wide_measure = colnames(df)[grep(pattern, colnames(df))]
    df = df %>% 
      dplyr::select(all_of(c(myid, get_wide_measure)))
    
    if (format == "wide") {
      df = df %>%
        tidyr::gather(measure, score, contains(pattern)) %>%
        dplyr::mutate(visit_number = str_extract(measure, regex("(\\d+)")))
    }
    else {
      # Check the input long format has the required columns
      ###TODO
    }
    
    # Cahnge all column names to lower case to avoid issues
    setnames(df, tolower(names(df)))
    df = data.table::setnames(
      df %>%
        mutate(measure = tolower(measure)),
      tolower(names(df))
    )
    myid = tolower(myid)
    
    
    moca_exclude = c("face", "velvet", "church", "red", "daisy",
                     "diag_adj", "screen_adj", "cat_screen", "cat_diag", 
                     "total")
    
    # Get all columns that do not have the total values
    df <- df %>%
      dplyr::select(!ends_with("total")) %>%
      dplyr::filter(!str_detect(.data[["measure"]], paste(moca_exclude, collapse = "|"))) %>%
      distinct(across(c(all_of(myid), measure)), .keep_all = T)
    
    measure_name = paste0(tolower(pattern), "_total")
    # Get the total value of the clinical assessment
    df_final <- df %>%
      dplyr::group_by(across( c({{ myid }}, "visit_number") )) %>% 
      #dplyr::mutate(moca_total= sum(score)) %>%
      dplyr::mutate("{measure_name}" := sum(score)) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::select(-c(measure, score)) %>%
      ungroup()
    
    df_final
  })
  
}


#########################################################################################
# Util to automate the performance of ttest and wilcoxon tests across group of vars  ###
#########################################################################################

my_statstest <- function(df, variable, by, env = parent.frame()){
  variable <- rlang::sym(variable)
  by <- rlang::sym(by)
  exp1 <- rlang::expr( !! variable ~ !! by)
  print(exp1) 
  
  ttest = t.test(formula = eval(exp1), data = df)[c("data.name",
                                                    "statistic",
                                                    "p.value",
                                                    "method")]
  ttest_df = data.frame(id=names(ttest),do.call(rbind,ttest))
  ttest_df$VAR = gsub(" .*", "", ttest_df$t[1])
  
  wilcox = wilcox.test(formula = eval(exp1), data = df)[c("data.name",
                                                          "statistic",
                                                          "p.value",
                                                          "method")]
  wilcox_df = data.frame(id=names(wilcox),do.call(rbind,wilcox))
  wilcox_df$VAR = gsub(" .*", "", wilcox_df$W[1])
  
  list(ttest = ttest_df, wilcoxtest = wilcox_df)
}

myvar = c("updrs_iii__total", "updrsIII_bl", "moca_total", "moca_bl")
grouper = "event_dyskinesia"
res = lapply(myvar, function(var) {
  my_statstest(tpd_visitnumber, variable = var, 
               by = grouper)
})






########################################################################
###### Utils to perform stepwise regression analyses  ##################
########################################################################

catchBadPredictor = function(p, fit0) {
  out <- tryCatch(
    {
      # Fitting the model
      model = stats::update(fit0, formula = stats::as.formula(paste(". ~ . +", p)) )
    },
    error=function(cond) {
      message("Failed with predictor: ", p)
      message("\n Original error: ", cond)
    })
  
  return(out)
}

stepper <- function(fit, upper) {
  preds <- setdiff(upper, all.vars(formula(fit)[[3]]))
  preds_out <- preds
  preds_in  <- character(0)
  fit0 <- fit
  
  n = 0
  while(TRUE) {
    n = n + 1
    fits <-
      c(list(fit0),
        lapply(preds_out, function(p) {catchBadPredictor(p, fit0 = fit0)}),
        lapply(preds_in,  function(p) {catchBadPredictor(p, fit0 = fit0)})
      )
    
    catchbadmods = fits[-which(sapply(fits, is.null))]
    if (!is_empty(catchbadmods)) {
      fitsnonull = fits[-which(sapply(fits, is.null))]
    } else {
      fitsnonull = fits
    }
    tab = sapply(fitsnonull, AIC)
    minqic <- which.min(tab)
    fit0 <- fitsnonull[[minqic]]
    preds_in  <- setdiff(all.vars(as.list(stats::formula(fit0))[[3]]), all.vars(formula(fit)[[3]]))
    preds_out <- dplyr::setdiff(preds, preds_in)
    
    cat("\n In Iteration ", n, "predictor", preds_in[n], "has been chosen \n")
    if (minqic == 1) {break}
  }
  fit0
}




########################################################################################
###### Function to automate the imputation MDS-UPDRSIII missing data  ##################
########################################################################################

# Function to automate the imputation of missing data
impute_UPDRSIII <- function(df, ID = "ID", score_column = "score", visit_column = "visit_number",
                            UPDRSIII_column = "UPDRSIII_measure", UPDRS_type = "total",
                            missing_threshold = c(6,5,1), keep_vars = "Yes") {
  # libs loading
  library(hash)      # To create key-value pairs
  library(tidyverse) # For tidy evaluation and wrangling
  library(glue)      # To dynamically name variables
  
  # Safety checks
  if(!"score" %in% colnames(df)) stop("Warning: score column not found")
  #if(!"visit_number" %in% colnames(df)) stop("Warining: visit_number column not found")
  if(length(colnames) > 4) stop("Consider using tidyr::gather to get your data
                               in long format as input for this function")
  
  # vector with UPDRS_limb terms
  limb_terms <- c("III_3a", "III_3b", "III_3c", "III_3d", "III_3e", 
                  "III_15a", "III_15b", "III_16a","III_16b", "III_17a", "III_17b",
                  "III_17c", "III_17d", "III_17e", "III_4a", "III_4b",
                  "III_5a", "III_5b", "III_6a", "III_6b", "III_7a", "III_7b",
                  "III_8a", "III_8b", "III_18",
                  "3_3a", "3_3b", "3_3c", "3_3d", "3_3e", 
                  "3_15a", "3_15b", "3_16a","3_16b", "3_17a", "3_17b",
                  "3_17c", "3_17d", "3_17e", "3_4a", "3_4b",
                  "3_5a", "3_5b", "3_6a", "3_6b", "3_7a", "3_7b",
                  "3_8a", "3_8b", "3_18",
                  "NP3RIGN", "NP3RIGRU", "NP3RIGLU", "NP3RIGRL", "NP3RIGLL",
                  "NP3PTRMR", "NP3PTRML", "NP3KTRMR", "NP3KTRML", "NP3RTARU", "NP3RTALU",
                  "NP3RTARL", "NP3RTALL", "NP3RTALJ", "NP3FTAPR", "NP3FTAPL",
                  "NP3HMOVR", "NP3HMOVL", "NP3PRSPR", "NP3PRSPL", 
                  "NP3TTAPR", "NP3TTAPL", "NP3LGAGR", "NP3LGAGL",
                  "NP3RTCON",
                  "RigidityNeck", "RigidityRightUpper", "RigidityRightLower", "RigidityLeftUpper", "RigidityLeftLower",
                  "PosturalTremorHandsRight", "PosturalTremorHandsLeft",
                  "KineticTremorHandsRight", "KineticTremorHandsLeft",
                  "RestTremorAmplitudeRightUpper", "RestTremorAmplitudeRightLower", "RestTremorAmplitudeLeftUpper", 
                  "RestTremorAmplitudeLeftLower", "RestTremorAmplitudeLipJaw",
                  "FingerTappingRight","FingerTappingLeft",
                  "HandMovementsRight", "HandMovementsLeft", "PSHandsRight", "PSHandsLeft",
                  "ToeTappingRight", "ToeTappingLeft", "LegAgilityRight", "LegAgilityLeft",
                  "ConstancyOfRestTremor",
                  "MOTOREXAM3NECK", "MOTOREXAM3RUE", "MOTOREXAM3LUE", "MOTOREXAM3RLE", "MOTOREXAM3LLE",
                  "MOTOREXAM15R", "MOTOREXAM15L", "MOTOREXAM16R", "MOTOREXAM16L", 
                  "MOTOREXAM17LIPJ", "MOTOREXAM17RUE", "MOTOREXAM17LUE", "MOTOREXAM17RLE", "MOTOREXAM17LLE",
                  "MOTOREXAM4R", "MOTOREXAM4L", "MOTOREXAM5R", "MOTOREXAM5L", "MOTOREXAM6R", "MOTOREXAM6L",
                  "MOTOREXAM7R", "MOTOREXAM7L", "MOTOREXAM8R", "MOTOREXAM8L", "MOTOREXAM18")
  
  # Vector with UPDRS_axial terms
  axial_terms <- c("III_1", "III_2", "III_9", "III_10", "III_11", "III_12", "III_13", "III_14", 
                   "3_1", "3_2", "3_9", "3_10", "3_11", "3_12", "3_13", "3_14",
                   "NP3SPCH", "NP3FACXP", "NP3RISNG", "NP3GAIT", "NP3FRZGT", "NP3PSTBL", 
                   "NP3POSTR", "NP3BRADY",
                   "Speech", "FacialExpression", "ArisingFromChair", "Gait", "FreezingOfGait",
                   "PosturalStability", "Posture", "GSM",
                   "MOTOREXAM1", "MOTOREXAM2", "MOTOREXAM9", "MOTOREXAM10", "MOTOREXAM11", "MOTOREXAM12", "MOTOREXAM13", "MOTOREXAM14")
  
  
  # Get all vars to merge later
  df_allvars <- df %>% select(!c(.data[[score_column]], .data[[UPDRSIII_column]])) %>%
    distinct(.data[[ID]], .data[[visit_column]], .keep_all = TRUE)
  
  
  # Removing UPDRSIII_total measure if present
  df <- df %>%
    dplyr::select(!ends_with(c("UPDRS_III", "total")))
  
  final_df = df %>% dplyr::select(.data[[ID]], .data[[visit_column]]) %>% 
    distinct()
  
  for (type_index in 1:length(UPDRS_type)) {
    
    # To be used when naming the imputed columns
    term = UPDRS_type[type_index]
    if (term == "total") {
      mydf <- df %>% 
        dplyr::select(.data[[ID]], .data[[visit_column]], .data[[score_column]], contains(c("UPDRS")))
      
      mydf <- mydf %>%
        dplyr::group_by(.data[[ID]], .data[[visit_column]]) %>% 
        dplyr::mutate(missing_total = sum(is.na(.data[[score_column]])),
                      UPDRSIIItotal_imputed = ifelse(missing_total > missing_threshold[1], NA,
                                                     (sum(score, na.rm = TRUE)) / (33-missing_total) * 33 )) 
      
      # UPDRSIII_limb imputation
    } else if (term == "limb") {
      
      mydf <- df %>%
        dplyr::select(.data[[ID]], .data[[score_column]], .data[[visit_column]], contains(c("UPDRS"))) %>%
        dplyr::filter(grepl(paste(limb_terms, collapse = "|"), .data[[UPDRSIII_column]]))
      
      mydf <- mydf %>%
        dplyr::group_by(.data[[ID]], .data[[visit_column]]) %>% 
        dplyr::mutate(missing_total = sum(is.na(.data[[score_column]])),
                      UPDRSIIIlimb_imputed = ifelse(missing_total > missing_threshold[2], NA,
                                                    (sum(score, na.rm = TRUE)) / (25-missing_total) * 25))
      # UPDRSIII_axial imputation
    } else if (term == "axial") {
      
      mydf <- df %>%
        dplyr::select(.data[[ID]], .data[[score_column]], .data[[visit_column]], contains(c("UPDRS"))) %>%
        dplyr::filter(!grepl(paste(limb_terms, collapse = "|"), .data[[UPDRSIII_column]]))
      
      mydf <- mydf %>%
        dplyr::group_by(.data[[ID]], .data[[visit_column]]) %>% 
        dplyr::mutate(missing_total = sum(is.na(.data[[score_column]])),
                      UPDRSIIIaxial_imputed = ifelse(missing_total > missing_threshold[3], NA,
                                                     (sum(score, na.rm = TRUE)) / (8-missing_total) * 8)) 
      
    } else {
      stop("Only UPDRS_total, UPDRS_limb, or UPDRS_axial supported")
    }
    
    measure_name = paste0("UPDRSIII_measure_", term)
    mydf <- mydf %>%
      dplyr::arrange(.data[[ID]], .data[[visit_column]]) %>%
      dplyr::filter(row_number() == 1) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate("{measure_name}" := paste0("V", .data[[visit_column]], "_UPDRS_III_", term)) %>%
      dplyr::select(-c(UPDRSIII_measure, score, missing_total))
    
    # final_df <- final_df %>% 
    #   dplyr::left_join(mydf, by = c(.data[[ID]], .data[[visit_column]]))
    
    final_df <- final_df %>% 
      dplyr::left_join(mydf, by = c(ID, visit_column))
    
  }
  
  if (keep_vars == "Yes") {
    final_df = df_allvars %>% inner_join(final_df)
  }
  
  return(final_df)
}


qc_metal = function(df, N) {
  
  rndn_snp = df[1,]$MarkerName
  snpFormat = checkFormat(rndn_snp)
  
  if (!snpFormat == "rs") {
    snp_splitted = as.data.frame(stringr::str_split_fixed(df$MarkerName, pattern = ":",  n = 2))
    names(snp_splitted) = c("CHR", "BP")
    df = cbind(df, snp_splitted)
    df = df %>% 
      dplyr::relocate(c(CHR,BP), .after=MarkerName) %>%
      dplyr::mutate(BP = as.integer(BP))
  }
  
  
  data_filtered <- df %>% 
    filter((HetDf >= N -2) & (TotalSampleSize > 1000)) 
  
  data_filtered_sorted <- data_filtered %>%
    arrange(`P-value`)
  
  #Filter out SNPs with HetPVal < 0.05 (Cochran's Q-test for heterogeneity)
  #Also filter out SNPs with HetISq > 80
  data_filtered_sorted_het <-data_filtered_sorted %>%
    filter(HetPVal > 0.05) %>%
    filter(HetISq < 80)
  
  #Check MAF variability - remove variants with MAF variability > 15%
  data_filtered_sorted_het_MAF <- data_filtered_sorted_het %>%
    mutate(MAF_variability = MaxFreq - MinFreq) %>%
    filter(MAF_variability <= 0.15)
  
  fwrite(data_filtered_sorted_het_MAF, paste0("QC_metaout", ".tbl"), quote = F, row.names = F, col.names = T, sep = "\t")
  
  #Export for FUMA
  export_FUMA <- data_filtered_sorted_het_MAF %>%
    select(MarkerName, CHR, BP, Allele1, Allele2, pval = `P-value`, Effect, StdErr, TotalSampleSize)
  
  fwrite(export_FUMA, "metaanalysis_FUMA.txt.gz", quote = F, row.names = F, col.names = T, sep = "\t")
  data_filtered_sorted_het_MAF
}