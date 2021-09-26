## ----message=FALSE, warning=FALSE, results="hide"-------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(ashr) # Also requires installation of RMosek, which needs a (free) licence. See the ashr Github page for help
library(mashr) # NB: This has multiple dependencies and was tricky to install. Read the Github page, and good luck!
library(glue)
library(kableExtra)
library(gridExtra)


# Get the list of SNPs (or chunks of 100% SNP clumps) that are in approx. LD with one another
get_mashr_data <- function(){
  
  SNPs_in_LD <- read.table("data/derived/SNPs_in_LD", header = FALSE)
  
  data_for_mashr <- read_csv("data/derived/all_univariate_GEMMA_results.csv") %>% 
    filter(str_detect(SNPs, ",") | SNPs %in% SNPs_in_LD$V1) %>% 
    select(SNPs, starts_with("beta"), starts_with("SE")) 
  
  single_snps <- data_for_mashr %>% filter(!str_detect(SNPs, ","))
  multiple_snps <- data_for_mashr %>% filter(str_detect(SNPs, ","))
  multiple_snps <- multiple_snps[enframe(strsplit(multiple_snps$SNPs, split = ", ")) %>% 
                                   unnest(value) %>% 
                                   filter(value %in% SNPs_in_LD$V1) %>% 
                                   pull(name), ]
  bind_rows(single_snps, multiple_snps) %>% arrange(SNPs)
}


## ----echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------
# setwd("/Users/lholman/Rprojects/fitnessGWAS")
# setwd("/data/projects/punim0243/DGRP_mashr")


## ----echo = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
n_snps <- nrow(read_csv("data/derived/all_univariate_GEMMA_results.csv"))
n_snps <- prettyNum(n_snps, big.mark=",", scientific=FALSE)

loci_tested <- nrow(get_mashr_data())
loci_tested <- prettyNum(loci_tested, big.mark=",", scientific=FALSE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
run_mashr <- function(beta_and_se, mashr_mode, ED_p_cutoff = NULL){
  
  print(str(beta_and_se))
  
  mashr_setup <- function(beta_and_se){
    betas <- beta_and_se %>% select(starts_with("beta")) %>% as.matrix()
    SEs <- beta_and_se %>% select(starts_with("SE")) %>% as.matrix()
    rownames(betas) <- beta_and_se$SNP
    rownames(SEs) <- beta_and_se$SNP
    mash_set_data(betas, SEs)
  }
  
  mash_data <- mashr_setup(beta_and_se)
  
  # Setting mashr_mode == "ED" makes mashr choose the covariance matrices for us, using the
  #  software Extreme Deconvolution. This software "reconstructs the error-deconvolved or 'underlying' 
  # distribution function common to all samples, even when the individual data points are samples from different distributions"
  # Following the mashr vignette, we initialise the algorithm in ED using the principal components of the strongest effects in the dataset
  # Reference for ED: https://arxiv.org/abs/0905.2979
  if(mashr_mode == "ED"){
    # Find the strongest effects in the data
    m.1by1 <- mash_1by1(mash_data) 
    strong <- get_significant_results(m.1by1, thresh = ED_p_cutoff)   
    # Obtain data-driven covariance matrices by running Extreme Deconvolution
    U.pca <- cov_pca(mash_data, npc = 4, subset = strong)
    U <- cov_ed(mash_data, U.pca, subset = strong)
  }
  
  # Otherwise, we define the covariance matrices ourselves (a long list of a priori interesting matrices are checked)
  if(mashr_mode == "canonical"){
    make_SA_matrix <- function(r) matrix(c(1,1,r,r,1,1,r,r,r,r,1,1,r,r,1,1), ncol=4)
    make_age_antag_matrix <- function(r) matrix(c(1,r,1,r,r,1,r,1,1,r,1,r,r,1,r,1), ncol=4)
    make_sex_specific <- function(mat, sex){
      if(sex == "F") {mat[, 3:4] <- 0;  mat[3:4, ] <- 0}
      if(sex == "M") {mat[, 1:2] <- 0; mat[1:2, ] <- 0}
      mat
    }
    make_age_specific <- function(mat, age){
      if(age == "early") {mat[, c(2,4)] <- 0; mat[c(2,4), ] <- 0}
      if(age == "late")  {mat[, c(1,3)] <- 0; mat[c(1,3), ] <- 0}
      mat
    }
    
    add_matrix <- function(mat, mat_list, name){
      mat_list[[length(mat_list) + 1]] <- mat
      names(mat_list)[length(mat_list)] <- name
      mat_list
    }
    id_matrix <- matrix(1, ncol=4, nrow=4)
    
    # Get the mashr default canonical covariance matrices: this includes the ones 
    # called "null", "uniform", and "same sign" in the list that precedes this code chunk
    U <- cov_canonical(mash_data)
    
    # And now our custom covariance matrices: 
    
    # Identical across ages, but sex-antagonistic
    U <- make_SA_matrix(-0.25) %>% add_matrix(U, "Sex_antag_0.25")   
    U <- make_SA_matrix(-0.5) %>% add_matrix(U, "Sex_antag_0.5")
    U <- make_SA_matrix(-0.75) %>% add_matrix(U, "Sex_antag_0.75")   
    U <- make_SA_matrix(-1) %>% add_matrix(U, "Sex_antag_1.0")       
    
    # Identical across sexes, but age-antagonistic
    U <- make_age_antag_matrix(-0.25) %>% add_matrix(U, "Age_antag_0.25")
    U <- make_age_antag_matrix(-0.5) %>% add_matrix(U, "Age_antag_0.5")
    U <- make_age_antag_matrix(-0.75) %>% add_matrix(U, "Age_antag_0.75")
    U <- make_age_antag_matrix(-1) %>% add_matrix(U, "Age_antag_1.0")
    
    # Sex-specific, identical effect in young and old
    U <- id_matrix %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_1")
    U <- id_matrix %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_1")
    
    # Age-specific, identical effect in males and females
    U <- id_matrix %>% make_age_specific("early") %>% add_matrix(U, "Early_life_specific_1")
    U <- id_matrix %>% make_age_specific("late")  %>% add_matrix(U, "Late_life_specific_1")
    
    # Positively correlated but variable effect across ages, and also sex-specific
    U <- make_age_antag_matrix(0.25) %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_0.25")
    U <- make_age_antag_matrix(0.5) %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_0.5")
    U <- make_age_antag_matrix(0.75) %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_0.75")
    U <- make_age_antag_matrix(0.25) %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_0.25")
    U <- make_age_antag_matrix(0.5) %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_0.5")
    U <- make_age_antag_matrix(0.75) %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_0.75")
    
    # Negatively correlated across ages, and also sex-specific
    U <- make_age_antag_matrix(-0.25) %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_age_antag_0.25")
    U <- make_age_antag_matrix(-0.5) %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_age_antag_0.5")
    U <- make_age_antag_matrix(-0.75) %>% make_sex_specific("F") %>% add_matrix(U, "Female_specific_age_antag_0.75")
    U <- make_age_antag_matrix(-0.25) %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_age_antag_0.25")
    U <- make_age_antag_matrix(-0.5) %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_age_antag_0.5")
    U <- make_age_antag_matrix(-0.75) %>% make_sex_specific("M") %>% add_matrix(U, "Male_specific_age_antag_0.75")
    
    # Positively correlated but variable effect across sexes, and also age-specific
    U <- make_SA_matrix(0.25) %>% make_age_specific("early") %>% add_matrix(U, "Early_life_specific_0.25")
    U <- make_SA_matrix(0.5) %>% make_age_specific("early") %>% add_matrix(U, "Early_life_specific_0.5")
    U <- make_SA_matrix(0.75) %>% make_age_specific("early") %>% add_matrix(U, "Early_life_specific_0.75")
    U <- make_SA_matrix(0.25) %>% make_age_specific("late") %>% add_matrix(U, "Late_life_specific_0.25")
    U <- make_SA_matrix(0.5) %>% make_age_specific("late") %>% add_matrix(U, "Late_life_specific_0.5")
    U <- make_SA_matrix(0.75) %>% make_age_specific("late") %>% add_matrix(U, "Late_life_specific_0.75")
    
    # Negatively correlated but variable effect across sexes, and also age-specific
    U <- make_SA_matrix(-0.25) %>% make_age_specific("early") %>% add_matrix(U, "Early_life_antag_0.25")
    U <- make_SA_matrix(-0.5) %>% make_age_specific("early") %>% add_matrix(U, "Early_life_antag_0.5")
    U <- make_SA_matrix(-0.75) %>% make_age_specific("early") %>% add_matrix(U, "Early_life_antag_0.75")
    U <- make_SA_matrix(-0.25) %>% make_age_specific("late") %>% add_matrix(U, "Late_life_antag_0.25")
    U <- make_SA_matrix(-0.5) %>% make_age_specific("late") %>% add_matrix(U, "Late_life_antag_0.5")
    U <- make_SA_matrix(-0.75) %>% make_age_specific("late") %>% add_matrix(U, "Late_life_antag_0.75")
  }
  
  return(mash(data = mash_data, Ulist = U)) # Run mashr
}


if(!file.exists("data/derived/mashr_results_canonical.rds")){
  
  data_for_mashr <- get_mashr_data()
  
  print("Starting the data-driven analysis")
  run_mashr(data_for_mashr, mashr_mode = "ED", ED_p_cutoff = 0.2) %>%
    saveRDS(file = "data/derived/mashr_results_ED.rds")
  
  print("Starting the canonical analysis")
  run_mashr(data_for_mashr, mashr_mode = "canonical") %>%
    saveRDS(file = "data/derived/mashr_results_canonical.rds")
} else {
  mashr_results_ED <- read_rds("data/derived/mashr_results_ED.rds")
  mashr_results_canonical <- read_rds("data/derived/mashr_results_canonical.rds")
}


## ----mashr_by_chromosome--------------------------------------------------------------------------------------------------------------------------------------------
mashr_one_chromosome <- function(chr){
  
  data_for_mashr <- get_mashr_data()
  
  focal_data <- data_for_mashr %>% 
    filter(grepl(glue("{chr}_"), SNPs)) %>%
    select(starts_with("beta"), starts_with("SE"))
  
  run_mashr(focal_data, mashr_mode = "canonical") %>%
    saveRDS(file = glue("data/derived/mashr_results_canonical_chr{chr}.rds"))
}

if(!file.exists("data/derived/mashr_results_canonical_chrX.rds")){
  print("Starting the chromosome-specific analysis")
  lapply(c("2L", "2R", "3L", "3R", "X"), mashr_one_chromosome)
} 


## ----get_mixture_assignments, message=FALSE-------------------------------------------------------------------------------------------------------------------------
# Get the mixture weights, as advised by mash authors here: https://github.com/stephenslab/mashr/issues/68
posterior_weights_cov <- mashr_results_canonical$posterior_weights 
colnames(posterior_weights_cov) <- sapply(
  str_split(colnames(posterior_weights_cov), '\\.'), 
  function(x) {
    if(length(x) == 1) return(x)
    else if(length(x) == 2) return(x[1])
    else if(length(x) == 3) return(paste(x[1], x[2], sep = "."))
  })
posterior_weights_cov <- t(rowsum(t(posterior_weights_cov), 
                                  colnames(posterior_weights_cov)))

data_for_mashr <- get_mashr_data()

# Make a neat dataframe
mixture_assignment_probabilities <- data.frame(
  SNP_clump = data_for_mashr$SNPs,
  posterior_weights_cov,
  stringsAsFactors = FALSE
) %>% as_tibble() %>%
  rename(P_equal_effects = equal_effects,
         P_female_specific = Female_specific_1,
         P_male_specific = Male_specific_1,
         P_null = null,
         P_sex_antag = Sex_antag_0.25)


## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------
## all_univariate_lmm_results <- read_csv("data/derived/all_univariate_GEMMA_results.csv") %>%
##   rename_at(vars(-SNPs), ~ str_c(., "_raw"))
## 
## mashr_snps <- data_for_mashr$SNPs
## 
## canonical_estimates <- get_pm(mashr_results_canonical) %>%
##   as_tibble() %>%
##   rename_all(~str_c(., "_mashr_canonical"))
## 
## ED_estimates <- get_pm(mashr_results_ED) %>%
##   as_tibble() %>%
##   rename_all(~str_c(., "_mashr_ED"))
## 
## lfsr_canonical <- get_lfsr(mashr_results_canonical) %>%
##   as_tibble() %>%
##   rename_all(~str_replace_all(., "beta", "LFSR")) %>%
##   rename_all(~str_c(., "_mashr_canonical"))
## 
## lfsr_ED <- get_lfsr(mashr_results_ED) %>%
##   as_tibble() %>%
##   rename_all(~str_replace_all(., "beta", "LFSR")) %>%
##   rename_all(~str_c(., "_mashr_ED"))
## 
## 
## all_mashr_results <- bind_cols(
##   tibble(SNPs = mashr_snps),
##   canonical_estimates,
##   ED_estimates,
##   lfsr_canonical,
##   lfsr_ED)
## 
## all_univariate_lmm_results <- left_join(all_univariate_lmm_results, all_mashr_results, by = "SNPs")
## 
## nested <- all_univariate_lmm_results %>% filter(str_detect(SNPs, ", "))
## split_snps <- strsplit(nested$SNPs, split = ", ")
## nested <- lapply(1:nrow(nested),
##                  function(i) {
##                    data.frame(SNP = split_snps[[i]],
##                               SNP_clump = nested$SNPs[i],
##                               nested[i,] %>% select(-SNPs), stringsAsFactors = FALSE)
##                    }) %>%
##   do.call("rbind", .) %>% as_tibble()
## rm(split_snps)
## 
## all_univariate_lmm_results <- all_univariate_lmm_results %>%
##   filter(!str_detect(SNPs, ", ")) %>%
##   rename(SNP_clump = SNPs) %>% mutate(SNP = SNP_clump) %>%
##   select(SNP, SNP_clump, everything()) %>%
##   bind_rows(nested) %>%
##   arrange(SNP)
## 
## # Merge in the mixture proportions
## all_univariate_lmm_results <-
##   all_univariate_lmm_results %>%
##   left_join(mixture_assignment_probabilities, by = "SNP_clump")
## 
## db <- DBI::dbConnect(RSQLite::SQLite(),
##                      "data/derived/annotations.sqlite3", create = FALSE)
## 
## # tbl(db, "univariate_lmm_results") %>% collect() %>% write_tsv("data/derived/mashr_results_early_2021/original_database_sheet.tsv.gz")
## 
## db %>% db_drop_table(table = "univariate_lmm_results")
## db %>% copy_to(all_univariate_lmm_results,
##                "univariate_lmm_results", temporary = FALSE)


## ----load_mashr_results, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------
mashr_results_canonical <- read_rds("data/derived/mashr_results_canonical.rds")
mashr_results_ED <- read_rds("data/derived/mashr_results_ED.rds")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
tibble(`Mashr version` = c("A. Data-driven covariance matrices",
                           "B. Cononical covariance matrices",
                           "Likelihood ratio (A / B)"),
       `Log likelihood` = c(get_loglik(mashr_results_ED),
                            get_loglik(mashr_results_canonical),
                            get_loglik(mashr_results_ED) / get_loglik(mashr_results_canonical))) %>%
  kable() %>%
  kable_styling()


## ----pairwise_sharing-----------------------------------------------------------------------------------------------------------------------------------------------
get_pairwise_sharing(mashr_results_ED, 
                     factor = 0, 
                     lfsr_thresh = 0.05) %>% 
  round(3)


## ----mashr_check_plot, fig.height=3.1, fig.width=7.3----------------------------------------------------------------------------------------------------------------
db <- DBI::dbConnect(RSQLite::SQLite(), 
                     "data/derived/annotations.sqlite3")

# Results for all 1,613,615 SNPs, even those that are in 100% LD with others (these are grouped up by the SNP_clump column)
all_snps <- tbl(db, "univariate_lmm_results")

# All SNPs and SNP groups that are in <100% LD with one another (n = 1,207,357)
SNP_clumps <- all_snps %>% select(-SNP) %>% distinct() %>% collect(n=Inf)

# Subsetting variable to get the approx-LD subset of SNPs
LD_subset <- !is.na(SNP_clumps$LFSR_female_early_mashr_ED)


hex_plot <- function(x, y, xlab, ylab){
  filter(SNP_clumps,LD_subset) %>% 
    ggplot(aes_string(x, y)) + 
    geom_abline(linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 3) +
    geom_hline(yintercept = 0, linetype = 3) +
    stat_binhex(bins = 200, colour = "#FFFFFF00") + 
    scale_fill_distiller(palette = "Purples") + 
    coord_cartesian(xlim = c(-1,1), ylim = c(-0.55, 0.3)) + 
    theme_minimal() + xlab(xlab) + ylab(ylab) +
    theme(legend.position = "none")
}
grid.arrange(
  hex_plot("beta_female_early_raw", 
           "beta_female_early_mashr_canonical", 
           "Raw estimate of SNP\neffect size from LMM", 
           "Corrected estimate from\nmashr (canonical)"),
  hex_plot("beta_female_early_raw", 
           "beta_female_early_mashr_ED", 
           "Raw estimate of SNP\neffect size from LMM", 
           "Corrected estimate from\nmashr (data-driven)"),
  hex_plot("beta_female_early_mashr_canonical", 
         "beta_female_early_mashr_ED", 
         "Corrected estimate from\nmashr (canonical)", 
         "Corrected estimate from\nmashr (data-driven)"),
  ncol = 3)

