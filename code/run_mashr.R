## ----message=FALSE, warning=FALSE, results="hide"--------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(ashr) # Also requires installation of RMosek, which needs a (free) licence. See the ashr Github page for help
library(mashr) # NB: This has multiple dependencies and was tricky to install. Read the Github page, and good luck!

# library(pander)
# library(reshape2)


## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
# setwd("/Users/lholman/Rprojects/fitnessGWAS")
setwd("/data/projects/punim0243/DGRP_mashr")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
run_mashr <- function(beta_and_se, mashr_mode, ED_p_cutoff = NULL){
  
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
    
    # Nwegatively correlated across ages, and also sex-specific
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

data_for_mashr <- read_csv("data/derived/all_univariate_GEMMA_results.csv") %>% 
  select(starts_with("beta"), starts_with("SE"))

run_mashr(data_for_mashr, mashr_mode = "ED", ED_p_cutoff = 0.2) %>%
  write_rds(path = "data/derived/mashr_results_ED.rds")

run_mashr(data_for_mashr, mashr_mode = "canonical") %>%
  write_rds(path = "data/derived/mashr_results_canonical.rds")


## ----eval = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
## # Get the mixture weights, as advised by mash authors here: https://github.com/stephenslab/mashr/issues/68
## posterior_weights_cov <- mashr_results_canonical$posterior_weights
## colnames(posterior_weights_cov) <- sapply(str_split(colnames(posterior_weights_cov), '\\.'),
##                                                  function(x) {
##                                                    if(length(x)==1) return(x)
##                                                    else if(length(x)==2) return(x[1])
##                                                    else if(length(x)==3) return(paste(x[1], x[2], sep = "."))
##                                                    })
## posterior_weights_cov <- t(rowsum(t(posterior_weights_cov), colnames(posterior_weights_cov)))
## 
## # Make a neat dataframe
## mixture_assignment_probabilities <- data.frame(
##   SNP_clump = read_csv("data/derived/all_univariate_GEMMA_results.csv")$SNPs,
##   posterior_weights_cov,
##   stringsAsFactors = FALSE
## ) %>% as_tibble() %>%
##   rename(P_equal_effects = equal_effects,
##          P_female_specific = Female_specific_1,
##          P_male_specific = Male_specific_1,
##          P_null = null,
##          P_sex_antag = Sex_antag_0.25)


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------
## all_univariate_lmm_results <- read_csv("data/derived/all_univariate_GEMMA_results.csv") %>%
##   rename_at(vars(-SNPs), ~ str_c(., "_raw"))
## mashr_results_canonical <- read_rds("data/derived/mashr_results_canonical.rds")
## mashr_results_ED <- read_rds("data/derived/mashr_results_ED.rds")
## 
## 
## canonical_estimates <- get_pm(mashr_results_canonical) %>% as_tibble() %>% rename_all(~str_c(., "_mashr_canonical"))
## ED_estimates        <- get_pm(mashr_results_ED) %>% as_tibble() %>% rename_all(~str_c(., "_mashr_ED"))
## 
## lfsr_canonical <- get_lfsr(mashr_results_canonical) %>% as_tibble() %>% rename_all(~str_replace_all(., "beta", "LFSR")) %>% rename_all(~str_c(., "_mashr_canonical"))
## lfsr_ED <- get_lfsr(mashr_results_ED) %>% as_tibble() %>% rename_all(~str_replace_all(., "beta", "LFSR")) %>% rename_all(~str_c(., "_mashr_ED"))
## 
## 
## all_univariate_lmm_results <- bind_cols(all_univariate_lmm_results, canonical_estimates, ED_estimates, lfsr_canonical, lfsr_ED)
## 
## nested <- all_univariate_lmm_results %>% filter(str_detect(SNPs, ", "))
## nested <- lapply(1:nrow(nested),
##                  function(i) {
##                    data.frame(SNP = strsplit(nested$SNPs[i], split = ", ")[[1]],
##                               SNP_clump = nested$SNPs[i],
##                               nested[i,] %>% select(-SNPs), stringsAsFactors = FALSE)
##                    }) %>%
##   do.call("rbind", .) %>% as_tibble()
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
## db <- DBI::dbConnect(RSQLite::SQLite(), "data/derived/annotations.sqlite3")
## 
## db %>% db_drop_table(table = "univariate_lmm_results")
## db %>% copy_to(all_univariate_lmm_results,
##                "univariate_lmm_results", temporary = FALSE)


## ----load_mashr_results, echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------
mashr_results_canonical <- read_rds("data/derived/mashr_results_canonical.rds")
mashr_results_ED <- read_rds("data/derived/mashr_results_ED.rds")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tibble(`Mashr version` = c("1. Data-driven covariance matrices",
                           "2. Cononical covariance matrices",
                           "Likelihood ratio (1 / 2)"),
       `Log likelihood` = c(get_loglik(mashr_results_ED),
                            get_loglik(mashr_results_canonical),
                            get_loglik(mashr_results_ED) / get_loglik(mashr_results_canonical))) %>%
  pander()


## ----plot_mixtrure_props---------------------------------------------------------------------------------------------------------------------------------------------------------
melt(sort(get_estimated_pi(mashr_results_canonical))) %>%
  rownames_to_column("Mixture_component") %>%
  filter(value > 0.01) %>% 
  spread(Mixture_component, value) %>%
  rename(`Equal effects on\neach fitness component` = equal_effects,
         `Female-specific effect` = Female_specific_1,
         `Male-specific effect` = Male_specific_1,
         `Sexually-antagonistic effect` = Sex_antag_0.25) %>%
  gather() %>%
  arrange(-value) %>%
  mutate(key = factor(key, rev(key))) %>%
  ggplot(aes(key, 100*value)) + 
  geom_bar(stat = "identity", fill = "#ff6f61", colour = "grey10") + 
  coord_flip() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 67)) + 
  theme_bw() + 
  theme(axis.ticks.y = element_blank()) +
  ylab("Estimated % SNPs") +
  xlab(NULL)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
melt(sort(get_estimated_pi(mashr_results_ED))) %>%
  rownames_to_column("Mixture_component") %>%
  filter(value > 0.01) %>%
  arrange(-value) %>%
  mutate(Mixture_component = factor(Mixture_component, rev(Mixture_component))) %>%
  ggplot(aes(Mixture_component, 100*value)) + 
  geom_bar(stat = "identity", fill = "#ff6f61", colour = "grey10") + 
  coord_flip() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 78)) + 
  theme_bw() + 
  theme(axis.ticks.y = element_blank()) +
  ylab("Estimated % SNPs") +
  xlab(NULL)

mashr_results_ED$fitted_g$Ulist[["ED_tPCA"]] %>% cov2cor()
mashr_results_ED$fitted_g$Ulist[["ED_PCA_1"]] %>% cov2cor()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
get_pairwise_sharing(mashr_results_ED, factor = 0, lfsr_thresh = 0.05) %>% round(3)

