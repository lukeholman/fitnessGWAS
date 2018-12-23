## ----message=FALSE, warning=FALSE, results="hide"------------------------
library(dplyr)
library(ashr) # Also requires installation of RMosek, which needs a (free) licence. See the ashr Github page for help
library(mashr) # NB: This has multiple dependencies and was tricky to install. Read the Github page, and good luck!
library(readr)
library(kableExtra)

kable_table <- function(df) {
  kable(df, "html") %>%
  kable_styling() %>%
  scroll_box(height = "300px")
}

# launch this from Terminal with: env R_MAX_VSIZE=700Gb R 
setwd("/Users/lholman/Rprojects/fitnessGWAS")

## ------------------------------------------------------------------------
load_gwas_results <- function(fitness.component){
  new_names <- c(x = "beta", y = "se", z = "p_wald")
  names(new_names) <- c(paste("beta_", fitness.component, sep = ""), 
                        paste("SE_", fitness.component, sep = ""), 
                        paste("pvalue_", fitness.component, sep = ""))
  
  paste("data/derived/output/", fitness.component, "_lmm.assoc.txt", sep = "") %>%
    read_tsv() %>%
    select(SNP, beta, se, p_wald) %>%
    rename(!!new_names)
}
all_univariate_lmm_results <- load_gwas_results("female_early") %>%
  full_join(load_gwas_results("female_late"), by = "SNP") %>%
  full_join(load_gwas_results("male_early"), by = "SNP") %>%
  full_join(load_gwas_results("male_late"), by = "SNP") 

print(paste("Number of SNPs that were analysed with GEMMA:", nrow(all_univariate_lmm_results)))

## ------------------------------------------------------------------------
all_univariate_lmm_results <- all_univariate_lmm_results %>%
  group_by(paste(beta_female_early, beta_female_late, beta_male_early, beta_male_late, # Group loci with identical GWAS results (these are in 100% LD with each other)
                 SE_female_early, SE_female_late, SE_male_early, SE_male_late)) %>%
  summarise(SNPs = paste0(SNP, collapse = ", "),    # If there are multiple SNPs, concatenate their names
            beta_female_early = unique(beta_female_early), 
            beta_female_late = unique(beta_female_late), 
            beta_male_early = unique(beta_male_early), 
            beta_male_late = unique(beta_male_late),
            SE_female_early = unique(SE_female_early), 
            SE_female_late = unique(SE_female_late), 
            SE_male_early = unique(SE_male_early), 
            SE_male_late = unique(SE_male_late),
            pvalue_female_early = unique(pvalue_female_early), 
            pvalue_female_late = unique(pvalue_female_late), 
            pvalue_male_early = unique(pvalue_male_early), 
            pvalue_male_late = unique(pvalue_male_late)) %>%
  ungroup() %>% select(-starts_with("paste")) %>% 
  arrange(pvalue_female_early + pvalue_female_late + pvalue_male_early + pvalue_male_late) # For the following table, arrange from most to least significant

## ------------------------------------------------------------------------
all_univariate_lmm_results %>% head(50) %>% kable_table()
print(paste("Number of SNPs (or SNP groups that are in <100% linkage disequilibrium with each other):", nrow(all_univariate_lmm_results)))

## ------------------------------------------------------------------------
run_mashr <- function(beta_and_se, mashr_mode){
  
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
    strong <- get_significant_results(m.1by1, 0.2)   
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

# Reorder by SNP name, and save the neatened results from the univariate GEMMA (also needed to interpret the mashr results later)
all_univariate_lmm_results <- all_univariate_lmm_results %>% arrange(SNPs)
write_csv(all_univariate_lmm_results, path = "data/derived/all_univariate_GEMMA_results.csv")
data_for_mashr <- all_univariate_lmm_results %>% select(starts_with("beta"), starts_with("SE"))
rm(all_univariate_lmm_results)

mashr_results_canonical <- run_mashr(data_for_mashr, mashr_mode = "canonical")
mashr_results_ED        <- run_mashr(data_for_mashr, mashr_mode = "ED")

write_rds(mashr_results_canonical, path = "data/derived/mashr_results_canonical.rds")
write_rds(mashr_results_ED, path = "data/derived/mashr_results_ED.rds")

## ------------------------------------------------------------------------
all_univariate_lmm_results <- read_csv("data/derived/all_univariate_GEMMA_results.csv") %>% 
  rename_at(vars(-SNPs), funs(str_c(., "_raw")))
mashr_results_canonical <- read_rds("data/derived/mashr_results_canonical.rds")
mashr_results_ED <- read_rds("data/derived/mashr_results_ED.rds")


canonical_estimates <- get_pm(mashr_results_canonical) %>% as_tibble() %>% rename_all(~str_c(., "_mashr_canonical"))
ED_estimates        <- get_pm(mashr_results_ED) %>% as_tibble() %>% rename_all(~str_c(., "_mashr_ED"))

lfsr_canonical <- get_lfsr(mashr_results_canonical) %>% as_tibble() %>% rename_all(~str_replace_all(., "beta", "LFSR")) %>% rename_all(~str_c(., "_mashr_canonical"))
lfsr_ED <- get_lfsr(mashr_results_ED) %>% as_tibble() %>% rename_all(~str_replace_all(., "beta", "LFSR")) %>% rename_all(~str_c(., "_mashr_ED"))


all_univariate_lmm_results <- bind_cols(all_univariate_lmm_results, canonical_estimates, ED_estimates, lfsr_canonical, lfsr_ED)

all_univariate_lmm_results %>% ggplot(aes(beta_female_early_raw, beta_female_early_mashr_canonical)) + geom_abline() + stat_binhex() + coord_fixed(xlim = c(-1,1), ylim = c(-1,1))
all_univariate_lmm_results %>% ggplot(aes(beta_female_early_raw, beta_female_early_mashr_ED)) + geom_abline() + stat_binhex() + coord_fixed(xlim = c(-1,1), ylim = c(-1,1))
all_univariate_lmm_results %>% ggplot(aes(beta_female_early_mashr_canonical, beta_female_early_mashr_ED)) + geom_abline() + stat_binhex() + coord_fixed(xlim = c(-1,1), ylim = c(-1,1))

all_univariate_lmm_results %>% ggplot(aes(-log10(pvalue_female_early_raw), -log10(LFSR_female_early_mashr_canonical))) + geom_abline() + stat_binhex() + coord_fixed(xlim = c(0, 8), ylim = c(0, 8))
all_univariate_lmm_results %>% ggplot(aes(-log10(pvalue_female_early_raw), -log10(LFSR_female_early_mashr_ED))) + geom_abline() + stat_binhex() + coord_fixed(xlim = c(0, 8), ylim = c(0, 8))
all_univariate_lmm_results %>% ggplot(aes(-log10(LFSR_female_early_mashr_canonical), -log10(LFSR_female_early_mashr_ED))) + geom_abline() + stat_binhex() + coord_fixed(xlim = c(0, 8), ylim = c(0, 8))

all_univariate_lmm_results %>%
  filter_at(vars(contains("pvalue")), any_vars(. < 10e-5)) %>%
  select(contains("LFSR") & contains("canonical"))

# Calcualte antagonism....
all_univariate_lmm_results %>%
  mutate(pp_female_early = ifelse(beta_female_early_mashr_ED > 0, LFSR_female_early_mashr_ED, 1 - LFSR_female_early_mashr_ED),
         pp_female_late  = ifelse(beta_female_late_mashr_ED > 0, LFSR_female_late_mashr_ED, 1 - LFSR_female_late_mashr_ED),
         pp_male_early   = ifelse(beta_male_early_mashr_ED > 0, LFSR_male_early_mashr_ED, 1 - LFSR_male_early_mashr_ED),
         pp_male_late    = ifelse(beta_male_late_mashr_ED > 0, LFSR_male_late_mashr_ED, 1 - LFSR_male_late_mashr_ED)) %>%
  select(starts_with("pp"))


## ------------------------------------------------------------------------
# Make a nice merged results file, and add it to the database...

