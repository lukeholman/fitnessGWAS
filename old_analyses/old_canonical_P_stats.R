## Statistical models: effects of site class, MAF and chromosome

The following models test whether site class, MAF (minor allele frequency) and chromosome predict how each locus affected our four phenotypes.

These analyses use the mixture assignment probabilities calculated by `mashr` (in the "canonical" analysis). For every locus, the `mashr` analysis calculated the probability that each locus was A) sexually antagonistic (meaning it had opposite effects on male and female fitness that were consistent across the two age categories), B) female-specific (meaning it affect female fitness but not male), C) male-specific (affecting male fitness only), and D) had equal effects on all four phenotypes (mean the allele that associated with higher fitness in males also was associated with higher fitness in females). Other types of loci were also considered in the model (e.g. null effect loci and ones with age-specific effects), but these were found to be quite rare, precluding analysis here.

Here, we use these mixture assignment probabilities as the response variable, in order to ask whether e.g. sexually antagonistic loci tend to appear in any particular type of site or chromosome, and whether they had a higher or lower MAF than average. Since these probabilities are bounded between zero and one, we use beta regression (a form of GLM designed to model variables in the range [0,1]). We fit a full model with site class, MAF and chromosome as predictor variables, as well as all the simpler nested models, and compare them using AICc model selection. We also present the results of the full model, and plot the parameter estimates from each of the four models below.

The dataset used in this model is

### Tables of statistical results (models of GWAS effects)

```{r warning=FALSE, message=FALSE}
library(betareg)
library(MuMIn)
library("lmtest")

LFSR_cutoff <- 0.05

dat <- univariate_lmm_results %>%
  filter(LFSR_female_early < LFSR_cutoff | LFSR_female_late < LFSR_cutoff |
           LFSR_male_early < LFSR_cutoff | LFSR_male_late < LFSR_cutoff) %>%
  mutate(chr = gsub("_", "", substr(SNP, 1, 2))) %>%
  select(SNP, SNP_clump, starts_with("P_"), MAF, site.class, chr) %>%
  distinct()

# Focus only on the commonest site classes:
dat <- dat %>%
  filter(site.class %in% c("Intron", "Intergenic", "Downstream", "Upstream",
                           "Synonymous coding", "Non-synonymous coding", "UTR 3-prime",
                           "Exon", "UTR 5-prime"))

# Remove chromosome 4 (too few sites)
dat <- dat %>% filter(chr != "4")

# If there are multiple site classes for a SNP, or multiple SNPs
# in the same 100% LD clump, pick one SNP and/or 1 site class at random
set.seed(1)
dat <- dat[sample(nrow(dat),nrow(dat)), ] %>%
  split(.$SNP_clump) %>%
  map_df(~ .x[1, ])

dat <- dat %>% arrange(site.class)
dat$site.class <- relevel(factor(dat$site.class), ref = "Synonymous coding")

n_loci <- prettyNum(nrow(dat), big.mark = ",", scientific = FALSE)

compare_mods <- function(MAF_and_Chromosome_and_site_class){

  MAF_and_Chromosome <- update(MAF_and_Chromosome_and_site_class, ~ . -site.class)
  MAF_and_siteclass <- update(MAF_and_Chromosome_and_site_class, ~ . -chr)
  Chromosome_and_siteclass <- update(MAF_and_Chromosome_and_site_class, ~ . -chr)

  MAF <- update(MAF_and_Chromosome_and_site_class, ~ . -site.class - chr)
  Chromosome <- update(MAF_and_Chromosome_and_site_class, ~ . -site.class - MAF)
  siteclass <- update(MAF_and_Chromosome_and_site_class, ~ .  - MAF -chr)
  Null_model <- update(MAF_and_Chromosome_and_site_class, ~ . -site.class - MAF -chr)

  AICc(MAF_and_Chromosome_and_site_class,
       MAF_and_Chromosome,
       MAF_and_siteclass,
       Chromosome_and_siteclass,
       MAF, Chromosome, siteclass,
       Null_model) %>%
    rownames_to_column("Model") %>%
    arrange(AICc) %>%
    mutate(delta = AICc - AICc[1],
           Model = str_replace_all(Model, "_and+_", " + "),
           Weight = round(exp(-0.5 * delta) / sum(exp(-0.5 * delta)), 3),
           delta = round(delta, 2)) %>%
    kable(digits = 2) %>% kable_styling(full_width = FALSE)
}
```


In each of the following analyses, the response variable is the mixture assignment probability to type $i$ for each of the `r n_loci` loci that affected at least one of the phenotypes significantly (defined as LFSR < `r LFSR_cutoff`), where $i$ is one of the four variant types shown in the above figure. In cases where multiple SNP or indel loci were in 100% linkage with one another, we picked a single locus at random and discarded the others.

Three predictors were available for each locus: the minor allele frequency (MAF) in the overall DGRP, the site class of the variant, and chromosome. Loci that were annotated with more than one site class (e.g. intron as well as exon, due to an overlap of genes) were assigned one of these site classes at random.

To evaluate the effects of three predictors on the response variable, we fit 8 nested models and compared them using AICc. We also present the `summary()` output for the top-ranked model according to AICc.

#### Probability of sexual antagonism {.tabset}

##### AICc table
```{r message=FALSE, warning=FALSE}
betareg(P_sex_antag ~ MAF + chr + site.class,
        data = dat) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
p1 <- summary(betareg(P_sex_antag ~ MAF + chr + site.class, data = dat), type = "deviance")
p1
```


#### Probability of being female-specific {.tabset}

##### AICc table
```{r message=FALSE, warning=FALSE}
betareg(P_female_specific ~ MAF + chr + site.class,
        data = dat) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
p2 <- summary(betareg(P_female_specific ~ MAF + chr + site.class, data = dat), type = "deviance")
p2
```


#### Probability of being male-specific {.tabset}

##### AICc table
```{r message=FALSE, warning=FALSE}
betareg(P_male_specific ~ MAF + chr + site.class,
        data = dat) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
p3 <- summary(betareg(P_male_specific ~ MAF + chr + site.class, data = dat), type = "deviance")
p3
```


#### Probability of equal effects {.tabset}

##### AICc table
```{r message=FALSE, warning=FALSE}
betareg(P_equal_effects ~ MAF + chr + site.class,
        data = dat) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
p4 <- summary(betareg(P_equal_effects ~ MAF + chr + site.class, data = dat), type = "deviance")
p4
```

### Tables of statistical results (models of TWAS effects)

```{r message=FALSE}
twas <- readRDS("data/derived/TWAS/TWAS_mixture_assignment_probabilities.rds") %>%
  left_join(tbl(db, "genes") %>% select(FBID, chromosome) %>% collect(), by = "FBID") %>%
  left_join(read_csv("data/derived/gene_expression_by_sex.csv"), by = "FBID") %>%
  filter(chromosome %in% c("2L", "2R", "3L", "3R", "X")) %>%
  mutate(h2 = (female_narrow_heritability + male_narrow_heritability) / 2)


compare_mods <- function(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability){

  Chromosome_and_Expression_level_and_Heritability <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                                             ~ . -male_bias_in_expression)
  Sex_bias_and_Chromosome_and_Expression_level <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                                         ~ . -h2)
  Sex_bias_and_Expression_level_and_Heritability <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                                           ~ . -chromosome)
  Sex_bias_and_Chromosome_and_Heritability <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                                     ~ . -AveExpr)

  Chromosome_and_Expression_level <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability, ~ . -male_bias_in_expression -h2)
  Sex_bias_and_Expression_level <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability, ~ . -chromosome -h2)
  Sex_bias_and_Chromosome <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability, ~ . -AveExpr -h2)

  Heritability_and_Expression_level <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                              ~ . -male_bias_in_expression -chromosome)
  Heritability_and_Sex_bias <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                      ~ . -chromosome -AveExpr)
  Heritability_and_Chromosome <- update(Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
                                        ~ . -male_bias_in_expression -AveExpr)

  Chromosome <- update(Chromosome_and_Expression_level, ~ . -AveExpr)
  Sex_bias <- update(Sex_bias_and_Expression_level, ~ . -AveExpr)
  Expression_level <- update(Sex_bias_and_Expression_level, ~ .  - male_bias_in_expression)
  Heritability <- update(Heritability_and_Sex_bias, ~ .  - male_bias_in_expression)
  Null_model <- update(Chromosome, ~ . -chromosome)

  AICc(
    Sex_bias_and_Chromosome_and_Expression_level_and_Heritability,
    Chromosome_and_Expression_level_and_Heritability,
    Sex_bias_and_Chromosome_and_Expression_level,
    Sex_bias_and_Expression_level_and_Heritability,
    Sex_bias_and_Chromosome_and_Heritability,
    Chromosome_and_Expression_level,
    Sex_bias_and_Expression_level,
    Sex_bias_and_Chromosome,
    Heritability_and_Expression_level,
    Heritability_and_Sex_bias,
    Heritability_and_Chromosome,
    Chromosome, Sex_bias, Expression_level, Heritability,
    Null_model) %>%
    rownames_to_column("Model") %>%
    arrange(AICc) %>%
    mutate(delta = AICc - AICc[1],
           Model = str_replace_all(Model, "_and+_", " + "),
           Weight = round(exp(-0.5 * delta) / sum(exp(-0.5 * delta)), 3),
           delta = round(delta, 2))
}
```


#### Probability of sexual antagonism {.tabset}

##### AICc table
```{r}
betareg(P_sex_antag ~ male_bias_in_expression + chromosome + AveExpr + h2,
        data = twas) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
q1 <- summary(betareg(P_sex_antag ~ male_bias_in_expression + AveExpr + chromosome + h2, data = twas), type = "deviance")
q1
```

#### Probability of being female-specific {.tabset}

##### AICc table
```{r}
betareg(P_female_specific ~ male_bias_in_expression + chromosome + AveExpr + h2,
        data = twas) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
q2 <- summary(betareg(P_female_specific ~ male_bias_in_expression + AveExpr + chromosome + h2, data = twas), type = "deviance")
q2
```

#### Probability of being female-specific {.tabset}

##### AICc table
```{r}
betareg(P_male_specific ~ male_bias_in_expression + chromosome + AveExpr + h2,
        data = twas) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
q3 <- summary(betareg(P_male_specific ~ male_bias_in_expression + AveExpr + chromosome + h2, data = twas), type = "deviance")
q3
```

#### Probability of equal effects {.tabset}

##### AICc table
```{r}
betareg(P_equal_effects ~ male_bias_in_expression + chromosome + AveExpr + chromosome + h2,
        data = twas) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
q4 <- summary(betareg(P_equal_effects ~ male_bias_in_expression + AveExpr + chromosome + h2, data = twas), type = "deviance")
q4
```

#### Probability of null effect {.tabset}

##### AICc table
```{r}
betareg(P_null ~ male_bias_in_expression + chromosome + AveExpr + chromosome + h2,
        data = twas) %>% compare_mods()
```

##### Full model
```{r message=FALSE, warning=FALSE}
q5 <- summary(betareg(P_null ~ male_bias_in_expression + AveExpr + chromosome + h2, data = twas), type = "deviance")
q5
```

### Plots showing the statistical results

```{r make_stats_figs}
gwas_stats_fig <- bind_rows(p1$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(sexually antagonistic)"),
                            p2$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(female-specific)"),
                            p3$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(male-specific)"),
                            p4$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(equal effects)")) %>%
  mutate(Parameter = str_replace_all(Parameter, "site[.]class", "Site class: "),
         Parameter = str_replace_all(Parameter, "chr", "Chromosome: "),
         sig = ifelse(`Pr(>|z|)` < 0.05, "yes", "no")) %>%
  filter(Parameter != "(Intercept)") %>%
  arrange(Parameter) %>%
  mutate(Parameter = factor(Parameter, rev(unique(Parameter)))) %>%
  ggplot(aes(Parameter, Estimate, colour = sig)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`*1.96, ymax = Estimate + `Std. Error`*1.96), width = 0) +
  geom_point() +
  coord_flip() +
  facet_wrap(~type) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        text = element_text(family = "Raleway", size = 12)) +
  scale_color_manual(values = c("black", "tomato")) +
  ylab("Estimate \u00B1 95% CIs")

twas_stats_fig <- bind_rows(q1$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(sexually antagonistic)"),
                            q2$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(female-specific)"),
                            q3$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(male-specific)"),
                            q4$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(equal effects)"),
                            q5$coefficients$mean %>%
                              as.data.frame() %>%
                              rownames_to_column("Parameter") %>%
                              mutate(type = "P(null effect)")) %>%
  mutate(Parameter = str_replace_all(Parameter, "AveExpr", "Average expression level"),
         Parameter = str_replace_all(Parameter, "male_bias_in_expression", "Male bias in expression"),
         Parameter = str_replace_all(Parameter, "chromosome", "Chromosome: "),
         sig = ifelse(`Pr(>|z|)` < 0.05, "yes", "no"),
         type = factor(type, c("P(null effect)", "P(equal effects)", "P(sexually antagonistic)",
                               "P(female-specific)", "P(male-specific)"))) %>%
  filter(Parameter != "(Intercept)") %>%
  arrange(Parameter) %>%
  mutate(Parameter = factor(Parameter, rev(unique(Parameter)))) %>%
  ggplot(aes(Parameter, Estimate, colour = sig)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`*1.96, ymax = Estimate + `Std. Error`*1.96), width = 0) +
  geom_point() +
  coord_flip() +
  facet_wrap(~type) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        text = element_text(family = "Raleway", size = 12)) +
  scale_color_manual(values = c("black", "tomato")) +
  ylab("Estimate \u00B1 95% CIs")


ggsave("figures/GWAS_stats_figure.pdf", gwas_stats_fig, width = 6, height = 6)
ggsave("figures/TWAS_stats_figure.pdf", twas_stats_fig, width = 7, height = 6)
```

#### Effect sizes for the GWAS analysis
```{r fig.showtext=T, fig.height=6, fig.width=6}
gwas_stats_fig
```

**Figure XX**: The figure shows estimated parameters from four separate beta regression models seeking predictors of the mixture assignment probabilities calculated using `mashr` on the GWAS results. Positive values indicate that the predictor is associated with an elevated probability across loci; estimates that are non-zero with 95% confidence are shown in red. The clearest result is the strong, positive relationship between minor allele frequency (MAF) and the probability that the locus was assigned to the 'sexually antagonistic' mixture component. There was also a strong, negative relationship between MAF and the chance that the locus was assigned to the 'equal effects' mixture component (i.e. that the locus had concordant effects on male and female fitness). Loci with a relatively high probability of being female-specific tended to have _lower_ MAF, while loci with a relatively high probability of being male-specific tended to have _higher_ MAF. Chromosome 2L seems to be enriched for loci with male-specific effects, and de-enriched for loci that affect females relative to the other chromosomes. Finally, there were some associations between site class and mixture assignment probability, which (though sometimes statistically significant) were small and subtle.

#### Effect sizes for the TWAS analysis

```{r fig.showtext=T, fig.height=6, fig.width=6}
twas_stats_fig
```

**Figure XX**: The figure shows estimated parameters from five separate beta regression models seeking predictors of the mixture assignment probabilities calculated using `mashr` on the TWAS results. Positive values indicate that the predictor is associated with an elevated probability across transcripts; estimates that are non-zero with 95% confidence are shown in red. The clearest result is the strong, negative relationship between the heritability of the transcript's expression level (calculated by Huang et al. 2015) and the probability that the transcript was assigned to any mixture component other than 'sexually antagonistic'. Transcripts with strongly male-biased expression were more likely to have no relationship with fitness (i.e. to be assigned to the 'null' mixture component); sex bias was a significant but weak predictor for the other mixture components. Transcripts with high expression levels were more likely to have a sexually concordant relationship with fitness, and less likely to be sexually antagonistic. Transcripts from genes located on the X chromosome were slightly more likely to have a female-specific correlation with fitness, or no correlation.
