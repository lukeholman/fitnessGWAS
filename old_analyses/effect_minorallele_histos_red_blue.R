

Evolutionary theory makes various testable predictions about the effects of the major and minor alleles on phenotypes, especially phenotypes that are closely correlated with fitness. For example, the minor allele should generally be associated with the lower-fitness phenotype (and the major allele with the higher-fitness phenotype), due to selection removing the lower-fitness allele. Furthermore, it is interesting to ask whether the female-beneficial allele or the male-beneficial allele tends to be the major allele at sexually-antagonistic loci? Lastly, we might expect different results at loci where the minor allele is quite rare vs loci where the minor allele is only slightly rarer than the major allele; loci with comparatively low MAF are more likely to reflect new polymorphisms (since most mutations are neutral or deleterious, and it takes time for the minor allele to drift towards high frequencies), and are more likely to be under selection (since selection against the minor allele keeps it rare).

To address these questions graphically, we here plot the effects on male and female early-life fitness of each allele, that was assigned to the female-specific, male-specific, or sexually antagonistic mixture component by `mashr`, with an assignment probability in the highest 0.1% across all the loci analysed. Since there were `r n_loci` loci, each column plots the male and female effect sizes for `r fraction_n` loci.

Among the top 0.1% of loci with female-specific effects on fitness (first column), the minor allele was usually associated with reduced female fitness among loci where the minor allele frequency was below 0.2, as expected if alleles that harm female fitness are removed by selection. Among loci with MAF > 0.2, the minor allele beneficial to female fitness almost as often as the major allele was.

There was no similar finding for the top 0.1% of loci with loci with _male_-specific effects (second column): the minor allele was equally likely to increase or reduce fitness, regardless of MAF. Also, most of these top male-specific loci had quite high minor allele frequencies, perhaps suggesting their effects on fitness were weaker than for the loci with female-specific fitness effects.

For the top 0.1% of loci with sexually-antagonistic effects on fitness (third column), the minor allele frequencies again tended to be high (above 0.2), suggesting either weak overall purifying selection (due to counteracting effects of selection on males and females), or even balancing selection. The minor allele tended to be the one that reduced female fitness and elevated male fitness, though alleles with the reverse effect were present as well.

Finally, among loci inferred to be under sexually concordant selection, the effect sizes were strongest among loci with MAF < 0.2, as expected if selection prevents alleles with strong, detrimental effects on the fitness of both sexes from becoming common.


```{r fig.width = 8, fig.height = 5}
histo_data <- tbl(db, "univariate_lmm_results") %>%
  left_join(tbl(db, "variants") %>%
              select(SNP, MAF, site.class),
            by = "SNP") %>%
  filter(!is.na(LFSR_female_early_mashr_ED)) %>%
  collect(n = Inf) %>%
  mutate(class =
           case_when(
             P_sex_antag > quantile(P_sex_antag, probs = 0.999)  ~ "Top 0.1%\nsexually antagonistic loci",
             P_equal_effects > quantile(P_equal_effects, probs = 0.999)  ~ "Top 0.1%\nsexually concordant loci",
             P_female_specific > quantile(P_female_specific, probs = 0.999)  ~ "Top 0.1%\nfemale-specific loci",
             P_male_specific > quantile(P_male_specific, probs = 0.999)  ~ "Top 0.1%\nmale-specific loci",
           )) %>%
  mutate(MAF = ifelse(MAF < 0.2, "0.05 < MAF < 0.2", "0.2 < MAF < 0.5")) %>%
  filter(!is.na(class)) %>%
  select(beta_female_early_mashr_ED, beta_male_early_mashr_ED, class, MAF)

eff_size_histos <- ggplot(histo_data, aes(beta_female_early_mashr_ED)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(bins=30, fill="tomato", colour=NA, alpha = 0.5) +
  geom_histogram(aes(beta_male_early_mashr_ED), bins=30, fill="steelblue", colour=NA, alpha = 0.5) +
  geom_histogram(bins=30, fill=NA, colour="black") +
  geom_histogram(aes(beta_male_early_mashr_ED), bins=30,fill=NA,colour="black") +
  facet_grid(MAF ~ class, scales = "free_y") +
  xlab("Effect of the minor allele on early-life fitness\n(females:red, males:blue)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        text = element_text(family = "Raleway", size = 12)) +
  ylab("Number of loci")

ggsave("figures/eff_size_histos.pdf", width = 8, height = 5)
eff_size_histos
```
