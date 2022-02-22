
## Plot the estimated mixture proportions from `mashr`

```{r mixture_props, fig.showtext = TRUE, fig.height = 8, fig.width = 5, message=FALSE, warning=FALSE}
# Make plot for the GWAS
mashr_results_canonical <- read_rds("data/derived/mashr_results_canonical.rds")
mashr_2L <- readRDS("data/derived/mashr_results_canonical_chr2L.rds")
mashr_2R <- readRDS("data/derived/mashr_results_canonical_chr2R.rds")
mashr_3L <- readRDS("data/derived/mashr_results_canonical_chr3L.rds")
mashr_3R <- readRDS("data/derived/mashr_results_canonical_chr3R.rds")
mashr_X <- readRDS("data/derived/mashr_results_canonical_chrX.rds")

mix <- bind_rows(
  enframe(sort(get_estimated_pi(mashr_results_canonical))) %>%
    mutate(Chromosome = "All"),
  enframe(sort(get_estimated_pi(mashr_2L))) %>%
    mutate(Chromosome = "2L"),
  enframe(sort(get_estimated_pi(mashr_2R))) %>%
    mutate(Chromosome = "2R"),
  enframe(sort(get_estimated_pi(mashr_3L))) %>%
    mutate(Chromosome = "3L"),
  enframe(sort(get_estimated_pi(mashr_3R))) %>%
    mutate(Chromosome = "3R"),
  enframe(sort(get_estimated_pi(mashr_X))) %>%
    mutate(Chromosome = "X")) %>%
  rename(Mixture_component = name)

to_keep <- mix %>%
  group_by(Mixture_component) %>%
  summarise(value = max(value), .groups = "drop") %>%
  filter(value > 0.01) %>%
  pull(Mixture_component)

mix <- mix %>%
  filter(Mixture_component %in% to_keep) %>%
  spread(Mixture_component, value) %>%
  rename(`Sexually concordant effect` = equal_effects,
         `Female-specific effect` = Female_specific_1,
         `Male-specific effect` = Male_specific_1,
         `Sexually antagonistic effect` = Sex_antag_0.25,
         `No effect on fitness` = null) %>%
  gather(Mixture_component, value, -Chromosome) %>%
  arrange(-value)

chr_levels <- mix %>%
  filter(Mixture_component == "Sexually antagonistic effect") %>%
  arrange(value) %>% pull(Chromosome)

mix <- mix %>%
  mutate(Chromosome = factor(Chromosome, chr_levels),
         Mixture_component = factor(Mixture_component,
                                    c("Sexually antagonistic effect",
                                      "Sexually concordant effect",
                                      "Female-specific effect",
                                      "Male-specific effect",
                                      "No effect on fitness")))


p1 <- mix %>%
  ggplot(aes(Chromosome, 100 * value)) +
  geom_bar(stat = "identity",aes(fill = Chromosome), colour = "grey10",  size = 0.3) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 70)) +
  scale_x_discrete(expand = c(0.14, 0.14)) +
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 0.8),
        text = element_text(family = "Raleway", size = 12),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank()) +
  ylab("Estimated % of loci") +
  facet_wrap(~ Mixture_component, ncol = 1)


# Make plot for the TWAS
mashr_results_canonical <- readRDS("data/derived/TWAS/TWAS_canonical.rds")
mashr_2L <- readRDS("data/derived/TWAS/TWAS_canonical_2L.rds")
mashr_2R <- readRDS("data/derived/TWAS/TWAS_canonical_2R.rds")
mashr_3L <- readRDS("data/derived/TWAS/TWAS_canonical_3L.rds")
mashr_3R <- readRDS("data/derived/TWAS/TWAS_canonical_3R.rds")
mashr_X <- readRDS("data/derived/TWAS/TWAS_canonical_X.rds")

mix <- bind_rows(
  enframe(sort(get_estimated_pi(mashr_results_canonical))) %>%
    mutate(Chromosome = "All"),
  enframe(sort(get_estimated_pi(mashr_2L))) %>%
    mutate(Chromosome = "2L"),
  enframe(sort(get_estimated_pi(mashr_2R))) %>%
    mutate(Chromosome = "2R"),
  enframe(sort(get_estimated_pi(mashr_3L))) %>%
    mutate(Chromosome = "3L"),
  enframe(sort(get_estimated_pi(mashr_3R))) %>%
    mutate(Chromosome = "3R"),
  enframe(sort(get_estimated_pi(mashr_X))) %>%
    mutate(Chromosome = "X")) %>%
  rename(Mixture_component = name) %>%
  mutate(Mixture_component = str_remove_all(Mixture_component, "_0.25"),
         Mixture_component = str_remove_all(Mixture_component, "_0.5"),
         Mixture_component = str_remove_all(Mixture_component, "_0.75"),
         Mixture_component = str_remove_all(Mixture_component, "_1.0")) %>%
  group_by(Mixture_component, Chromosome) %>%
  summarise(value = sum(value), .groups = "drop")

to_keep <- mix %>%
  group_by(Mixture_component) %>%
  summarise(value = max(value), .groups = "drop") %>%
  filter(value > 0.01) %>%
  pull(Mixture_component)

mix <- mix %>%
  filter(Mixture_component %in% to_keep) %>%
  spread(Mixture_component, value) %>%
  rename(`Sexually concordant effect` = equal_effects,
         `Female-specific effect` = Female_specific_1,
         `Male-specific effect` = Male_specific_1,
         `Sexually antagonistic effect` = Sex_antag,
         `No effect on fitness` = null) %>%
  gather(Mixture_component, value, -Chromosome) %>%
  arrange(-value)

chr_levels <- rev(c("X", "3L", "3R", "All", "2R", "2L"))

mix <- mix %>%
  mutate(Chromosome = factor(Chromosome, chr_levels),
         Mixture_component = factor(Mixture_component,
                                    c("Sexually antagonistic effect",
                                      "Sexually concordant effect",
                                      "Female-specific effect",
                                      "Male-specific effect",
                                      "No effect on fitness")))

p2 <- mix %>%
  ggplot(aes(Chromosome, 100 * value)) +
  geom_bar(stat = "identity",aes(fill = Chromosome), colour = "grey10",  size = 0.3) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 65)) + #
  scale_x_discrete(expand = c(0.14, 0.14)) +
  theme(axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 0.8),
        text = element_text(family = "Raleway", size = 12),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank()) +
  ylab("Estimated % of transcripts") +
  xlab(" ") +
  facet_wrap(~ Mixture_component, ncol = 1)


# Save composite figure of the GWAS and TWAS mixture proportions
ggsave(filename = "figures/composite_mixture_figure.pdf",
       cowplot::plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12),
       width = 5, height = 8)

cowplot::plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
```


**Figure X**: Proportions of each type of locus, as estimated using the mixture model computed by `mashr` in the analysis using canonical covariance matrices. This analysis involved a number of pre-specified covariance matrices, each corresponding to a type of locus that we hypothesised to exist (shown in panel headings). The analysis fit some other matrix types not shown here, because the corresponding locus type was inferred to be rare/absent (these included neutral loci, estimated to comprise 0.5-1% of those tested, and age-antagonistic loci, none of which were detected). The analysis was run either using all 1,207,357 loci for which data were available (labeled 'All') or for all loci on each of the chromosomes (chromosomes 4 and Y had insufficient data). Notably, sexually-antagonistic loci were inferred to be especially common on the X chromosome, loci that affected males only were inferred to be rarer than those affecting females only, and chromosome 2R had many more female-specific than male-specific loci.
