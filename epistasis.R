# Have a look at the epistasis data to see which models differ the most

# data
d_epi_means <- read.table(paste0(DATA_PATH, "d_epi_mean.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "minEW", "maxEW", "q025EW",
                                        "q25EW", "q50EW", "q75EW", "q975EW", 
                                        "count", "freqAboveDB"))

d_epi_means_plt <- AddCombosToDF(d_epi_means %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier, filter to tau = 0.0125
d_epi_means_plt <- d_epi_means_plt %>%
  filter(tau == 0.0125, r %in% r_subsample)

# Average over nloci, r, and optPerc (no effect)
d_epi_means_plt %>%
  group_by(model) %>%
  summarise(minEW = mean(minEW),
            q025EW = mean(q025EW),
            q25EW = mean(q25EW),
            q50EW = mean(q50EW),
            q75EW = mean(q75EW),
            q975EW = mean(q975EW),
            maxEW = mean(maxEW),
            meanFreqAboveDB = mean(freqAboveDB),
            CIFreqAboveDB = CI(freqAboveDB),
            n = sum(count)) -> d_epi_means_plt2

ggplot(d_epi_means_plt2 %>% 
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = model, y = meanFreqAboveDB, fill = model, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.2) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", y = "Proportion of samples with\nnon-negligible fitness epistasis", 
       fill = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_y_continuous(limits = c(0, 0.03)) +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_ew_freq_db
plt_ew_freq_db

# Per molecular component epistasis
# data
d_epi_means_percomp <- read.table(paste0(DATA_PATH, "d_epi_mean_percomp.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "molComp", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "minEW", "maxEW", "q025EW",
                                        "q25EW", "q50EW", "q75EW", "q975EW",
                                        "count", "freqAboveDB"))

d_epi_means_pc_plt <- AddCombosToDF(d_epi_means_percomp %>% 
                                   mutate(modelindex = as.factor(modelindex)))

d_epi_means_pc_plt <- d_epi_means_pc_plt %>%
  filter(tau == 0.0125, r %in% r_subsample,
         model != "Add")

# Average over r, nloci and optPerc
d_epi_means_pc_plt %>%
  mutate(molCompSt = molComp %>%
           str_split("_") %>%
           purrr::map(~ sort(.x) %>% paste(collapse = "_"))
  ) %>%
  distinct(optPerc, modelindex, molCompSt, .keep_all = T) %>% # Keep only 3_4s, not 4_3s
  mutate(molComp = as_factor(molComp)) %>%
  ungroup() -> d_epi_means_pc_plt

d_epi_means_pc_plt %>%
  group_by(model, molComp) %>%
  summarise(minEW = mean(minEW),
            q025EW = mean(q025EW),
            q25EW = mean(q25EW),
            q50EW = mean(q50EW),
            q75EW = mean(q75EW),
            q975EW = mean(q975EW),
            maxEW = mean(maxEW),
            meanFreqAboveDB = mean(freqAboveDB),
            CIFreqAboveDB = se(freqAboveDB),
            n = sum(count)) -> d_epi_means_pc_plt_sum


molCompComparisons <- c(
  TeX("$\\alpha_Z$ vs $\\alpha_Z$", output = "character"),
  TeX("$\\alpha_Z$ vs $\\beta_Z$", output = "character"),
  TeX("$\\beta_Z$ vs $\\beta_Z$", output = "character"),
  TeX("$\\alpha_Z$ vs $K_Z$", output = "character"),
  TeX("$\\alpha_Z$ vs $K_{XZ}$", output = "character"),
  TeX("$\\beta_Z$ vs $\\K_Z$", output = "character"),
  TeX("$\\beta_Z$ vs $\\K_{XZ}$", output = "character"),
  TeX("$K_Z$ vs $K_Z$", output = "character"),
  TeX("$K_Z$ vs $K_{XZ}$", output = "character"),
  TeX("$K_{XZ}$ vs $K_{XZ}$", output = "character")
)

d_epi_means_pc_plt_sum <- d_epi_means_pc_plt_sum %>%
  mutate(molCompLabel = molCompComparisons[as.integer(molComp)])

# Handle layout: center facets
design <- c("
            AABBCCDD
            #EEFFGG#
            ##HHII##
            ###JJ###
            "
)

KXZ_comparisons <- molCompComparisons[c(5, 7, 10, 9)]

# Plot frequency of epistasis effects greater than the size of the drift barrier
ggplot(d_epi_means_pc_plt_sum %>% 
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = model, y = meanFreqAboveDB, fill = model, colour = model)) +
  facet_manual(.~molCompLabel, design = design, axes = "x",
               labeller = labeller(molCompLabel = label_parsed)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.5) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1)[2:3],
                      guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1)[2:3],
                    guide = "none") +
  labs(x = "Model", y = TeX("Samples with non-negligible fitness epistasis (%)"), 
       fill = "Model") +
  scale_x_discrete(labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plt_freq_pc_db
plt_freq_pc_db
ggsave("plt_ew_freq_pc_molComp.png", width = 15, height = 5, device = png)


# Only draw KXZ comparisons
ggplot(d_epi_means_pc_plt_sum %>% filter(molCompLabel %in% KXZ_comparisons) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = molCompLabel, y = meanFreqAboveDB, fill = model, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.5) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1)[2:3],
                    guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1)[2:3],
                      guide = "none") +
  labs(x = "Molecular component comparison", 
       y = "Proportion of samples with\nnon-negligible fitness epistasis", 
       fill = "Model") +
  scale_x_discrete(labels = parse(text = KXZ_comparisons)) +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_freq_pc_db_KXZ
plt_freq_pc_db_KXZ
ggsave("plt_ew_freq_pc_KXZ.png", width = 5, height = 5, device = png)


# Combined figure
## Shared Y axis
y.grob <- textGrob("Proportion of samples with non-negligible fitness epistasis", 
                   gp=gpar(fontsize=14), rot=90)

plt_ew <- plot_grid(plt_ew_freq_db + theme(plot.margin = margin(5.5, 150, 5.5, 150),
                                           axis.title.y = element_blank(),
                                             axis.title.x = element_blank()),
                    plt_freq_pc_db + theme(axis.title.y = element_blank()),
                    nrow = 2, rel_heights = c(0.25, 1),
                    labels = "AUTO")

g <- arrangeGrob(plt_ew, left = y.grob)

ggsave("plt_ew_molComp.png", g, width = 7, height = 10, device = png)

# KXZ only figure
plt_ew <- plot_grid(plt_ew_freq_db,
                    plt_freq_pc_db_KXZ,
                    ncol = 2,
                    labels = "AUTO")
plt_ew
ggsave("plt_ew_KXZ.png", width = 9, height = 4, device = png)


# Table of results
d_epi_means_pc_plt_sum %>% select(model, molComp, minEW, q50EW, maxEW, meanFreqAboveDB, n)
d_epi_means_plt2 %>% select(model, minEW, q50EW, maxEW, meanFreqAboveDB, n)
