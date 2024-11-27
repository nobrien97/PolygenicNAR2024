# Have a look at the epistasis data to see which models differ the most

# data
d_epi_means <- read.table(paste0(DATA_PATH, "d_epi_mean.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "count"))

d_epi_means_plt <- AddCombosToDF(d_epi_means %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier, filter to tau = 0.0125
d_epi_means_plt <- d_epi_means_plt %>%
  filter(tau == 0.0125, r %in% r_subsample)

d_epi_means_plt_sum <- d_epi_means_plt %>%
  group_by(model, r) %>%
  summarise(meanEWBar = mean(meanEW),
            CIEWBar = CI(meanEW),
            varEWBar = var(meanEW),
            n = n())

ggplot(d_epi_means_plt %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = model, y = meanEW, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_epi_means_plt_sum %>% 
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanEWBar, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", y = "Average fitness epistasis", colour = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave("plt_ew_sml.png", width = 4, height = 7, device = png)

# Per molecular component epistasis
# data
d_epi_means_percomp <- read.table(paste0(DATA_PATH, "d_epi_mean_percomp.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "molComp", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "count"))

d_epi_means_plt <- AddCombosToDF(d_epi_means_percomp %>% 
                                   mutate(modelindex = as.factor(modelindex)))

d_epi_means_plt <- d_epi_means_plt %>%
  filter(tau == 0.0125, r %in% r_subsample)

# Combine/average 3_4 and 4_3 etc.
d_epi_means_plt <- d_epi_means_plt %>%
  mutate(molCompSt = molComp %>%
           str_split("_") %>%
           purrr::map(~ sort(.x) %>% paste(collapse = "_"))
         ) %>%
  group_by(optPerc, modelindex, molCompSt) %>%
  mutate(meanEW = mean(meanEW, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(optPerc, modelindex, molCompSt, .keep_all = T) %>%
  mutate(molComp = as_factor(molComp))

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

d_epi_means_plt <- d_epi_means_plt %>%
  mutate(molCompLabel = molCompComparisons[as.integer(molComp)])

d_epi_means_plt_sum <- d_epi_means_plt %>%
  group_by(model, r, molCompLabel) %>%
  summarise(meanEWBar = mean(meanEW),
            CIEWBar = CI(meanEW),
            varEWBar = var(meanEW),
            n = n())

ggplot(d_epi_means_plt %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = model, y = meanEW, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Molecular component comparison" + molCompLabel, labeller = labeller(molCompLabel = label_parsed)) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_epi_means_plt_sum %>% 
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanEWBar, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", y = "Average fitness epistasis", colour = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))
ggsave("plt_ew_molComp.png", width = 20, height = 5, device = png)
