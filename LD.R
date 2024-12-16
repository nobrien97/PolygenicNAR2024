# Plot LD

# Load data: all pairs
d_ld <- read.table(paste0(DATA_PATH, "out_LD.csv"), header = F,
                   col.names = c("gen", "seed", "modelindex", 
                                 "meanD", "sdD", "nD",
                                 "nDP", "nDN", "nDHalf", 
                                 paste0("n", 1:20)))


# LD between loci with similar MAF
d_ld_freq <- data.table::fread(paste0(DATA_PATH, "out_LDf.csv"), header = F,
                               col.names = c("gen", "seed", "modelindex", "freqBin",
                                             "meanD", "sdD", "nD",
                                             "nDP", "nDN", "nDHalf",
                                             paste0("n", 1:20)))


d_ld <- d_ld %>% 
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex))


# inner join optPerc
d_ld <- left_join(d_ld, d_qg, by = c("gen", "seed", "modelindex"))


# Add on variables
d_ld <- AddCombosToDF(d_ld)

# Proportion of estimates with positive/negative D
d_ld <- d_ld %>%
  group_by(optPerc, model, nloci, tau, r) %>%
  mutate(propDP = nDP / nD,
         propDN = nDN / nD)

# average across replicates
d_ld_sum <- d_ld %>%
  group_by(optPerc, model, nloci, tau, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

# Reorganise D estimates for plotting
bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_dist <- d_ld_sum %>% select(optPerc, model, nloci, tau, r, 12:31) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_dist$col <- bins[as.numeric(str_extract(d_ld_dist$col, "[[0-9]]*(?=_)"))]

d_ld_dist_sd <- d_ld_sum %>% select(optPerc, model, nloci, tau, r, 40:59) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count_sd")

d_ld_dist_hist <- d_ld %>% select(gen, seed, optPerc, model, nloci, tau, r, 10:29) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, optPerc, model, nloci, tau, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Plot histogram
# Offset x axis since we group leftwise [x, y), offset is half the bin size, 0.025/2
ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 col_num = as.numeric(col) + 0.025/2,
                                 r_title = "Recombination rate (log10)",
                                 nloci_title = "Number of loci",
                                 tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1,
                tau == 0.0125), 
       aes(x = col_num, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_sml
plt_ld_sml
ggsave("plt_ld_sml.png", device = png, width = 9, height = 6)

#################################
# Frequency adjusted
d_ld_freq <- d_ld_freq %>% 
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex)) %>%
  distinct()

# inner join optPerc
d_ld_freq <- left_join(d_ld_freq, d_qg, by = c("gen", "seed", "modelindex"))

# Add on variables
d_ld_freq <- AddCombosToDF(d_ld_freq)

# average across replicates
d_ld_freq_sum <- d_ld_freq %>%
  group_by(optPerc, freqBin, model, nloci, tau, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

d_ld_freq_sum <- d_ld_freq_sum %>%
  mutate(freqBin = if_else(freqBin > 0.5, 1 - freqBin, freqBin))

# Reshape for plotting
bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_freq_dist <- d_ld_freq_sum %>% select(optPerc, freqBin, model, nloci, tau, r, 13:32) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_freq_dist$col <- bins[as.numeric(str_extract(d_ld_freq_dist$col, "[[0-9]]*(?=_)"))]

d_ld_freq_dist_sd <- d_ld_freq_sum %>% select(optPerc, freqBin, model, nloci, tau, r, 39:58) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count_sd")

# Outliers: histogram of all estimates
d_ld_freq_dist_hist <- d_ld_freq %>% filter(freqBin > 0.1) %>% 
  select(gen, seed, optPerc, model, nloci, tau, r, 11:30) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, optPerc, model, nloci, tau, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Number of loci doesn't matter
# Doesn't appear to be related to effect size variance
ggplot(d_ld_freq_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                      col_num = as.numeric(col) + 0.025/2,
                                      r_title = "Recombination rate (log10)",
                                      nloci_title = "Number of loci",
                                      tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1) %>%
         filter(tau == 0.0125), 
       aes(x = col_num, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.2) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_freq_sml
plt_ld_freq_sml
ggsave("plt_ld_freq_sml.png", device = png, width = 9, height = 4)

# Grid figure
leg <- get_legend(plt_ld_freq_sml)

plt_ld_sml_com <- plot_grid(plt_ld_sml + theme(legend.position = "none"),
                            plt_ld_freq_sml + theme(legend.position = "none"),
                            ncol = 1, labels = "AUTO")

plt_ld_sml_com
ggsave("plt_ld_com.png", device = png, bg = "white",
       width = 9, height = 8)
