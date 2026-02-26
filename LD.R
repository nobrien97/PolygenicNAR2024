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
# Additive under low recombination, has very few data points - highlight that
# Use custom transparent boxplot
alpha_lookup <- d_ld_dist_hist %>%
  filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1,
         tau == 0.0125) %>%
  mutate(poorData = if_else((model == "Add" & log10(r) > -5) | model != "Add",
                            1, 0.2)) %>%
  distinct(model, r, poorData) %>%
  arrange(model)

pal <- rep(paletteer_d("nationalparkcolors::Everglades", 3, direction = -1), each = 3)
pal_alpha <- alpha(pal, alpha_lookup$poorData)
pal_alpha <- colorspace::desaturate(pal_alpha, 1 - alpha_lookup$poorData)

ggplot(d_ld_dist_hist %>%
         mutate(poorData = if_else((model == "Add" & log10(r) > -5) | model != "Add",
                                   1, 0.2)
                )%>%
         mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 col_num = as.numeric(col) + 0.025/2,
                                 r_title = "Recombination rate (log10)",
                                 nloci_title = "Number of loci",
                                 tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add",
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1,
                tau == 0.0125),
       aes(x = col_num, y = prop, colour = interaction(r, model), group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(
               position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = interaction(r,model))
  ) +
  scale_colour_manual(values = pal_alpha, guide = "none") +
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
ggplot(d_ld_freq_dist_hist %>%
         mutate(poorData = if_else((model == "Add" & log10(r) > -5) | model != "Add",
                                   1, 0.2)
         )%>%
         mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                      col_num = as.numeric(col) + 0.025/2,
                                      r_title = "Recombination rate (log10)",
                                      nloci_title = "Number of loci",
                                      tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add",
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1) %>%
         filter(tau == 0.0125),
       aes(x = col_num, y = prop, colour = interaction(r, model),
           group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.2) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = interaction(r, model))
  ) +
  scale_colour_manual(values = pal_alpha,
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_freq_sml
plt_ld_freq_sml
ggsave("plt_ld_freq_sml.png", device = png, width = 9, height = 4)

# Grid figure: Fig. 4
leg <- get_legend(plt_ld_freq_sml)

plt_ld_sml_com <- plot_grid(plt_ld_sml + theme(legend.position = "none"),
                            plt_ld_freq_sml + theme(legend.position = "none"),
                            ncol = 1, labels = "AUTO")

plt_ld_sml_com
ggsave("plt_ld_com.png", device = png, bg = "white", dpi = 350,
       width = 9, height = 8)


# LD among molecular components
d_ld_molcomp <- data.table::fread(paste0(DATA_PATH, "out_LDm.csv"), header = F,
                               col.names = c("gen", "seed", "modelindex", "mutType",
                                             "meanD", "sdD", "nD",
                                             "nDP", "nDN", "nDHalf",
                                             paste0("n", 1:20)))

d_ld_molcomp <- d_ld_molcomp %>%
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex))


# inner join optPerc
d_ld_molcomp <- left_join(d_ld_molcomp, d_qg, by = c("gen", "seed", "modelindex"))


# Add on variables
d_ld_molcomp <- AddCombosToDF(d_ld_molcomp)

# Combine redundant mutTypes 5_6 vs 6_5
d_ld_molcomp <- d_ld_molcomp %>%
  separate(mutType, into = c("mutTypeA", "mutTypeB")) %>%
  mutate(mutTypeA = as.numeric(mutTypeA),
         mutTypeB = as.numeric(mutTypeB)) %>%
  rowwise() %>%
  mutate(mutType = paste0(min(mutTypeA, mutTypeB), "_", max(mutTypeA, mutTypeB))) %>%
  select(-c(mutTypeA, mutTypeB))

# Proportion of estimates with positive/negative D
d_ld_molcomp <- d_ld_molcomp %>%
  group_by(optPerc, mutType, model, nloci, tau, r) %>%
  mutate(propDP = nDP / nD,
         propDN = nDN / nD)

# average across replicates
d_ld_molcomp_sum <- d_ld_molcomp %>%
  group_by(optPerc, mutType, model, nloci, tau, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

# Reorganise D estimates for plotting
bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_molcomp_dist <- d_ld_molcomp_sum %>% select(optPerc, mutType, model, nloci, tau, r, 12:31) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_molcomp_dist$col <- bins[as.numeric(str_extract(d_ld_molcomp_dist$col, "[[0-9]]*(?=_)"))]

d_ld_dist_molcomp_hist <- d_ld_molcomp %>% select(gen, seed, mutType, optPerc, model, nloci, tau, r, 10:29) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, mutType, optPerc, model, nloci, tau, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Plot histogram
# Offset x axis since we group leftwise [x, y), offset is half the bin size, 0.025/2
mutType_comp_labeller <-
  c("3_3" = TeX("$\\alpha_Z \\, vs \\, \\alpha_Z$", output = "character"),
    "3_4" = TeX("$\\alpha_Z \\, vs \\, \\beta_Z$", output = "character"),
    "3_5" = TeX("$\\alpha_Z \\, vs \\, K_Z$", output = "character"),
    "3_6" = TeX("$\\alpha_Z \\, vs \\, K_{XZ}$", output = "character"),
    "4_4" = TeX("$\\beta_Z \\, vs \\, \\beta_Z$", output = "character"),
    "4_5" = TeX("$\\beta_Z \\, vs \\, K_Z$", output = "character"),
    "4_6" = TeX("$\\beta_Z \\, vs \\, K_{XZ}$", output = "character"),
    "5_5" = TeX("$K_Z \\, vs \\, K_Z$", output = "character"),
    "5_6" = TeX("$K_Z \\, vs \\, K_{XZ}$", output = "character"),
    "6_6" = TeX("$K_{XZ} \\, vs \\, K_{XZ}$", output = "character"))

ggplot(d_ld_dist_molcomp_hist %>% filter(mutType %in% c("3_3", "3_4", "3_5", "3_6", "4_4")) %>%
         mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                col_num = as.numeric(col) + 0.025/2,
                r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add",
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1,
                tau == 0.0125, model == "K+"),
       aes(x = col_num, y = prop, colour = model,
           group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ "Molecular component comparison" + mutType,
               labeller = labeller(mutType = as_labeller(mutType_comp_labeller,
                                                         default = label_parsed))) +
  geom_boxplot(
    position = position_identity(), outlier.shape = 1,
    outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = 1)[2],
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_molcomp_sml_pt1
plt_ld_molcomp_sml_pt1

ggplot(d_ld_dist_molcomp_hist %>% filter(mutType %in% c("4_5", "4_6", "5_5", "5_6", "6_6")) %>%
         mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                col_num = as.numeric(col) + 0.025/2,
                r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add",
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(log10(r) == -10 | log10(r) == -5 | log10(r) == -1,
                tau == 0.0125, model == "K+"),
       aes(x = col_num, y = prop, colour = model,
           group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ "Molecular component comparison" + mutType,
               labeller = labeller(mutType = as_labeller(mutType_comp_labeller,
                                                         default = label_parsed))) +
  geom_boxplot(
    position = position_identity(), outlier.shape = 1,
    outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = 1)[2],
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_molcomp_sml_pt2
plt_ld_molcomp_sml_pt2

plot_grid(plt_ld_molcomp_sml_pt1,
          plt_ld_molcomp_sml_pt2,
          nrow = 2)

ggsave("plt_ld_molcomp_sml.png", device = png, width = 9, height = 8)

d_ld_molcomp_gls <- d_ld_molcomp %>%
  filter(model == "K", !is.na(sdD), sdD > 0,
         r %in% r_subsample, tau == 0.0125) %>%
  mutate(r = factor(r),
         mutType = factor(mutType))

gls.ld.molcomp.k <- gls(log(sdD) ~ mutType * r, d_ld_molcomp_gls,
                        weights = varIdent(form = ~ 1 | r))
summary(gls.ld.molcomp.k)
plot(gls.ld.molcomp.k)
anova(gls.ld.molcomp.k)

em.ld.molcomp <- emmeans(gls.ld.molcomp.k, ~ mutType * r)
pairs(em.ld.molcomp, simple = "mutType")
plot(em.ld.molcomp, comparisons = T)
joint_tests(em.ld.molcomp)

em_pred <- update(em.ld.molcomp)
em_pred_sum <- summary(em_pred, infer = c(TRUE, FALSE), adjust = "none", type = "lp")

ggplot(em_pred_sum %>% mutate(r = factor(log10(as.numeric(as.character(r))))),
       aes(x = emmean, y = mutType, colour = r)) +
  geom_segment(aes(x = lower.CL, xend = upper.CL, y = mutType, yend = mutType,
                   colour = r)) +
  geom_point() +
  scale_y_discrete(labels = parse(text = mutType_comp_labeller)) +
  scale_colour_manual(values = c("#41476BFF", "#9E6374FF", "#EFBC82FF")) +
  labs(x = TeX("Predicted standard deviation of D ($log_{e}$)"),
       y = "Molecular component comparison",
       colour = TeX("Recombination rate ($log_{10}$)")) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom") -> plt_ld_pred
plt_ld_pred
ggsave("plt_ld_molcomp_fx.png", device = png, width = 5, height = 5)

# Combined model for mutType comparison "contains K_Z"
d_ld_molcomp_gls_kz <- d_ld_molcomp_gls %>%
  mutate(mutTypeKZ = if_else(grepl("5", mutType), T, F))
gls.ld.molcomp.kz <- gls(log(sdD) ~ mutTypeKZ * r, d_ld_molcomp_gls_kz,
                         weights = varIdent(form = ~ 1 | r))
summary(gls.ld.molcomp.kz)
plot(gls.ld.molcomp.kz)
anova(gls.ld.molcomp.kz)

em.ld.molcomp.kz <- emmeans(gls.ld.molcomp.kz, ~ mutTypeKZ * r)
pairs(em.ld.molcomp.kz, simple = "mutTypeKZ")
plot(em.ld.molcomp.kz, comparisons = T)
joint_tests(em.ld.molcomp.kz)

summary(em.ld.molcomp.kz, infer = T, type = "response")

em_pred_kz <- update(em.ld.molcomp.kz)
em_pred_kz_sum <- summary(em_pred_kz, infer = c(TRUE, FALSE), adjust = "none", type = "lp")

ggplot(em_pred_kz_sum %>% mutate(r = factor(log10(as.numeric(as.character(r))))),
       aes(x = emmean, y = mutTypeKZ, colour = r)) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL, y = mutTypeKZ,
                    colour = r), position = position_dodge(0.3, orientation = "y"),
                width = 0.2) +
  geom_point(position = position_dodge(0.3, orientation = "y")) +
  scale_y_discrete(labels =
                     parse(text = c(TeX("Does not contain $K_Z", output = "character"),
                                    TeX("Contains $K_Z", output = "character")))) +
  scale_colour_manual(values = c("#41476BFF", "#9E6374FF", "#EFBC82FF")) +
  labs(x = TeX("Predicted standard deviation of D ($log_{e}$)"),
       y = "D pairwise comparison identity",
       colour = TeX("Recombination rate ($log_{10}$)")) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom",
        axis.text.y = element_text(angle = 90, hjust = 0.5)) -> plt_ld_pred_sd
plt_ld_pred_sd
ggsave("plt_ld_molcomp_fx.png", device = png, width = 5, height = 5)

# Mean D
gls.ld.molcomp.k.mean <- gls(meanD ~ mutType * r, d_ld_molcomp_gls,
                        weights = varIdent(form = ~ 1 | r))
summary(gls.ld.molcomp.k.mean)
plot(gls.ld.molcomp.k.mean)
anova(gls.ld.molcomp.k.mean)

em.ld.molcomp.k.mean <- emmeans(gls.ld.molcomp.k.mean, ~ mutType * r)
pairs(em.ld.molcomp.k.mean, simple = "mutType")
plot(em.ld.molcomp.k.mean, comparisons = T)
joint_tests(em.ld.molcomp)


gls.ld.molcomp.kz.mean <- gls(meanD ~ mutTypeKZ * r, d_ld_molcomp_gls_kz,
                         weights = varIdent(form = ~ 1 | r))
summary(gls.ld.molcomp.kz.mean)
plot(gls.ld.molcomp.kz.mean)
anova(gls.ld.molcomp.kz.mean)

em.ld.molcomp.kz.mean <- emmeans(gls.ld.molcomp.kz.mean, ~ mutTypeKZ * r)
pairs(em.ld.molcomp.kz.mean, simple = "mutTypeKZ", reverse = T)
plot(em.ld.molcomp.kz.mean, comparisons = T)
joint_tests(em.ld.molcomp.kz.mean)

# Comparisons
pwpm(em.ld.molcomp.kz.mean)

summary(em.ld.molcomp.kz.mean, infer = T, type = "response")
em_pred_kz_mean <- update(em.ld.molcomp.kz.mean)
em_pred_kz_mean_sum <- summary(em_pred_kz_mean, infer = c(TRUE, FALSE), adjust = "none", type = "lp")

ggplot(em_pred_kz_mean_sum %>% mutate(r = factor(log10(as.numeric(as.character(r))))),
       aes(x = emmean, y = mutTypeKZ, colour = r)) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL, y = mutTypeKZ,
                   colour = r), position = position_dodge(0.3, orientation = "y"),
                width = 0.2) +
  geom_point(position = position_dodge(0.3, orientation = "y")) +
  scale_y_discrete(labels =
                     parse(text = c(TeX("Does not contain $K_Z", output = "character"),
                                    TeX("Contains $K_Z", output = "character")))) +
  scale_colour_manual(values = c("#41476BFF", "#9E6374FF", "#EFBC82FF")) +
  labs(x = TeX("Predicted mean of D"),
       y = "D pairwise comparison identity",
       colour = TeX("Recombination rate ($log_{10}$)")) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom",
        axis.text.y = element_text(angle = 90, hjust = 0.5)) -> plt_ld_pred_mean
plt_ld_pred_mean

leg <- get_legend(plt_ld_pred_mean)

plot_grid(plt_ld_pred_mean + theme(legend.position = "none"),
         plt_ld_pred_sd + theme(legend.position = "none"),
         align = "h", axis = "bt", rel_widths = c(1,0.9),
         labels = "AUTO",
         nrow = 1) -> plt_ld_pred

plot_grid(plt_ld_pred,
          leg,
          nrow = 2,
          rel_heights = c(1, 0.1))

ggsave("plt_ld_molcomp.png", device = png, bg = "white", dpi = 350,
       width = 8, height = 4)

# Molecular component has small effect overall, recombination rate is the large effect
# without r, R2 is 0.017, with r, it's 0.681
# however, under high recombination, LD between involving K_Z loci tend to have more linkage
# if they are effectively neutral, they would just hitch hike everywhere, but why not expect that
# in the other recombination rates?
# Mean doesn't really change either


