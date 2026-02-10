# Adaptive walk figures

# Add predictors to quant gen data
d_qg <- AddCombosToDF(d_qg)

# Proportion of each model that adapted
d_prop_adapted <- d_qg %>% group_by(model, nloci, tau, r) %>%
  filter(gen == 59950) %>%
  summarise(n = n(),
                   nAdapted = sum(isAdapted),
                   pAdapted = mean(isAdapted)
  )

# Average over nloci
# How many adapted in low recombination scenarios?
d_prop_adapted <- d_qg %>%
  mutate(r_cut = cut(r, breaks = c(0, 1e-7, 1),
                     labels = c("Low", "High"))) %>%
  group_by(model, tau, r_cut) %>%
  filter(gen == 59950) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )
d_prop_adapted

# Output to table
stargazer(d_prop_adapted)

# Summarise phenotype trajectories
d_adapted_sum <- d_qg %>%
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

# Small fx only: Fig. 1
ggplot(d_adapted_sum %>% filter(tau == 0.0125, r %in% r_subsample),
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~nloci) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean,
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)",
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~ ., name = "Number of QTLs",
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines"))
ggsave("plt_adapt_smlfx.png", width = 12, height = 5, device = png, dpi = 350)

#######
# Mean Z (K+) as a function of mean Z (K-) and mean Z(Additive) (Supp fig)
d_adapted_permod_sum <- d_adapted_sum %>%
  filter(tau == 0.0125, r %in% r_subsample, nloci == 1024) %>%
  select(gen, model, nloci, tau, r, meanPhenomean, SEPhenomean) %>%
  pivot_wider(names_from = "model",
              values_from = c("meanPhenomean", "SEPhenomean")
              )

ggplot(d_adapted_permod_sum,
       aes(x = meanPhenomean_ODE, y = meanPhenomean_K, colour = factor(log10(r))),
       group = as.factor(gen)) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_manual(values = c("#41476BFF", "#9E6374FF", "#EFBC82FF")) +
  labs(x = "Mean phenotype (K-)", y = "Mean phenotype (K+)",
       colour = "Recombination rate (log10)") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) -> plt_meanpheno_KvsODE

ggplot(d_adapted_permod_sum,
       aes(x = meanPhenomean_Add, y = meanPhenomean_K, colour = factor(log10(r))),
       group = as.factor(gen)) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_manual(values = c("#41476BFF", "#9E6374FF", "#EFBC82FF")) +
  labs(x = "Mean phenotype (Additive)", y = "Mean phenotype (K+)",
       colour = "Recombination rate (log10)") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) -> plt_meanpheno_KvsAdd

leg_meanpheno <- get_legend(plt_meanpheno_KvsODE)

plt_meanpheno_permodel <-
  plot_grid(plt_meanpheno_KvsAdd + theme(legend.position = "none"),
          plt_meanpheno_KvsODE + theme(legend.position = "none"),
          nrow = 1,
          labels = "AUTO")

plot_grid(plt_meanpheno_permodel,
            leg_meanpheno,
            nrow = 2,
          rel_heights = c(1, 0.1))
ggsave("plt_adapt_permod.png", width = 9, height = 4.5, bg = "white",
       device = png, dpi = 350)


# Burn-in evolution of molecular components vs phenotype (Fig. SY)
d_qg_burnin <- d_qg %>%
  filter(gen <= 50000, tau == 0.0125, model == "K", r %in% r_subsample) %>%
  select(gen, seed, r, phenomean, aZ, bZ, KZ, KXZ) %>%
  pivot_longer(cols = 4:8, names_to = "trait", values_to = "meanTrait")

ggplot(d_qg_burnin,
       aes(x = gen, y = meanTrait, colour = as.factor(log10(r)),
           group = factor(seed))) +
  facet_grid(as.factor(log10(r)) ~ trait) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("#41476BFF", "#9E6374FF", "#EFBC82FF")) +
  labs(x = "Generation", y = "Mean phenotype/molecular component",
       colour = "Recombination rate") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines"))



# Time to adaptation
d_adaptTime <- d_qg %>% filter(gen >= 49500) %>%
  mutate(gen = gen - 50000,
                isAdapted = between(phenomean, 1.9, 2.1),
                isAdapting = between(phenomean, 1.0, 1.9)) %>%
  group_by(seed, model, nloci, tau, r) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1),
                initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1)) %>%
  distinct(seed, model, nloci, tau, r, .keep_all = T) %>%
  ungroup()

# Distribution of adaptation times for larger effect sizes
ggplot(d_adaptTime %>% filter(adaptTime > -1,
                              r %in% r_subsample,
                              tau > 0.0125) %>%
         mutate(r_title = "Recombination rate (log10)",
                tau_title = "Mutational effect size variance"),
       aes(x = adaptTime, fill = model), alpha = 0.4) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_x_continuous(labels = scales::comma) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  labs(x = "Time to adaptation (Generations)", y = "Density", fill = "Model")
ggsave("sfig_timetoadaptation_largetau.png", device = png, width = 10, height = 6)

# Average adaptation time in the larger effect size models
mean(d_adaptTime[d_adaptTime$tau > 0.0125,]$adaptTime)
CI(d_adaptTime[d_adaptTime$tau > 0.0125,]$adaptTime)

# Mean phenotype figures for large effect models
{
  ggplot(d_adapted_sum %>% filter(tau == 0.0125),
         aes(x = gen, y = meanPhenomean, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    geom_hline(yintercept = 2, linetype = "dashed") +
    ggtitle("Tau = 0.0125") +
    geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean,
                    ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
                alpha = 0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)",
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Number of QTLs",
                                           breaks = NULL, labels = NULL)) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
    labs(x = "Generations post-optimum shift", y = "Mean phenotype",
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_smlFX
  ggsave("adapt_grid_smlFX.png", adapt_grid_smlFX, width = 14, height = 10, device = png)

  ggplot(d_adapted_sum %>% filter(tau == 0.125),
         aes(x = gen, y = meanPhenomean, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    geom_hline(yintercept = 2, linetype = "dashed") +
    ggtitle("Tau = 0.125") +
    geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean,
                    ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
                alpha = 0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)",
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Number of QTLs",
                                           breaks = NULL, labels = NULL)) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
    labs(x = "Generations post-optimum shift", y = "Mean phenotype",
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_medFX
  ggsave("adapt_grid_medFX.png", adapt_grid_medFX, width = 14, height = 10, device = png)

  ggplot(d_adapted_sum %>% filter(tau == 1.25),
         aes(x = gen, y = meanPhenomean, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    geom_hline(yintercept = 2, linetype = "dashed") +
    ggtitle("Tau = 1.25") +
    geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean,
                    ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
                alpha = 0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)",
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Number of QTLs",
                                           breaks = NULL, labels = NULL)) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
    labs(x = "Generations post-optimum shift", y = "Mean phenotype",
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_lrgFX
  ggsave("adapt_grid_lrgFX.png", adapt_grid_lrgFX, width = 14, height = 10, device = png)

  adapt_grid <- plot_grid(adapt_grid_smlFX, adapt_grid_medFX, adapt_grid_lrgFX,
                          labels = "AUTO", nrow = 3)

  ggsave("adapt_grid.png", adapt_grid, width = 14, height = 30, device = png)
}

# Adaptive walk following randomised starting points
d_qg_rs <- data.table::fread(paste0(DATA_PATH, "slim_qg_randomisedStarts.csv"), header = F,
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    rep("numeric", times = 12)),
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"),
                          fill = T)

# Add on optPerc
d_qg_rs %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg_rs

d_qg_rs$optPerc <- d_qg_rs$phenomean - 1
d_qg_rs$optPerc <- cut(d_qg_rs$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg_rs_optPerc <- d_qg_rs %>% select(gen, seed, modelindex, optPerc, isAdapted) %>% filter(gen >= 49500)

d_qg_rs <- AddCombosToDF(d_qg_rs)

# Proportion of each model that adapted
d_prop_adapted_rs <- d_qg_rs %>% group_by(model, nloci, tau, r) %>%
  filter(gen == 59950) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )

# Average over nloci
# How many adapted in low recombination scenarios?
d_prop_adapted_rs <- d_qg_rs %>%
  mutate(r_cut = cut(r, breaks = c(0, 1e-7, 1),
                     labels = c("Low", "High"))) %>%
  group_by(model, tau, r_cut) %>%
  filter(gen == 59950) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )
d_prop_adapted_rs

d_adapted_rs_sum <- d_qg_rs %>%
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

# Small fx only
ggplot(d_adapted_rs_sum %>% filter(tau == 0.0125, r %in% r_subsample),
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~.) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean,
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)",
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                    labels = c("K+", "K-"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines"))
ggsave("plt_adapt_rs.png", width = 6, height = 5, device = png)
