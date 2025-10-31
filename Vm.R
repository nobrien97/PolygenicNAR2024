# Mutational variance
DATA_PATH <- "/mnt/c/GitHub/PolygenicNAR2024/data/"
# Read initial data: pre-adjusted simulations (nloci = 1024, tau = 0.0125)
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar_percomp <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_percomp.csv"), header = F,
                                      col.names = c("replicate", "seed", "modelindex", "cov_aZ", "cov_bZ",
                                                    "cov_KZ", "cov_KXZ"))

d_mutvar <- d_mutvar %>%
  mutate(replicate = replicate %/% 2 + 1, # convert replicate id from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))


d_mutvar_percomp <- d_mutvar_percomp %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_percomp <- AddCombosToDF(d_mutvar_percomp)

d_mutvar_percomp <- d_mutvar_percomp %>%
  pivot_longer(cols = starts_with("cov_"), names_to = "molComp", values_to = "cov") %>%
  # Clean data
  filter(cov < 1e1, cov > -1e1)



# Plot Vm in adjusted tau runs (nloci = 1024, tau = 0.0125*scalingFactor)
# scaling factor given by sensitivity analysis
d_mutvar_adjtau <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_adjtau.csv"), header = F,
                                     col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar_adjtau <- d_mutvar_adjtau %>%
  mutate(replicate = replicate %/% 2 + 1, # convert replicate id from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_adjtau$normalised <- T
d_mutvar$normalised <- F


# Join with regular
d_mutvar2 <- full_join(d_mutvar, d_mutvar_adjtau)
d_mutvar2 <- AddCombosToDF(d_mutvar2)

d_mutvar2 <- d_mutvar2 %>%
  mutate(scaled = if_else(normalised, TeX("Scaled $\\tau$", output = 'character'), 
                          TeX("Unscaled $\\tau$", output = 'character')),
         scaled = factor(scaled, levels = c("'Unscaled '*tau", "'Scaled '*tau")))

d_mutvar_sum <- d_mutvar2 %>%
  group_by(model, scaled) %>%
  summarise(meanVar = mean(log10(variance)),
            CIVar = CI(log10(variance)))

# Repeat for covariance
d_mutvar_percomp2 <- full_join(d_mutvar_percomp, d_mutvar_percomp_adjtau, 
                               by = c("replicate", "seed", "modelindex",
                                      "normalised", "cov", "molComp", "model",
                                      "nloci", "tau", "r"))

# Filter to the tau we tested
d_mutvar_percomp2 <- d_mutvar_percomp2 %>%
  filter(tau == 0.0125)

d_mutvar_percomp2 <- d_mutvar_percomp2 %>%
  mutate(scaled = if_else(normalised, TeX("Scaled $\\tau$", output = 'character'), 
                          TeX("Unscaled $\\tau$", output = 'character')),
         scaled = factor(scaled, levels = c("'Unscaled '*tau", "'Scaled '*tau")),
         scaledCOV = scale(cov)) %>%
  group_by(modelindex, seed, replicate, scaled) %>%
  mutate(contribution = abs(cov) / sum(abs(cov))) %>%
  ungroup()

d_mutvar_percomp_sum <- d_mutvar_percomp2 %>%
  group_by(model, molComp, scaled) %>%
  summarise(meanCOV = mean(cov),
            CICOV = CI(cov),
            meanAbsCov = mean(abs(cov)),
            meanScaledCOV = mean(scaledCOV),
            CIUScaledCov = CI(scaledCOV),
            meanCont = mean(contribution))


# Plot variance
ggplot(d_mutvar2 %>%
         mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE")), 
       aes(x = interaction(model, scaled), y = log10(variance), colour = model)) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_sum %>% ungroup() %>%
               mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE")), 
             aes(x = interaction(model, scaled), y = meanVar),
             shape = 3, size = 2, colour = "black", inherit.aes = F,
             stroke = 1) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1)) +
  scale_x_discrete(guide = "axis_nested") +
  labs(x = "Model", 
       y = "Mutational variance (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none") -> plt_vm
plt_vm
ggsave("plt_vm_scaled.png", device = png, width = 6, height = 5)  


# Covariance
ggplot(d_mutvar_percomp2 %>% filter(cov != 0.0) %>%
         mutate(model = fct_recode(model, "K+" = "K", "K-" = "ODE")), 
       aes(x = molComp, y = abs(cov), colour = model)) +
  facet_nested(. ~ scaled, labeller = labeller(scaled = label_parsed)) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_percomp_sum %>% ungroup() %>% filter(meanAbsCov > 0.0) %>%
               mutate(model = fct_recode(model, "K+" = "K", "K-" = "ODE")),
             aes(x = molComp, y = meanAbsCov, group = model),
             position = position_dodge(0.8),
             fill = "white", stroke = 1,
             shape = 21, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, 
                                           direction = 1)) +
  #coord_cartesian(ylim = c(-0.00001, 0.00001)) +
  scale_x_discrete(labels = c(TeX("$\\alpha_Z$"),
                              TeX("$\\beta_Z$"),
                              TeX("$K_Z$"),
                              TeX("$K_{XZ}$"))) +
  scale_y_log10() +
  labs(x = "Molecular component", 
       y = "Absolute mutational\ncovariance (log10)",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_covm
plt_covm
ggsave("plt_covm_scaled_new.png", device = png, width = 9, height = 5)  

# Plot adaptive trajectory

# load trait evolution data
d_qg_adjTau <- data.table::fread(paste0(DATA_PATH, "slim_qg_adjTau.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

d_qg_adjTau %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg_adjTau

# Attach additive replicates as well
d_qg_add <- AddCombosToDF(d_qg) %>% filter(model == "Add", r %in% r_subsample,
                                nloci == 1024, tau == 0.0125) 
# join
d_adapted_adjTau <- full_join(d_qg_adjTau, d_qg_add)

d_adapted_adjTau <- AddCombosToDF(d_adapted_adjTau)

d_adapted_adjTau_sum <- d_adapted_adjTau %>% 
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

ggplot(d_adapted_adjTau_sum,
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
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        panel.spacing = unit(0.75, "lines")) -> plt_adapt_adjTau
plt_adapt_adjTau
ggsave("plt_adapt_mutScale.png", width = 6, height = 5, device = png)

# Combine figures
# Plot as a grid
layout <-
  "
AAC
BBC
"
plt_vm + plt_covm + plt_adjtau_pheno -> plt_combined
plt_combined + plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') -> plt_vm_final

ggsave("plt_vm_final.png", plt_vm_final, width = 10, height = 8, device = png, bg = "white")

