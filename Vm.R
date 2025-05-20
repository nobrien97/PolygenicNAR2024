# Mutational variance

# Read initial data: pre-adjusted simulations (nloci = 1024, tau = 0.0125)
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar <- d_mutvar %>%
  mutate(replicate = replicate %/% 2 + 1, # convert replicate id from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

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
  mutate(scaled = if_else(normalised, "Scaled tau", "Unscaled tau"),
         scaled = factor(scaled, levels = c("Unscaled tau", "Scaled tau")))

d_mutvar_sum <- d_mutvar2 %>%
  group_by(model, scaled) %>%
  summarise(meanVar = mean(log10(variance)),
            CIVar = CI(log10(variance)))

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

leg <- get_legend(plt_adapt_adjTau)

# Combine figures
plt_vm_final <- plot_grid(plt_vm + theme(legend.position = "none"),
                          plt_adapt_adjTau + theme(legend.position = "none"),
                          nrow = 2,
                          labels = "AUTO")

plot_grid(plt_vm_final, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_vm_final
plt_vm_final
ggsave("plt_vm_final.png", plt_vm_final, width = 8, height = 10, device = png, bg = "white")

