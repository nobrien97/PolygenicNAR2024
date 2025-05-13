# Mutational variance

d_h2 <- read.csv(paste0(DATA_PATH, "d_VA_Z.csv"))
d_h2 <- d_h2[,-1]


# Read initial data: pre-adjusted simulations (nloci = 1024, tau = 0.0125)
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", "variance"))


# Plot Vm in adjusted tau runs (nloci = 1024, tau = 0.0125*scalingFactor)
# scaling factor given by sensitivity analysis
d_mutvar_adjtau <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_adjtau.csv"), header = F,
                                     col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar_adjtau$normalised <- T
d_mutvar$normalised <- F

d_mutvar_adjtau <- d_mutvar_adjtau %>%
  mutate(replicate = replicate %/% 2 + 1, # convert replicate id from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_adjtau <- AddCombosToDF(d_mutvar_adjtau)

# Join with regular
d_mutvar2 <- full_join(d_mutvar, d_mutvar_adjtau, by = c("replicate", "seed", "modelindex",
                                                         "normalised", "variance", "model",
                                                         "nloci", "tau", "r"))

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
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1)) +
  scale_x_discrete(guide = "axis_nested") +
  labs(x = "Model", 
       y = "Mutational variance (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none")
ggsave("plt_vm_scaled.png", device = png, width = 6, height = 5)  

