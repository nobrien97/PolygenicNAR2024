# Analysis of effect sizes for each molecular component

d_fx <- data.table::fread(paste0(DATA_PATH, "d_fx.csv"), header = F,
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    "factor", "character", "numeric"),
                          col.names = c("gen", "seed", "modelindex",
                                        "mutType", "mutID", "s"),
                          fill = T)

# Join with phenotypic data
d_fx <- d_fx %>% distinct()
d_fx <- left_join(d_fx, d_qg_optPerc, by = c("gen", "seed", "modelindex"))
d_fx <- AddCombosToDF(d_fx)

# Filter to groups we're looking at
d_fx <- d_fx %>%
  filter(isAdapted, tau == 0.0125, r %in% r_subsample)

# plot distribution of fitness effects
ggplot(d_fx %>%
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = s)) +
  facet_nested("Model" + model ~ "Mutation Type" + mutType) +
  geom_boxplot(position = position_dodge(0.9)) +
  labs(x = "Selection coefficient (s)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Number of neutral mutations per molecular component/model
d_fx %>%
  mutate(mutType = factor(mutType, levels = c("2", levels(mutType)))) %>%
  group_by(model, mutType) %>%
  summarise(propNeutral = sum(abs(s) < 1e-10) / n())
# # A tibble: 7 × 3
# # Groups:   model [3]
# model mutType propNeutral
# <chr> <fct>         <dbl>
# Add   3       0.000000464
# K     3       0.000235
# K     4       0.000185
# K     5       0.958
# K     6       0.0000882
# ODE   3       0.000176
# ODE   4       0.000179


# Plot the number of beneficial mutations per molecular component
d_fx_ben <- d_fx %>% filter(s > 0)

mutTypes_vec <- c(TeX("Additive", output = "character"),
                  TeX("$\\alpha_Z$", output = "character"),
                  TeX("$\\beta_Z$", output = "character"),
                  TeX("$K_Z$", output = "character"),
                  TeX("$K_{XZ}$", output = "character"))


# Proportion of mutations that are beneficial in each model
d_fx_propBen <- d_fx %>%
  mutate(mutType = factor(mutType, levels = c("2", levels(mutType)))) %>%
  group_by(seed, model, mutType, nloci, tau, r) %>%
  mutate(isBen = (s > 0)) %>%
  summarise(propBen = sum(isBen)/n()) %>%
  ungroup()

# Set muttype for additive to a unique id
d_fx_propBen[d_fx_propBen$model == "Add", "mutType"] <- "2"

# Summarise
d_fx_propBen_sum <- d_fx_propBen %>%
  group_by(model, mutType, r) %>%
  summarise(meanPropBen = mean(propBen),
            CIPropBen = CI(propBen))

# Proportion of beneficial mutations per molecular component - Fig 8
ggplot(d_fx_propBen %>%
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = propBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_propBen_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanPropBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9),
             stroke = 1) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component",
       y = "Proportion of beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_propben_muts_wholewalk
plt_propben_muts_wholewalk
ggsave("plt_propben_muts_wholewalk.png", plt_propben_muts_wholewalk, dpi = 350,
       width = 10, height = 4, device = png)

# beta regression
# adjust for 0/1 values - need to inflate (Smithson + Verkuilen 2006)
# Set up combinations
# Use marginaleffects since it supports betareg
library(marginaleffects)

d_fx_propBen <- transform(d_fx_propBen, modelMutType = factor(interaction(model, mutType)))
d_fx_propBen$propBen_adj <- (d_fx_propBen$propBen * (nrow(d_fx_propBen)-1) + 0.5)/nrow(d_fx_propBen)
d_fx_propBen$model <- as.factor(d_fx_propBen$model)

# Bimodality driven by mutation type alone
plot(density(d_fx_propBen[d_fx_propBen$mutType != 5,]$propBen))
plot(density(d_fx_propBen[d_fx_propBen$mutType == 5,]$propBen))

# beta regression for everything
br.benmut <- betareg(propBen_adj ~ modelMutType, d_fx_propBen)

summary(br.benmut)
plot(br.benmut)

saveRDS(br.benmut, "betareg_benmut.RDS")
br.benmut <- readRDS(paste0(DATA_PATH, "betareg_benmut.RDS"))

#anova(br.benmut)
em.benmut <- emmeans(br.benmut, ~ modelMutType)

pairs(em.benmut, simple = "modelMutType")
plot(em.benmut, comparisons = T)
emmip(em.benmut,  ~ modelMutType)

# Differences between models (apart from KZ) range from <1% to ~5%
# KZ is about 50% change, almost never beneficial



# Calculate beta regression separately for each model:
#   Different molecular components in each
# In additive only one mutation type, so calc mean and CI
mean.benmut.add <- d_fx_propBen %>% filter(model == "Add") %>%
  summarise(meanPropBen = mean(propBen_adj),
            SEPropBen = se(propBen_adj),
            CIPropBen = CI(propBen_adj),
            upperCL = meanPropBen + CIPropBen,
            lowerCL = meanPropBen - CIPropBen)

br.benmut.km <- betareg(propBen_adj ~ mutType, data = d_fx_propBen, subset = model == "ODE")
br.benmut.kp <- betareg(propBen_adj ~ mutType, data = d_fx_propBen, subset = model == "K")

em.benmut.km <- emmeans(br.benmut.km, ~ mutType)
em.benmut.kp <- emmeans(br.benmut.kp, ~ mutType)

# Were beneficial mutations more likely in some molecular components than others?
joint_tests(em.benmut.km)
joint_tests(em.benmut.kp)

# Which ones in particular?
xtable(em.benmut.km)
xtable(em.benmut.kp)
# And the additive model
mean.benmut.add

# What about effect size: if K_Z has very large effects, generating that variation
# could be important
d_fx_ben$mutType <- as.character(d_fx_ben$mutType)
d_fx_ben[d_fx_ben$model == "Add", "mutType"] <- "2"
d_fx_ben$mutType <- factor(d_fx_ben$mutType)

d_fx_ben_sum <- d_fx_ben %>%
  group_by(model, mutType) %>%
  summarise(meanBen = mean(s),
            CIBen = CI(s))

# Average effect size of mutations per model/molecular component - Fig. 9
ggplot(d_fx_ben %>%
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = s, colour = model)) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_ben_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9),
             stroke = 1) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component",
       y = "Average fitness effect\nof beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_ben_muts_s
plt_ben_muts_s
ggsave("plt_ben_muts_s.png", plt_ben_muts_s, dpi = 350,
       width = 6, height = 4, device = png)


# GLS - fitness effect of beneficial mutations

summary(gls.s.km <- gls(s ~ mutType, (d_fx_ben %>% filter(model == "ODE")),
                        weights = varIdent(form = ~ 1 | mutType)))

plot(gls.s.km)

# No difference between alpha and beta in effect size, so can just calc the mean
# like in the additive model across both alpha and beta
mean.s.km <- d_fx_ben %>% filter(model == "ODE") %>%
  summarise(meanS = mean(s),
            SES = se(s),
            CIS = CI(s),
            upperCL = meanS + CIS,
            lowerCL = meanS - CIS)

# Per model gls
d_fx_ben_km <- d_fx_ben %>% select(-optPerc) %>% filter(model == "ODE") %>% distinct()
summary(gls.s.km <- gls(s ~ mutType, d_fx_ben_km,
                        weights = varIdent(form = ~ 1 | mutType)))

d_fx_ben_k <- d_fx_ben %>% select(-optPerc) %>% filter(model == "K") %>% distinct()
summary(gls.s.kp <- gls(s ~ mutType, d_fx_ben_k,
                        weights =
                          varIdent(form = ~ 1 | mutType)))

plot(gls.s.kp)

# Marginal means
em.s.kp <- emmeans(gls.s.kp, ~mutType, mode = "appx-satterthwaite")

# Model/mutType anova
anova(gls.s.kp)
anova(gls.s.km)

em.s.kp

# Additive effect
mean.s.add <- d_fx_ben %>% filter(model == "Add") %>%
  summarise(meanS = mean(s),
            SES = se(s),
            CIS = CI(s),
            lowerCL = meanS - CIS,
            upperCL = meanS + CIS)

# format into table
xtable(mean.s.add, digits = 6)
xtable(mean.s.km, digits = 6)
xtable(em.s.kp, digits = 9)

# Frequency data
d_freqs <- data.table::fread(paste0(DATA_PATH, "d_mutfreq.csv"), header = F,
                          sep = ",", colClasses = c("factor", "factor", "factor",
                                                    "factor", "character", "numeric",
                                                    "integer", "integer"),
                          col.names = c("optPerc", "seed", "modelindex",
                                        "mutType", "mutID", "freq",
                                        "count", "fixGen"),
                          fill = T)

d_freqs <- AddCombosToDF(d_freqs)

# Remove duplicates
d_freqs <- d_freqs %>%
  distinct()

# Look at relevant treatments
d_freqs <- d_freqs %>%
  filter(tau == 0.0125, r %in% r_subsample)

# Set muttype for additive to a unique id
d_freqs <- d_freqs %>%
  mutate(mutType = factor(mutType, levels = c("2", levels(mutType))))
d_freqs[d_freqs$model == "Add", "mutType"] <- "2"

# Look at number of fixations
d_fix <- d_freqs[!is.na(d_freqs$fixGen),]

d_fix_n <- d_fix %>%
  filter(optPerc == "[0.9, Inf)") %>%
  group_by(seed, model, r, nloci, mutType) %>%
  summarise(nFixations = n(),
            meanFixGen = mean(fixGen),
            CIFixGen = CI(fixGen))

gls.numfix <- gls(nFixations ~ model, d_fix_n,
    weights = varIdent(form = ~ 1 | model))
summary(gls.numfix)
plot(gls.numfix)
anova(gls.numfix)

ggplot(d_fix_n %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = mutType, y = nFixations,
           colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_boxplot() +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component",
       y = "Number of fixations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Among models - no difference between components
d_fix_n <- d_fix %>%
  filter(optPerc == "[0.9, Inf)") %>%
  group_by(seed, model, r, nloci) %>%
  summarise(nFixations = n(),
            meanFixGen = mean(fixGen),
            CIFixGen = CI(fixGen))

ggplot(d_fix_n %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = model, y = nFixations,
           colour = model)) +
  #facet_nested(r_title + log10(r) ~ .) +
  geom_boxplot() +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  labs(x = "Model",
       y = "Number of fixations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "none") -> plt_numfix
plt_numfix
ggsave("plt_numfix.png", plt_numfix, dpi = 350,
       width = 6, height = 4, device = png)

# Age of fixations
ggplot(d_fix_n %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = model, y = meanFixGen,
           colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_boxplot() +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  labs(x = "Model",
       y = "Mean fixation generation",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "none")

# Mean frequency of segregating alleles
d_seg <- d_freqs %>% filter(is.na(fixGen))

d_seg_sum <- d_seg %>%
  group_by(optPerc, model, r, mutType) %>%
  summarise(meanFreq = mean(freq),
            CIFreq = CI(freq))

mutTypes_labeller <- as_labeller(c("2" = mutTypes_vec[1],
                                   "3" = mutTypes_vec[2],
                                   "4" = mutTypes_vec[3],
                                   "5" = mutTypes_vec[4],
                                   "6" = mutTypes_vec[5]),
                                 default = label_parsed)

ggplot(d_seg_sum %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = optPerc, y = meanFreq,
           colour = model, group = interaction(model, mutType))) +
  facet_nested(r_title + log10(r) ~ "Molecular component" + mutType,
               labeller = labeller(mutType = mutTypes_labeller)) +
  geom_line() +
  geom_ribbon(aes(ymin = meanFreq - CIFreq, ymax = meanFreq + CIFreq, fill = model),
              colour = NA, alpha = 0.2, show.legend = F) +
  scale_x_discrete(labels = c("0", "25", "50", "75", "100")) +
  labs(x = "Progress to the optimum (%)",
       y = "Mean allele frequency",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_meanfreq
plt_meanfreq
ggsave("plt_meanfreq.png", plt_meanfreq, dpi = 350,
       width = 12, height = 4.5, device = png)

# Combine
plot_grid(plt_meanfreq,
          plt_numfix,
          labels = "AUTO", align = "h",
          axis = "bt",
          rel_widths = c(2, 1))
ggsave("plt_alleles.png", dpi = 350, bg = "white",
       width = 12, height = 5, device = png)
