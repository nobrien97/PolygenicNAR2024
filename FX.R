# Analysis of effect sizes for each molecular component

d_fx <- data.table::fread(paste0(DATA_PATH, "d_fx_new.csv"), header = F, 
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
ggplot(d_fx %>% filter(mutType != "5") %>%
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = s)) +
  facet_nested("Model" + model ~ "Mutation Type" + mutType) +
  geom_density() +
  labs(x = "Selection coefficient (s)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


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
  group_by(optPerc, seed, model, mutType, nloci, tau, r) %>%
  mutate(isBen = (s > 0)) %>%
  summarise(propBen = sum(isBen)/n()) %>%
  ungroup()

# Set muttype for additive to a unique id
d_fx_propBen[d_fx_propBen$model == "Add", "mutType"] <- "2"

d_fx_propBen_sum <- d_fx_propBen %>%
  group_by(model, mutType, r) %>%
  summarise(meanPropBen = mean(propBen),
            CIPropBen = CI(propBen))

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
             shape = 3, size = 2, position = position_dodge(0.9)) +
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
ggsave("plt_propben_muts_wholewalk.png", plt_propben_muts_wholewalk, 
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
br.benmut <- betareg(propBen_adj ~ mutType, d_fx_propBen)

summary(br.benmut)
plot(br.benmut)

saveRDS(br.benmut, "betareg_benmut.RDS")
br.benmut <- readRDS(paste0(DATA_PATH, "betareg_benmut.RDS"))

#anova(br.benmut)
#em.benmut <- emmeans(br.benmut, ~ modelMutType)


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

joint_tests(em.benmut.km)
joint_tests(em.benmut.kp)

xtable(em.benmut.km)
xtable(em.benmut.kp)

pairs(em.benmut, simple = "modelMutType")
plot(em.benmut, comparisons = T)
emmip(em.benmut,  ~ modelMutType)

# Differences between models (apart from KZ) range from <1% to ~5%
# KZ is about 50% change, almost never beneficial

# What about effect size: if K_Z has very large effects, generating that variation
# could be important
d_fx_ben$mutType <- as.character(d_fx_ben$mutType)
d_fx_ben[d_fx_ben$model == "Add", "mutType"] <- "2"
d_fx_ben$mutType <- factor(d_fx_ben$mutType)

d_fx_ben_sum <- d_fx_ben %>%
  group_by(model, mutType) %>%
  summarise(meanBen = mean(s),
            CIBen = CI(s))

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
             shape = 3, size = 2, position = position_dodge(0.9)) +
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
ggsave("plt_ben_muts_s.png", plt_ben_muts_s, width = 6, height = 4, device = png)

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

summary(gls.s.kp <- gls(s ~ mutType, (d_fx_ben %>% filter(model == "K")), 
                        weights = 
                          varIdent(form = ~ 1 | mutType)))

plot(gls.s.kp)

em.s.kp <- emmeans(gls.s.kp, ~ mutType)

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