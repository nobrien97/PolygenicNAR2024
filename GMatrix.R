# G matrix and additive variance figures/analysis

# Load VA for Z
d_h2 <- read.csv(paste0(DATA_PATH, "d_VA_Z.csv"))
d_h2 <- d_h2[,-1]

# summarise
d_h2_sum <- d_h2 %>% 
  filter(r %in% r_subsample) %>%
  group_by(optPerc, model, tau, r, method, isAdapted) %>%
  summarise(meanH2Z = mean(h2_Z, na.rm = T),
                   seH2Z = se(h2_Z, na.rm = T),
                   meanVAZ = mean(VA_Z, na.rm = T),
                   seVAZ = se(VA_Z, na.rm = T))
d_h2_sum$model <- as.factor(d_h2_sum$model)

# Additive variance
# Small effects as separate figure
ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
ggsave("plt_va_sml.png", device = png, bg = "white",
       width = 560*4, height = (980*4)/3, units = "px")

# Larger effects
ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_sml

ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 0.7)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_med

ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 1.25),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 1.25),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_lrg

leg <- get_legend(plt_add_va_lrg)

plt_add_va <- plot_grid(plt_add_va_sml + theme(legend.position = "none"),
                        plt_add_va_med + theme(legend.position = "none"),
                        plt_add_va_lrg + theme(legend.position = "none"),
                        ncol = 1, labels = "AUTO")

plt_add_va <- plot_grid(plt_add_va,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_add_va
ggsave("plt_va.png", device = png, bg = "white",
       width = 560*4, height = 980*4, units = "px")

# Infinitesimal expects zero change in additive variance due to selection
# So see how much it changes between timepoints
# Scale by the total variance as well -> a large effect model will produce a lot
# of variance, so the differences are more likely to be greater
d_h2 %>%
  group_by(model, seed, tau, r, nloci, isAdapted, method) %>%
  filter(n() > 1) %>%
  summarise(totalDeltaVA = sum(diff(VA_Z))/sum(VA_Z)) -> d_h2_deltaVA

d_h2_deltaVA %>%
  group_by(model, tau, r, method, isAdapted) %>%
  summarise(meanDeltaVA = mean(totalDeltaVA, na.rm = T),
            seDeltaVA = se(totalDeltaVA, na.rm = T)) -> d_h2_deltaVA_sum

ggplot(d_h2_deltaVA %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = model, y = totalDeltaVA, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_deltaVA_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = model, y = meanDeltaVA, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Model", 
       y = TeX("Change in additive variance $(\\Delta V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           3, direction = -1),
                      guide = "none",
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggsave("plt_deltaVA.png", device = png, width = 9, height = 4)

###########################################################
# Molecular G comparisons on data w/o Z included
# molecular components scaled
d_h2 <- read.csv(paste0(DATA_PATH, "d_VA_noZ.csv"))

d_h2_sum <- d_h2 %>% 
  filter(r %in% r_subsample) %>%
  group_by(optPerc, model, tau, r, method, isAdapted) %>%
  summarise(meanH2a = mean(h2_a, na.rm = T),
                   seH2a = se(h2_a, na.rm = T),
                   meanVAa = mean(VA_a, na.rm = T),
                   seVAa = se(VA_a, na.rm = T))
d_h2_sum$model <- as.factor(d_h2_sum$model)

# G matrix analysis
d_h2$optPerc <- factor(d_h2$optPerc)

d_h2 %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2

# Separate into model indices
# each sublist is replicates of a model index
sourceCpp("./getCovarianceMatrices.cpp")
lapply(split_h2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices
lapply(split_h2, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex

# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat <- unlist(cov_matrices, recursive = F)

# get ids from the matrix
cov_matrix_modelindex <- GetMatrixIDs(split_h2)


# Distance between G matrices
# Analysis across all timepoints, timepoint doesn't affect the tree structure
sourceCpp("./distanceFunctions.cpp")

dist_matrix <- distanceMatrix(h2_mat)
colnames(dist_matrix) <- paste("Matrix", 1:nrow(dist_matrix))
rownames(dist_matrix) <- colnames(dist_matrix)

# Clustering based on distance
hc <- hclust(as.dist(dist_matrix), method="average")
plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 2 seems to be the best
# elbow plot
fviz_nbclust(dist_matrix, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc, main = "Power Euclidean distances between molecular G matrices", labels = F)
rect.hclust(hc, 2, border = 2)

clus <- cutree(hc, 2)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)

id <- rbindlist(cov_matrix_modelindex, 
                fill = T)
id$label <- as.character(1:nrow(id))
id$modelindex <- as.factor(id$modelindex)
id <- AddCombosToDF(id)
id$nloci_group <- "[4, 64)"
id$nloci_group[id$nloci >= 64 & id$nloci < 1024] <- "[64, 256]"
id$nloci_group[id$nloci == 1024] <- "[1024]"
id$nloci_group <- factor(id$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id$clus, id$model)
names(dimnames(tab)) <- c("cluster", "model")
tab <- as.data.frame(tab)

# Model describes the clustering
glm.clus <- glm(Freq~model,family=poisson(),data=tab)
summary(glm.clus)
report::report(glm.clus)

id %>% ungroup() %>%
  group_by(model, clus) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  mutate(prop = n/sum(n)) -> cluster_percs

phylo <- full_join(as.phylo(phylo), id, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(size = 2) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(shape=16, size = 5,
                                                   linetype = 0))) -> tree_full

tree_full
ggsave("plt_tree_gmatrix_full_noZ.png", device = png, bg = "white",
       width = 7/2, height = 9/2)


# Evolvability metrics
# First convert to nearest positive definite matrix
h2_pd <- lapply(h2_mat, function(x) {
  if (!is.positive.definite(x)) {as.matrix(nearPD(x)$mat)}
})


d_ecr <- CalcECRA(h2_pd, 
                  id, noZ = T)
d_ecr <- AddCombosToDF(d_ecr)

# Remove Z -> not included in this model
d_ecr <- d_ecr %>% select(-cev_Z)

# Outliers
lofscores <- lofactor(d_ecr$cev, 20)
threshold <- 1.5
outliers <- lofscores > threshold

plot(density(lofscores[lofscores < 4]))

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score")
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)
boxplot(d_ecr[!outliers,]$cev)

# filter out outliers
d_ecr <- d_ecr[!outliers,]

# Need to calculate cev means separately for the K- and K+ models
# K- shouldn't mean over cev_KZ and KXZ

d_ecr_sum <- d_ecr %>%
  group_by(optPerc, model, r) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean conditional evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = res, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = res_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean respondability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = aut, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = aut_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean autonomy",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = ev, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = ev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                      plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO", label_size = 12)

plt_evol <- plot_grid(plt_evol,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol_noZ.png", device = png, bg = "white",
       width = 10, height = 7)

# Compare K+ and K- among recombination rates 
# (nloci doesn't affect, neither does optPerc)
# Variance differs between groups, use gls to account for unequal variance
library(nlme)
summary(gls.cev <- gls(cev ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.cev)
report(gls.cev)

anova(gls.cev)

# Marginal means
library(emmeans)
library(xtable)
em.cev <- emmeans(gls.cev, ~ model * r)
pairs(em.cev, simple = "model")
pairs(em.cev, simple = "r")
plot(em.cev, comparisons = T)

xtable(em.cev)


# Format to a nice table
library(stargazer)
stargazer(gls.cev)

# Repeat for autonomy, respondability, evolvability
summary(gls.aut <- gls(aut ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))

anova(gls.aut)
plot(gls.aut)
em.aut <- emmeans(gls.aut, ~ model * r)
pairs(em.aut, simple = "model")
pairs(em.aut, simple = "r")
plot(em.aut) # one of the comparisons has negative length because error is so small (ODE 0.1)

xtable(em.aut)

stargazer(gls.aut)

summary(gls.res <- gls(res ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.res)
anova(gls.res)

em.res <- emmeans(gls.res, ~ model * r)
pairs(em.res, simple = "model")
pairs(em.res, simple = "r")
plot(em.res, comparisons = T)
stargazer(gls.res)

xtable(em.res)


summary(gls.ev <- gls(ev ~ model * as.factor(r), d_ecr, 
                      weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.ev)
anova(gls.ev)
em.ev <- emmeans(gls.ev, ~ model * r, type = "response")
pairs(em.ev, simple = "model")
pairs(em.ev, simple = "r")
plot(em.ev, comparisons = T)

stargazer(gls.ev)
xtable(em.ev)

##############################
# Comparing the shape of G matrices with PCA similarity

# Bootstrap PCA similarity factor test:
# sample two matrices at random within group, get PCA similarity
# compare mean similarity between groups
# data input: dataframe with ids and a column with the matrix

krz_in <- id %>%
  mutate(g = h2_pd,
                group = interaction(model, r))

# Remove null matrices (no nearest matrix found)
krz_in <- krz_in[!sapply(krz_in$g,is.null)];


# Bootstrap in ten parts for RAM reasons
# This is slow: uncomment to run, otherwise read in precalculated data
# Generate seeds
# # newseed <- sample(1:.Machine$integer.max, 10)
# # 1360932387 1900268993  991875895 1523108407  197897667  199526283 2070940443
# # 128221287 1383031956  970870370
# newseed <- c(1360932387, 1900268993,  991875895, 1523108407, 197897667, 199526283, 
#              2070940443, 128221287, 1383031956, 970870370)
# 
# bootPCASim <- vector(mode = "list", length = 10)
# 
# for (i in seq_along(newseed)) {
#   # Set seed
#   set.seed(newseed[i])
#   
#   # Run replicate
#   res <- mcreplicate::mc_replicate(1000, bootKrzCorFn(krz_in, "group", T))
#   bootPCASim[[i]] <- unnest(as.data.frame(t(res)), cols = everything())
# }
# 
# # Output list into combined df
# bootPCASim2 <- bind_rows(bootPCASim)
# 
# bootPCASim <- bootPCASim2 %>%
#   separate(group1, c("model1", "r1"), "\\.",
#            extra = "merge") %>%
#   separate(group2, c("model2", "r2"), "\\.",
#            extra = "merge") %>%
#   mutate(r1 = log10(as.numeric(r1)),
#          r2 = log10(as.numeric(r2)),
#          model1 = factor(model1, levels = c("ODE", "K")),
#          model2 = factor(model2, levels = c("ODE", "K"))) %>%
#   rename(PCASim = krzCor)

#saveRDS(bootPCASim, "d_bootPCASim.RDS")
bootPCASim <- readRDS(paste0(DATA_PATH, "d_bootPCASim.RDS"))

# Plot
# Split into three by model
bootPCASim <- bootPCASim %>%
  mutate(modelCombo = ifelse(model1 != model2, "Mix",
                             as.character(model1)),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)))

# recomb by modelCombo
bootPCASim_sum <- bootPCASim %>%
  group_by(r1, r2, modelCombo) %>%
  summarise(meanPCASim = mean(PCASim),
                   ciPCASim = CI(PCASim))


ggplot(bootPCASim_sum, aes(
  x = as.factor(r1), y = as.factor(r2)
)) +
  facet_nested(. ~ "Model comparison" + modelCombo,
               labeller = labeller(modelCombo = as_labeller(c("K" = "K+ vs K+",
                                                              "ODE" = "K- vs K-",
                                                              "Mix" = "K+ vs K-")))) + 
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 12), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_r_modelCombo_noZ.png", device = png, width = 7, height = 5)

# beta regression
# Distributions
ggplot(bootPCASim, 
       aes(x = modelCombo, y = PCASim)) +
  geom_quasirandom(dodge.width = 0.9)

# Definitely different between models, none are normally distributed

# Floating point error: clamp to 1
bootPCASim <- bootPCASim %>%
  mutate(PCASim = raster::clamp(PCASim, 0, 1))

# Run regression: this is slow! Uncomment to run, otherwise load in saved object
# br.pcasim <- betareg(PCASim ~ modelCombo * as.factor(rCombo), bootPCASim)

# Save output
# saveRDS(br.pcasim, "betareg_pcaSim_big.RDS")
br.pcasim <- readRDS(paste0(DATA_PATH, "betareg_pcaSim_big.RDS"))
summary(br.pcasim)
plot(br.pcasim)

car::Anova(br.pcasim, type = 3, test.statistic = "F")

em.pcasim <- emmeans(br.pcasim, ~modelCombo * rCombo)
pairs(em.pcasim, simple = "modelCombo")
pairs(em.pcasim, simple = "rCombo")
plot(em.pcasim, comparisons = T)
pwpp(em.pcasim, by = "modelCombo", type = "response")
emmip(br.pcasim,  ~ modelCombo | rCombo)
joint_tests(em.pcasim)
xtable(em.pcasim)
