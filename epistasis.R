# Have a look at the epistasis data to see which models differ the most

# data
d_epi_means <- read.table(paste0(DATA_PATH, "d_epi_mean.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "count"))


# Get pairwise differences between in mean epistasis between models at start
# vs end of simulation
d_epi_means_sbst <- d_epi_means %>% 
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)") %>% 
  distinct()

d_epi_means_plt <- AddCombosToDF(d_epi_means %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier, filter to tau = 0.0125
d_epi_means_plt <- d_epi_means_plt %>%
  filter(modelindex != 396, tau == 0.0125, r %in% r_subsample)

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