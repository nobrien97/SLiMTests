source("wrangle_data.R")
# Maladapted figures
# phenomean
ggplot(d_maladapted_sum %>% filter(gen > 49000) %>% mutate(gen = gen - 50000), 
       aes(x = gen, y = meanPheno, colour = modelindex)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno, 
                  fill = modelindex), alpha = 0.2, colour = NA
  ) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  scale_fill_paletteer_d("ggsci::nrc_npg", guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean population phenotype",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenomean
plt_phenomean
ggsave("phenomean_maladapted.png", plt_phenomean, png)

# fixed effects
## Additive
ggplot(d_fix_maladapted %>% filter(modelindex == 1), 
       aes(x = value)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Phenotypic effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_phenotype_mal_additive
dfe_phenotype_mal_additive

ggsave("dfe_phenotype_mal_additive.png", dfe_phenotype_mal_additive, png)


ggplot(d_fix_mal_add %>% filter(modelindex == 1), 
       aes(x = avFit)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Fitness effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_fitness_mal_additive
dfe_fitness_mal_additive

ggsave("dfe_fitness_mal_additive.png", dfe_fitness_mal_additive, png)

## NAR
ggplot(d_fix_mal_nar, 
       aes(x = avFX, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Phenotypic effect", y = "Density",
       colour = "Molecular component") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5), 
        legend.position = "bottom") -> dfe_phenotype_mal_nar
dfe_phenotype_mal_nar
ggsave("dfe_phenotype_mal_nar.png", dfe_phenotype_mal_nar, png)

ggplot(d_fix_mal_nar, 
       aes(x = avFit, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Fitness effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5),
        legend.position = "none") -> dfe_fitness_mal_nar
dfe_fitness_mal_nar
ggsave("dfe_fitness_mal_nar.png", dfe_fitness_mal_nar, png)

# Combine the above figures
leg <- get_legend(dfe_phenotype_mal_nar)

fig3_mal <- plot_grid(dfe_phenotype_mal_additive,
                      dfe_phenotype_mal_nar + theme(legend.position = "none"),
                      dfe_fitness_mal_additive,
                      dfe_fitness_mal_nar + theme(legend.position = "none"),
                      ncol = 2, nrow = 2, labels = "AUTO")

fig3_mal <- plot_grid(fig3_mal, get_legend(dfe_phenotype_mal_nar), 
                      ncol = 1, rel_heights = c(1, 0.1))
fig3_mal
ggsave("dfe_combined_mal.png", fig3_mal, png, bg = "white")

# Adaptive walks - these are really maladapted, not even close
d_qg_maladapting <- d_maladapted %>% filter(gen >= 50000) %>%
  group_by(seed, modelindex) %>%
  filter(any(gen >= 59800 & w <= 0.95)) %>% 
  ungroup()

d_qg_maladapting$gen_width <- scales::rescale(d_qg_maladapting$gen, to = c(0.1, 1))

seed <- sample(0:.Machine$integer.max, 1)
set.seed(seed)
#set.seed(1269262162)
sampled_seed <- sample(d_qg_maladapting[d_qg_maladapting$modelindex == 2,]$seed, 2)
walk <- plotPairwiseScatter(d_qg_maladapting %>% filter(modelindex == 2, seed %in% sampled_seed), 
                            "aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
suppressWarnings(walk)
ggsave("example_maladaptiveWalk.png", walk, png)


# Ranked fixed effects
d_rank_av_nar_mal <- d_fix_mal_ranked %>%
  group_by(rank) %>%
  summarise(CIFit = CI(avFit),
            meanFit = mean(avFit),
            CIPheno = CI(phenomean),
            meanPheno = mean(phenomean))

ggplot(d_rank_av_nar_mal %>% filter(rank > 0), aes(x = rank, y = meanFit)) +
  geom_point(position = position_dodge(1)) +
  geom_line(position = position_dodge(1)) +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4, position = position_dodge(1)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Mutation order", y = "Average effect on fitness", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16))

d_rank_av_add_mal <- d_fix_ranked_add_mal %>%
  group_by(rank) %>%
  summarise(CIFit = CI(avFit),
            meanFit = mean(avFit),
            CIPheno = CI(phenomean),
            meanPheno = mean(phenomean))

ggplot(d_rank_av_mal_add %>% filter(rank > 0), aes(x = rank, y = meanFit)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4) +
  labs(x = "Mutation order", y = "Average effect on fitness") +
  theme_bw() +
  theme(text = element_text(size = 16))

# diminishing returns of each step in the walk
ggplot(d_rank_av_nar_mal, aes(x = rank, y = meanPheno)) +
  geom_point(position = position_dodge(1)) +
  geom_line(position = position_dodge(1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4, position = position_dodge(1)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Mutation step", y = "Mean population phenotype", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_rank_av_add_mal, aes(x = rank, y = meanPheno)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4) +
  labs(x = "Mutation step", y = "Mean population phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16))
