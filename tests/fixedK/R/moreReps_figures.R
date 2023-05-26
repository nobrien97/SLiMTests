library(tidyverse)
library(paletteer)
library(latex2exp)
library(GGally)
library(cowplot)
library(ggnewscale)
library(ggalt)
library(ggarrow)

setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("wrangle_data.R")

# Fig 2 - phenotype mean
ggplot(d_adapted_sum %>% filter(gen > 49000) %>% mutate(gen = gen - 50000), 
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
ggsave("phenomean.png", plt_phenomean, png)

d_qg %>% group_by(modelindex) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted))


# Fig 3 - effect sizes

## Additive
ggplot(d_fix_adapted %>% filter(modelindex == 1), 
       aes(x = value)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Phenotypic effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_phenotype_additive
dfe_phenotype_additive

ggsave("dfe_phenotype_additive.png", dfe_phenotype_additive, png)


ggplot(d_fix_add %>% filter(modelindex == 1), 
       aes(x = avFit)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Fitness effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_fitness_additive
dfe_fitness_additive

ggsave("dfe_fitness_additive.png", dfe_fitness_additive, png)



ggplot(d_fix_nar, 
       aes(x = avFX, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Phenotypic effect", y = "Density",
       colour = "Molecular component") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5), 
        legend.position = "bottom") -> dfe_phenotype_nar
dfe_phenotype_nar
ggsave("dfe_phenotype_nar.png", dfe_phenotype_nar, png)

ggplot(d_fix_nar, 
       aes(x = avFit, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Fitness effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5),
        legend.position = "none") -> dfe_fitness_nar
dfe_fitness_nar
ggsave("dfe_fitness_nar.png", dfe_fitness_nar, png)

# Combine the above figures
leg <- get_legend(dfe_phenotype_nar)

fig3 <- plot_grid(dfe_phenotype_additive,
          dfe_phenotype_nar + theme(legend.position = "none"),
          dfe_fitness_additive,
          dfe_fitness_nar + theme(legend.position = "none"),
          ncol = 2, nrow = 2, labels = "AUTO")

fig3 <- plot_grid(fig3, get_legend(dfe_phenotype_nar), ncol = 1, rel_heights = c(1, 0.1))
fig3
ggsave("dfe_combined.png", fig3, png, bg = "white")

# Fig 4: Plot adaptive walks
cc2 <- paletteer_c("grDevices::Blues", 50, -1)

plotPairwiseScatter <- function(dat, x, y, labels) {
  GRID_RES <- 200
  d_grid <- expand.grid(seq(from = min(dat[[x]]), 
                                   to = max(dat[[x]]), length.out = GRID_RES), 
                        seq(from = min(dat[[y]]), 
                            to = max(dat[[y]]), length.out = GRID_RES), 1, 1)
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
  ggplot(dat %>% filter(modelindex == 2), aes(x = .data[[x]], y = .data[[y]], 
                                              group = seed)) +
    geom_raster(data = d_landscape, mapping = 
                     aes(x = .data[[x]], y = .data[[y]], fill = fitness)) +
    geom_point() +
    geom_encircle(s_shape = 0.2, expand = 0.2, colour = "black", mapping = aes(group = seed)) +
    scale_fill_gradientn(colors = c(cc[1], cc), limits = c(minFit, 1),
                         values = wValues) +
    labs(fill = "Log fitness (w)") +
    
    new_scale_colour() +
    geom_arrow_segment(data = dat %>% filter(seed %in% sampled_seed),
                       mapping = aes(x = lag(.data[[x]]), y = lag(.data[[y]]),
                                     xend = .data[[x]], yend = .data[[y]],
                                     group = seed, linewidth_head = gen_width,
                                     linewidth_fins = gen_width * 0.8,
                                     colour = gen_width),
                       arrow_head = arrow_head_line()) +
    scale_colour_gradientn(colors = cc2, labels = scales::comma(c(0, 0.25*10000, 0.5*10000,
                                                    0.75*10000, 10000))) +
    scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2], colour = "Generation") +
    theme_bw() + 
    theme(legend.position = "bottom", text = element_text(size = 14)) +
    guides(colour=guide_colourbar(barwidth=10),
           fill = guide_colourbar(barwidth = 10))
  )
}

d_qg_adapting <- d_adapted %>% filter(gen >= 50000)
d_qg_adapting$gen_width <- scales::rescale(d_qg_adapting$gen, to = c(0.1, 1))

seed <- sample(0:.Machine$integer.max, 1)
set.seed(seed)
#set.seed(1269262162)
sampled_seed <- sample(d_qg_adapting[d_qg_adapting$modelindex == 2,]$seed, 2)
walk <- plotPairwiseScatter(d_qg_adapting %>% filter(modelindex == 2, seed %in% sampled_seed), 
                    "aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
suppressWarnings(walk)
ggsave("example_adaptiveWalk.png", walk, png)

# Ranked fixed effects
d_rank_av_nar <- d_fix_ranked %>%
  group_by(rank) %>%
  summarise(CIFit = CI(avFit),
            meanFit = mean(avFit),
            CIPheno = CI(phenomean),
            meanPheno = mean(phenomean))

ggplot(d_rank_av_nar %>% filter(rank > 0), aes(x = rank, y = meanFit)) +
  geom_point(position = position_dodge(1)) +
  geom_line(position = position_dodge(1)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, by = 0.005)) +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4, position = position_dodge(1)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Mean fitness effect", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize
plt_adaptivestepsize

d_rank_av_add <- d_fix_ranked_add %>%
  group_by(rank) %>%
  summarise(CIFit = CI(avFit),
            meanFit = mean(avFit),
            CIPheno = CI(phenomean),
            meanPheno = mean(phenomean))

ggplot(d_rank_av_add %>% filter(rank > 0), aes(x = rank, y = meanFit)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, by = 0.005)) +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4) +
  labs(x = "Adaptive step", y = "Mean fitness effect") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize_add
plt_adaptivestepsize_add

# diminishing returns of each step in the walk
ggplot(d_rank_av_nar, aes(x = rank, y = meanPheno)) +
  geom_point(position = position_dodge(1)) +
  geom_line(position = position_dodge(1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.9, 2.1), breaks = seq(1, 2, by = 0.25)) +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4, position = position_dodge(1)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Mean phenotype", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno

ggplot(d_rank_av_add, aes(x = rank, y = meanPheno)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.9, 2.1), breaks = seq(1, 2, by = 0.25)) +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4) +
  labs(x = "Adaptive step", y = "Mean phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno_add

plot_grid(plt_adaptivestepsize_add,
          plt_adaptivestepsize,
          plt_phenostepPheno_add,
          plt_phenostepPheno, 
          nrow = 2, labels = "AUTO") -> plt_steps

ggsave("stepsize_combined.png", plt_steps, png, bg = "white")

# Box plots
ggplot(d_fix_ranked %>% filter(rank > 0), aes(x = as.factor(rank), y = avFit)) +
  geom_boxplot() +
  #scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, by = 0.005)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Fitness effect", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize_bp
plt_adaptivestepsize_bp

ggplot(d_fix_ranked_add %>% filter(rank > 0), aes(x = as.factor(rank), y = avFit)) +
  geom_boxplot() +
  #scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, by = 0.005)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Fitness effect", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize_add_bp
plt_adaptivestepsize_add_bp

# diminishing returns of each step in the walk
ggplot(d_fix_ranked, aes(x = as.factor(rank), y = phenomean)) +
  geom_boxplot() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.64, 2.6), breaks = seq(1, 2.5, by = 0.5)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Mean phenotype", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno_bp
plt_phenostepPheno_bp

ggplot(d_fix_ranked_add, aes(x = as.factor(rank), y = phenomean)) +
  geom_boxplot() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.64, 2.6), breaks = seq(1, 2.5, by = 0.5)) +
  labs(x = "Adaptive step", y = "Mean phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno_add_bp
plt_phenostepPheno_add_bp

plot_grid(plt_adaptivestepsize_add_bp,
          plt_adaptivestepsize_bp,
          plt_phenostepPheno_add_bp,
          plt_phenostepPheno_bp, 
          nrow = 2, labels = "AUTO") -> plt_steps_bp
plt_steps_bp
ggsave("stepsize_combined_boxplot.png", plt_steps_bp, png, bg = "white")

# Cumulative fitness effects for a few replicates
d_com_nar_sample %>% filter(gen > 49000) %>%
  group_by(gen, seed) %>%
  mutate(sumFit = sum(avFit)) -> d_com_nar_sample

d_com_add_sample %>% filter(gen > 49000) %>%
  group_by(gen, seed) %>%
  mutate(sumFit = sum(avFit)) -> d_com_add_sample

d_com_add_sample %>%
  ungroup(gen) %>%
  distinct(gen, .keep_all = T) %>%
  mutate(cumSumFit = cumsum(sumFit)) -> d_add_sample_sum

ggplot(d_com_sample_sum %>% distinct(), aes(x = gen, y = cumSumFit, colour = seed)) +
  geom_line() +
  theme_bw()

ggplot(d_add_sample_sum %>% distinct(), aes(x = gen, y = cumSumFit, colour = seed)) +
  geom_line() +
  theme_bw()
