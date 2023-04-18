# plot molecular components vs phenotype w/ samples overlaid 
# as points connected by arrows (for time)
library(tidyverse)
library(latex2exp)
library(lattice)
library(latticeExtra)
library(paletteer)
library(reshape2)
library(scales)

library(GGally)
library(ggarrow)
library(ggnewscale)

# Load data
d_ode_phasemeasures <- readRDS("d_phaseDiagramMeasures.RDS")

d_ode_phasemeasures$id <- as.factor(paste(d_ode_phasemeasures$seed, d_ode_phasemeasures$modelindex, sep = "_"))

# Get range of each mol trait
molTraitRange <- list(aZ = quantile(d_ode_phasemeasures$aZ, c(0.01, 0.99)),
                      bZ = quantile(d_ode_phasemeasures$bZ, c(0.01, 0.99)),
                      KZ = quantile(d_ode_phasemeasures$KZ, c(0.01, 0.99)),
                      KXZ = quantile(d_ode_phasemeasures$KXZ, c(0.01, 0.99)))

# Filter data
# d_ode_phasemeasures <- d_ode_phasemeasures %>% 
#   filter(aZ >= molTraitRange$aZ[1] & aZ <= molTraitRange$aZ[2],
#          bZ >= molTraitRange$bZ[1] & bZ <= molTraitRange$bZ[2],
#          KZ >= molTraitRange$KZ[1] & KZ <= molTraitRange$KZ[2],
#          KXZ >= molTraitRange$KXZ[1] & KXZ <= molTraitRange$KXZ[2])

d_ode_phasemeasures$pheno_rescaled <- ifelse(d_ode_phasemeasures$phenomean > 3, 4, 
                                             d_ode_phasemeasures$phenomean)


# Line width for time - scale gens
d_ode_phasemeasures$gen_width <- rescale(d_ode_phasemeasures$gen, to = c(0.001, 1))

sampled_id <- sample(d_ode_phasemeasures$id, 1)

# Choose one with a huge KZ
sampled_id <- d_ode_phasemeasures[order(d_ode_phasemeasures$KZ, decreasing = T),]$id[1]

cc <- paletteer_c("ggthemes::Orange-Blue Diverging", 3)
cc <- c(cc[1], cc[3], cc[2])
cc2 <- paletteer_c("grDevices::Emrld", 50, -1)

d_ode_phasemeasures2 <- d_ode_phasemeasures
d_ode_phasemeasures <- d_ode_phasemeasures2
d_ode_phasemeasures <- d_ode_phasemeasures %>% filter(id %in% sampled_id)

# one with a low nloci and sigma
d_ode_phasemeasures <- d_ode_phasemeasures %>% filter(id == "35_30")

# # plot qg
# d_qg$id <- as.factor(paste(d_qg$seed, d_qg$modelindex, sep = "_"))
# d_qg$pheno_rescaled <- ifelse(d_qg$phenomean > 3, 4, d_qg$phenomean)
# d_qg$gen_width <- rescale(d_qg$gen, to = c(0.001, 1))
# 
# d_ode_phasemeasures <- d_qg %>% filter(id %in% sampled_id)


plotPairwiseScatter <- function(x, y, labels) {
  ggplot(d_ode_phasemeasures %>% mutate(KZ = log10(KZ)), aes(x = .data[[x]], y = .data[[y]], colour = pheno_rescaled)) +
    geom_point() +
    scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                           limits = c(0, 4),
                           labels = c(seq(0, 2, 1), "3"),
                           breaks = seq(0, 3, 1)
    ) +
    labs(colour = "Phenotype (Z)") +
    
    new_scale_colour() +
    geom_arrow_segment(data = d_ode_phasemeasures %>% mutate(KZ = log10(KZ)) %>% filter(id %in% sampled_id), 
                       mapping = aes(x = lag(.data[[x]]), y = lag(.data[[y]]), 
                                     xend = .data[[x]], yend = .data[[y]], 
                                     group = id, linewidth_head = gen_width, 
                                     linewidth_fins = gen_width * 0.8,
                                     colour = gen_width), 
                       arrow_head = arrow_head_line()) +
    scale_colour_gradientn(colors = cc2, labels = c(0, 0.25*2500, 0.5*2500, 0.75*2500, 2500)) +
    scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2], colour = "Generation") +
    theme_bw() + 
    theme(legend.position = "bottom") +
    guides(colour=guide_colourbar(barwidth=20))
}

aZbZScatter <- plotPairwiseScatter("aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
aZKZScatter <- plotPairwiseScatter("aZ", "KZ", c(TeX("$\\alpha_Z$"), TeX("$K_Z$")))
aZKXZScatter <- plotPairwiseScatter("aZ", "KXZ", c(TeX("$\\alpha_Z$"), TeX("$K_{XZ}$")))
bZKZScatter <- plotPairwiseScatter("bZ", "KZ", c(TeX("$\\beta_Z$"), TeX("$K_Z$")))
bZKXZScatter <- plotPairwiseScatter("bZ", "KXZ", c(TeX("$\\beta_Z$"), TeX("$K_{XZ}$")))
KZKXZScatter <- plotPairwiseScatter("KZ", "KXZ", c(TeX("$K_Z$"), TeX("$K_{XZ}$")))

plot_list <- list(ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = aZ)), 
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = aZ, y = bZ)),
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = aZ, y = KZ)),
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = aZ, y = KXZ)),
                  
                  aZbZScatter, 
                  ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = bZ)), 
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = bZ, y = KZ)),
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = bZ, y = KXZ)),
                  
                  aZKZScatter,
                  bZKZScatter,
                  ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = KZ)), 
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = KZ, y = KXZ)),
                  
                  aZKXZScatter,
                  bZKXZScatter,
                  KZKXZScatter,
                  ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = KXZ)))

xlabs <- c(TeX("$\\alpha_Z$", output = "character"), TeX("$\\beta_Z$", output = "character"), 
           TeX("$log_{10}(K_Z)$", output = "character"), TeX("$K_{XZ}$", output = "character"))

ggmatrix(plot_list, nrow = 4, ncol = 4, xAxisLabels = xlabs, yAxisLabels = xlabs,
         progress = T, byrow = T, labeller = "label_parsed", 
         legend = grab_legend(plot_list[[5]])) + theme_bw() + 
  theme(legend.position = "bottom") -> pair_mat

ggsave("molTrait_landscape_singleWalk_largesigma_smallnloci.png", pair_mat, width = 11, height = 8)

ggpairs(d_ode_phasemeasures, columns = c(6:9))

# Why is KZ so big?
d_com <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/d_combined_after.RDS")
d_com$id <- as.factor(paste(d_com$seed, d_com$modelindex, sep = "_"))
d_com_sampled <- d_com %>% filter(id %in% sampled_id)
View(d_com_sampled %>% filter(mutType == 5))
View(d_com_sampled %>% filter(mutType == 5, id == "13_32", gen == 49500))
d_com_sampled %>% filter(mutType == 6, id == "13_32", gen == 49500) %>% summarise(value = sum(value))

# How often is KZ bigger than the other values?
d_ode_phasemeasures %>% 
  mutate(KZ_bigger = (2 * KZ > aZ) & (2 * KZ > bZ), (2 * KZ > KXZ)) %>%
  group_by(gen) %>%
  summarise(mean(KZ_bigger)) -> KZ_bigger

# Lets look at that over burn-in
d_qg <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/checkpoint/d_qg.RDS")
d_qg %>%
  group_by(seed, nloci, sigma) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  mutate(phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  ungroup() -> d_qg

d_qg %>% mutate(KZ_bigger = (2 * KZ > aZ) & (2 * KZ > bZ), (2 * KZ > KXZ)) %>%
  group_by(gen) %>%
  summarise(mean(KZ_bigger)) -> KZ_bigger

plot(KZ_bigger$gen, KZ_bigger$`mean(KZ_bigger)`, type = "b")

d_pheno <- d_ode_phasemeasures %>% mutate(KZ = log10(KZ)) %>%
  pivot_longer(cols = c(phenomean, aZ, bZ, KZ, KXZ), names_to = "trait", values_to = "value")

library(cowplot)

molTrait_names <- c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), 
                    TeX("$K_{XZ}$"), TeX("$log_{10}(K_Z)$"), TeX("Z"))

ggplot(d_pheno %>% filter(gen > 49000, id %in% sampled_id) %>% mutate(gen = gen - 50000),
       aes(x = gen, y = value, colour = trait)) +
  geom_line() +
  scale_colour_paletteer_d("ggsci::nrc_npg", 1, labels = molTrait_names) +
  geom_point() +
  labs(x = "Generations post-optimum shift", y = "Mean trait/\ncomponent value",
       colour = "Trait/component") +
  theme_bw() + 
  theme(text = element_text(size = 16), legend.position = "bottom") -> singlewalk2_pheno

plot_grid(singlewalk2_pheno + coord_cartesian(ylim = c(0, 3)) + theme(legend.position = "none"), 
          singlewalk2_pheno + 
            guides(colour = "none"),
          nrow = 2) -> pheno_grid

plot_grid(pheno_grid, get_legend(singlewalk2_pheno), ncol = 1, rel_heights = c(1, 0.1))


ggsave("molTrait_phenowalk2.png", width = 7, height = 6, bg = "white")
# Evidence for KXZ neutral evolution between 0 - 1
