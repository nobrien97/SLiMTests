library(tidyverse)
library(paletteer)
library(latex2exp)
library(GGally)
library(cowplot)
library(ggnewscale)
library(ggalt)
library(ggarrow)
library(gghalves)
library(ggstance)
library(ggridges)

setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("wrangle_data.R")

# Fig 1: phenomean and adaptive walk
# A - phenomean ridgeline plot
d_adapted_walk <- d_adapted %>% filter(gen > 49000)
breaks <- seq(min(d_adapted_walk$gen - 50000), max(d_adapted_walk$gen - 50000), by = 1000)

d_adapted_walk$gen_group <- breaks[findInterval(d_adapted_walk$gen - 50000, breaks, rightmost.closed = TRUE)]

ggplot(d_adapted_walk,
       aes(y = as.factor(gen_group), x = phenomean, fill = modelindex)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(y = "Generations post-optimum shift", x = "Phenotype mean", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_phenomean_dist
plt_phenomean_dist

# B: phenotype at each step
ggplot(d_fix_ranked_combined,
       aes(y = as.factor(rank), x = phenomean, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Phenotype mean", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_adaptivewalk_pheno_dist
plt_adaptivewalk_pheno_dist

leg <- get_legend(plt_phenomean_dist)

plot_grid(plt_phenomean_dist + theme(legend.position = "none"), 
          plt_adaptivewalk_pheno_dist, 
          ncol = 2,
          labels = "AUTO") -> plt_fig1

plot_grid(plt_fig1, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_fig1
ggsave("fig1.png", plt_fig1, device = png, bg = "white")

# supp fig 1 - adaptive step timing
ggplot(d_fix_ranked_combined %>% filter(rank > 0),
       aes(y = as.factor(rank), x = gen - 50000, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Generations post-optimum shift", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_adaptivestepgen_dist
plt_adaptivestepgen_dist
ggsave("s_fig1.png", plt_adaptivestepgen_dist, device = png)

# text: populations adapted
d_qg %>% group_by(modelindex) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted))

# fig 2 - distribution of fixations and step size over time
d_fix_ranked_combined$model <- if_else(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")

ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_distfixed
plt_distfixed


ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = s, y = as.factor(rank), fill = model)) +
  geom_density_ridges(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_distfixed_time
plt_distfixed_time

plot_grid(plt_distfixed, 
          plt_distfixed_time, 
          ncol = 2,
          labels = "AUTO") -> plt_fig2

plot_grid(plt_fig2, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_fig2
ggsave("fig2.png", plt_fig2, device = png, bg = "white")

# fig 3: space of possible mutations
source("mutationScreenExp.R")

# A: all possible mutations at each step
ggplot(mutExp_combined, aes(y = as.factor(rank), x = s, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_effectsizerandom_time
plt_effectsizerandom_time

# B: Proportion of mutations that are beneficial
ggplot(mutExp_sum_combined, aes(x = as.factor(rank), y = percBeneficial, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = percBeneficial - CIperc, 
                              ymax = percBeneficial + CIperc),
                width = 0.2) +
  geom_line(aes(group = model)) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Adaptive step", y = "Proportion of\nbeneficial mutations (s > 0)",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_propbeneficial
plt_propbeneficial

# C: Waiting time to a beneficial mutation
ggplot(mutExp_sum_combined %>% mutate(waitingTime = 1/(10000 * (9.1528*10^-6) * percBeneficial),
                                      CIWaitingTime_lower = 1/(10000 * (9.1528*10^-6) * (percBeneficial - CIperc)),
                                      CIWaitingTime_upper = 1/(10000 * (9.1528*10^-6) * (percBeneficial + CIperc))),
       aes(x = as.factor(rank), y = waitingTime, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = CIWaitingTime_lower, 
                              ymax = CIWaitingTime_upper),
                width = 0.2) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  geom_line(aes(group=model)) +
  labs(x = "Adaptive step", y = "Expected waiting time\nto beneficial mutation", colour = "Model") +
  theme_bw() + 
  theme(text = element_text(size = 12), legend.position = "none") -> plt_waitingtime
plt_waitingtime

plot_grid(plt_effectsizerandom_time, 
          plt_propbeneficial, 
          plt_waitingtime,
          nrow = 3,
          labels = "AUTO") -> plt_fig3

plot_grid(plt_fig3, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_fig3
plt_fig3
ggsave("fig3.png", plt_fig3, width = 4, height = 10, device = png, bg = "white")

# supp fig 2 - distribution of percentage beneficial
ggplot(mutExp_perc_combined, aes(y = as.factor(rank), x = percBeneficial, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Proportion of beneficial mutations (s > 0)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_percBeneficial_time
plt_percBeneficial_time

ggsave("s_fig2.png", plt_percBeneficial_time, device = png)

# supp fig 3 EVD fit and stats

# supp fig 4 - distribution of new beneficial mutations
ggplot(mutExp_combined %>% filter(s > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_evt1
plt_evt1

ggsave("s_fig4.png", plt_evt1, device = png)

# fig 4 - fitness landscape and aZbZ ratio
# A - fitness landscape
plotaZbZLandscape <- function(minVal, maxVal) {
  GRID_RES <- 400
  d_grid <- expand.grid(seq(from = minVal, 
                            to = maxVal, length.out = GRID_RES), 
                        seq(from = minVal, 
                            to = maxVal, length.out = GRID_RES), 1, 1)
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZ, y = bZ, fill = fitness)) +
      geom_tile() +
      #geom_abline(slope = 1/1.27) +
      # geom_point(data = d_landscape %>% mutate(aZbZ = aZ/bZ) 
      #            %>% filter(aZbZ > 1.25, aZbZ < 1.35), size = 0.1, shape = 4) +
      scale_fill_gradientn(colors = c(cc[1], cc), 
                           limits = c(minFit, 1),
                           values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\alpha_Z$"), y = TeX("$\\beta_Z$"), 
           fill = "Fitness (w)", size = TeX("$\\frac{\\alpha_Z}{\\beta_Z} ratio$")) +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        fill = guide_colourbar(barwidth = 20, title.vjust = 0.87)) # i love magic numbers
  )
}
plotaZbZLandscape(0, 3) -> plt_aZbZ_landscape
plt_aZbZ_landscape

# B - GPW map of aZbZ
plotRatioLandscape <- function(minRatio, maxRatio) {
  GRID_RES <- 50000
  rational <- MASS:::.rat(seq(minRatio, maxRatio, length.out = GRID_RES),
                          max.denominator = 20)$rat
  d_grid <- data.frame(aZ = rational[,1], 
                       bZ = rational[,2],
                       KZ = 1, 
                       KXZ = 1) %>% distinct()
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8) %>%
    filter(pheno > 0)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape  %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZbZ, y = pheno, colour = fitness)) +
      geom_point() +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                             limits = c(ifelse(minFit < 0.8, 0.8, minFit), 1),
                             values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\frac{\\alpha_Z}{\\beta_Z} ratio$"), y = "Phenotype", 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        colour = guide_colourbar(barwidth = 20, title.vjust = 0.87)) # i love magic numbers
  )
}

plotRatioLandscape(0.5, 3) -> plt_aZbZratio
plt_aZbZratio

# C - difference in evolution among alpha and beta
ggplot(d_molCompDiff,
       aes(x = molCompDiff)) +
  geom_density() +
  labs(x = "Molecular component\ncontribution", y = "Density") +
  theme_bw() + 
  theme(text = element_text(size = 16)) -> plt_molCompDiff
plt_molCompDiff

leg <- get_legend(plt_aZbZ_landscape)

plot_grid(plt_aZbZ_landscape + theme(legend.position = "none"), 
          plt_aZbZratio + theme(legend.position = "none"), 
          plt_molCompDiff,
          nrow = 3,
          labels = "AUTO") -> plt_fig4

plot_grid(plt_fig4, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_fig4
plt_fig4
ggsave("fig4.png", plt_fig4, width = 4, height = 12, device = png, bg = "white")
