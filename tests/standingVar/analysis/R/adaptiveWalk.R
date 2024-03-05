# Figures of adaptive walk

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- read.table(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                   sep = ",", colClasses = c("integer", "factor", "factor", 
                                             rep("numeric", times = 12)), 
                   col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                 "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                 "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                   fill = T)

# Add predictors
d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg

View(d_qg %>% group_by(model, nloci, tau, r) %>%
       filter(gen == 59950) %>%
       summarise(n = n(),
                 pAdapted = mean(isAdapted),
                 CIAdapted = CI(isAdapted)))

d_adapted_sum <- d_qg %>% 
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))
{
ggplot(d_adapted_sum %>% filter(tau == 0.0125),
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~nloci) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 0.0125") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean, 
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  scale_fill_paletteer_d("nationalparkcolors::Badlands",
                         guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_smlFX

ggplot(d_adapted_sum %>% filter(tau == 0.125),
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~nloci) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 0.125") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean, 
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  scale_fill_paletteer_d("nationalparkcolors::Badlands",
                         guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_medFX

ggplot(d_adapted_sum %>% filter(tau == 1.25),
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~nloci) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 1.25") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean, 
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  scale_fill_paletteer_d("nationalparkcolors::Badlands",
                         guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_lrgFX

adapt_grid <- plot_grid(adapt_grid_smlFX, adapt_grid_medFX, adapt_grid_lrgFX,
          nrow = 3)

ggsave("adapt_grid.png", adapt_grid, width = 14, height = 30, device = png)
  }

{
  ggplot(d_adapted_sum %>% filter(tau == 0.0125),
         aes(x = gen, y = meanPhenovar, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    ggtitle("Tau = 0.0125") +
    geom_ribbon(aes(ymin = meanPhenovar - sdPhenovar, 
                    ymax = meanPhenovar + sdPhenovar, fill = model), colour = NA,
                alpha = 0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma, 
                       sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    scale_fill_paletteer_d("nationalparkcolors::Badlands",
                           guide = "none") +
    labs(x = "Generations post-optimum shift", 
         y = "Mean within-population\nphenotypic variance", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adaptvar_grid_smlFX
  
  ggplot(d_adapted_sum %>% filter(tau == 0.125),
         aes(x = gen, y = meanPhenovar, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    ggtitle("Tau = 0.125") +
    geom_ribbon(aes(ymin = meanPhenovar - sdPhenovar, 
                    ymax = meanPhenovar + sdPhenovar, fill = model), colour = NA,
                alpha = 0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma, 
                       sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    scale_fill_paletteer_d("nationalparkcolors::Badlands",
                           guide = "none") +
    labs(x = "Generations post-optimum shift", 
         y = "Mean within-population\nphenotypic variance", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adaptvar_grid_medFX
  
  ggplot(d_adapted_sum %>% filter(tau == 1.25),
         aes(x = gen, y = meanPhenovar, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    coord_cartesian(ylim = c(0, 0.5)) +
    ggtitle("Tau = 1.25") +
    geom_ribbon(aes(ymin = meanPhenovar - sdPhenovar, 
                    ymax = meanPhenovar + sdPhenovar, fill = model), colour = NA,
                alpha = 0.2) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma, 
                       sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    scale_fill_paletteer_d("nationalparkcolors::Badlands",
                           guide = "none") +
    labs(x = "Generations post-optimum shift", 
         y = "Mean within-population\nphenotypic variance", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adaptvar_grid_lrgFX
  
  adaptvar_grid <- plot_grid(adaptvar_grid_smlFX, adaptvar_grid_medFX, adaptvar_grid_lrgFX,
                          nrow = 3)
  
  ggsave("adapt_grid.png", adapt_grid, width = 14, height = 30, device = png)
}