# Figures of adaptive walk

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- fread::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
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
                 nAdapted = sum(isAdapted),
                 pAdapted = mean(isAdapted)
                 )
     )

d_adapted_sum <- d_qg %>% 
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))
# Mean phenotype figures
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

# Within-population phenotype variance figures
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
  
  ggsave("adaptvar_grid.png", adaptvar_grid, width = 14, height = 30, device = png)
}

# Across seed phenotype plot
{
  ggplot(d_qg %>% filter(tau == 0.0125, isAdapted, gen >= 49500),
         aes(x = gen - 50000, y = phenomean, colour = model),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    geom_hline(yintercept = 2, linetype = "dashed") +
    ggtitle("Tau = 0.0125") +
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
    theme(legend.position = "bottom", text = element_text(size = 14))
}

# dpdt data
GENS_BETWEEN_SAMPLES <- 50
d_dpdt <- d_qg %>% filter(isAdapted, gen >= 49500) %>%
  mutate(dPdt = abs(deltaPheno - lag(deltaPheno)) / GENS_BETWEEN_SAMPLES) %>%
  select(gen, seed, modelindex, dPdt) %>%
  group_by(gen, modelindex) %>%
  summarise(meanDPDT = mean(dPdt),
            sdDPDT = sd(dPdt))
d_dpdt$modelindex <- as.factor(d_dpdt$modelindex)

# d_dpdt <- data.table::fread(paste0(DATA_PATH, "mutationStats/d_dpdt.csv"), header = F, 
#                   sep = ",", colClasses = c("character", "integer", "numeric", 
#                                             "numeric"), 
#                   col.names = c("optPerc", "modelindex", "meanDPDT", "sdDPDT"), 
#                   fill = T)
# d_dpdt$optPerc <- factor(d_dpdt$optPerc, 
#                          levels = c("[-Inf,0.25)", 
#                                     "[0.25,0.5)",
#                                     "[0.5,0.75)",
#                                     "[0.75,0.9)",
#                                     "[0.9, Inf)"))
d_dpdt <- AddCombosToDF(d_dpdt) 
d_dpdt$gen <- d_dpdt$gen - 50000

# dPdt plots
{
  ggplot(d_dpdt %>% filter(tau == 0.0125), 
         aes(x = gen, y = meanDPDT, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    ggtitle("Tau = 0.0125") +
    geom_ribbon(aes(ymin = meanDPDT - sdDPDT,
                    ymax = meanDPDT + sdDPDT, fill = model), colour = NA,
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
         y = "Mean per-generation change in phenotype", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14),
          panel.spacing.x = unit(2, "lines")) -> dpdt_grid_smlFX

  ggplot(d_dpdt %>% filter(tau == 0.125), 
         aes(x = gen, y = meanDPDT, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    ggtitle("Tau = 0.125") +
    geom_ribbon(aes(ymin = meanDPDT - sdDPDT, 
                    ymax = meanDPDT + sdDPDT, fill = model), colour = NA,
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
         y = "Mean per-generation change in phenotype", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14),
          panel.spacing.x = unit(2, "lines")) -> dpdt_grid_medFX
  
  ggplot(d_dpdt %>% filter(tau == 1.25), 
         aes(x = gen, y = meanDPDT, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    ggtitle("Tau = 1.25") +
    geom_ribbon(aes(ymin = meanDPDT - sdDPDT, 
                    ymax = meanDPDT + sdDPDT, fill = model), colour = NA,
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
         y = "Mean per-generation change in phenotype", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14),
          panel.spacing.x = unit(2, "lines")) -> dpdt_grid_lrgFX
  
  dpdt_grid <- plot_grid(dpdt_grid_smlFX, dpdt_grid_medFX, dpdt_grid_lrgFX,
                             nrow = 3)
  
  ggsave("dpdt_grid.png", dpdt_grid, width = 14, height = 30, device = png)
  
# Scale dpdt plots to better look at non-peak trajectories
dpdt_grid_scaled <- plot_grid(dpdt_grid_smlFX + coord_cartesian(ylim = c(0, 0.0025)), 
                              dpdt_grid_medFX + coord_cartesian(ylim = c(0, 0.0025)), 
                              dpdt_grid_lrgFX + coord_cartesian(ylim = c(0, 0.0025)),
                       nrow = 3)

ggsave("dpdt_grid_scaled.png", dpdt_grid_scaled, width = 14, height = 30, device = png)
}


