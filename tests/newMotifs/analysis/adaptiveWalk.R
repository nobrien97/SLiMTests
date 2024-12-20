# Adaptive walk figures
library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
library(stargazer)

setwd("/mnt/c/GitHub/SLiMTests/tests/newMotifs/analysis")
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/pilot/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/analysis/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

# Cowplot 1.1.3 bug: won't get legend, this fixes
get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


d_combos <- read.table("../R/combos.csv", header = F,
                       col.names = c("model", "r"))

# load trait evolution data
d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 24)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", 
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var",
                                        "trait3_var", "trait4_var", "dist", "w", 
                                        "deltaPheno", "deltaw", "mc1_mean", 
                                        "mc2_mean", "mc3_mean", "mc4_mean", "mc5_mean",
                                        "mc6_mean", "mc7_mean", "mc8_mean", "mc9_mean",
                                        "mc10_mean", "mc11_mean"), 
                          fill = T)

# Add predictors
d_qg <- AddCombosToDF(d_qg) 

# Initial optimum distance
INIT_DIST <- sqrt(-2 * log(0.95))

# Optimum: fitness > 98%

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & w > 0.98)) %>%
  ungroup() -> d_qg

# Proportion of each model that adapted
d_prop_adapted <- d_qg %>% group_by(model, r) %>%
  filter(gen == 60000) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )

# Average over nloci
# How many adapted in low recombination scenarios?
d_prop_adapted <- d_qg %>% 
  group_by(model, tau) %>%
  filter(gen == 59950) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )
d_prop_adapted

# Output to table
stargazer(d_prop_adapted)

# Summarise phenotype trajectories
d_qg_sum <- d_qg %>% 
  filter(gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r) %>%
  summarise(meanDist = mean(dist),
            SEDist = se(dist),
            meanFitness = mean(w),
            SEFitness = sd(w))

# plot
ggplot(d_qg_sum,
       aes(x = gen, y = meanDist, colour = model)) +
  facet_grid(log10(r)~.) +
  geom_line() +
  #geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanDist - SEDist, 
                  ymax = meanDist + SEDist, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean distance from the optimum", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_adapt_smlfx.png", width = 12, height = 5, device = png)

# Time to adaptation
d_adaptTime <- d_qg %>% filter(gen >= 49500) %>%
  mutate(gen = gen - 50000,
         isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.0, 1.9)) %>%
  group_by(seed, model, nloci, tau, r) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1),
         initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1)) %>%
  distinct(seed, model, nloci, tau, r, .keep_all = T) %>%
  ungroup()

# Distribution of adaptation times for larger effect sizes
ggplot(d_adaptTime %>% filter(adaptTime > -1, 
                              r %in% r_subsample,
                              tau > 0.0125) %>%
         mutate(r_title = "Recombination rate (log10)",
                tau_title = "Mutational effect size variance"), 
       aes(x = adaptTime, fill = model), alpha = 0.4) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_x_continuous(labels = scales::comma) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  labs(x = "Time to adaptation (Generations)", y = "Density", fill = "Model")
ggsave("sfig_timetoadaptation_largetau.png", device = png, width = 10, height = 6)

# Average adaptation time in the larger effect size models
mean(d_adaptTime[d_adaptTime$tau > 0.0125,]$adaptTime)
CI(d_adaptTime[d_adaptTime$tau > 0.0125,]$adaptTime)

# Mean phenotype figures for large effect models
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
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
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
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
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
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
    labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
         colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> adapt_grid_lrgFX
  
  adapt_grid <- plot_grid(adapt_grid_smlFX, adapt_grid_medFX, adapt_grid_lrgFX,
                          nrow = 3)
  
  ggsave("adapt_grid.png", adapt_grid, width = 14, height = 30, device = png)
}