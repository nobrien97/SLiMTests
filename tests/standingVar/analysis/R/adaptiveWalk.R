# Figures of adaptive walk
library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
library(stargazer)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/analysis/R")
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
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


d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# d_combos <- read.table("~/tests/standingVar/R/combos.csv", header = F,
#                        col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
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

# We have many recombination rates: choose a few
r_subsample <- c(1e-10, 1e-5, 1e-1)


d_prop_adapted <- d_qg %>% group_by(model, nloci, tau, r) %>%
       filter(gen == 59950) %>%
       dplyr::summarise(n = n(),
                 nAdapted = sum(isAdapted),
                 pAdapted = mean(isAdapted)
                 )

# Average over nloci
# How many adapted in low recombination scenarios?
d_prop_adapted <- d_qg %>% 
  mutate(r_cut = cut(r, breaks = c(0, 1e-7, 1),
                      labels = c("Low", "High"))) %>%
  group_by(model, tau, r_cut) %>%
  filter(gen == 59950) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )
d_prop_adapted

stargazer(d_prop_adapted)
  

# Plot adapted populations - almost all populations adapted, except the small effect
# additive case with low recombination. High recombination allows the additive to 
# recover
ggplot(d_prop_adapted %>% mutate(r_title = "Recombination rate (log10)",
                                 nloci_title = "Number of loci",
                                 tau_title = "Mutational effect size variance"),
       aes(x = as.factor(tau), y = pAdapted, colour = model)) +
  facet_nested(r_title + log10(r) ~ nloci_title + nloci) +
  geom_point(position = position_dodge(0.7), size = 2) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Mutational effect size variance", 
       y = "Probability of adaptation", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom")
ggsave("pAdapted.png", width = 12, height = 10, device = png)


ggplot(d_prop_adapted %>% filter(tau == 0.0125) %>% 
         mutate(r_title = "Recombination rate (log10)",
                                 nloci_title = "Number of loci",
                                 tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")),
       aes(x = model, y = pAdapted, colour = model)) +
  facet_nested(r_title + log10(r) ~ nloci_title + nloci) +
  geom_point(position = position_dodge(0.7), size = 2) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", 
       y = "Probability of adaptation", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom")
ggsave("pAdapted_sml.png", width = 10, height = 10, device = png)


d_adapted_sum <- d_qg %>% 
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

# Small fx only
ggplot(d_adapted_sum %>% filter(tau == 0.0125, r %in% r_subsample),
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~nloci) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
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
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_adapt_smlfx.png", width = 12, height = 5, device = png)

# Time to adaptation
d_adaptTime <- d_qg %>% filter(gen >= 49500) %>%
  dplyr::mutate(gen = gen - 50000,
         isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.0, 1.9)) %>%
  group_by(seed, model, nloci, tau, r) %>%
  dplyr::mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1),
         initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1)) %>%
  distinct(seed, model, nloci, tau, r, .keep_all = T) %>%
  ungroup()

# Distribution of adaptation times
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

mean(d_adaptTime[d_adaptTime$tau > 0.0125,]$adaptTime)
CI(d_adaptTime[d_adaptTime$tau > 0.0125,]$adaptTime)



d_adaptTime %>%
filter(adaptTime != -1, initRespTime != -1) %>%
  group_by(model, nloci, tau, r) %>%
  dplyr::summarise(meanAdaptTime = mean(adaptTime), 
            CIAdaptTime = CI(adaptTime),
            meanInitRespTime = mean(initRespTime),
            CIInitRespTime = CI(initRespTime)) -> d_meanAdaptTime

ggplot(d_meanAdaptTime, 
       aes(x = model, y = meanAdaptTime, colour = as.factor(nloci))) + 
  facet_grid(tau~log10(r)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanAdaptTime - CIAdaptTime, ymax = meanAdaptTime + CIAdaptTime)) +
  labs(y = "Time to adaptation\n(generations)", x = "Model", 
       colour = "Number of loci") +
  theme_bw()

d_adaptTime %>%
  filter(adaptTime != -1, initRespTime != -1) %>%
  group_by(model, tau, r) %>%
  dplyr::summarise(meanAdaptTime = mean(adaptTime), 
            CIAdaptTime = CI(adaptTime),
            meanInitRespTime = mean(initRespTime),
            CIInitRespTime = CI(initRespTime)) -> d_meanAdaptTime

View(d_meanAdaptTime %>% filter(tau > 0.0125))

# Improvement with increasing recombination



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
  
  # variance: are the spikes in K+ tau = 1.25 consistent or random among models?
  ggplot(d_qg %>% AddCombosToDF(.) %>% filter(gen >= 49500, tau == 1.25, model == "K"),
         aes(x = gen - 50000, y = phenovar, colour = as.factor(seed)),
         group = as.factor(seed)) +
    facet_grid(log10(r)~nloci) +
    geom_line() +
    ggtitle("Tau = 1.25") +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(labels = scales::comma, 
                       sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_discrete(guide = "none") +
    labs(x = "Generations post-optimum shift", 
         y = "Mean within-population\nphenotypic variance", 
         colour = "Replicate") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> phenovar_K

  phenovar_K
  phenovar_K + coord_cartesian(ylim = c(0, 5))
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


# allele frequency spectrum
d_SFS <- data.table::fread(paste0(DATA_PATH, "mutationStats/d_SFS.csv"), header = F, 
                          sep = ",", colClasses = c("factor", "factor", "numeric", 
                                                    "factor",
                                                    rep("numeric", times = 3)), 
                          col.names = c("optPerc", "modelindex", "mutType", "freqBin",
                                        "countFreqBin", "meanValue", "sdValue"), 
                          fill = T)

d_SFS <- AddCombosToDF(d_SFS)

freqBins <- seq(from = 0.1, to = 1, by = 0.1)
d_SFS <- d_SFS %>%
  mutate(freqBin = freqBins[as.numeric(freqBin)])

# change in SFS over the walk
d_deltaSFS <- d_SFS %>%
  group_by(mutType, freqBin, model, nloci, tau, r) %>%
  summarise(deltaCount = sum(diff(countFreqBin)))

ggplot(d_deltaSFS %>% 
         mutate(freqBin = freqBin - 0.1) %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE"),
                mutType = as.factor(mutType),
                mutType = fct_recode(mutType, "$\\alpha_Z$" = "3",
                                     "$\\beta_Z$" = "4",
                                     "K_Z" = "5",
                                     "K_{XZ}" = "6")) %>%
         filter(tau == 0.0125, 
                nloci == 1024, r %in% r_subsample) %>%
         uncount(countFreqBin),
       aes(x = freqBin, y = deltaCount, fill = model)) +
  facet_nested(log10(r) ~ model + mutType) +
  geom_col(position = position_nudge(x = 0.05)) +
  #stat_binline(bins = 10, binwidth = 0.1, position = position_nudge(x = 0.05), 
  #             scale = 0.95) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_y_discrete(labels = c("25%", "50%", "75%", "100%")) +
  labs(x = "Allele frequency", y = "Progress to the optimum", 
       fill = "Model") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))


ggplot(d_SFS %>% 
         mutate(freqBin = freqBin - 0.1) %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE"),
                mutType = as.factor(mutType),
                mutType = fct_recode(mutType, "$\\alpha_Z$" = "3",
                                     "$\\beta_Z$" = "4",
                                     "K_Z" = "5",
                                     "K_{XZ}" = "6")) %>%
         filter(tau == 0.0125, 
                nloci == 1024, r %in% r_subsample) %>%
         uncount(countFreqBin),
       aes(x = freqBin, y = optPerc, fill = model)) +
  facet_nested(log10(r) ~ model + mutType) +
  stat_binline(bins = 10, binwidth = 0.1, position = position_nudge(x = 0.05), 
                      scale = 0.95) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_y_discrete(labels = c("25%", "50%", "75%", "100%")) +
  labs(x = "Allele frequency", y = "Progress to the optimum", 
       fill = "Model") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12))
ggsave("plt_sfs.png", device = png, width = 10, height = 4)

# Heterozygosity
# add header row with maximum number of columns
head_name <- c("gen","seed",c("gen", "seed", "modelindex", paste0("Ho_", 1:1024)))
write.csv(het_header, "het_header")


# load in observed heterozygosity data
d_het <- data.table::fread(paste0(DATA_PATH, "slim_locusHo_hdr.csv"), 
                    header = T, sep = ",", 
                    colClasses = c("integer", "factor", "factor", 
                                   rep("numeric", times = 1024)), 
                    col.names = c("gen", "seed", "modelindex", paste0("Ho_", 1:1024)), 
                    fill = T)

d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))
d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, optPerc) %>% filter(gen >= 49500)

d_het <- left_join(d_het, d_qg_optPerc, by = c("gen", "seed", "modelindex"))

d_het <- AddCombosToDF(d_het)

# Summarise, mean H_O over time
d_het %>%
  group_by(optPerc, seed, model, nloci, tau, r) %>%
  summarise(meanHo = mean(cbind(select(starts_with("Ho"))), na.rm = T),
            CIHo = CI(cbind(select(starts_with("Ho"))), na.rm = T)) -> d_Ho_sum

# Mean heterozygosity across all contributing loci
d_het %>%
  rowwise() %>%
  mutate(meanHo = mean(c_across(4:1027), na.rm = T),
         CIHo = CI(c_across(4:1027), na.rm = T)) %>%
  select(1:3, meanHo, CIHo) -> d_Ho_sum

# Define a custom function to compute mean and confidence interval
compute_stats <- function(row) {
  mean_val <- mean(row, na.rm = TRUE)
  ci_val <- CI(row, na.rm = TRUE)
  return(list(meanHo = mean_val, CIHo = ci_val))
}

setDT(d_het)

# Apply the function to each row and store the results in two new columns
result <- d_het[, .(meanHo = mean(unlist(.SD), na.rm = TRUE), 
                    CIHo = CI(unlist(.SD), na.rm = TRUE)), 
                .SDcols = 4:1027, by = 1:nrow(d_het)]

# Combine the results with the first 3 columns
d_Ho_sum <- cbind(d_het[, 1:3], result[, -1, with = FALSE])

# Optionally convert back to a dataframe
d_Ho_sum <- as.data.frame(d_Ho_sum)

write.table(d_Ho_sum, "d_Ho_sum.csv", sep = ",", col.names = F, row.names = F)

d_Ho_sum <- data.table::fread(paste0(DATA_PATH, "d_Ho_sum.csv"), 
                           header = F, sep = ",", 
                           colClasses = c("integer", "factor", "factor", 
                                          "numeric", "numeric"), 
                           col.names = c("gen", "seed", "modelindex", "meanHo", "CIH"), 
                           fill = T)

d_Ho_sum <- left_join(d_Ho_sum, d_qg_optPerc, by = c("gen", "seed", "modelindex"))
d_Ho_sum <- AddCombosToDF(d_Ho_sum)


# Distribution
ggplot(d_Ho_sum %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = optPerc, y = meanHo, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Heterozygosity $(H)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_Ho_sml

ggplot(d_qg %>% 
         mutate(gen = gen - 50000) %>%
         filter(gen > -1000) %>%
         group_by(gen, model, tau, r, nloci) %>%
         summarise(CIH = CI(meanH),
                   meanH = mean(meanH)) %>%
         filter(tau == 0.0125) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
         aes(x = gen, y = meanH, colour = model)) +
  facet_nested(r_title + log10(r) ~ nloci_title + nloci) +
  geom_line() +
  geom_ribbon(aes(ymin = meanH - CIH, ymax = meanH + CIH, fill = model), 
              alpha = 0.2, colour = NA) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post optimum shift", 
       y = TeX("Heterozygosity $(H)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_het_sml

ggplot(d_qg %>% 
         mutate(gen = gen - 50000) %>%
         filter(gen > -1000) %>%
         group_by(gen, model, tau, r, nloci) %>%
         summarise(CIH = CI(meanH),
                   meanH = mean(meanH)) %>%
         filter(tau == 0.125) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = gen, y = meanH, colour = model)) +
  facet_nested(r_title + log10(r) ~ nloci_title + nloci) +
  geom_line() +
  geom_ribbon(aes(ymin = meanH - CIH, ymax = meanH + CIH, fill = model), 
              alpha = 0.2, colour = NA) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post optimum shift", 
       y = TeX("Heterozygosity $(H)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-"), guide = "none") +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_het_med

ggplot(d_qg %>% 
         mutate(gen = gen - 50000) %>%
         filter(gen > -1000) %>%
         group_by(gen, model, tau, r, nloci) %>%
         summarise(CIH = CI(meanH),
                   meanH = mean(meanH)) %>%
         filter(tau == 1.25) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = gen, y = meanH, colour = model)) +
  facet_nested(r_title + log10(r) ~ nloci_title + nloci) +
  geom_line() +
  geom_ribbon(aes(ymin = meanH - CIH, ymax = meanH + CIH, fill = model), 
              alpha = 0.2, colour = NA) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post optimum shift", 
       y = TeX("Heterozygosity $(H)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-"), guide = "none") +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_het_lrg

leg <- get_legend(plt_het_lrg)

plt_het <- plot_grid(plt_het_sml + theme(legend.position = "none"),
                     plt_het_med + theme(legend.position = "none"),
                     plt_het_lrg + theme(legend.position = "none"),
                        ncol = 1, labels = "AUTO")

plt_het <- plot_grid(plt_het,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_het
ggsave("plt_het.png", device = png, bg = "white",
       width = 560*6, height = 980*6, units = "px")
