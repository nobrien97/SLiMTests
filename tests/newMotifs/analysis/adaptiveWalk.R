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
library(ggridges)
library(GGally)

setwd("/mnt/c/GitHub/SLiMTests/tests/newMotifs/analysis")
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"#"/mnt/d/SLiMTests/tests/newMotifs/pilot/"
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
                                                    rep("numeric", times = 29)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", 
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var",
                                        "trait3_var", "trait4_var", "mahal_dist", "dist1",
                                        "dist2", "dist3", "dist4", "w", "var_w",
                                        "deltaPheno", "deltaw", "mc1_mean", 
                                        "mc2_mean", "mc3_mean", "mc4_mean", "mc5_mean",
                                        "mc6_mean", "mc7_mean", "mc8_mean", "mc9_mean",
                                        "mc10_mean", "mc11_mean"), 
                          fill = T)

# Add predictors
d_qg <- AddCombosToDF(d_qg) 

# Initial optimum distance
INIT_DIST <- sqrt(-2 * log(0.90))

# Optimum: fitness > 95%

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & w > 0.95)) %>%
  ungroup() -> d_qg

# Proportion of each model that adapted
d_prop_adapted <- d_qg %>% group_by(model, r) %>%
  filter(gen == 60000) %>%
  summarise(n = n(),
            nAdapted = sum(isAdapted),
            pAdapted = mean(isAdapted)
  )

# Output to table
stargazer(d_prop_adapted %>% mutate(r = as.character(r)) %>% as.data.frame(.), type = "html", summary = F)

# Summarise phenotype trajectories
d_qg_sum <- d_qg %>% 
  #filter(gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r) %>%
  summarise(meanDist = mean(mahal_dist),
            SEDist = se(mahal_dist),
            meanFitness = mean(w),
            SEFitness = se(w),
            meanFitnessVar = mean(var_w),
            SEFitnessVar = se(var_w))

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
  labs(x = "Generations post-optimum shift", y = "Mean Mahalanobis distance\nfrom the optimum", 
       colour = "Model") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
    override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_adapt_mdist.png", width = 12, height = 5, device = png)

# Mean fitness
ggplot(d_qg_sum,
       aes(x = gen, y = meanFitness, colour = model)) +
  facet_grid(log10(r)~.) +
  geom_line() +
  #geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanFitness - SEFitness, 
                  ymax = meanFitness + SEFitness, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean Fitness", 
       colour = "Model") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_randomisedStarts_adapt_w.png", width = 12, height = 5, device = png)

# Separate runs
d_qg$simID <- interaction(d_qg$seed, d_qg$modelindex)

sampled_examples <- d_qg %>%
  group_by(model, r, isAdapted) %>%
  slice_sample(n = 3) %>%
  ungroup() %>%
  select(simID) %>% unlist(.)

ggplot(d_qg %>% filter(simID %in% sampled_examples) %>%
         filter(gen > 40000) %>%
         mutate(gen = gen - 50000),
       aes(x = gen, y = w, colour = model, group = simID)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Did the population adapt?" + isAdapted) +
  geom_line() +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR")) +
  labs(x = "Generations post-optimum shift", y = "Mean population fitness",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines"),
        legend.position = "bottom") 
ggsave("plt_randomisedStarts_adapt_w_eg.png", width = 12, height = 5, device = png)


# Variance in fitness
ggplot(d_qg_sum,
       aes(x = gen, y = meanFitnessVar, colour = model)) +
  facet_grid(log10(r)~.) +
  geom_line() +
  #geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanFitnessVar - SEFitnessVar, 
                  ymax = meanFitnessVar + SEFitnessVar, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean variance in fitness", 
       colour = "Model") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 

ggsave("plt_randomisedStarts_adapt_wvar.png", width = 12, height = 5, device = png)

# Plot individual simulations among adapted populations
ggplot(d_qg %>% filter(isAdapted == T),
       aes(x = gen, y = w, colour = model, group = interaction(modelindex,seed))) +
  facet_grid(log10(r)~model) +
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR")) +
  labs(x = "Generations post-optimum shift", y = "Mean fitness", 
       colour = "Model") +
  theme_bw() +
  # guides(colour = guide_legend(position = "bottom",
  #                              override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12), legend.position = "none",
        panel.spacing = unit(0.75, "lines")) 

ggsave("plt_randomisedStarts_adapt_ind_w.png", width = 12, height = 5, device = png)

# What direction is selection going for each adapted simulation?
d_opt <- data.table::fread(paste0(DATA_PATH, "slim_opt.csv"), header = F, 
                          sep = ",", colClasses = c("factor", "factor", 
                                                    rep("numeric", times = 29)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", 
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var",
                                        "trait3_var", "trait4_var", "mahal_dist", "dist1",
                                        "dist2", "dist3", "dist4", "w", "var_w",
                                        "deltaPheno", "deltaw", "mc1_mean", 
                                        "mc2_mean", "mc3_mean", "mc4_mean", "mc5_mean",
                                        "mc6_mean", "mc7_mean", "mc8_mean", "mc9_mean",
                                        "mc10_mean", "mc11_mean"), 
                          fill = T)

# Each trait separately
d_qg_sum <- d_qg %>% 
  #filter(gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r) %>%
  summarise(meanDist1 = mean(dist1),
            SEDist1 = se(dist1),
            meanDist2 = mean(dist2),
            SEDist2 = se(dist2),
            meanDist3 = mean(dist3),
            SEDist3 = se(dist3),
            meanDist4 = mean(dist4),
            SEDist4 = se(dist4))


# trait sigmas for NAR, PAR, FFLC1, FFLI1, FFBH in that order
traitsigma <- matrix(c(0.111064, 0.0785957, -1, -1,
                       0.0217844, 0.00179048, -1, -1,
                       0.290934, 0.0217844, 0.0364384, -1,
                       0.101508, 0.0732563, 1.15457, -1,
                       0.335821, 0.0390923, 0.225045, 0.00138527),
                     nrow = 5, byrow = T)

rownames(traitsigma) <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")

d_qg_sum <- d_qg %>%
  group_by(model) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("dist"), names_to = "trait_num", names_prefix = "dist",
               values_to = "dist") %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r, trait_num) %>%
  mutate(dist = dist / traitsigma[model, as.integer(trait_num)]) %>%          # Adjust distance by selection strength
  summarise(meanDist = mean(dist),
            SEDist = se(dist))

# Variance in fitness
ggplot(d_qg_sum,
       aes(x = gen, y = meanDist, colour = model)) +
  facet_wrap(log10(r)~trait_num, scales = "free") +
  geom_line() +
  #geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanDist - SEDist, 
                  ymax = meanDist + SEDist, fill = model), colour = NA,
              alpha = 0.2) +
  # scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
  #                                        breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("FFBH", "FFL-C1", "FFL-I1", "NAR", "PAR"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean distance from trait optimum", 
       colour = "Model") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_adapt_dist_pertrait.png", width = 12, height = 5, device = png)

# Correlations among traits
ggpairs(d_qg %>% filter(model == "NAR"), 
        columns = 5:6, 
        columnLabels = c("Response time", "Steady state concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_trait_corr_nar.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "PAR"), 
        columns = 5:6, 
        columnLabels = c("Response time", "Steady state concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_trait_corr_par.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "FFLC1"), 
        columns = 5:7, 
        columnLabels = c("Response time", "Response delay", "Steady state concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_trait_corr_fflc1.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "FFLI1"), 
        columns = 5:7, 
        columnLabels = c("Time to half-max\nconcentration", "Maximum concentration", 
                         "Time above\nhalf-max concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_trait_corr_ffli1.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "FFBH"), 
        columns = 5:8, 
        columnLabels = c("Time to half-max\nconcentration", "Maximum concentration", 
                        "Response time to\nsecond steady state", "Second steady\nstate concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_trait_corr_ffbh.png", width = 12, height = 5, device = png)


# Correlations among distances
ggpairs(d_qg %>% filter(model == "NAR"), 
        columns = 14:15, 
        columnLabels = c("Response time", "Steady state concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_dist_corr_nar.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "PAR"), 
        columns = 14:15, 
        columnLabels = c("Response time", "Steady state concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_dist_corr_par.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "FFLC1"), 
        columns = 14:16, 
        columnLabels = c("Response time", "Response delay", "Steady state concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_dist_corr_fflc1.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "FFLI1"), 
        columns = 14:16, 
        columnLabels = c("Time to half-max\nconcentration", "Maximum concentration", 
                         "Time above\nhalf-max concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_dist_corr_ffli1.png", width = 12, height = 5, device = png)

ggpairs(d_qg %>% filter(model == "FFBH"), 
        columns = 14:17, 
        columnLabels = c("Time to half-max\nconcentration", "Maximum concentration", 
                         "Response time to\nsecond steady state", "Second steady\nstate concentration"),
        lower = list(continuous = wrap("points", shape = 1)))
ggsave("plt_dist_corr_ffbh.png", width = 12, height = 5, device = png)



# Fitness distribution (assuming normal)

# are the different traits correlated?
d_qg %>%
  filter(gen == 60000) %>%
  group_by(model, r) %>%
  summarise(cor(dist1, dist2))

cor(d_qg$dist1, d_qg_dist2)

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