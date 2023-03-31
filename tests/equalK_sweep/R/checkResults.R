# Compare the QG results with those from newTestCross to make sure they are the same
library(tidyverse)
library(gghighlight)
library(latex2exp)
library(ggh4x)
library(cowplot)
library(grid)
library(gridExtra)
library(factoextra)
library(FactoMineR)
library(gganimate)
library(akima)
library(paletteer)
library(lattice)
library(latticeExtra)


se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

# interpolation for heatmaps
# https://stackoverflow.com/a/37533314
# https://stackoverflow.com/a/65873728
dpinterp <- function(df, x, y, z) {
  x <- unlist(df[,x])
  y <- unlist(df[,y])
  z <- unlist(df[,z])
  interp_df <- akima::interp(x=x, y=y, z=z, duplicate = "median")
  interp2xyz(interp_df, data.frame=TRUE)
}


# Running on HPC for RAM reasons
path <- "/g/data/ht96/nb9894/equalK_sweep/"
local_path <- "/mnt/d/SLiMTests/tests/equalK_sweep/data/"
setwd(path)



cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")

cc_10cols <- c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#012749")

d_new <- readRDS(paste0(path, "checkpoint/d_qg.RDS"))

d_new <- readRDS(paste0(local_path, "checkpoint/d_qg.RDS"))

d_new %>% mutate(id = as_factor(paste(seed, modelindex, sep = "_"))) -> d_new

d_new$nloci_cat <- cut(d_new$nloci,
                       breaks=c(-Inf, 10, 50, 250, Inf))

d_new$sigma_cat <- cut(d_new$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))


d_new %>%
  group_by(seed, nloci, sigma) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  mutate(phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  ungroup() -> d_new_adapted


d_new %>%
  group_by(seed, nloci, sigma) %>%
  filter(all(phenomean < 1.9 | phenomean > 2.1)) %>%
  ungroup() -> d_new_maladapted


d_indPheno <- readRDS("./checkpoint/d_indPheno.RDS")

seed <- sample(1:.Machine$integer.max, 1)
set.seed(seed)
# 563655642
set.seed(563655642)

sampled_seeds <- d_new_adapted %>% filter(gen > 49500, phenomean < 5) %>%
  group_by(nloci, sigma) %>% 
  select(nloci, sigma, seed, modelindex, id) %>%
  sample_n(1)

d_qg_sum <- d_new %>% group_by(gen, nloci_cat, sigma_cat) %>%
  summarise(He = mean(meanH),
            ciHe = qnorm(0.975) * se(meanH),
            meanPheno = mean(phenomean),
            ciPheno = qnorm(0.975) * se(phenomean),
            meanDist = mean(dist),
            ciDist = qnorm(0.975) * se(dist))


pheno_time <- ggplot(d_qg_sum %>% 
                       filter(gen > 49000) %>% mutate(gen = gen - 50000),
       aes(x = gen, y = meanPheno, colour = sigma_cat)) +
  facet_grid(nloci_cat~.) + 
  geom_line() +
  scale_colour_manual(values = cc_ibm) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = meanPheno - ciPheno, ymax = meanPheno + ciPheno, fill = sigma_cat),
                  colour = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", colour = "Mutational\neffect\nvariance") +
  theme_bw() +
  theme(text = element_text(size=16))

pheno_time
ggsave("pheno_time.png", pheno_time, width = 8, height = 8)

if (!dir.exists("./pheno_anim"))
{
  dir.create("./pheno_anim")
}

pheno_time <- ggplot(d_new %>% filter(gen > 49000) %>% mutate(gen = gen - 50000),
       aes(x = nloci, y = sigma, colour = phenomean)) +
  geom_jitter(size = 10) +
  scale_colour_distiller(palette = "Spectral") +
  labs(x = "Number of loci", y = "Mutational effect variance", 
    colour = "Phenotypic\nmean") +
  theme_bw() +
  theme(text = element_text(size=36), legend.key.size = unit(2, 'cm'))


pheno_time <- pheno_time + transition_states(gen) +
  labs(title = "Generations post-optimum shift: {closest_state}")

animate(pheno_time, nframes = 42, duration = 10, width = 1920, height = 1920, 
  renderer = file_renderer(dir = "./pheno_anim"))

# Single population
sample_num <- 1
pheno_time_single <- ggplot(d_new %>% filter(seed %in% sampled_seeds$seed[sample_num],
                                             modelindex %in% sampled_seeds$modelindex[sample_num]) %>% 
                       mutate(gen = gen - 50000),
                     aes(x = gen, y = phenomean)) +
  geom_line() +
  labs(x = "Generations post-optimum shift", y = "Trait value") +
  theme_bw() +
  theme(text = element_text(size=20))
pheno_time_single



# time to optimum
d_new %>% filter(gen >= 49500) %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  group_by(seed, modelindex) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1),
         initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1)) %>%
  distinct(seed, modelindex, .keep_all = T) %>%
  ungroup() %>%
  filter(adaptTime != -1, initRespTime != -1) %>%
  group_by(nloci, sigma) %>%
  summarise(adaptTime = mean(adaptTime), 
            CIAdaptTime = CI(adaptTime),
            initRespTime = mean(initRespTime),
            CIInitRespTime = CI(initRespTime)) -> d_new_adapttime


grid_adapt_time <- levelplot(adaptTime ~ nloci*sigma, d_new_adapttime %>% 
                               mutate(adaptTime = adaptTime - 50000) %>% 
                 select(nloci, sigma, adaptTime), 
                 col.regions = paletteer_c("viridis::viridis", 100, -1),
                 panel = panel.levelplot.points, 
                 xlab = list(label = "Number of loci", cex = 1.2), 
                 ylab = list(label = "Mutational effect variance", cex = 1.2), 
                 colorkey = list(space = "bottom",
                                 title = "Time to adaptation (generations)",
                                 labels = list(cex = 1.2)),
                 par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                 scales = list(tck = c(1, 0),
                               cex = 1.2,
                               x = list(at = c(200, 400, 600, 800, 1000)),
                               y = list(at = c(0.0, 0.5, 1.0, 1.5, 2.0))), 
                 pretty = T) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_adapt_time <- 
  contourplot(
    adaptTime ~ nloci*sigma, d_new_adapttime %>% 
      mutate(adaptTime = adaptTime - 50000) %>% 
      select(nloci, sigma, adaptTime), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_adaptTime <- grid_adapt_time + contours_adapt_time
plt_adaptTime

grid_resp_time <- levelplot(initRespTime ~ nloci*sigma, d_new_adapttime %>% 
                               mutate(initRespTime = initRespTime - 50000) %>% 
                               select(nloci, sigma, initRespTime), 
                             col.regions = paletteer_c("viridis::viridis", 100, -1),
                             panel = panel.levelplot.points, 
                             xlab = list(label = "Number of loci", cex = 1.2), 
                            ylab = list(label = "Mutational effect variance", cex = 1.2), 
                             colorkey = list(space = "bottom",
                                             title = "Time to initial response (generations)",
                                             labels = list(cex = 1.2)),
                             par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                             scales = list(tck = c(1, 0),
                                           cex = 1.2,
                                           x = list(at = c(200, 400, 600, 800, 1000)),
                                           y = list(at = c(0.0, 0.5, 1.0, 1.5, 2.0))), 
                            pretty = T) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_resp_time <- 
  contourplot(
    initRespTime ~ nloci*sigma, d_new_adapttime %>% 
      mutate(initRespTime = initRespTime - 50000) %>% 
      select(nloci, sigma, initRespTime), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_respTime <- grid_resp_time + contours_resp_time
plt_respTime

plot_grid(plt_respTime, plt_adaptTime, ncol = 2, labels = "AUTO")
ggsave("phenoRespTimeLattice.png", width = 15, height = 8, bg = "white")

# single population
d_new %>% filter(gen >= 49500, seed %in% sampled_seeds$seed[sample_num],
                 modelindex %in% sampled_seeds$modelindex[sample_num]) %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]) - 50000, -1),
         initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]) - 50000, -1)) -> d_adapttime_single


# size effects
d_com <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/d_combined_after.RDS")

# single population

d_com %>%
  filter(gen >= 49500, seed %in% sampled_seeds$seed[sample_num],
         modelindex %in% sampled_seeds$modelindex[sample_num]) %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  mutate(gen = gen - 50000, 
         molTrait = recode_factor(mutType,
                                  `1`="Neutral", `2`="Del", `3`="$$\\alpha_Z$$", `4`="$$\\beta_Z$$", 
                                  `5`="$$K_Z$$", `6`="$$K_{XZ}$$"),
         adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1),
         initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1))-> d_com_single


mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$"),
  TeX("$K_Z$"),
  TeX("$K_{XZ}$")
)

plt_meanVal_single <- levelplot(phenomean ~ Freq*value, d_com_single %>% 
                               select(value, Freq, phenomean), 
                             col.regions = paletteer_c("viridis::viridis", 100, -1),
                             panel = panel.levelplot.points, 
                             xlab = list(label = "Frequency", cex = 1.2), 
                             ylab = list(label = "Effect size", cex = 1.2), 
                             colorkey = list(space = "bottom",
                                             title = "Phenotype mean",
                                             labels = list(cex = 1.2)),
                             par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                             scales = list(tck = c(1, 0),
                                           cex = 1.2),
                                           #x = list(at = c(200, 400, 600, 800, 1000)),
                                           #y = list(at = c(-0.5, -0.25, 0.0, 0.25, 0.5))), 
                             pretty = T) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_adapt_time <- 
  contourplot(
    phenomean ~ Freq*value, d_com_single %>% 
      select(value, Freq, phenomean), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_adaptTime <- plt_meanVal_single + contours_adapt_time
plt_adaptTime


plt_meanVal_single <- ggplot(d_com_single, 
                             aes(x = gen, y = value, color = molTrait)) +
  geom_boxplot() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  labs(x = "Generations post-optimum shift", y = "Mutational effect", color = "Molecular component") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal_single

ggsave("meanVal.png", plt_meanVal, width = 14, height = 10, bg = "white")


plt_meanFreq <- ggplot(d_meanfreq, 
                       aes(x = gen, y = meanFreq, color = molTrait)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanFreq - ci95Freq), ymax = (meanFreq + ci95Freq), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  
  labs(x = "Generations post-optimum shift", y = "Mean frequency", color = "Molecular component") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal

ggsave("meanFreq.png", plt_meanFreq, width = 14, height = 10, bg = "white")


d_com %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  group_by(seed, nloci, sigma, mutType) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1)) %>%
  # get mean freq/value of mutations segregating in models at the adaptation time
  filter(gen == adaptTime) %>%
  mutate(meanFreq = mean(Freq), meanValue = mean(abs(value))) %>%
  ungroup() %>%
  filter(adaptTime != -1) %>%
  group_by(nloci, sigma, mutType) %>%
  summarise(CIAdaptTime = CI(adaptTime),
            adaptTime = mean(adaptTime),
            meanFreq = mean(Freq), 
            CIMeanFreq = CI(Freq),
            meanValue = mean(value),
            CIMeanValue = CI(value)) %>% ungroup() -> d_com_adapttime

grid_adapt_time <- levelplot(adaptTime ~ nloci*meanValue, d_com_adapttime %>% 
                               mutate(adaptTime = adaptTime - 50000) %>% 
                               select(nloci, meanValue, meanFreq, adaptTime), 
                             col.regions = paletteer_c("viridis::viridis", 100, -1),
                             panel = panel.levelplot.points, 
                             xlab = list(label = "Number of loci", cex = 1.2), 
                             ylab = list(label = "Mean effect size", cex = 1.2), 
                             colorkey = list(space = "bottom",
                                             title = "Time to adaptation (generations)",
                                             labels = list(cex = 1.2)),
                             par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                             scales = list(tck = c(1, 0),
                                           cex = 1.2,
                                           x = list(at = c(200, 400, 600, 800, 1000)),
                                           y = list(at = c(-0.5, -0.25, 0.0, 0.25, 0.5))), 
                             pretty = T) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_adapt_time <- 
  contourplot(
    adaptTime ~ nloci*meanValue, d_com_adapttime %>% 
      mutate(adaptTime = adaptTime - 50000) %>% 
      select(nloci, meanValue, meanFreq, adaptTime), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_adaptTime <- grid_adapt_time + contours_adapt_time
plt_adaptTime


d_com %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  group_by(seed, nloci, sigma, mutType) %>%
  mutate(initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1)) %>%
  # get mean freq/value of mutations segregating in models at the adaptation time
  filter(gen == initRespTime) %>%
  mutate(meanFreq = mean(Freq), meanValue = mean(abs(value))) %>%
  ungroup() %>%
  filter(initRespTime != -1) %>%
  group_by(nloci, sigma, mutType) %>%
  summarise(CIInitRespTime = CI(initRespTime),
            initRespTime = mean(initRespTime),
            meanFreq = mean(Freq),
            CIMeanFreq = CI(Freq),
            meanValue = mean(value),
            CIMeanalue = CI(value)) %>% ungroup() -> d_com_adapttime

grid_resp_time <- levelplot(initRespTime ~ nloci*meanValue, d_com_adapttime %>% 
                               mutate(initRespTime = initRespTime - 50000) %>% 
                               select(nloci, meanValue, meanFreq, initRespTime), 
                             col.regions = paletteer_c("viridis::viridis", 100, -1),
                             panel = panel.levelplot.points, 
                             xlab = list(label = "Number of loci", cex = 1.2), 
                             ylab = list(label = "Mean effect size", cex = 1.2), 
                             colorkey = list(space = "bottom",
                                             title = "Time to initial response (generations)",
                                             labels = list(cex = 1.2)),
                             par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                             scales = list(tck = c(1, 0),
                                           cex = 1.2,
                                           x = list(at = c(200, 400, 600, 800, 1000))), 
                             pretty = T, ylim = c(-0.86, 0.36)) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_resp_time <- 
  contourplot(
    initRespTime ~ nloci*meanValue, d_com_adapttime %>% 
      mutate(initRespTime = initRespTime - 50000) %>% 
      select(nloci, meanValue, meanFreq, initRespTime), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_respTime <- grid_resp_time + contours_resp_time
plt_respTime

plot_grid(plt_respTime, plt_adaptTime, ncol = 2, labels = "AUTO")
ggsave("respTimeFX.png", width = 15, height = 8, bg = "white")

# frequency by effect size

d_com %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  group_by(seed, nloci, sigma, mutType) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]), -1)) %>%
  # get mean freq/value of mutations segregating in models at the adaptation time
  filter(gen == adaptTime) %>%
  mutate(meanFreq = mean(Freq), meanValue = mean(abs(value))) %>%
  ungroup() %>%
  filter(adaptTime != -1) %>%
  group_by(nloci, sigma, mutType) %>%
  summarise(CIAdaptTime = CI(adaptTime),
            adaptTime = mean(adaptTime),
            meanFreq = mean(Freq), 
            CIMeanFreq = CI(Freq),
            meanValue = mean(value),
            CIMeanValue = CI(value)) %>% ungroup() -> d_com_adapttime

grid_adapt_time <- levelplot(adaptTime ~ meanFreq*meanValue, d_com_adapttime %>% 
                               mutate(adaptTime = adaptTime - 50000) %>% 
                               select(meanValue, meanFreq, adaptTime), 
                             col.regions = paletteer_c("viridis::viridis", 100, -1),
                             panel = panel.levelplot.points, 
                             xlab = list(label = "Mean frequency", cex = 1.2), 
                             ylab = list(label = "Mean effect size", cex = 1.2), 
                             colorkey = list(space = "bottom",
                                             title = "Time to adaptation (generations)",
                                             labels = list(cex = 1.2)),
                             par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                             scales = list(tck = c(1, 0),
                                           cex = 1.2,
                                           y = list(at = c(-0.5, -0.25, 0.0, 0.25, 0.5))), 
                             pretty = T, xlim = c(0.05, 0.47)) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_adapt_time <- 
  contourplot(
    adaptTime ~ meanFreq*meanValue, d_com_adapttime %>% 
      mutate(adaptTime = adaptTime - 50000) %>% 
      select(meanValue, meanFreq, adaptTime), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_adaptTime <- grid_adapt_time + contours_adapt_time
plt_adaptTime

write_csv(d_com_adapttime %>% mutate(adaptTime = adaptTime - 50000), "d_adaptedTime.csv")


d_com %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  group_by(seed, nloci, sigma, mutType) %>%
  mutate(initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]), -1)) %>%
  # get mean freq/value of mutations segregating in models at the adaptation time
  filter(gen == initRespTime) %>%
  mutate(meanFreq = mean(Freq), meanValue = mean(abs(value))) %>%
  ungroup() %>%
  filter(initRespTime != -1) %>%
  group_by(nloci, sigma, mutType) %>%
  summarise(CIInitRespTime = CI(initRespTime),
            initRespTime = mean(initRespTime),
            meanFreq = mean(Freq),
            CIMeanFreq = CI(Freq),
            meanValue = mean(value),
            CIMeanValue = CI(value)) %>% ungroup() -> d_com_adapttime

grid_resp_time <- levelplot(initRespTime ~ meanFreq*meanValue, d_com_adapttime %>% 
                              mutate(initRespTime = initRespTime - 50000) %>% 
                              select(meanValue, meanFreq, initRespTime), 
                            col.regions = paletteer_c("viridis::viridis", 100, -1),
                            panel = panel.levelplot.points, 
                            xlab = list(label = "Mean frequency", cex = 1.2), 
                            ylab = list(label = "Mean effect size", cex = 1.2), 
                            colorkey = list(space = "bottom",
                                            title = "Time to initial response (generations)",
                                            labels = list(cex = 1.2)),
                            par.settings = list(layout.heights = list(xlab.key.padding = 4)),
                            scales = list(tck = c(1, 0),
                                          cex = 1.2,
                                          x = list(at = c(0.1, 0.2, 0.3, 0.4))), 
                            pretty = T, xlim = c(0.05, 0.47)) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_resp_time <- 
  contourplot(
    initRespTime ~ meanFreq*meanValue, d_com_adapttime %>% 
      mutate(initRespTime = initRespTime - 50000) %>% 
      select(meanValue, meanFreq, initRespTime), 
    panel=panel.2dsmoother, col = "white",
    labels = list(col = "white",
                  cex = 1.2))

plt_respTime <- grid_resp_time + contours_resp_time
plt_respTime

write_csv(d_com_adapttime %>% mutate(initRespTime = initRespTime - 50000), "d_initRespTime.csv")


plot_grid(plt_respTime, plt_adaptTime, ncol = 2, labels = "AUTO")
ggsave("respTimeFreqFX.png", width = 15, height = 8, bg = "white")




# Mean allele frequency
d_com <- readRDS("d_combined_after.RDS")

d_com$nloci_cat <- cut(d_com$nloci,
                       breaks=c(10, 50, 250, Inf))

d_com$sigma_cat <- cut(d_com$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))



d_com %>%
  mutate(molTrait = recode_factor(mutType, 
                                  `1`="Neutral", `2`="Del", `3`="$$\\alpha_Z$$", `4`="$$\\beta_Z$$", 
                                  `5`="$$K_Z$$", `6`="$$K_{XZ}$$")) -> d_com

d_com %>%
  group_by(seed, nloci_cat, sigma_cat) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  mutate(phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  ungroup() -> d_com_adapted


d_com_adapted %>%
  mutate(gen = gen - 50000,
         value = value) %>%
  group_by(gen, molTrait, nloci_cat, sigma_cat) %>%
  summarise(meanValue = mean(value),
            ci95Value = CI(value),
            meanFreq = mean(Freq),
            ci95Freq = CI(Freq)) -> d_meanfreq

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$"),
  TeX("$K_Z$"),
  TeX("$K_{XZ}$")
)



plt_meanVal <- ggplot(d_meanfreq, 
                      aes(x = gen, y = meanValue, color = molTrait)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanValue - ci95Value), ymax = (meanValue + ci95Value), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                           breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                           breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean mutational effect", color = "Molecular component") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal

ggsave("meanVal.png", plt_meanVal, width = 14, height = 10, bg = "white")


plt_meanFreq <- ggplot(d_meanfreq, 
                      aes(x = gen, y = meanFreq, color = molTrait)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanFreq - ci95Freq), ymax = (meanFreq + ci95Freq), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                           breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                           breaks = NULL, labels = NULL)) +

  labs(x = "Generations post-optimum shift", y = "Mean frequency", color = "Molecular component") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal

ggsave("meanFreq.png", plt_meanFreq, width = 14, height = 10, bg = "white")


# SFS
d_com_adapted %>%
  mutate(Freq_bin = cut(Freq,
                        breaks = seq(from = 0, to = 1, by = 0.05))) %>% # categorise freq
  group_by(gen, seed, molTrait, nloci_cat, sigma_cat) %>%
  mutate(binCount = n()) %>%
  group_by(gen, seed, molTrait, nloci_cat, sigma_cat, Freq_bin) %>%
  mutate(freqBinCount = n()/binCount) %>% # Get the number in each bin - normalise!
  ungroup() -> d_freqs

d_freqs %>%
  group_by(gen, molTrait, nloci_cat, sigma_cat, Freq_bin) %>%
  summarise(meanFreqBinCount = mean(freqBinCount),
            seFreqBinCount = se(freqBinCount),
            ci95FreqBinCount = qnorm(0.975)*seFreqBinCount) %>%
  ungroup() %>%
  complete(gen, molTrait, nloci_cat, sigma_cat, Freq_bin, 
           fill = list(meanFreqBinCount = 0)) %>%
  mutate(Freq = as.numeric(Freq_bin) * 0.05, # convert freq back to continuous
         gen = as.integer(gen - 50000)) -> d_freqs 


plt_freqs <- ggplot(d_freqs %>% filter(gen == 1950) , 
                    aes(x = Freq, y = meanFreqBinCount, color = molTrait)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  geom_errorbar(aes(ymin = (meanFreqBinCount - ci95FreqBinCount), ymax = (meanFreqBinCount + ci95FreqBinCount), 
                    color = molTrait)) +
  scale_y_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                           breaks = NULL, labels = NULL)) +
  scale_x_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                           breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency", y = "Mean proportion of mutations", color = "Molecular component") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))

plt_freqs
ggsave("sfs_end.png", plt_freqs, width = 14, height = 10, bg = "white")



plt_freqs <- ggplot(d_freqs %>% filter(gen < 2000), 
       aes(x = Freq, y = meanFreqBinCount, color = molTrait)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line() +
  geom_errorbar(aes(ymin = (meanFreqBinCount - ci95FreqBinCount), ymax = (meanFreqBinCount + ci95FreqBinCount), 
                    color = molTrait)) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_y_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency", y = "Mean proportion of mutations", color = "Molecular component") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom")

plt_freqs

# gganimate stuff
plt_freqs <- plt_freqs + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_freqs, nframes = 42, duration = 10, width = 1920, height = 1920, 
  renderer = file_renderer(dir = "./sfs_anim", overwrite = T))



# PCA mol traits

d_indPheno <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/d_indPheno.RDS")
d_indPheno %>%
  filter(aZ < 10, bZ < 10, KZ < 10, KXZ < 10) %>%
  mutate(isAdapted = between(phenotype, 1.9, 2.1)) %>%
  mutate(seed = as_factor(seed)) -> d_isAdapted

d_isAdapted$nloci_cat <- cut(d_isAdapted$nloci,
                       breaks=c(-Inf, 10, 50, 250, Inf))

d_isAdapted$sigma_cat <- cut(d_isAdapted$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))



res.pca <- PCA(d_isAdapted %>% select(aZ, bZ, KZ, KXZ), scale.unit = T, graph = F)
fviz <- fviz_eig(res.pca, addlabels = TRUE)
ggsave("fviz_moltraits.png", fviz, bg = "white")
var <- get_pca_var(res.pca)
head(var$contrib)

library(mmtable2)

var$contrib %>% mmtable(cells = value)

# https://tem11010.github.io/Plotting-PCAs/
d_isAdapted$pc1 <- res.pca$ind$coord[, 1]
d_isAdapted$pc2 <- res.pca$ind$coord[, 2]

pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars 

ggplot(d_isAdapted %>%
         mutate(gen = gen - 50000), aes(x = pc1, y = pc2, colour = isAdapted)) +
  #facet_grid2(nloci~sigma, scales = "free", independent = T) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(shape = 1, size = 2) +
  scale_colour_manual(values = c(cc_ibm[3], cc_ibm[1])) +
  labs(x = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), y = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), colour = "Population adapted") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_isAdapted_pca

ggsave("moltrait_pca_adapted.png", plt_isAdapted_pca, width = 10, height = 10)
ggsave("moltrait_pca_adapted_facet.png", plt_isAdapted_pca + 
         facet_grid(nloci_cat~sigma_cat) +
         scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                breaks = NULL, labels = NULL)) +
         scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                breaks = NULL, labels = NULL)), 
       width = 10, height = 10)

plt_isAdapted_pca <- plt_isAdapted_pca + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_isAdapted_pca, nframes = 41, duration = 10, width = 1920, height = 1920, renderer = ffmpeg_renderer())
anim_save("moltrait_pca.mp4", last_animation())

plt_isAdapted_pca <- plt_isAdapted_pca + 
  facet_grid(nloci~sigma) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_freqs, nframes = 42, duration = 10, width = 1920, height = 1920, 
  renderer = file_renderer(dir = "./moltraitpca_anim", overwrite = T))


animate(plt_isAdapted_pca, nframes = 5, duration = 5, width = 720, height = 720, renderer = ffmpeg_renderer())
anim_save("moltrait_facet_pca.mp4", last_animation())





d_new_adapted %>% filter(gen >= 49500, phenomean < 5) %>%
  mutate(gen = gen - 50000) %>% 
ggplot(aes(x = gen, y = phenomean, group = seed, colour = seed)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  gghighlight(id %in% sampled_seeds$id, calculate_per_facet = T, use_direct_label = F) +
  scale_colour_manual(values = rep(cc_ibm[3], 6), guide = "none") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adapted
plt_adapted

ggsave("pheno_adapted.png", plt_adapted, width = 8, height = 8)


d_indPheno <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/equalK/slim_sampled_moltrait.csv", col_names = F)
d_indPheno %>% pivot_longer(cols = 4:ncol(d_indPheno), names_to = NULL, values_to = "phenotype") -> d_indPheno

d_indPheno %>% mutate(type = rep(c("phenotype", "aZ", "bZ", "KZ", "KXZ"), 
                                 length.out = nrow(d_indPheno))) -> d_indPheno

d_indPheno %>% 
  group_by(type) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = type, values_from = phenotype) -> d_indPheno

colnames(d_indPheno) <- c("gen", "seed", "modelindex", "ind", "phenotype", "aZ", "bZ", "KZ", "KXZ")

d_indPheno %>%
  group_by(gen, seed, modelindex) %>%
  mutate(ind = 1:n()) %>%
  ungroup() %>%
  mutate(ind = as_factor(ind)) -> d_indPheno


d_indPheno %>% mutate(ind = as_factor(paste(seed, modelindex, ind, sep = "_")),
                      id = as_factor(paste(seed, modelindex, sep = "_")),
                      nloci = d_combos$nloci[.$modelindex],
                      sigma = d_combos$locisigma[.$modelindex]) -> d_indPheno

saveRDS(d_indPheno, "d_indPheno.RDS")

d_indPheno <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/equalK/d_indPheno.RDS")

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$"),
  TeX("$K_Z$"),
  TeX("$K_{XZ}$")
)

ggplot(d_indPheno %>% filter(gen >= 49500, id %in% sampled_seeds$id) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = phenotype, group = as.factor(gen))) +
  facet_grid(nloci~sigma) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 10)) +
  #gghighlight(ind_mod %in% chosen_inds, calculate_per_facet = T) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Phenotypic value") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_indpheno
plt_indpheno

ggsave("pheno_inds_adapted_total.png", plt_indpheno, width = 8, height = 8)

# Molecular traits
ggplot(d_indPheno %>% filter(gen >= 49500, id %in% sampled_seeds$id) %>%
         pivot_longer(cols = c(aZ, bZ, KZ, KXZ), names_to = "molTrait", values_to = "molTraitVal") %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = molTraitVal, colour = molTrait, group = interaction(as.factor(gen), molTrait))) +
  facet_grid(nloci~sigma, scales = "free") +
  geom_boxplot(outlier.size = 0.5) +
  #coord_cartesian(ylim = c(0, 11)) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Molecular trait value", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indMol
plt_indMol

ggsave("molTrait_inds_adapted.png", plt_indMol, width = 8, height = 8)
ggsave("molTrait_inds_adapted_zoom.png", plt_indMol + coord_cartesian(ylim = c(0, 11)), width = 8, height = 8)


# Heterozygosity
ggplot(d_new_adapted %>% filter(gen >= 49500, id %in% sampled_seeds$id) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = meanH)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean population heterozygosity") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_meanH
plt_meanH

ggsave("meanH_adapted.png", plt_meanH, width = 8, height = 8)

# Mutations for all individuals in the sampled models (red highlight in pheno figure)
ggplot(d_ind_com %>% filter(gen >= 49500, mutType != 1) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = log(value), colour = as.factor(mutType), 
           group = interaction(as.factor(gen), as.factor(mutType)))) +
  facet_grid(nloci~sigma) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  #geom_point() +
  geom_boxplot(outlier.size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Molecular effect size", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts
plt_indmuts

ggsave("pheno_indMuts_adapted_allInds.png", plt_indmuts, width = 8, height = 8)

# Mutations in one individual only - the maximum outlier
d_ind_com %>% 
  filter(gen >= 49500, mutType != 1) %>%
  group_by(gen, modelindex, mutType) %>%
  mutate(valMed = median(log(value)), 
         devVal = (log(value) - valMed)) %>%
  filter(devVal == max(devVal)) -> d_ind_com_filter


ggplot(d_ind_com %>% filter(gen >= 49500, mutType != 1, ind_mod %in% chosen_inds, sigma == 0.1) %>% 
         mutate(gen = gen - 50000, sigma = fct_drop(sigma)), 
       aes(x = gen, y = value, colour = as.factor(mutType))) +
  facet_grid(nloci~sigma, drop = F) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  geom_point(size = 0.5) +
  #geom_boxplot(outlier.size = 0.1) +
  #scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
  #                                       breaks = NULL, labels = NULL)) +
  #scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
  #                                       breaks = NULL, labels = NULL)) +
  labs(x = NULL, y = "Molecular effect size", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts_s01
plt_indmuts_s01

ggplot(d_ind_com %>% filter(gen >= 49500, mutType != 1, ind_mod %in% chosen_inds, sigma == 1) %>% 
         mutate(gen = gen - 50000, sigma = fct_drop(sigma)), 
       aes(x = gen, y = value, colour = as.factor(mutType))) +
  facet_grid(nloci~sigma, drop = F) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  geom_point(size = 0.5) +
  #geom_boxplot(outlier.size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  #scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
  #                                       breaks = NULL, labels = NULL)) +
  labs(x = NULL, y = NULL, colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts_s1
plt_indmuts_s1

indMuts_leg <- get_legend(plt_indmuts_s01)
y.grob <- textGrob("Number of QTLs", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Mutational effect size variance", 
                     gp = gpar(fontsize = 16))
x_lab <- textGrob("Generations post-optimum shift", 
                  gp = gpar(fontsize = 16))

plt_indmuts <- plot_grid(plt_indmuts_s01 + theme(strip.text.y = element_blank(), legend.position = "none"), 
                         plt_indmuts_s1 + theme(axis.title.y.left = element_blank(), legend.position = "none"), 
                         ncol = 2)
plt_indmuts <- plot_grid(top.grob, plt_indmuts, x_lab, indMuts_leg, ncol = 1, rel_heights = c(0.03, 1, 0.03, 0.1))
plt_indmuts

ggsave("pheno_indMuts_adapted_selectedInds.png", plt_indmuts, width = 8, height = 8, bg = "white")


# Number adapted
d_new %>%
  group_by(seed, nloci, sigma) %>%
  mutate(isAdapted = any(gen >= 51800 & between(phenomean, 1.9, 2.1)),
         phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  distinct(.keep_all = T) %>%
  select(seed, nloci, sigma, isAdapted) %>%
  ungroup(seed) %>%
  summarise(pAdapted = mean(isAdapted)) %>%
  ungroup() -> d_isAdapted

ggplot(d_isAdapted, aes(x = nloci, y = pAdapted, fill = sigma)) +
  geom_col(position = position_dodge(0.9)) +
  labs(x = "Number of loci", y = "Probability of adaptation", fill = "Mutational effect size variance") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_pAdapted

ggsave("probAdapted.png", plt_pAdapted, width = 5, height = 5)

# Probability of adaptation by molecular trait
d_indPheno %>%
  filter(aZ < 10, bZ < 10, KZ < 10, KXZ < 10) %>%
  mutate(isAdapted = between(phenotype, 1.9, 2.1)) %>%
  mutate(seed = as_factor(seed)) -> d_isAdapted



res.pca <- PCA(d_isAdapted %>% select(aZ, bZ, KZ, KXZ), scale.unit = T, graph = F)
fviz_eig(res.pca, addlabels = TRUE)
ggsave("fviz.png", bg = "white")
var <- get_pca_var(res.pca)
head(var$contrib)

library(mmtable2)

var$contrib %>% mmtable(cells = value)

# https://tem11010.github.io/Plotting-PCAs/
d_isAdapted$pc1 <- res.pca$ind$coord[, 1]
d_isAdapted$pc2 <- res.pca$ind$coord[, 2]

pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars 

ggplot(d_isAdapted %>%
         mutate(gen = gen - 50000), aes(x = pc1, y = pc2, colour = isAdapted)) +
  #facet_grid2(nloci~sigma, scales = "free", independent = T) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(shape = 1, size = 2) +
  scale_colour_manual(values = c(cc_ibm[3], cc_ibm[1])) +
  labs(x = "PC 1 (45.3%)", y = "PC 2 (27%)", colour = "Population adapted") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_isAdapted_pca

ggsave("moltrait_pca_adapted.png", plt_isAdapted_pca, width = 10, height = 10)
ggsave("moltrait_pca_adapted_facet_equalK.png", plt_isAdapted_pca + 
         facet_grid(nloci~sigma) +
         scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                breaks = NULL, labels = NULL)) +
         scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                breaks = NULL, labels = NULL)), 
       width = 10, height = 10)

plt_isAdapted_pca <- plt_isAdapted_pca + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_isAdapted_pca, nframes = 41, duration = 10, width = 720, height = 720, renderer = ffmpeg_renderer())
anim_save("moltrait_pca.mp4", last_animation())

plt_isAdapted_pca <- plt_isAdapted_pca + 
  facet_grid(nloci~sigma) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_isAdapted_pca, nframes = 5, duration = 5, width = 720, height = 720, renderer = ffmpeg_renderer())
anim_save("moltrait_facet_pca.mp4", last_animation())


# Eigenvector of pheno/moltrait/freq/value
d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/equalK/d_com_prefiltered.RDS")
d_com %>% ungroup() %>% mutate(isAdapted = between(phenomean, 1.9, 2.1),
                 alleleAge = gen - originGen,
                 ) -> d_isAdapted


res.pca <- PCA(d_isAdapted %>% 
                 select(S, beta, Freq, value), 
               scale.unit = T, graph = F)
fviz_eig(res.pca, addlabels = TRUE)
ggsave("scree_phenofreq.png", bg = "white")
var <- get_pca_var(res.pca)
var$contrib


# https://tem11010.github.io/Plotting-PCAs/
d_isAdapted$pc1 <- res.pca$ind$coord[, 1]
d_isAdapted$pc2 <- res.pca$ind$coord[, 2]
pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars 

ggplot(d_isAdapted %>%
         mutate(gen = gen - 50000), aes(x = pc1, y = pc2, colour = isAdapted)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(shape = 1, size = 2) +
  scale_colour_manual(values = c(cc_ibm[3], cc_ibm[1])) +
  labs(x = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), y = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), colour = "Is adapted?") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_isAdapted_pca

ggsave("phenofreq_pca_adapted.png", plt_isAdapted_pca, width = 10, height = 10)
ggsave("phenofreq_pca_adapted_facet_equalK.png", plt_isAdapted_pca + 
         facet_grid(nloci~sigma) +
         scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                breaks = NULL, labels = NULL)) +
         scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                breaks = NULL, labels = NULL)), 
       width = 10, height = 10)

anim <- plt_isAdapted_pca + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(anim, nframes = 41, duration = 10, width = 720, height = 720, renderer = ffmpeg_renderer())
anim_save("freqvalue_pca.mp4", last_animation())

anim <- plt_isAdapted_pca + 
  facet_grid(nloci~sigma) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(anim, nframes = 41, duration = 10, width = 720, height = 720, renderer = ffmpeg_renderer())
anim_save("freqvalue_facet_pca.mp4", last_animation())


# just freq vs effect surface
ggplot(d_isAdapted %>%
         mutate(gen = gen - 50000), aes(x = Freq, y = scale(value), colour = isAdapted)) +
  #geom_bin_2d() +
  geom_point(shape = 1, size = 2) +
  scale_colour_manual(values = c(cc_ibm[3], cc_ibm[1])) +
  labs(x = "Allele frequency", y = "Allelic effect size", colour = "Is adapted?") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_isAdapted_valfreq

ggsave("phenofreq_adapted.png", plt_isAdapted_valfreq, width = 10, height = 10)
ggsave("phenofreq_adapted_facet_equalK.png", plt_isAdapted_valfreq + 
         facet_grid(nloci~sigma) +
         scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                breaks = NULL, labels = NULL)) +
         scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                breaks = NULL, labels = NULL)), 
       width = 10, height = 10)
