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
library(ggpmisc)
library(deSolve)
library(DescTools)

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
  geom_vline(xintercept = 2, linetype = "dashed") +
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
  geom_vline(xintercept = 2, linetype = "dashed") +
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

# supp fig 1 - adaptive step timing for populations not yet at the optimum
ggplot(d_fix_ranked_combined %>% filter(rank > 0, phenomean < 1.9),
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

# text: timing of fixations
d_fix_ranked_combined %>% filter(rank > 0, phenomean < 1.9) %>%
  mutate(gen = gen - 50000) %>%
  group_by(model, rank) %>%
  summarise(meanGen = mean(gen), CIGen = CI(gen)) -> d_meanGenTiming

ggplot(d_meanGenTiming, 
       aes(x = as.factor(rank), y = meanGen, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = meanGen - CIGen, 
                              ymax = meanGen + CIGen),
                width = 0.2) +
  geom_line(aes(group = model)) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Adaptive step", y = "Mean adaptive step fixation time",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_meanGen
plt_meanGen

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

# Supp fig: dist of beneficial fixations at each step - zoom in of Fig. 2B
ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = s, y = as.factor(rank), fill = model)) +
  geom_density_ridges(alpha = 0.4) + 
  coord_cartesian(xlim = c(0, 0.1)) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
ggsave("s_fig_driftbarrier.png", device = png)

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

# Combined across all time points
ggplot(mutExp_combined, aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs( x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_effectsizerandom
plt_effectsizerandom

# Find modes
d <- density(mutExp$s)

modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

mutExp_modes <- modes(d)
mutExp_modes[order(mutExp_modes$y, decreasing = T),]

d_add <- density(mutExp_add$s)
mutExp_add_modes <- modes(d_add)
mutExp_add_modes[order(mutExp_add_modes$y, decreasing = T),]


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
ggplot(mutExp_sum_combined %>% filter(rank < 4) %>% mutate(waitingTime = 1/(10000 * (9.1528*10^-6) * percBeneficial),
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
        fill = guide_colourbar(barwidth = 10, title.vjust = 0.87)) # i love magic numbers
  )
}
plotaZbZLandscape(0, 3) -> plt_aZbZ_landscape
plt_aZbZ_landscape

# B - GPW map of aZbZ
ODEs_FBA <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * (t > Xstart && t <= Xstop) * (X^Hilln)/(KXZ^Hilln + X^Hilln) * (KZ^Hilln)/(KZ^Hilln + Z^Hilln) - aZ*Z
    dZnoFB <- aZ * (t > Xstart && t <= Xstop) - aZ*ZnoFB
    return(list(c(dZ, dZnoFB)))
  })
}

plotDynamics_FBA <- function(Xstart = 1,
                             Xstop = 6,
                             tmax = 10,
                             dt = 0.1,
                             pars = list(aZ = 1,
                                         bZ = 1,
                                         KXZ = 1,
                                         KZ = 1)) {
  params <- c(Xstart = Xstart, Xstop = Xstop, 
              aZ = pars$aZ, bZ = pars$bZ, KXZ = pars$KXZ, KZ = pars$KZ,
              Hilln = 8)
  iniState <- c(Z=0, ZnoFB = 0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_FBA, params) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0),
           aZ = as.factor(pars$aZ), bZ = as.factor(pars$bZ), 
           KXZ = as.factor(pars$KXZ), KZ = as.factor(pars$KZ)) %>%
    dplyr::select(time, aZ, bZ, KXZ, KZ, X, Z, ZnoFB)
  return(AUC(solution$time, solution$Z, absolutearea = T))
}


genRatioLandscapeData <- function(minRatio, maxRatio) {
  GRID_RES <- 1000
  rational <- MASS:::.rat(seq(minRatio, maxRatio, length.out = GRID_RES),
                          max.denominator = 2000)$rat
  d_grid <- data.frame(aZ = rational[,1], 
                       bZ = rational[,2],
                       KZ = 1, 
                       KXZ = 1) %>% distinct()
  # write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  # 
  # return(runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8) %>%
  #          filter(pheno > 0))
  
  result <- data.frame(pheno = numeric(nrow(d_grid)), 
                       aZ = numeric(nrow(d_grid)), 
                       bZ = numeric(nrow(d_grid)), 
                       KZ = numeric(nrow(d_grid)),  
                       KXZ = numeric(nrow(d_grid)))
  
  for (i in seq_len(nrow(d_grid))) {
    pheno <- plotDynamics_FBA(pars = as.list(d_grid[i,]))
    result[i,] <- c(pheno, d_grid[i,])
  }
  result$fitness <- calcAddFitness(result$pheno, 2, 0.05)
  return(result)
}
plotRatioLandscape <- function(minRatio, maxRatio) {
  d_landscape <- genRatioLandscapeData(minRatio, maxRatio)
  
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

plt_aZbZratio +
  stat_poly_line(colour = "#AAAAAA", linetype = "dashed") +
  stat_poly_eq(use_label(c("eq", "adj.R2", "R2.CI", "p.value"), sep = "*\"; \"*"), 
               label.x = "right", colour = "#000000") -> plt_aZbZratio

plt_aZbZratio

# Optimum phenotype - what aZbZ ratio gives exactly 2: from the above figure,
# it's around 1.25
opt_pheno_ratio <- genRatioLandscapeData(1.24, 1.26)
opt_pheno_ratio$aZbZ <- opt_pheno_ratio$aZ / opt_pheno_ratio$bZ
opt_pheno_ratio[match(max(opt_pheno_ratio$fitness), opt_pheno_ratio$fitness),]

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

# percentage of models with more than 1 step that used only alpha, only beta, or both
d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(evoBybZ = all(value_aZ == 0, na.rm = T),
         evoByaZ = all(value_bZ == 0, na.rm = T)) %>%
  ungroup() %>%
  distinct(seed, .keep_all = T) %>% 
  summarise(percEvoByaZ = mean(evoByaZ),
            percEvoBybZ = mean(evoBybZ),
            percEvoByBoth = 1 - (percEvoByaZ + percEvoBybZ),
            countEvoByaZ = sum(evoByaZ),
            countEvoBybZ = sum(evoBybZ),
            countEvoByBoth = n() - (countEvoByaZ + countEvoBybZ))

# Fig 5 - Phenotype fixed vs segregating effects
# A - correlation between phenotype of fixations only with mean phenotype
group_means <- d_fix_ranked_combined %>% filter(rank > 0) %>%
  group_by(model) %>%
  summarise(ratio = mean(AA_pheno/phenomean),
            CIRatio = CI(AA_pheno/phenomean))

# set seed for geom_jitter
set.seed(seed)
ggplot(d_fix_ranked_combined %>% filter(rank > 0), 
       aes(x = as.factor(model), y = AA_pheno/phenomean)) +
  geom_jitter(size = 0.5, shape = 1, alpha = 0.3) +
  geom_point(data = group_means, aes(y = ratio, colour = model), size = 2) +
  geom_errorbar(data = group_means, aes(y = ratio, ymin = ratio - CIRatio,
                                        ymax = ratio + CIRatio,
                                        colour = model), width = 0.1) +
  scale_colour_paletteer_d("ggsci::nrc_npg", guide = NULL) +
  scale_x_discrete(labels = c("Additive", "NAR")) +
  labs(x = "Model", y = "Fixed effect/mean phenotype ratio") +
  theme_bw() + 
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_fig5
plt_fig5

# Segregating variation: frequency > 25%
d_seg_ranked_combined <- rbind(d_seg_ranked_add, d_seg_ranked)
ggplot(d_seg_ranked_combined %>% filter(Freq >= 0.25),
       aes(x = AA_pheno/aa_pheno, y = Freq, colour = modelindex)) +
  geom_jitter()

ggsave("fig5.png", plt_fig5, device = png)

# Balancing selection example
ggplot(d_com_adapted %>% filter((modelindex == 1 & seed == 1448101263 & mutID == 4624809) | 
                                  (modelindex == 2 & seed == 2270695859 & mutID == 4607309), 
                                gen >= 49000) %>% 
         distinct() %>%
         mutate(gen = gen - 50000, 
                model = if_else(modelindex == 1, "Additive", "NAR")), 
       aes(x = gen, y = Freq, colour = model)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Generations post-optimum shift", y = "Allele frequency", 
       colour = "Model") +
  theme(text = element_text(size = 16), legend.position = "bottom")
ggsave("sfig_balsel.png", device = png)

# View the mutations' effects
View(d_com_adapted %>% filter((modelindex == 1 & seed == 1448101263 & mutID == 4624809) | 
                                (modelindex == 2 & seed == 2270695859 & mutID == 4607309), 
                              gen >= 49000))
