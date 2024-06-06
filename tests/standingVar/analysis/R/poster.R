# Poster: Do genetic networks modulate the effectiveness of recombination
# in adaptation? 

# Fig 1: adaptive walk - small mutational effects
## Differences between models - why does K+ adapt so much faster?
## Doesn't care about recombination?

# Fig 2: additive variance
## K+ has a lot of VA relative to the others
## Overwhelms the lack of recombination
## Why does so much variance persist?

# Fig 3: epistasis
## K+ has a lot of positive epistasis
## many gene combinations are synergistic in fitness effect
## Recombination can break apart these effects
# creating beneficial allele combinations?

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/R")

# Fig 1
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
r_subsample <- c(1e-10, 1e-1)

source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

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

d_adapted_sum <- d_qg %>% 
  filter(r %in% r_subsample) %>%
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, tau, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

# Mean phenotype figures
  ggplot(d_adapted_sum %>% filter(tau == 0.0125) %>%
           mutate(r_title = "Recombination rate (log10)"),
         aes(x = gen, y = meanPhenomean, colour = model),
         group = as.factor(seed)) +
    facet_nested(r_title + log10(r) ~.) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean, 
                    ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
                alpha = 0.2) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                        labels = c("Additive", "K+", "K-")) +
    scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
         colour = "Model") +
    theme_bw() +
    guides(colour = guide_legend(override.aes=list(linewidth = 6))) +
    theme(legend.position = "bottom", panel.spacing.y = unit(1, "line"),
          text = element_text(size = 14)) -> fig1
  fig1
ggsave("fig1_poster.png", width = 5, height = 5, device = png)


# Fig 2
d_h2_mkr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mrr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mrr_pt1 <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr_pt1.csv"), header = F,
                                  col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                                "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                                "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                                "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                                "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                                "h2_KXZ"))

d_h2_mkr_pt1 <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr_pt1.csv"), header = F,
                                  col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                                "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                                "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                                "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                                "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                                "h2_KXZ"))

d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, optPerc) %>% filter(gen >= 49500)


# Combine
d_h2_mkr <- rbind(d_h2_mkr, d_h2_mkr_pt1) %>% distinct()
d_h2_mrr <- rbind(d_h2_mrr, d_h2_mrr_pt1) %>% distinct()

d_h2_mkr$method <- "mkr"
d_h2_mrr$method <- "mrr"

d_h2 <- rbind(d_h2_mkr, d_h2_mrr)

# Clean data
d_h2 <- d_h2 %>%
  distinct(gen, seed, modelindex, method, .keep_all = T) %>%
  mutate(modelindex = as.factor(modelindex),
         seed = as.factor(seed)) %>%
  drop_na(VA_Z) %>% distinct()

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex"))

d_h2 <- AddCombosToDF(d_h2)

d_h2 <- d_h2 %>% filter(method == "mkr")

boxplot(d_h2$VA_Z)
d_h2_all <- d_h2

# Detect outliers: LOF
library(DMwR2)

lofscores <- lofactor(scale(d_h2$VA_Z), 10)
threshold <- 1.6
outliers <- lofscores > threshold
d_h2 <- d_h2[!outliers,]

# summarise
d_h2_sum <- d_h2 %>%
  filter(r %in% r_subsample) %>%
  group_by(optPerc, model, tau, r, method) %>%
  summarise(meanH2Z = mean(h2_Z, na.rm = T),
            seH2Z = se(h2_Z, na.rm = T),
            meanVAZ = mean(VA_Z, na.rm = T),
            seVAZ = se(VA_Z, na.rm = T))

d_h2_sum %>% filter(method == "mkr", r %in% r_subsample, tau == 0.0125) %>%
  group_by(optPerc, model) %>%
  summarise(diffVA = meanVAZ[2]/meanVAZ[1])
  

# Additive variance
# Again, nloci not important
ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.8)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14), panel.spacing.y = unit(1, "line"),
        legend.position = "bottom") -> fig2
fig2
ggsave("fig2_poster.png", width = 5, height = 5, device = png)

# Fig 3
# Fitness epistasis
# data
d_epi_means <- read.table(paste0(DATA_PATH, "epistasisDensity/d_epi_mean.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "count"))


# Get pairwise differences between in mean epistasis between models 
d_epi_means <- d_epi_means %>% 
  filter(optPerc == "[-Inf,0.25)") %>%
  distinct()

d_epi_means_plt <- AddCombosToDF(d_epi_means %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier
d_epi_means_plt <- d_epi_means_plt %>%
  filter(modelindex != 396, tau == 0.0125, r %in% r_subsample)

ggplot(d_epi_means_plt, aes(x = model, y = meanEW, colour = model)) +
  geom_boxplot(width = 0.5, linewidth = 0.4) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", y = "Average fitness epistasis", colour = "Model") +
  scale_x_discrete(labels = c("K-", "K+", "Additive"), limits = rev) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14)) -> fig3
fig3
ggsave("fig3_poster.png", width = 5, height = 2, device = png)


d_epi_exp <- data.frame(
  genotype = rep(c("ab", "Ab", "aB", "AB"), times = 2),
  fitness = c(0.1, 0.3, 0.3, 0.6, 0.1, 0.3, 0.3, 1),
  model = rep(c("Control", "Positive Interactions"), each = 4)
)

pal <- paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)

ggplot(d_epi_exp %>% mutate(VA = "VA"),
       aes(x = genotype, y = fitness, colour = model)) +
  geom_point(position = position_dodge(0.4), size = 2.5) +
  annotate("segment", x = 4.3, xend = 4.3, y = 0.1, yend = 0.6, 
           arrow = arrow(length = unit(0.02, "npc"), 
                         angle = 90, ends = "both"),
           color = pal[1], linewidth = 1) +
  annotate("segment", x = 4.5, xend = 4.5, y = 0.1, yend = 1, 
           arrow = arrow(length = unit(0.02, "npc"), 
                         angle = 90, ends = "both"),
           color = pal[2], linewidth = 1) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1)) +
  labs(x = "Genotype", y = "Fitness", colour = "Model", size = "VA") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2)) -> fig4
fig4
ggsave("fig4_poster.png", fig4, width = 4, height = 4, device = png)


# Fitness adjusted LD
d_ld_freq <- data.table::fread(paste0(DATA_PATH, "calcLD/fitness_adjusted/out_LDf.csv"), header = F,
                               col.names = c("gen", "seed", "modelindex", "freqBin",
                                             "meanD", "sdD", "nD",
                                             "nDP", "nDN", "nDHalf",
                                             paste0("n", 1:20)))

d_ld_freq <- d_ld_freq %>% 
  filter(seed > 0.01) %>%
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex)) %>%
  distinct()

# inner join optPerc
d_ld_freq <- left_join(d_ld_freq, d_qg, by = c("gen", "seed", "modelindex"))

# Add on variables
d_ld_freq <- AddCombosToDF(d_ld_freq)

# Fold freqbin
d_ld_freq <- d_ld_freq %>%
  mutate(freqBin = if_else(freqBin > 0.5, 1 - freqBin, freqBin))

d_ld_freq_dist_hist <- d_ld_freq %>% select(gen, seed, optPerc, model, nloci, tau, r, 11:30) %>%
  ungroup() %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, optPerc, model, nloci, tau, r) %>%
  distinct() %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

model_labels <- c(
  "Add"="Additive",
  "K"="K+",
  "ODE"="K-"
)

ggplot(d_ld_freq_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                      r_title = "Recombination rate (log10)",
                                      nloci_title = "Number of loci",
                                      tau_title = "Mutational effect size variance") %>%
         filter(nloci == 1024, tau == 0.0125, r %in% r_subsample), 
       aes(x = col, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model, 
               labeller = labeller(model = as_labeller(model_labels))) +
  geom_boxplot(position = position_identity()) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"),
                      guide = "none") +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom")
