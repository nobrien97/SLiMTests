# Load in summary data
library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/R")

source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"

d_ld <- read.table(paste0(DATA_PATH, "calcLD/sumLD_d.csv"), header = F,
                          col.names = c("gen", "seed", "modelindex", "meanD", "sdD",
                                        "meanDZeros", "sdDZeros", "nD",
                                        "nDP", "nDN", "nDHalf", 
                                        paste0("n", 1:21)))

d_ld_freq <- data.table::fread(paste0(DATA_PATH, "calcLD/sumLD_df.csv"), header = F,
                   col.names = c("gen", "seed", "modelindex", "freqBin", 
                                 "meanD", "sdD",
                                 "meanDZeros", "sdDZeros", "nD",
                                 "nDP", "nDN", "nDHalf", 
                                 paste0("n", 1:20)))

d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                     sep = ",", colClasses = c("integer", "factor", "factor", 
                                               rep("numeric", times = 12)), 
                     col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                   "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                   "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                     fill = T)

# Add on optPerc
d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg <- d_qg %>% select(gen, seed, modelindex, optPerc) %>% filter(gen >= 49500)

# First: convert gen to optPerc: will be a bit complicated since they don't all
# have 5 samples - need to match against d_qg

# Remove rows with invalid seeds
d_ld <- d_ld %>% 
  filter(seed > 0.01) %>%
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex))

# inner join optPerc
d_ld <- left_join(d_ld, d_qg, by = c("gen", "seed", "modelindex"))

# Add on variables
d_ld <- AddCombosToDF(d_ld)

# average across replicates
d_ld_sum <- d_ld %>%
  group_by(optPerc, model, nloci, tau, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

# plot average distributions

bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_dist <- d_ld_sum %>% select(optPerc, model, nloci, tau, r, 15:34) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_dist$col <- bins[as.numeric(str_extract(d_ld_dist$col, "[[0-9]]*(?=_)"))]

d_ld_dist_sd <- d_ld_sum %>% select(optPerc, model, nloci, tau, r, 43:62) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count_sd")

d_ld_dist$count_sd <- d_ld_dist_sd$count_sd

# Outliers: histogram of all estimates
d_ld_dist_hist <- d_ld %>% select(gen, seed, optPerc, model, nloci, tau, r, 12:32) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, optPerc, model, nloci, tau, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 r_title = "Recombination rate (log10)",
                                 nloci_title = "Number of loci") %>%
         filter(between(col, -0.1, 0.1), as.numeric(optPerc) < 3), 
       aes(x = col, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ nloci_title + nloci) +
  geom_boxplot(position = position_identity()) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
    #position = position_dodge(width = 0.9)
  ) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  #scale_x_continuous(labels = bins[7:15]) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom")


# Mean counts - after adaptation, mean plots
{
  ggplot(d_ld_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 0.0125),
         aes(x = meanD_mean, y = nD_mean, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    #geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd)) +
    #geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd)) +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = "Log Mean Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> plt_mean_LD_sml
  
  ggplot(d_ld_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 0.125),
         aes(x = meanD_mean, y = nD_mean, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = "Log Mean Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> plt_mean_LD_med
  
  ggplot(d_ld_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 1.25),
         aes(x = meanD_mean, y = nD_mean, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = "Log Mean Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> plt_mean_LD_lrg
  
  mean_LD_grid <- plot_grid(plt_mean_LD_sml, plt_mean_LD_med, plt_mean_LD_lrg,
                          nrow = 3)
  ggsave("mean_LD_grid_end.png", mean_LD_grid, width = 14, height = 30, device = png)
}

# Mean counts - before adaptation, mean plots
{
  ggplot(d_ld_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(as.numeric(optPerc) < 3, tau == 0.0125),
         aes(x = meanD_mean, y = nD_mean, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    #geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd)) +
    #geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd)) +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = "Log Mean Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> plt_mean_LD_sml
  
  ggplot(d_ld_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(as.numeric(optPerc) < 3, tau == 0.125),
         aes(x = meanD_mean, y = nD_mean, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = "Log Mean Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> plt_mean_LD_med
  
  ggplot(d_ld_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(as.numeric(optPerc) < 3, tau == 1.25),
         aes(x = meanD_mean, y = nD_mean, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = "Log Mean Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> plt_mean_LD_lrg
  
  mean_LD_grid <- plot_grid(plt_mean_LD_sml, plt_mean_LD_med, plt_mean_LD_lrg,
                            nrow = 3)
  ggsave("mean_LD_grid_beforeshift.png", mean_LD_grid, width = 14, height = 30, device = png)
  
}


# Mean counts - after adaptation, distribution
{
ggplot(d_ld_dist %>% mutate(count = log10(count), count_sd = log10(count_sd)) %>% 
         filter(optPerc == "(0.75, Inf]", tau == 0.0125),
       aes(x = col, y = count, colour = model)) +
  facet_grid(log10(r)~nloci) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
  ggtitle("Tau = 0.0125") +
  labs(x = "D", y = TeX("Mean Count (log_10)"), colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> D_smlFX

ggplot(d_ld_dist %>% mutate(count = log10(count), count_sd = log10(count_sd)) %>% 
         filter(optPerc == "(0.75, Inf]", tau == 0.125),
       aes(x = col, y = count, colour = model)) +
  facet_grid(log10(r)~nloci) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
  ggtitle("Tau = 0.125") +
  labs(x = "D", y = TeX("Mean Count (log_10)"), colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> D_medFX

ggplot(d_ld_dist %>% mutate(count = log10(count), count_sd = log10(count_sd)) %>%
         filter(optPerc == "(0.75, Inf]", tau == 1.25),
       aes(x = col, y = count, colour = model)) +
  facet_grid(log10(r)~nloci) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("Additive", "K+", "K-")) +
  geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
  ggtitle("Tau = 1.25") +
  labs(x = "D", y = TeX("Mean Count (log_10)"), colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> D_lrgFX

D_grid <- plot_grid(D_smlFX, D_medFX, D_lrgFX,
                        nrow = 3)

ggsave("LD_grid_d_end.png", D_grid, width = 14, height = 30, device = png)
}


# Mean counts - before adaptation, distribution
{
  ggplot(d_ld_dist %>% mutate(count = log10(count), count_sd = log10(count_sd)) %>% 
           filter(optPerc == "(-Inf,0.25]", tau == 0.0125),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_smlFX
  
  ggplot(d_ld_dist %>% mutate(count = log10(count), count_sd = log10(count_sd)) %>% 
           filter(optPerc == "(-Inf,0.25]", tau == 0.125),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_medFX
  
  ggplot(d_ld_dist %>% mutate(count = log10(count), count_sd = log10(count_sd)) %>%
           filter(optPerc == "(-Inf,0.25]", tau == 1.25),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_lrgFX
  
  D_grid <- plot_grid(D_smlFX, D_medFX, D_lrgFX,
                      nrow = 3)
  
  ggsave("LD_grid_d_beforeshift.png", D_grid, width = 14, height = 30, device = png)
}

# 


##### 
# Frequency adjusted
# Remove rows with invalid seeds
d_ld_freq <- d_ld_freq %>% 
  filter(seed > 0.01) %>%
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex))

# inner join optPerc
d_ld_freq <- left_join(d_ld_freq, d_qg, by = c("gen", "seed", "modelindex"))

# Add on variables
d_ld_freq <- AddCombosToDF(d_ld_freq)

# average across replicates
d_ld_freq_sum <- d_ld_freq %>%
  group_by(optPerc, freqBin, model, nloci, tau, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

# plot average distributions

bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_freq_dist <- d_ld_freq_sum %>% select(optPerc, model, nloci, tau, r, 15:34) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_freq_dist$col <- bins[as.numeric(str_extract(d_ld_freq_dist$col, "[[0-9]]*(?=_)"))]

d_ld_freq_dist_sd <- d_ld_freq_sum %>% select(optPerc, model, nloci, tau, r, 43:62) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count_sd")

d_ld_freq_dist$count_sd <- d_ld_freq_dist_sd$count_sd

# Normalise mean counts by max elements
MAX_ELEMENTS <- 1024 * 1024
d_ld_freq_sum <- d_ld_freq_sum %>%
  mutate(nD_maxel_prop = nD_mean / MAX_ELEMENTS)

# Mean counts - after adaptation, sum data: check variance of means
{
  ggplot(d_ld_freq_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 0.0125),
         aes(x = meanD_mean, y = nD_mean, colour = model, size = as.factor(freqBin))) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_size_manual(values = c(seq(1, 5, length.out = 5), 
                                 seq(5, 1, length.out = 5))) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    # geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd), linewidth = 0.5,
    #               linetype = "dashed") +
    # geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd),
    #               linewidth = 0.5, linetype = "dashed") +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), 
         colour = "Model", size = "Allele Frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_mean_freq_sml
  
  ggplot(d_ld_freq_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 0.125),
         aes(x = meanD_mean, y = nD_mean, colour = model, size = as.factor(freqBin))) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_size_manual(values = c(seq(1, 5, length.out = 5), 
                                 seq(5, 1, length.out = 5))) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    # geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd), linewidth = 0.5,
    #               linetype = "dashed") +
    # geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd),
    #               linewidth = 0.5, linetype = "dashed") +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), 
         colour = "Model", size = "Allele Frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_mean_freq_med
  
  ggplot(d_ld_freq_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 1.25),
         aes(x = meanD_mean, y = nD_mean, colour = model, size = as.factor(freqBin))) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_size_manual(values = c(seq(1, 5, length.out = 5), 
                                 seq(5, 1, length.out = 5))) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    # geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd), linewidth = 0.5,
    #               linetype = "dashed") +
    # geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd),
    #               linewidth = 0.5, linetype = "dashed") +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), 
         colour = "Model", size = "Allele Frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_mean_freq_lrg
  
  D_freq_grid <- plot_grid(D_mean_freq_sml, D_mean_freq_med, D_mean_freq_lrg,
                      nrow = 3)
  
  ggsave("LD_freq_grid_d_end.png", D_freq_grid, width = 14, height = 30, device = png)
  
}

# Mean counts - before adaptation, sum data: check variance of means
{
  ggplot(d_ld_freq_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(as.numeric(optPerc) < 3, tau == 0.0125),
         aes(x = meanD_mean, y = nD_mean, colour = model, size = as.factor(freqBin))) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_size_manual(values = c(seq(1, 5, length.out = 5), 
                                 seq(5, 1, length.out = 5))) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    # geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd), linewidth = 0.5,
    #               linetype = "dashed") +
    # geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd),
    #               linewidth = 0.5, linetype = "dashed") +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), 
         colour = "Model", size = "Allele Frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_mean_freq_sml
  
  ggplot(d_ld_freq_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(as.numeric(optPerc) < 3, tau == 0.125),
         aes(x = meanD_mean, y = nD_mean, colour = model, size = as.factor(freqBin))) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_size_manual(values = c(seq(1, 5, length.out = 5), 
                                 seq(5, 1, length.out = 5))) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    # geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd), linewidth = 0.5,
    #               linetype = "dashed") +
    # geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd),
    #               linewidth = 0.5, linetype = "dashed") +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), 
         colour = "Model", size = "Allele Frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_mean_freq_med
  
  ggplot(d_ld_freq_sum %>% 
           mutate(nD_mean = log(nD_mean), nD_sd = log(nD_sd)) %>%
           filter(as.numeric(optPerc) < 3, tau == 1.25),
         aes(x = meanD_mean, y = nD_mean, colour = model, size = as.factor(freqBin))) +
    facet_grid(log10(r)~nloci) +
    geom_point(shape = 1) +
    coord_cartesian(xlim = c(-0.25, 0.25)) +
    scale_size_manual(values = c(seq(1, 5, length.out = 5), 
                                 seq(5, 1, length.out = 5))) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    # geom_errorbar(aes(ymin = nD_mean - nD_sd, ymax = nD_mean + nD_sd), linewidth = 0.5,
    #               linetype = "dashed") +
    # geom_errorbar(aes(xmin = meanD_mean - meanD_sd, xmax = meanD_mean + meanD_sd),
    #               linewidth = 0.5, linetype = "dashed") +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = TeX("Mean Count (log_10)"), 
         colour = "Model", size = "Allele Frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_mean_freq_lrg
  
  D_freq_grid <- plot_grid(D_mean_freq_sml, D_mean_freq_med, D_mean_freq_lrg,
                           nrow = 3)
  
  ggsave("LD_freq_grid_d_beforeshift.png", D_freq_grid, width = 14, height = 30, device = png)
  
}

# Mean counts - after adaptation
{
  ggplot(d_ld_freq_dist %>% 
           mutate(freqBelow50 = (freqBin < 0.5)) %>%
           filter(optPerc == "(0.75, Inf]", tau == 0.0125),
         aes(x = col, y = count, colour = model, linetype = freqBelow50)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = "Count", colour = "Model", linetype = "Allele frequency") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_smlFX
  
  ggplot(d_ld_dist %>% filter(optPerc == "(0.75, Inf]", tau == 0.125),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_medFX
  
  ggplot(d_ld_dist %>% filter(optPerc == "(0.75, Inf]", tau == 1.25),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_lrgFX
  
  D_grid <- plot_grid(D_smlFX, D_medFX, D_lrgFX,
                      nrow = 3)
  
  ggsave("LD_grid_d.png", D_grid, width = 14, height = 30, device = png)
}

# Mean counts - before adaptation
{
  ggplot(d_ld_dist %>% filter(as.numeric(optPerc) < 3, tau == 0.0125),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 0.0125") +
    labs(x = "D", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_smlFX
  
  ggplot(d_ld_dist %>% filter(as.numeric(optPerc) < 3, tau == 0.125),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 0.125") +
    labs(x = "D", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_medFX
  
  ggplot(d_ld_dist %>% filter(as.numeric(optPerc) < 3, tau == 1.25),
         aes(x = col, y = count, colour = model)) +
    facet_grid(log10(r)~nloci) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                           breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                           breaks = NULL, labels = NULL)) +
    scale_colour_paletteer_d("nationalparkcolors::Badlands",
                             labels = c("Additive", "K+", "K-")) +
    geom_errorbar(aes(ymin = count - count_sd, ymax = count + count_sd)) +
    ggtitle("Tau = 1.25") +
    labs(x = "D", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_lrgFX
  
  D_grid <- plot_grid(D_smlFX, D_medFX, D_lrgFX,
                      nrow = 3)
  
  ggsave("LD_grid_d_beforeshift.png", D_grid, width = 14, height = 30, device = png)
}