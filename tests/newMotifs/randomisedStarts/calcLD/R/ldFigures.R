# Load in summary data
library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(ggbeeswarm)
library(cowplot)
library(gganimate)
setwd("/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcLD/R")

source("/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/helperFunctionsAndSetup.R")

d_combos <- read.table("../../../R/combos.csv", header = F,
                       col.names = c("model", "r"))

model_names <- c("NAR", "PAR", "FFLC1", 
                 "FFLI1", "FFBH")
model_labels <- c("NAR", "PAR", "FFL-C1", "FFL-I1", "FFBH")


DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"

# Fitness adjusted frequencies
d_ld <- read.table(paste0(DATA_PATH, "calcLD/out_LD.csv"), header = F,
                   colClasses = c("integer", "factor", "factor",
                                  rep("numeric", times = 26)),
                   col.names = c("gen", "seed", "modelindex", 
                                 "meanD", "sdD", "nD",
                                 "nDP", "nDN", "nDHalf", 
                                 paste0("n", 1:20)))


# Fitness adjusted
d_ld_freq <- data.table::fread(paste0(DATA_PATH, "calcLD/out_LDf.csv"), header = F,
                               colClasses = c("integer", "factor", "factor",
                                            rep("numeric", times = 27)),
                   col.names = c("gen", "seed", "modelindex", "freqBin",
                                 "meanD", "sdD", "nD",
                                 "nDP", "nDN", "nDHalf",
                                 paste0("n", 1:20)))


d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 29)), 
                          col.names = c("gen", "seed", "modelindex", "meanH",
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var", 
                                        "trait3_var", "trait4_var", "dist", 
                                        "dist1", "dist2", "dist3", "dist4", "mean_w",
                                        "var_w", "deltaPheno", "deltaW", 
                                        "meanMC1", "meanMC2", "meanMC3", "meanMC4", 
                                        "meanMC5", "meanMC6", "meanMC7", "meanMC8", 
                                        "meanMC9", "meanMC10", "meanMC11"), 
                          fill = T)
d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.95)) %>%
  mutate(model = factor(model, levels = model_names)) %>%
  ungroup() -> d_qg

d_qg <- d_qg %>% filter(gen >= 49500)


# inner join optPerc
d_ld <- left_join(d_ld, d_qg, by = c("gen", "seed", "modelindex"))

# Proportion of estimates with positive/negative D
d_ld <- d_ld %>%
  group_by(gen, isAdapted, model, r) %>%
  mutate(propDP = nDP / nD,
         propDN = nDN / nD)

# Clean data
d_ld <- d_ld %>% filter(between(meanD, -0.25, 0.25))

# average across replicates
d_ld_sum <- d_ld %>%
  group_by(gen, isAdapted, model, r) %>%
  summarise_at(vars(-seed,-modelindex), list(mean = mean, se = se), na.rm = T)

# plot average
ggplot(d_ld_sum %>% mutate(gen = (gen - 50000) / 1000), 
       aes(x = interaction(gen, model), y = meanD_mean, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ 
                 "Did the population adapt?" + isAdapted) +
  geom_point() +
  scale_x_discrete(guide = "axis_nested") +
  geom_errorbar(aes(ymin = meanD_mean - meanD_se, ymax = meanD_mean + meanD_se)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  labs(x = TeX("Generations post-optimum shift ($x 10^3$) / Model"), y = "Mean LD (D)") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_avg_ld
plt_avg_ld
ggsave("plt_ld_avg.png", device = png, width = 12, height = 6, bg = "white")

# plot sd as well
ggplot(d_ld_sum %>% mutate(gen = (gen - 50000) / 1000), 
       aes(x = interaction(gen, model), y = sdD_mean, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ 
                 "Did the population adapt?" + isAdapted) +
  geom_point() +
  scale_x_discrete(guide = "axis_nested") +
  geom_errorbar(aes(ymin = sdD_mean - sdD_se, ymax = sdD_mean + sdD_se)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  labs(x = TeX("Generations post-optimum shift ($x 10^3$) / Model"), y = "Standard deviation of LD (D)") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_sd_ld
plt_sd_ld
ggsave("plt_ld_sd.png", device = png, width = 12, height = 6, bg = "white")

# Plot together
leg <- get_legend(plt_sd_ld)

plt_ld_sum <- plot_grid(plt_avg_ld + theme(legend.position = "none"), 
                       plt_sd_ld + theme(legend.position = "none"),
                       nrow = 2, labels = "AUTO", label_size = 12)

plt_ld_sum <- plot_grid(plt_ld_sum,
                       leg, nrow = 2, rel_heights = c(1, 0.05))
plt_ld_sum
ggsave("plt_ld_sum.png", device = png, bg = "white",
       width = 16, height = 10)

  
# plot average distributions

bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_dist <- d_ld_sum %>% dplyr::select(gen, model, r, 12:31) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_dist$col <- bins[as.numeric(str_extract(d_ld_dist$col, "[[0-9]]*(?=_)"))]

# Outliers: histogram of all estimates
d_ld_dist_hist <- d_ld %>% dplyr::select(gen, seed, model, r, 10:29) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  ungroup() %>%
  mutate(time_level = as.numeric(factor(gen)), # convert to 1, 2, 3 etc.
         y_offset = time_level * 2) %>% 
  group_by(gen, seed, model, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Plot histogram
# Offset x axis since we group leftwise [x, y), offset is half the bin size, 0.025/2
ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 col_num = as.numeric(col) + 0.025/2,
                                 r_title = "Recombination rate (log10)") %>%
         filter(!isAdapted), 
       aes(x = col_num, y = prop + y_offset, colour = model, 
           group = interaction(gen, col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = interaction(gen, model), colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  scale_y_continuous(
    name = "Proportion (ridgeline-stacked)",
    breaks = unique(d_ld_dist_hist$y_offset),
    labels = scales::comma(unique(d_ld_dist_hist$gen) - 50000)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld
plt_ld
ggsave("plt_ld.png", device = png, width = 9, height = 6)

# Animation
ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 col_num = as.numeric(col) + 0.025/2,
                                 r_title = "Recombination rate (log10)") %>%
         filter(!isAdapted), 
       aes(x = col_num, y = prop, colour = model, 
           group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  transition_states(gen, wrap = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(title = "Generations post-optimum shift: {closest_state}", 
       x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_anim_maladapt
anim_maladapted <- animate(plt_ld_anim_maladapt, nframes = 10, duration = 10,
                        width = 1080, height = 720,
                        renderer = ffmpeg_renderer())
anim_save("plt_ld_maladapted_anim.mp4", anim_maladapted)

ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 col_num = as.numeric(col) + 0.025/2,
                                 r_title = "Recombination rate (log10)") %>%
         filter(isAdapted), 
       aes(x = col_num, y = prop, colour = model, 
           group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  transition_states(gen, wrap = FALSE) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(title = "Generations post-optimum shift: {closest_state}", 
       x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_anim_adapt

anim_adapted <- animate(plt_ld_anim_adapt, nframes = 10, duration = 10,
                        width = 1080, height = 720,
                renderer = ffmpeg_renderer())
anim_save("plt_ld_adapted_anim.mp4", anim_adapted)


#################################
# Frequency adjusted
d_ld_freq <- d_ld_freq %>% 
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex)) %>%
  distinct()

# inner join optPerc
d_ld_freq <- left_join(d_ld_freq, d_qg, by = c("gen", "seed", "modelindex"))

# average across replicates
d_ld_freq_sum <- d_ld_freq %>%
  group_by(timePoint, freqBin, model, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

d_ld_freq_sum <- d_ld_freq_sum %>%
  mutate(freqBin = if_else(freqBin > 0.5, 1 - freqBin, freqBin))

# Reshape for plotting
bins <- seq(-0.25, 0.25, length.out = 21)
d_ld_freq_dist <- d_ld_freq_sum %>% select(timePoint, freqBin, model, r, 13:32) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_freq_dist$col <- bins[as.numeric(str_extract(d_ld_freq_dist$col, "[[0-9]]*(?=_)"))]

# Outliers: histogram of all estimates
d_ld_freq_dist_hist <- d_ld_freq %>% filter(freqBin > 0.1) %>% 
  select(gen, seed, timePoint, model, r, 11:30) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, timePoint, model, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Number of loci doesn't matter
# Doesn't appear to be related to effect size variance
ggplot(d_ld_freq_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                      col_num = as.numeric(col) + 0.025/2,
                                      r_title = "Recombination rate (log10)"), 
       aes(x = col_num, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + log10(r) ~ model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.2) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      guide = "none") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_ld_freq
plt_ld_freq
ggsave("plt_ld_freq.png", device = png, width = 9, height = 4)

# Grid figure
plt_ld_com <- plot_grid(plt_ld + theme(legend.position = "none"),
                            plt_ld_freq + theme(legend.position = "none"),
                            ncol = 1, labels = "AUTO")

plt_ld_com
ggsave("plt_ld_com.png", device = png, bg = "white",
       width = 9, height = 8)
