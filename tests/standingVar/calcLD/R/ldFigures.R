# Load in summary data

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)

source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"

d_ld <- read.table(paste0(DATA_PATH, "calcLD/sumLD.csv"), header = F,
                          col.names = c("gen", "seed", "modelindex", "meanD", "sdD",
                                        "meanDZeros", "sdDZeros", "nD",
                                        "nDP", "nDN", "nDHalf", 
                                        paste0("n", 1:21)))

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
  summarise_at(vars(-seed,-gen,-modelindex), list(mean, sd))

# plot average distributions

bins <- seq(-1, 1, by = 0.1)
d_ld_dist <- d_ld_sum %>% select(optPerc, model, nloci, tau, r, 14:34) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count")
d_ld_dist$col <- bins[as.numeric(str_extract(d_ld_dist$col, "[[0-9]]*(?=_)"))]

d_ld_dist_sd <- d_ld_sum %>% select(optPerc, model, nloci, tau, r, 43:63) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count_sd")

d_ld_dist$count_sd <- d_ld_dist_sd$count_sd


# Mean counts - after adaptation
{
ggplot(d_ld_dist %>% filter(optPerc == "(0.75, Inf]", tau == 0.0125),
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
  labs(x = "D'", y = "Count", colour = "Model") +
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
  labs(x = "D'", y = "Count", colour = "Model") +
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
  labs(x = "D'", y = "Count", colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) -> D_lrgFX

D_grid <- plot_grid(D_smlFX, D_medFX, D_lrgFX,
                        nrow = 3)

ggsave("LD_grid.png", D_grid, width = 14, height = 30, device = png)
}

# Mean counts - before adaptation
{
  ggplot(d_ld_dist %>% filter(optPerc == "(-Inf,0.25]", tau == 0.0125),
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
    labs(x = "D'", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_smlFX
  
  ggplot(d_ld_dist %>% filter(optPerc == "(-Inf,0.25]", tau == 0.125),
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
    labs(x = "D'", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_medFX
  
  ggplot(d_ld_dist %>% filter(optPerc == "(-Inf,0.25]", tau == 1.25),
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
    labs(x = "D'", y = "Count", colour = "Model") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 14)) -> D_lrgFX
  
  D_grid <- plot_grid(D_smlFX, D_medFX, D_lrgFX,
                      nrow = 3)
  
  ggsave("LD_grid_beforeshift.png", D_grid, width = 14, height = 30, device = png)
}