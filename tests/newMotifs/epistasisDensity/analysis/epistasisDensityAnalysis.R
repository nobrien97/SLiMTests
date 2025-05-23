# Have a look at the epistasis data to see which models differ the most

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
library(legendry)

# functions
source("/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/helperFunctionsAndSetup.R")

# combos
d_combos <- read.table("../../R/combos.csv", header = F,
                            col.names = c("model", "r"))

model_labels <- c("NAR", "PAR", "FFL-C1", "FFL-I1", "FFBH")
model_levels <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")


DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"


# data
d_epi_means_percomp <- read.table(paste0(DATA_PATH, "epistasisDensity/d_epi_mean.csv"), header = F, sep = ",",
                          col.names = c("timePoint", "modelindex", "isAdapted", "molComp", 
                                        "meanEW", "sdEW", "minEW", "maxEW", "q025EW",
                                        "q25EW", "q50EW", "q75EW", "q975EW", 
                                        "count", "freqAboveDB"))

d_epi_means <- read.table(paste0(DATA_PATH, "epistasisDensity/d_epi_nomolcomp_mean.csv"), header = F, sep = ",",
                          col.names = c("timePoint", "modelindex", "isAdapted",
                                        "meanEW", "sdEW", "minEW", "maxEW", "q025EW",
                                        "q25EW", "q50EW", "q75EW", "q975EW", 
                                        "count", "freqAboveDB"))

d_epi_means <- d_epi_means %>%
  mutate(isAdapted = as.logical(isAdapted)) 

d_epi_means_plt <- AddCombosToDF(d_epi_means %>% 
                                   mutate(modelindex = as.factor(modelindex))) %>%
  mutate(model = factor(model, levels = model_levels))

# Average over r, no effect
d_epi_means_plt %>%
  group_by(model, timePoint, isAdapted) %>%
  summarise(meanEW = mean(meanEW),
            CIEW = mean((sdEW / sqrt(count)) * qnorm(0.975)),
            minEW = mean(minEW),
            q025EW = mean(q025EW),
            q25EW = mean(q25EW),
            q50EW = mean(q50EW),
            q75EW = mean(q75EW),
            q975EW = mean(q975EW),
            maxEW = mean(maxEW),
            meanFreqAboveDB = mean(freqAboveDB),
            CIFreqAboveDB = CI(freqAboveDB),
            n = sum(count)) -> d_epi_means_plt2

ggplot(d_epi_means_plt2 %>% mutate(timePoint = (timePoint - 50000) / 1000) %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = interaction(timePoint, model), y = meanFreqAboveDB, fill = model, colour = model)) +
  facet_nested(.~"Did the population adapt?" + isAdapted) +
  geom_point() +
  scale_x_discrete(guide = "axis_nested") +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.2) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                         5, direction = 1),
                    labels = model_labels, guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels, guide = "none") +
  labs(x = "Model", y = "Proportion of samples with\nnon-negligible fitness epistasis", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_ew_freq_db_time
plt_ew_freq_db_time
ggsave("plt_ew_freq_time.png", width = 12, height = 4, device = png)

# Without time either
d_epi_means_plt %>%
  group_by(model, isAdapted) %>%
  summarise(meanEW = mean(meanEW),
            CIEW = mean((sdEW / sqrt(count)) * qnorm(0.975)),
            minEW = mean(minEW),
            q025EW = mean(q025EW),
            q25EW = mean(q25EW),
            q50EW = mean(q50EW),
            q75EW = mean(q75EW),
            q975EW = mean(q975EW),
            maxEW = mean(maxEW),
            meanFreqAboveDB = mean(freqAboveDB),
            CIFreqAboveDB = CI(freqAboveDB),
            n = sum(count)) -> d_epi_means_plt2

ggplot(d_epi_means_plt2 %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = meanFreqAboveDB, fill = model, colour = model)) +
  facet_nested(.~"Did the population adapt?" + isAdapted) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.2) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                         5, direction = 1),
                    labels = model_labels, guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels, guide = "none") +
  labs(x = "Model", y = "Proportion of samples with\nnon-negligible fitness epistasis", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_ew_freq_db
plt_ew_freq_db
ggsave("plt_ew_freq.png", width = 6, height = 4, device = png)


# Per molecular component figures
d_epi_means_percomp <- d_epi_means_percomp %>%
  mutate(isAdapted = as.logical(isAdapted)) 

d_epi_means_pc_plt <- AddCombosToDF(d_epi_means_percomp %>% 
                                   mutate(modelindex = as.factor(modelindex))) %>%
  mutate(model = factor(model, levels = model_levels))

d_epi_means_pc_plt <- ReconcileMutTypeComparisonNames(d_epi_means_pc_plt)

# Average across time
d_epi_means_pc_plt %>%
  group_by(model, isAdapted, molCompNamed) %>%
  summarise(minEW = mean(minEW),
            q025EW = mean(q025EW),
            q25EW = mean(q25EW),
            q50EW = mean(q50EW),
            q75EW = mean(q75EW),
            q975EW = mean(q975EW),
            maxEW = mean(maxEW),
            meanFreqAboveDB = mean(freqAboveDB),
            CIFreqAboveDB = se(freqAboveDB),
            n = sum(count)) -> d_epi_means_pc_plt_sum

# Look at all comparisons of molecular components for adapted vs maladapted
ggplot(d_epi_means_pc_plt_sum %>% filter(!isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = meanFreqAboveDB, fill = model, colour = model)) +
  facet_wrap(isAdapted~molCompNamed, 
               labeller = labeller(molCompNamed = label_parsed)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.5) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                         5, direction = 1),
                    labels = model_labels, guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels, guide = "none") +
  labs(x = "Model", y = TeX("Samples with non-negligible fitness epistasis (%)"), 
       fill = "Model") +
  #scale_x_discrete(labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plt_ew_freq_pc_db
plt_ew_freq_pc_db
ggsave("plt_ew_freq_pc_molComp_maladapted.png", width = 15, height = 15, device = png)

ggplot(d_epi_means_pc_plt_sum %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = meanFreqAboveDB, fill = model, colour = model)) +
  facet_wrap(isAdapted~molCompNamed, 
             labeller = labeller(molCompNamed = label_parsed)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.5) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                         5, direction = 1),
                    labels = model_labels, guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels, guide = "none") +
  labs(x = "Model", y = TeX("Samples with non-negligible fitness epistasis (%)"), 
       fill = "Model") +
  #scale_x_discrete(labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plt_ew_freq_pc_db
plt_ew_freq_pc_db
ggsave("plt_ew_freq_pc_molComp_adapted.png", width = 15, height = 15, device = png)


# Average across all comparisons - identify components that are consistently epistatic
d_epi_means_pc_plt %>%
  # Filter out any frequencies with counts < 10
  filter(count >= 10) %>%
  pivot_longer(cols = c(molCompNamed1, molCompNamed2), names_to = "pos", 
               values_to = "molComp") %>%
  group_by(model, isAdapted, molComp) %>%
  summarise(n = sum(count),
            minEW = mean(minEW),
            q025EW = mean(q025EW),
            q25EW = mean(q25EW),
            q50EW = mean(q50EW),
            q75EW = mean(q75EW),
            q975EW = mean(q975EW),
            maxEW = mean(maxEW),
            meanFreqAboveDB = mean(freqAboveDB),
            CIFreqAboveDB = CI(freqAboveDB),
            .groups = "drop") -> d_epi_means_pc_plt_sum

# Plot 
ggplot(d_epi_means_pc_plt_sum %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = meanFreqAboveDB, fill = model, colour = model)) +
  facet_nested("Molecular component" + molComp~"Has the population adapted?" + isAdapted, 
             labeller = labeller(molComp = label_parsed)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFreqAboveDB - CIFreqAboveDB,
                    ymax = meanFreqAboveDB + CIFreqAboveDB), width = 0.5) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                         5, direction = 1),
                    labels = model_labels, guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels, guide = "none") +
  labs(x = "Model", y = TeX("Samples with non-negligible fitness epistasis (%)"), 
       fill = "Model") +
  #scale_x_discrete(labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plt_ew_freq_pc_db
plt_ew_freq_pc_db
ggsave("plt_ew_freq_pc_av_molComp.png", width = 8, height = 9, device = png)

# How many times does epistasis change sign?
d_epi_sign <- read.table(paste0(DATA_PATH, "epistasisDensity/d_epi_sign_nomolcomp.csv"), header = F, sep = ",",
                          col.names = c("seed", "modelindex",
                                        "n", "nChangesEW", 
                                        "nChangesEW_s"),
                         colClasses = c("factor", "factor", 
                                        rep("integer", times = 3)))

# Is it adapted? Attach info from qg
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
  mutate(model = factor(model, levels = model_levels)) %>%
  ungroup() -> d_qg

d_qg <- d_qg %>% filter(gen >= 49500)

d_epi_sign <- left_join(d_epi_sign, 
                        d_qg %>% filter(gen == 60000), 
                         by = c("seed", "modelindex"))


d_epi_sign_mean <- d_epi_sign %>%
  group_by(model, r, isAdapted) %>%
  summarise(meanEWChanges = mean(nChangesEW),
            CIEWChanges = CI(nChangesEW),
            meanEWsChanges = mean(nChangesEW_s),
            CIEWsChanges = CI(nChangesEW_s))

ggplot(d_epi_sign_mean %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = meanEWsChanges, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~"Did the population adapt?" + isAdapted) +
  geom_point() +
  geom_errorbar(aes(ymin = meanEWsChanges - CIEWsChanges,
                    ymax = meanEWsChanges + CIEWsChanges), 
                position = position_dodge(0.9)) +
  scale_x_discrete(guide = "axis_nested") + 
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels) +
  labs(x = "Model", 
       y = "Mean number of sign changes in fitness epistasis", colour = "Model") +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_ew_s_sign  
plt_ew_s_sign
ggsave("plt_ew_s_sign.png", width = 6, height = 6, device = png)

# Combine these figures to compare across model combinations
plot_grid(plt_ew_freq_db,
          plt_ew_s_sign + theme(legend.position = "none"),
          labels = "AUTO",
          ncol = 2, rel_widths = c(0.9, 1))

ggsave("plt_ew.png", width = 12, height = 5, device = png)

# Table of results
print(xtable(d_epi_means_pc_plt_sum %>% 
        mutate(minEW = round(minEW, digits = 4),
               maxEW = round(maxEW, digits = 4),
          range = paste0(minEW, " - ", maxEW)) %>%
        select(model, molComp, q50EW, range, meanFreqAboveDB, n)),
      include.rownames = F
)
d_epi_means_plt2 %>% select(model, minEW, q50EW, maxEW, meanFreqAboveDB, n)
d_epi_sign_mean


print(xtable(d_epi_means_plt2 %>% 
               mutate(minEW = round(minEW, digits = 4),
                      maxEW = round(maxEW, digits = 4),
                      range = paste0(minEW, " - ", maxEW)) %>%
               select(model, q50EW, range, meanFreqAboveDB, n)),
      include.rownames = F
)

