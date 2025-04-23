library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)

DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

# read SFS
d_sfs <- data.table::fread(paste0(DATA_PATH, "calcMutationStats/d_SFS.csv"), header = F, 
                     sep = ",", colClasses = c("factor", "factor", "factor",
                                               "factor", "numeric", "numeric", 
                                               "numeric"), 
                     col.names = c("timePoint", "modelindex", "mutType", "freqBin",
                                   "countFreqBin", "meanValue", "sdValue"), 
                     fill = T)

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


# Combine counts across mutation types, adjust for the number of adapted seeds
# per model
d_sfs_combined <- d_sfs %>%
  group_by(timePoint, modelindex, freqBin) %>%
  summarise(countFreqBin = sum(countFreqBin))

d_sfs <- AddCombosToDF(d_sfs)
ggplot(d_sfs, 
       aes(x = model, y = meanDPDT, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ .) +
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  labs(x = "Model", 
       y = "Mean per-generation change in phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"))
