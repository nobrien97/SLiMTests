library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

# read SFS
d_sfs <- data.table::fread(paste0(DATA_PATH, "mutationStats/d_SFS.csv"), header = F, 
                     sep = ",", colClasses = c("factor", "factor", "factor",
                                               "factor", "numeric", "numeric", 
                                               "numeric"), 
                     col.names = c("optPerc", "modelindex", "mutType", "freqBin",
                                   "countFreqBin", "meanValue", "sdValue"), 
                     fill = T)

# Read qg data
d_qg <- fread::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                     sep = ",", colClasses = c("integer", "factor", "factor", 
                                               rep("numeric", times = 12)), 
                     col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                   "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                   "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                     fill = T)



# Combine counts across mutation types, adjust for the number of adapted seeds
# per model
d_sfs_combined <- d_sfs %>%
  group_by(optPerc, modelindex, freqBin) %>%
  summarise(countFreqBin = sum(countFreqBin))

d_sfs <- AddCombosToDF(d_sfs)
ggplot(d_sfs %>% filter(tau == 0.0125, optPerc == "(0.75, Inf]"), 
       aes(x = optPerc, y = meanDPDT, colour = model)) +
  facet_grid(log10(r)~nloci) +
  geom_line() +
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
  labs(x = "Generations post-optimum shift", 
       y = "Mean per-generation change in phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        panel.spacing.x = unit(2, "lines")) -> dpdt_grid_lrgFX