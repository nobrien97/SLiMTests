# Figures of adaptive walk

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- read.table(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
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

View(d_qg %>% group_by(model, nloci, tau, r) %>%
       filter(gen == 59950) %>%
       summarise(n = n(),
                 pAdapted = mean(isAdapted),
                 CIAdapted = CI(isAdapted)))
