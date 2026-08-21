library(tidyverse)
library(data.table)

DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/randomisedStartsM/"
R_PATH <- "~/tests/newMotifs/paper1/randomisedStartsM/calcMutationStats/R/"
COMBO_PATH <- "~/tests/newMotifs/R/"

# Use the right mutate/summarise functions
mutate <- dplyr::mutate
summarise <- dplyr::summarise
select <- dplyr::select

d_combos <- read.table(paste0(COMBO_PATH, "combos.csv"), header = F,
                       col.names = c("model", "r"))
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

model_levels <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")

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


d_fx <- data.table::fread(paste0(DATA_PATH, "calcMutationStats/d_fx_new.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    "factor", "character", "numeric"), 
                          col.names = c("gen", "seed", "modelindex", 
                                        "mutType", "mutID", "s"), 
                          fill = T)

# Join with phenotypic data
d_fx <- d_fx %>% distinct()
d_fx <- left_join(d_fx, d_qg, by = c("gen", "seed", "modelindex"))

# Drop any NAs and Infs
d_fx <- d_fx %>% drop_na(s)

driftBarrier <- 1 / (2 * 5000)

# Filter out lower recombination rates
d_fx <- d_fx %>% filter(log10(r) == -1)

# More extreme ruggedness = more neutral mutations, more canalisation
# Also fewer deleterious mutations
# So look at proportion of deleterious mutations and those above barrier
# Run on HPC - memory usage
d_fx_sum <- d_fx %>%
  drop_na(s) %>%
  group_by(gen, seed, model, isAdapted) %>%
  mutate(nMuts = n(),
         propNeutral = sum((abs(s) < driftBarrier)) / nMuts,
         propDel = sum(s < -driftBarrier) / nMuts,
         propBen = 1 - (propNeutral + propDel)) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("prop"),
               names_to = "mutClass", values_to = "prop") %>%
  group_by(gen, model, isAdapted, mutClass) %>%
  summarise(meanProp = mean(prop),
            CIProp = CI(prop))

saveRDS(d_fx_sum, paste0(DATA_PATH, "calcMutationStats/d_fx_sum.RDS"))
