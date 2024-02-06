# Calculates epistasis for a subset of mutations
library(dplyr)
library(tibble)

# Get command line arguments
## 1: model (ODE, K, or Additive)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

# Path to write output
WRITE_PATH <- paste0("/scratch/ht96/nb9894/standingVar/calcEpistasis/out_e_", model, ".csv")

# Load mutation data from database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(dataPath, "standingVarMuts.db"))

# Get rows for this model only
d_muts <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model, gen >= 49500) %>%
  select(!c(mutID:pos, chi))
d_muts <- d_muts %>% collect()

d_adapted <- tbl(con, "slim_qg") %>%
  filter(gen >= 49500, modelindex == model) %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  filter(any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup()
d_adapted <- d_adapted %>% collect()

# Set data types for incorrect columns
d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_muts$mutType <- as.factor(d_muts$mutType)
d_muts$fixGen <- as.numeric(d_muts$fixGen)
d_muts <- d_muts %>% 
  rename(value = effect)