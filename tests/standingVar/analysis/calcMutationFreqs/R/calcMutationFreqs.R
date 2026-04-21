# Calculates statistics involving mutation data
# Runs per model across all seeds of that model
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(data.table)

# Get command line arguments
## 1: model (ODE, K, or Additive)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

# Paths
R_PATH <- "~/tests/standingVar/calcMutationFreqs/R/"
GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/"

WRITE_PATH <- "/scratch/ht96/nb9894/standingVar/calcMutationFreqs/"
FREQ_FILE <- paste0(WRITE_PATH, "d_mutfreq_", model, ".csv")

# Load combo information
d_combos <- read.table(paste0(R_PATH, "combos.csv"), header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# Load mutation data from database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "standingVarMuts.db"))

# Quantitative data
d_adapted <- tbl(con, "slim_qg") %>%
  filter(gen >= 49500, modelindex == model) %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  filter(any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup()
d_adapted <- d_adapted %>% collect()

# If the model has 0 adapted runs, end early
if (nrow(d_adapted) == 0) {
  print(paste("Model ", run, " never adapted across all seeds, closing R."))
  q(save = "no")
  
}

# Mutation data
adapted_seeds <- unique(d_adapted$seed)
d_muts_adapted <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model, gen >= 49500, seed %in% adapted_seeds) %>%
  select(gen, seed, modelindex, mutType, mutID, freq, count, fixGen)
d_muts_adapted <- d_muts_adapted %>% collect()

# dP/dt
d_dpdt <- d_adapted %>%
  filter(gen >= 50000) %>%
  mutate(optPerc = (phenomean - 1))    # percent to optimum

# Determine when we first reach 25%, 50%, 75%, 90% of the optimum 
# (90% being our cutoff for adaptation)
d_dpdt$optPerc <- cut(d_dpdt$optPerc, c(-Inf, 0.25, 0.5, 0.75, 0.9, Inf),
                      right = F)

# filter by optPerc to select timepoints where populations 
# first reached 50% adapted etc.
d_adapted_optPerc <- d_dpdt %>%
  group_by(optPerc, seed, modelindex) %>%
  filter(row_number() == 1)

# Filter mutations by optPerc generations
# inner join the mutation w/ quantitative data + add model info
d_com_adapted <- inner_join(d_adapted_optPerc, d_muts_adapted, 
                            by = c("gen", "seed", "modelindex"))
d_com_adapted <- AddCombosToDF(d_com_adapted)

d_com_adapted <- d_com_adapted %>%
    select(optPerc, seed, modelindex, mutType, mutID, freq, count, fixGen)

# write
data.table::fwrite(d_com_adapted, 
                   FREQ_FILE, sep = ",", 
                   col.names = F, row.names = F)
