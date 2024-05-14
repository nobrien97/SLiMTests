# R Script to calculate LD from allele frequencies and effect sizes
library(tidyverse)

# Get command line arguments
## 1: run
## 2: chunk
args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])

# Path to write output
WRITE_PATH <- paste0("/scratch/ht96/nb9894/standingVar/calcLD/")
GDATA_PATH <- paste0("/g/data/ht96/nb9894/standingVar/")
FILE_LD <- paste0(WRITE_PATH, "out_LD_", run, ".csv")
FILE_FX <- paste0(WRITE_PATH, "out_fx_", run, ".csv")

# Load functions for loading relatedness/haplotype matrices
source("~/tests/standingVar/calcLD/LDHelpFns.R")

# Load in correct line
d_freqs <- scan(paste0("slim_sharedmutfreqs.csv"), skip = run - 1, 
                    nlines = 1, sep = ",")

model_info <- d_freqs[1:3]
d_freqs <- d_freqs[-(1:3)]

# Get mutational effects
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "standingVarMuts.db"))

# Quantitative data
d_qg <- tbl(con, "slim_qg") %>%
  filter(modelindex == model_info[3], gen == model_info[1], seed == model_info[2]) %>%
  distinct() %>%
  ungroup()
d_qg <- d_qg %>% collect()

d_muts <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model_info[3], gen == model_info[1], seed == model_info[2])
d_muts <- d_muts %>% collect()

d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_muts$mutType <- as.factor(d_muts$mutType)
d_muts$fixGen <- as.numeric(d_muts$fixGen)
d_muts <- d_muts %>% 
  rename(value = effect)
