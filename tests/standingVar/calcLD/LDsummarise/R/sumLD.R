# Summarise LD information

library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(data.table)

# Get command line arguments
## 1: model (ODE, K, or Additive)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

R_PATH <- "~/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))
GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/calcLD/"

# TODO: check model off by one error
ld_pos <- scan(paste0("slim_ld_pos.csv"), skip = model, sep = ",")

model_info <- ld_pos[1:3]
ld_pos <- ld_pos[-(1:3)]

ld_val <- scan(paste0("slim_ld_val.csv"), skip = model, sep = ",")

# Expand to matrix

