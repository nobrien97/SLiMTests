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

GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/calcLD/"

FILE_NAME <- paste0("/scratch/ht96/nb9894/standingVar/calcLD/sumLD", model, ".csv")

# TODO: check model off by one error
# ld_pos <- scan(paste0("slim_ld_pos.csv"), skip = model, sep = ",")


ld_val <- scan(paste0(GDATA_PATH, "slim_ld_val.csv"), skip = model, sep = ",")
model_info <- ld_val[1:3]
ld_val <- ld_val[-(1:3)]

# Maximum number of non-zero element
max_elements <- 1024*1024

# Calculate means (don't need pos?)
# means/sd with and without the zeroes added (just change denominator to include zeroes)
result <- tibble() %>%
    mutate(
        gen = model_info[1],
        seed = model_info[2],
        modelindex = model_info[3],
        meanD = mean(ld_val),
        sdD = sd(ld_val),
        nD = length(ld_val),
        propNonZeroD = length(ld_val)/max_elements,
        meanDZeros = sum(ld_val) / (max_elements),
        sdDZeros = sqrt( ( sum((ld_val - meanDZeros)^2) ) / max_elements )
    )

# Write output
write.table(result, FILE_NAME, row.names = F, col.names = F)
