# Summarise LD information
library(dplyr)
library(data.table)

# Get command line arguments
## 1: model (ODE, K, or Additive)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/calcLD/neutral/"

FILE_NAME <- paste0("/scratch/ht96/nb9894/standingVar/calcLD/LDsummarise/sumLD_d", model, ".csv")

# TODO: check model off by one error
ld_val <- unlist(fread(paste0(GDATA_PATH, "slim_ld_val_d.csv"), 
               skip = model, nrows = 1, sep = ","))

model_info <- ld_val[1:3]
ld_val <- ld_val[-(1:3)]

ld_pos <- unlist(fread(paste0(GDATA_PATH, "slim_ld_pos_d.csv"), 
                skip = model, nrows = 1, sep = ","))[-(1:3)]

# Expand to matrix form
n <- 1024
max_elements <- n*n

decompressLD <- function(pos, val, n = 1024L) {
  # Converts a compressed haplotype matrix to a sparse form
  result <- matrix(integer(n*n), nrow = n, ncol = n)
  result[pos+1] <- val
  return(result)
}

LDmat <- decompressLD(ld_pos, ld_val, n)

# Just get the upper triangle w/o diagonal or zero elements
LDmat <- LDmat[upper.tri(LDmat) & abs(LDmat) > 0]

# Maximum number of non-zero elements

# Calculate means (don't need pos?)
# means/sd with and without the zeroes added (just change denominator to include zeroes)
result <- tibble(
        gen = model_info[1],
        seed = model_info[2],
        modelindex = model_info[3],
        meanD = mean(LDmat),
        sdD = sd(LDmat),
        meanDZeros = sum(LDmat) / (max_elements),
        sdDZeros = sqrt( ( sum((LDmat - meanDZeros)^2) ) / max_elements ),
        nD = length(LDmat),
        nDP = length(LDmat[LDmat > 0.05]),
        nDN = length(LDmat[LDmat < -0.05]),
        nDHalf = length(LDmat[abs(LDmat) > 0.05])
    )

# Add counts of 10% groups for a histogram with 21 bins
labels <- paste0("n", 1:21)
bins <- seq(-1, 1.1, by = 0.1)

LDbins <- cut(LDmat, breaks = bins, right = F)
bin_labels <- levels(LDbins)

for (i in seq_along(bin_labels)) {
    result[,labels[i]] <- length(LDbins[LDbins == bin_labels[i]])
}

# Write output
write.table(result, FILE_NAME, row.names = F, col.names = F)
