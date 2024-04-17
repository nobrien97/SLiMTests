# Summarise LD information
library(dplyr)
library(data.table)

# Get command line arguments
## 1: model (ODE, K, or Additive)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/calcLD/"

FILE_NAME <- paste0("/scratch/ht96/nb9894/standingVar/calcLD/LDsummarise/sumLD", model, ".csv")

# TODO: check model off by one error
ld_val <- unlist(fread(paste0(GDATA_PATH, "slim_ld_freq.csv"), 
                       skip = model, nrows = 1, sep = ","))

model_info <- ld_val[1:3]
ld_val <- ld_val[-(1:3)]

# Separate ld data into LD-frequency pairs

ld_frame <- data.frame(D = ld_val[seq(1, length(ld_val), by = 2)],
                       freqBin = ld_val[seq(2, length(ld_val), by = 2)])


# Calculate means (don't need pos?)
# means/sd with and without the zeroes added (just change denominator to include zeroes)
max_elements <- 1024*(1024-1)/2

result <- ld_frame %>%
  mutate(freqBin = factor(freqBin)) %>%
  group_by(freqBin) %>%
  summarise(meanD = mean(D),
            sdD = sd(D),
            meanDZeros = sum(D) / (max_elements),
            sdDZeros = sqrt( ( sum((D - meanDZeros)^2) ) / max_elements ),
            nD = length(D),
            nDP = length(D[D > 0.05]),
            nDN = length(D[D < -0.05]),
            nDHalf = length(D[abs(D) > 0.05]))

result$gen <- model_info[1]
result$seed <- model_info[2]
result$modelindex <- model_info[3]

result <- result %>% relocate(gen, seed, modelindex)

# Add counts of 10% groups for a histogram with 21 bins
labels <- paste0("n", 1:21)
bins <- seq(-0.25, 0.25, length.out = 21)

LDbins <- cut(ld_frame$D, breaks = bins, right = F)
bin_labels <- levels(LDbins)

for (i in seq_along(bin_labels)) {
  result[,labels[i]] <- length(LDbins[LDbins == bin_labels[i]])
}

# Write output
write.table(result, FILE_NAME, row.names = F, col.names = F)
