# R Script to calculate LD from allele frequencies and effect sizes
library(tidyverse)
library(data.table)

# Get command line arguments
## 1: run
## 2: chunk
args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])

# Path to write output
WRITE_PATH <- paste0("/scratch/ht96/nb9894/standingVar/sanityChecks/calcLDR/")
GDATA_PATH <- paste0("/g/data/ht96/nb9894/standingVar/sanityChecks/")
#GDATA_PATH <- paste0("/mnt/d/SLiMTests/tests/standingVar/calcLD/")
FILE_LD <- paste0(WRITE_PATH, "out_LD_", run, ".csv")
FILE_LD_F <- paste0(WRITE_PATH, "out_LDf_", run, ".csv")
FILE_LD_TABLE <- paste0(WRITE_PATH, "out_LD_raw_", run, ".csv")

FNS_PATH <- "~/tests/standingVar/calcMutationStats/R/"
#FNS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"

source(paste0(FNS_PATH, "helperFunctionsAndSetup.R"))

# Load functions for loading relatedness/haplotype matrices
source("~/tests/standingVar/sanityChecks/calcLD/R/LDHelpFns.R")
#source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/R/LDHelpFns.R")

# Load in correct line
d_freqs <- scan(paste0(GDATA_PATH, "slim_sharedmutfreqs.csv"), skip = run, 
                    nlines = 1, sep = ",")

# d_freqs <- scan(paste0(GDATA_PATH, "frequencies_test/slim_sharedmutfreqs.csv"), skip = run, 
#                 nlines = 1, sep = ",")

model_info <- d_freqs[1:3]
d_freqs <- d_freqs[-(1:3)]

# Model info into separate variables for dbplyr
run_gen <- model_info[1]
run_seed <- model_info[2]
run_modelindex <- model_info[3]

# Slight mismatch after loading population: adjust 
if (run_gen == 50001) {
  run_gen <- 50000
}

# Get number of entries
n <- length(d_freqs) / 7

if (n != as.integer(n)) {
  stop("Invalid data: length of dataset should be a multiple of 7")
  q(save = "n")
}

# Separate data
mutAID <- as.integer(d_freqs[seq(from = 1, to = n*7, by = 7)])
mutBID <- as.integer(d_freqs[seq(from = 2, to = n*7, by = 7)])
mutpA <- d_freqs[seq(from = 3, to = n*7, by = 7)]
mutpB <- d_freqs[seq(from = 4, to = n*7, by = 7)]
mutpab <- d_freqs[seq(from = 5, to = n*7, by = 7)]
mutpaB <- d_freqs[seq(from = 6, to = n*7, by = 7)]
mutpAb <- d_freqs[seq(from = 7, to = n*7, by = 7)]
mutpAB <- 1 - (mutpab + mutpAb + mutpaB)
mutIDs <- unique(c(mutAID, mutBID))

mut_freqs <- list(pAB = mutpAB,
                  pAb = mutpAb,
                  paB = mutpaB,
                  pab = mutpab)

# Get mutational effects
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "standingVarMuts.db"))
# Load mutation data: matching mutIDs, or fixations for the parental genotype
d_muts <- tbl(con, "slim_muts") %>% 
  filter(gen == run_gen, seed == run_seed, modelindex == run_modelindex, 
         mutID %in% mutIDs | fixGen != "NA")
d_muts <- d_muts %>% collect() %>% distinct()
# write_csv(d_muts, "d_muts_test.csv")
# d_muts <- read_csv(paste0(GDATA_PATH, "frequencies_test/d_muts_test.csv"))

d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_muts$mutType <- as.factor(d_muts$mutType)
d_muts$fixGen <- as.numeric(d_muts$fixGen)
d_muts <- d_muts %>% 
  rename(value = effect)

# If there are no segregating mutations, exit early
if (nrow(d_muts %>% filter(is.na(fixGen))) < 2) {
  stop(paste0("<2 segregating mutations for run", run, "!"))
  q(save = "n")
}

# Read in combos
d_combos <- read.table("~/tests/standingVar/sanityChecks/R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "w", "model"))

d_muts <- AddCombosToDF(d_muts)

# Calculate parental genotype, intermediates, and derived
## In example parental effect is 0
d_rank <- PairwiseFitnessRank(d_muts %>% filter(!is.na(fixGen)), 
                             d_muts %>% filter(is.na(fixGen)),
                             mutAID, mutBID)

relabeled_freqs <- RelabelGenotypeFrequencies(d_rank, mut_freqs)

# With the relabelled frequencies, calculate D
D <- CalcLD(relabeled_freqs)

# Check if the comparison is valid
mutA_valid <- match(mutAID, d_rank$mutIDA, nomatch = 0) > 0
mutB_valid <- match(mutBID, d_rank$mutIDB, nomatch = 0) > 0

d_LD <- data.frame(gen = rep(model_info[1], times = length(D)),
                   seed = rep(model_info[2], times = length(D)),
                   modelindex = rep(model_info[3], times = length(D)),
                   mutID_A = d_rank$mutIDA,
                   mutID_B = d_rank$mutIDB,
                   freqDiff = mutpA[mutA_valid & mutB_valid] - mutpB[mutA_valid & mutB_valid],
                   freqBin = round(mutpA[mutA_valid & mutB_valid], 1),
                   freq_intermediateFit = d_rank$wparAb * d_rank$wparaB,
                   freq_extremeFit = d_rank$wparAB * d_rank$wparab,
                   mutType_AB = rep("3_3", times = length(D)),
                   D = D)

if (d_muts$model[1] != "Add") {
  d_LD$mutType_AB <- d_rank$mutType_ab
}

# Summarise: overall
sum_LD <- d_LD %>%
  summarise(meanD = mean(D),
            sdD = sd(D),
            nD = length(D),
            nDP = length(D[D > 0.05]),
            nDN = length(D[D < -0.05]),
            nDHalf = length(D[abs(D) > 0.05]))

sum_LD$gen <- model_info[1]
sum_LD$seed <- model_info[2]
sum_LD$modelindex <- model_info[3]

sum_LD <- sum_LD %>% relocate(gen, seed, modelindex)

# Add counts of 10% groups for a histogram with 21 bins
labels <- paste0("n", 1:21)
bins <- seq(-0.25, 0.25, length.out = 21)

LDbins <- cut(d_LD$D, breaks = bins, right = F)
bin_labels <- levels(LDbins)

for (i in seq_along(bin_labels)) {
  sum_LD[,labels[i]] <- length(LDbins[LDbins == bin_labels[i]])
}


# Summarise: by frequency
sum_LD_f <- d_LD %>%
  filter(abs(freqDiff) <= 0.1) %>%
  mutate(freqBin = factor(freqBin)) %>%
  group_by(freqBin) %>%
  summarise(meanD = mean(D),
            sdD = sd(D),
            nD = length(D),
            nDP = length(D[D > 0.05]),
            nDN = length(D[D < -0.05]),
            nDHalf = length(D[abs(D) > 0.05]))

sum_LD_f$gen <- model_info[1]
sum_LD_f$seed <- model_info[2]
sum_LD_f$modelindex <- model_info[3]

sum_LD_f <- sum_LD_f %>% relocate(gen, seed, modelindex)

for (i in seq_along(bin_labels)) {
  sum_LD_f[,labels[i]] <- length(LDbins[LDbins == bin_labels[i]])
}

# Write output
write.table(sum_LD, FILE_LD, row.names = F, col.names = F)
write.table(sum_LD_f, FILE_LD_F, row.names = F, col.names = F)
write.table(d_LD, FILE_LD_TABLE, row.names = F, col.names = F)
