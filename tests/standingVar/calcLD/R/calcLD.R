# R Script to calculate LD from allele frequencies and effect sizes
library(tidyverse)
library(data.table)

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

FNS_PATH <- "~/tests/standingVar/calcMutationStats/R/"
#FNS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"

source(paste0(FNS_PATH, "helperFunctionsAndSetup.R"))

# Load functions for loading relatedness/haplotype matrices
source("~/tests/standingVar/calcLD/R/LDHelpFns.R")
source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/R/LDHelpFns.R")

# Load in correct line
d_freqs <- scan(paste0(WRITE_PATH, "slim_sharedmutfreqs.csv"), skip = run - 1, 
                    nlines = 1, sep = ",")

d_freqs <- scan(paste0("~/Desktop/slim_sharedmutfreqs4536158911153045504_1.csv"), skip = 0,
                nlines = 1, sep = ",")

model_info <- d_freqs[1:3]
d_freqs <- d_freqs[-(1:3)]

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

# Quantitative data
d_qg <- tbl(con, "slim_qg") %>%
  filter(modelindex == model_info[3], gen == model_info[1], seed == model_info[2]) %>%
  distinct() %>%
  ungroup()
d_qg <- d_qg %>% collect()

# Load mutation data: matching mutIDs, or fixations for the parental genotype
d_muts <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model_info[3], gen == model_info[1], seed == model_info[2],
         mutID %in% mutIDs || freq == 1)
d_muts <- d_muts %>% collect()

d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_muts$mutType <- as.factor(d_muts$mutType)
d_muts$fixGen <- as.numeric(d_muts$fixGen)
d_muts <- d_muts %>% 
  rename(value = effect)

# Read in combos
d_combos <- read.table("~/tests/standingVar/R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))
#testing
d_combos <- read.table("/mnt/c/GitHub/SLiMTests/tests/standingVar/R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

## Generate some effects for testing, otherwise load from d_muts ##
d_muts <- data.frame(gen = rep(model_info[1], times = length(mutIDs)),
                     seed = rep(model_info[2], times = length(mutIDs)),
                     modelindex = rep(model_info[3], times = length(mutIDs)),
                     mutID = mutIDs,
                     mutType = 3,
                     value = rnorm(length(mutIDs)),
                     fixGen = NA)
d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_muts$mutType <- as.factor(d_muts$mutType)
d_muts$fixGen <- as.numeric(d_muts$fixGen)


d_muts <- AddCombosToDF(d_muts)
d_muts$model <- "ODE"
# From the mutations, calculate the fitness of genotypes to determine which is AB/ab

# Calculate parental genotype, intermediates, and derived
## In example parental effect is 0
d_rank <- PairwiseFitnessRank(d_muts %>% filter(!is.na(fixGen)), 
                             d_muts %>% filter(is.na(fixGen)),
                             mutAID, mutBID)

relabeled_freqs <- RelabelGenotypeFrequencies(d_rank, mut_freqs)

# With the relabelled frequencies, calculate D
D <- CalcLD(relabeled_freqs)

d_LD <- data.frame(gen = rep(model_info[1], times = length(D)),
                   seed = rep(model_info[2], times = length(D)),
                   modelindex = rep(model_info[3], times = length(D)),
                   mutID_A = d_rank$mutIDA,
                   mutID_B = d_rank$mutIDB,
                   freq_A = mutpA,
                   freq_B = mutpB,
                   mutType_AB = d_rank$mutType_ab,
                   D = D,
                   fixGen = NA)

# Summarise
result <- d_LD %>%
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
