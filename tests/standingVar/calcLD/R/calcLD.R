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

FNS_PATH <- "~/tests/standingVar/calcMutationStats/R/"
#FNS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"

source(paste0(FNS_PATH, "helperFunctionsAndSetup.R"))

# Load functions for loading relatedness/haplotype matrices
source("~/tests/standingVar/calcLD/LDHelpFns.R")

# Load in correct line
d_freqs <- scan(paste0(WRITE_PATH, "slim_sharedmutfreqs.csv"), skip = run - 1, 
                    nlines = 1, sep = ",")

d_freqs <- scan(paste0("~/Desktop/slim_sharedmutfreqs4190440546994487296_1.csv"), skip = 0,
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
d_muts$model <- "Add"
# From the mutations, calculate the fitness of genotypes to determine which is AB/ab

# Calculate parental genotype, intermediates, and derived
## In example parental effect is 0
d_rank <- PairwiseFitnessRank(d_muts %>% filter(!is.na(fixGen)), 
                             d_muts %>% filter(is.na(fixGen)),
                             mutAID, mutBID)

# Rerank genotypes from parental-based to fitness-based
genotype_names <- substr(colnames(d_rank)[6:9], 5, 6)
ab_id <- genotype_names[apply(d_rank[, 6:9], 1, which.min)]
AB_id <- genotype_names[apply(d_rank[, 6:9], 1, which.max)]

# Intermediates can be randomly assigned, doesn't matter too much
Ab_id <- genotype_names[apply(d_rank[, 6:9], 1, function(x) {(1:4)[-c(which.max(x), which.min(x))][1]})]
aB_id <- genotype_names[apply(d_rank[, 6:9], 1, function(x) {(1:4)[-c(which.max(x), which.min(x))][2]})]

# Now with the correct reranking, we need to reassign our frequencies
pAB <- sym(paste0("mutp", AB_id))
