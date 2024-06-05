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
R_PATH <- "~/tests/standingVar/sanityChecks/calcMutationStats/R/"
source("~/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")
GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/sanityChecks/"

WRITE_PATH <- "/scratch/ht96/nb9894/standingVar/sanityChecks/calcMutationStats/"
EPISTASIS_FILE <- paste0(WRITE_PATH, "d_epistasis_", model, ".csv")
EPISTASIS_WEIGHTED_FILE <- paste0(WRITE_PATH, "d_freqweight_epistasis_", model, ".csv")
EFFECTS_FILE <- paste0(WRITE_PATH, "d_fx_", model, ".csv")
DPDT_FILE <- paste0(WRITE_PATH, "d_dpdt_", model, ".csv")
SFS_FILE <- paste0(WRITE_PATH, "d_SFS_", model, ".csv")

AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           nloci = d_combos$nloci[as.numeric(levels(modelindex))[modelindex]],
           tau = d_combos$tau[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]],
           width = d_combos$width[as.numeric(levels(modelindex))[modelindex]])
}

# Load combo information
d_combos <- read.table(paste0(R_PATH, "combos.csv"), header = F,
                       col.names = c("nloci", "tau", "r", "width", "model"))

# Load mutation data from database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "standingVarMuts.db"))

d_qg <- data.table::fread(paste0(GDATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

# Quantitative data
d_qg %>%
  distinct() %>%
  filter(modelindex == model) %>%
  group_by(seed, modelindex) %>%
  filter(gen >= 49500) %>% distinct() %>%
  ungroup() -> d_adapted


# If the model has 0 adapted runs, end early
if (nrow(d_adapted) == 0) {
  print(paste("Model ", run, " never adapted across all seeds, closing R."))
  q(save = "no")
  
}

# Mutation data
adapted_seeds <- unique(d_adapted$seed)
d_muts_adapted <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model, gen >= 49500, seed %in% adapted_seeds) %>%
  select(!c(pos, chi))
d_muts_adapted <- d_muts_adapted %>% collect()


# Set data types for incorrect columns
d_muts_adapted$seed <- as.factor(d_muts_adapted$seed)
d_muts_adapted$modelindex <- as.factor(d_muts_adapted$modelindex)
d_muts_adapted$mutType <- as.factor(d_muts_adapted$mutType)
d_muts_adapted$fixGen <- as.numeric(d_muts_adapted$fixGen)
d_muts_adapted <- d_muts_adapted %>% 
  rename(value = effect)

d_adapted$seed <- as.factor(d_adapted$seed)
d_adapted$modelindex <- as.factor(d_adapted$modelindex)

# dP/dt
sampleRate <- 50 # sample every 50 generations, so divide deltaP by 50
d_adapted %>%
  mutate(dPdT = deltaPheno / sampleRate) -> d_adapted

d_dpdt <- d_adapted %>%
  filter(gen >= 50000) %>%
  mutate(optPerc = (phenomean - 1))    # percent to optimum

# Determine when we first reach 25%, 50%, 75%, 90% of the optimum 
# (90% being our cutoff for adaptation)
d_dpdt$optPerc <- cut(d_dpdt$optPerc, c(-Inf, 0.25, 0.5, 0.75, 0.9, Inf),
                      right = F)

# Mean change within each of these groups (from 25% to 50%, from 50% to 75% etc.)
d_dpdt %>%
  group_by(optPerc, modelindex) %>%
  summarise(meandPdT = mean(dPdT),
            sddPdT = sd(dPdT)) -> d_dpdt_sum

# write
data.table::fwrite(d_dpdt_sum, 
                   DPDT_FILE, sep = ",", 
                   col.names = F, row.names = F)

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

d_fixed_adapted <- d_com_adapted %>% filter(!is.na(fixGen))

# SFS
d_SFS <- CalcSFS(d_com_adapted)

d_SFS %>%
  group_by(optPerc, modelindex, mutType, freqBin) %>%
  summarise(countFreqBin = n(),
            meanValue = mean(value),
            sdValue = sd(value)) -> d_SFS_sum


data.table::fwrite(d_SFS_sum, 
                   SFS_FILE, sep = ",", 
                   col.names = F, row.names = F)


# Calculate phenotype effects
d_phenofx <- CalcPhenotypeEffects(d_com_adapted %>% filter(is.na(fixGen)),
                                  d_fixed_adapted)

d_phenofx <- d_phenofx %>%
  select(gen, seed, modelindex, mutType, mutID, s)


d_epistasis <- PairwiseEpistasis(d_fixed_adapted,
                                 d_com_adapted %>% 
                                   filter(is.na(fixGen)) %>%
                                   select(gen, seed, modelindex, mutType, freq, value),
                                 m = 48, n = 1000, F, F)

d_epistasis_freqweight <- PairwiseEpistasis(d_fixed_adapted,
                                                  d_com_adapted %>% 
                                                    filter(is.na(fixGen)) %>%
                                                    select(gen, seed, modelindex, mutType, freq, value),
                                                  m = 48, n = 1000, F, T)

# write to file
data.table::fwrite(d_epistasis, 
  EPISTASIS_FILE, sep = ",", col.names = F, row.names = F)

data.table::fwrite(d_epistasis_freqweight,
  EPISTASIS_WEIGHTED_FILE, sep = ",", col.names = F, row.names = F)

data.table::fwrite(d_phenofx,
  EFFECTS_FILE, sep = ",", col.names = F, row.names = F)
