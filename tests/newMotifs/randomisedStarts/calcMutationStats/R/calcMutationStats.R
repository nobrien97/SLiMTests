# Calculates statistics involving mutation data
# Runs per model across all seeds of that model
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(data.table)

# Get command line arguments
## 1: model (modelindex)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

# Paths
R_PATH <- "~/tests/newMotifs/randomisedStarts/calcMutationStats/R/"
COMBO_PATH <- "~/tests/newMotifs/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))
GDATA_PATH <- "/g/data/ht96/nb9894/newMotifs/randomisedStarts/"

WRITE_PATH <- "/scratch/ht96/nb9894/newMotifs/randomisedStarts/calcMutationStats/"
EPISTASIS_FILE <- paste0(WRITE_PATH, "d_epistasis_", model, ".csv")
EPISTASIS_WEIGHTED_FILE <- paste0(WRITE_PATH, "d_freqweight_epistasis_", model, ".csv")
EFFECTS_FILE <- paste0(WRITE_PATH, "d_fx_", model, ".csv")
DPDT_FILE <- paste0(WRITE_PATH, "d_dpdt_", model, ".csv")
SFS_FILE <- paste0(WRITE_PATH, "d_SFS_", model, ".csv")

# Load combo information
d_combos <- read.table(paste0(COMBO_PATH, "combos.csv"), header = F,
                       col.names = c("model", "r"))

# Load mutation data from database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "newMotifsMuts.db"))

# Quantitative data
d_qg <- tbl(con, "slim_qg") %>%
  filter(gen >= 49500, modelindex == model) %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.95)) %>%
  ungroup()
d_qg <- d_qg %>% collect()

# Measure at start and end of adaptive walk
d_qg <- d_qg %>%
filter(gen == 50000 | gen == 60000) %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))

# Mutation data
d_muts <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model, gen >= 49500) %>%
  select(!c(pos))
d_muts <- d_muts %>% collect()


# Set data types for incorrect columns
d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_muts$mutType <- as.factor(d_muts$mutType)
d_muts$fixGen <- as.numeric(d_muts$fixGen)
d_muts <- d_muts %>% 
  rename(value = effect)

d_qg$seed <- as.factor(d_qg$seed)
d_qg$modelindex <- as.factor(d_qg$modelindex)

# dP/dt
sampleRate <- 50 # sample every 50 generations, so divide deltaP by 50
d_qg %>%
  mutate(dPdT = deltaPheno / sampleRate) -> d_dpdt

# Mean change within each of these groups
d_dpdt %>%
  group_by(timePoint, modelindex) %>%
  summarise(meandPdT = mean(dPdT),
            sddPdT = sd(dPdT)) -> d_dpdt_sum

# write
data.table::fwrite(d_dpdt_sum, 
                   DPDT_FILE, sep = ",", 
                   col.names = F, row.names = F)

# filter by optPerc to select timepoints where populations 
# first reached 50% adapted etc.

# Filter mutations by optPerc generations
# inner join the mutation w/ quantitative data + add model info
d_com <- inner_join(d_dpdt, d_muts, 
                            by = c("gen", "seed", "modelindex"))
d_com <- AddCombosToDF(d_com)

d_fixed <- d_com %>% filter(!is.na(fixGen))

# SFS
d_SFS <- CalcSFS(d_com)

d_SFS %>%
  group_by(timePoint, modelindex, mutType, freqBin) %>%
  summarise(countFreqBin = n(),
            meanValue = mean(value),
            sdValue = sd(value)) -> d_SFS_sum


data.table::fwrite(d_SFS_sum, 
                   SFS_FILE, sep = ",", 
                   col.names = F, row.names = F)


# Calculate phenotype effects

## Get optima and selection strength
d_opt <- data.table::fread(paste0(GDATA_PATH, "slim_opt.csv"), header = F, 
                          sep = ",", colClasses = c("factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("seed", "modelindex", "trait1_opt", "trait2_opt",
                                        "trait3_opt", "trait4_opt", "trait1_sig", "trait2_sig",
                                        "trait3_sig", "trait4_sig", "trait1_dir", "trait2_dir",
                                        "trait3_dir", "trait4_dir"), 
                          fill = T)

# Get the optimum traits and width of the fitness function for this simulation
# We are running this on multiple seeds which have their own optima
d_opt <- d_opt %>% filter(modelindex == model)

d_phenofx <- CalcPhenotypeEffects(d_com %>% filter(is.na(fixGen)),
                                  d_fixed, d_opt)

d_phenofx <- d_phenofx %>%
  select(gen, seed, modelindex, mutType, mutID, s)

# m needs to be a multiple of the number of molecular components - set to twice the number of components
mValue <- GetNMutTypes(d_com$model[1]) * 2

d_epistasis <- PairwiseEpistasis(d_fixed,
                                 d_com %>% 
                                   filter(is.na(fixGen)) %>%
                                   select(gen, seed, modelindex, mutType, freq, value),
                                 m = mValue, n = 100, F, F)

d_epistasis_freqweight <- PairwiseEpistasis(d_fixed,
                                                  d_com %>% 
                                                    filter(is.na(fixGen)) %>%
                                                    select(gen, seed, modelindex, mutType, freq, value),
                                                  m = mValue, n = 100, F, T)

# write to file
data.table::fwrite(d_epistasis, 
  EPISTASIS_FILE, sep = ",", col.names = F, row.names = F)

data.table::fwrite(d_epistasis_freqweight,
  EPISTASIS_WEIGHTED_FILE, sep = ",", col.names = F, row.names = F)

data.table::fwrite(d_phenofx,
  EFFECTS_FILE, sep = ",", col.names = F, row.names = F)
