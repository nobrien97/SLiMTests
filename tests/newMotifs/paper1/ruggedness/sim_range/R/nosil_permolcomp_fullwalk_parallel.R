library(tidyverse)
library(deSolve)
library(mvtnorm)
library(broom)

DATA_PATH <- "/home/564/nb9894/tests/newMotifs/paper1/ruggedness/sim_range/R/" 
# DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/" 
SAVE_PATH <- "/scratch/ht96/nb9894/newMotifs/paper1/ruggedness/sim_range/"

setwd(DATA_PATH)

# Load functions
source("./fitnesslandscapefunctions.R")

# input arguments: index to load parameters
args <- commandArgs(trailingOnly = T)
par_idx <- as.numeric(args[1])


# Nosil method
# Generate Latin hypercube starting points

NUM_BACKGROUNDS <- 10
NUM_STEPS <- 2
REPS_PER_RUN <- 10

# Load parameter combinations and seeds
pars_all <- readRDS(paste0(DATA_PATH, "pars.RDS"))
seeds_all <- readRDS(paste0(DATA_PATH, "seeds.RDS"))

# Read in parallel/orthogonal/randomised directions
parallel_opt_dir <- read_csv(paste0(DATA_PATH, "parallel_traitdir.csv"), col_names = F)
orth_opt_dir <- read_csv(paste0(DATA_PATH, "orth_traitdir.csv"), col_names = F)

# output
d_ruggedness <- list()

# Load in range of values to use for each motif
d_molcomp_maxvals <- readRDS("~/tests/newMotifs/paper1/ruggedness/sim_range/R/d_molcomp_maxvals.RDS")

# Iterate over models
for (current_model in models) {
  active_comps <- molComp_names[[current_model]] 
  
  molcomp_min_values <- d_molcomp_maxvals %>% ungroup() %>%
    filter(model == current_model) %>% select(min_value) %>% unlist()
  names(molcomp_min_values) <- active_comps 
  
  molcomp_max_values <- d_molcomp_maxvals %>% ungroup() %>%
    filter(model == current_model) %>% select(max_value) %>% unlist()
  names(molcomp_max_values) <- active_comps 
  
  nComps <- length(active_comps)

  # 10 backgrounds evaluated per run
  # 10 replicates per run, each run will return a dataframe with nComps * 10 * 10 rows in it
  # for 1000 total files to combine
  ROWS_PER_RUN <- nComps * NUM_BACKGROUNDS * REPS_PER_RUN 

  # range of input rows to evaluate this run
  par_idx_range <- (ROWS_PER_RUN * (par_idx - 1) + 1):(ROWS_PER_RUN * par_idx)

# Read in parameters:
# Data frame in blocks of nComps * NUM_BACKGROUNDS
# each block is one replicate mutation applied in 10 backgrounds in each different molecular components
# 10000 total replicates applied in each backgrounds and mol comp

  # Filter parameters to the current model
  pars <- pars_all %>% filter(model == current_model) %>% select(-model) %>%
    select(all_of(active_comps))
  
  pars <- pars[par_idx_range,]
  
  # Read in seeds and choose the correct model
  seeds <- seeds_all[[current_model]]
  seed <- as.integer(seeds[par_idx_range])
  
  # Sample seed for optimum (different between replicates, same between backgrounds)
  opt_seed <- sample(1:.Machine$integer.max, 1)

  # randomly sample an optimum
  set.seed(opt_seed)

  optMolComps <- as.data.frame(t(runif(nComps, molcomp_min_values, molcomp_max_values)))
  colnames(optMolComps) <- colnames(pars)
  startSolution <- SolveModel(optMolComps, current_model)
  startTraits <- GetTraitValues(startSolution, current_model, optMolComps)
  sigma <- CalcSelectionSigmas(startTraits, 0.05, 0.1, 0.1)
  par_dir_model <- unlist(parallel_opt_dir[which(models == current_model),c(1:length(startTraits), ncol(parallel_opt_dir))])
  orth_dir_model <- unlist(orth_opt_dir[which(models == current_model),c(1:length(startTraits), ncol(orth_opt_dir))])

  opt_rand <- CalcOptima(startTraits, sigma, 0.95)
  opt_par <- CalcNewOptimumAlongVector(startTraits, sigma, 0.95, par_dir_model)
  opt_orth <- CalcNewOptimumAlongVector(startTraits, sigma, 0.95, orth_dir_model)

  
  RugRes_rand <- CalculateRuggednessLandscaper(pars, current_model, "Randomised", opt_rand, sigma,
                                        n = NUM_STEPS,
                                        nCores = 1,
                                        seed = seed)

  RugRes_par <- CalculateRuggednessLandscaper(pars, current_model, "Parallel", opt_par, sigma,
                                        n = NUM_STEPS,
                                        nCores = 1,
                                        seed = seed)
                                      
  RugRes_orth <- CalculateRuggednessLandscaper(pars, current_model, "Orthogonal", opt_orth, sigma,
                                        n = NUM_STEPS,
                                        nCores = 1,
                                        seed = seed)

  RugRes <- rbind(RugRes_rand, RugRes_par, RugRes_orth)
  
  # Set identifiers
  RugRes$molComp <- rep(active_comps[(par_idx_range - 1) %% nComps + 1], each = NUM_STEPS+1)
  RugRes$bkg <- rep(c(rep(rep(1:NUM_BACKGROUNDS, each = nComps), times = REPS_PER_RUN)), each = NUM_STEPS+1)
  
  output_index <- match(current_model, models)
  d_ruggedness[[output_index]] <- distinct(RugRes)
}

d_ruggedness <- data.table::rbindlist(d_ruggedness, fill = T)

# Rearrange columns
d_ruggedness <- d_ruggedness %>%
  select(step, model, dataset, fitness, startW, endW, netChangeW, sumChangeW, numFitnessHoles, 
  nSteps, aX, KZX, aY, bY, KY, KZ, KXZ,
  aZ, bZ, Hilln, XMult, base,
  molComp, bkg)

write_csv(d_ruggedness, paste0(SAVE_PATH, "d_ruggedness_", par_idx, ".csv"), col_names = F)
