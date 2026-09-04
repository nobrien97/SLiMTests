library(tidyverse)
library(deSolve)
library(mvtnorm)
library(broom)

DATA_PATH <- "/home/564/nb9894/tests/newMotifs/fitnessLandscape/ruggedness/R/" 
# DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/" 
SAVE_PATH <- "/scratch/ht96/nb9894/newMotifs/fitnessLandscape/ruggedness/"

setwd(DATA_PATH)

# Load functions
source("./fitnesslandscapefunctions.R")

CalculateRuggednessParallel <- function(g, model, dataset, optima, sigma, n = 10,
                                        width = 0.004,
                                        nCores,
                                        seed,
                                        path) {
  # g = genotypes (molecular components). Replicate starting points for the walk
  # w = fitnesses of the starting points
  # n = number of steps in the walk
  # seed = replicate seed for the run
  
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  
  df_result <- foreach (row_index = seq_len(nrow(g)), .combine = rbind) %dopar% {
    require(tidyverse)
    require(deSolve)
    require(mvtnorm)
    
    setwd(path)
    source("./fitnesslandscapefunctions.R")
    
    nComps <- ncol(g)
    rollingGenotypes <- g[1:(n+1), ]
    rollingFitnesses <- numeric(n+1)
    
    # Set the seed for each walk
    set.seed(seed[row_index])
    # Sample n steps per genotype per a normal distribution with a given width
    # Assume width is split evenly across the components
    mutations <- rmvnorm(n, sigma = diag(nComps) * ( width / nComps ))
    mutations <- rbind(rep(0.0, nComps), mutations)
    
    # cumulative sum each column to add it to rollingGenotypes
    mutations <- apply(mutations, 2, cumsum)
    rollingGenotypes <- exp(log(g[rep(row_index, times = n+1),]) + mutations)
    for (j in seq_len(n+1)) {
      rollingFitnesses[j] <- CalcTraitAndFitness(rollingGenotypes[j,], 
                                                 model,
                                                 optima, 
                                                 sigma)
    }
    # Calculate results - add in original fitness
    # remove invalid fitnesses from bad solutions
    changeFitnesses <- rollingFitnesses[rollingFitnesses >= 0.0]
    
    result <- data.frame(step = 1:(n+1),
                         model = rep(model, times = n+1),
                         dataset = rep(dataset, times = n+1),
                         fitness = rollingFitnesses,
                         startW = rep(rollingFitnesses[1], times = n+1),
                         endW = rep(rollingFitnesses[n+1], times = n+1),
                         netChangeW = rep(changeFitnesses[length(changeFitnesses)] - changeFitnesses[1], times = n+1),
                         sumChangeW = rep(sum(abs(diff(changeFitnesses))), times = n+1),
                         numFitnessHoles = rep(sum(rollingFitnesses <= 0.0)), times = n+1)
 
    result[,comps] <- rollingGenotypes

    return(result)
  }
  
  stopCluster(cl)
  return(df_result)
}


# input arguments: index to load parameters
args <- commandArgs(trailingOnly = T)
par_idx <- as.numeric(args[1])


# Nosil method
# Generate Latin hypercube starting points
models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")

NUM_BACKGROUNDS <- 10
NUM_STEPS <- 10
REPS_PER_RUN <- 10
MAX_COMP_SIZE <- log(3)
nComps <- length(comps)

# 10 backgrounds evaluated per run
# 10 replicates per run, each run will return a dataframe with 12 * 10 * 10 = 12000 rows in it
# for 1000 total files to combine - 120,000,000 rows
ROWS_PER_RUN <- nComps * NUM_BACKGROUNDS * REPS_PER_RUN 

# range of input rows to evaluate this run
par_idx_range <- (ROWS_PER_RUN * (par_idx - 1) + 1):(ROWS_PER_RUN * par_idx)

# Read in parameters:
# Data frame in blocks of 120 (nComps * NUM_BACKGROUNDS)
# each block is one replicate mutation applied in 10 backgrounds in 12 different molecular components
# 10000 total replicates for 1200000 applications of that replicate in the backgrounds and mol comps
# 
pars <- readRDS(paste0(DATA_PATH, "pars.RDS"))
pars <- pars[par_idx_range,]

seeds <- readRDS(paste0(DATA_PATH, "seeds.RDS"))
seed <- seeds[par_idx_range]

# Read in parallel/orthogonal/randomised directions
parallel_opt_dir <- read_csv(paste0(DATA_PATH, "parallel_traitdir.csv"), col_names = F)
orth_opt_dir <- read_csv(paste0(DATA_PATH, "orth_traitdir.csv"), col_names = F)

OUTPUT_LENGTH <- ROWS_PER_RUN * NUM_STEPS * length(models) * 3

d_ruggedness <- list()

# Sample seed for optimum (different between replicates, same between backgrounds)
opt_seed <- sample(1:.Machine$integer.max, 1)

# Iterate over models
for (model in models) {
  # randomly sample an optimum
  set.seed(opt_seed)
  parsMasked <- ParsMask(pars, model)
  optMolComps <- as.data.frame(t(runif(ncol(parsMasked), 0, MAX_COMP_SIZE)))
  colnames(optMolComps) <- colnames(parsMasked)
  startSolution <- SolveModel(exp(optMolComps), model)
  startTraits <- GetTraitValues(startSolution, model, exp(optMolComps))  
  sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
  par_dir_model <- unlist(parallel_opt_dir[which(models == model),c(1:length(startTraits), ncol(parallel_opt_dir))])
  orth_dir_model <- unlist(orth_opt_dir[which(models == model),c(1:length(startTraits), ncol(orth_opt_dir))])

  opt_rand <- CalcOptima(startTraits, sigma, 0.95)
  opt_par <- CalcNewOptimumAlongVector(startTraits, sigma, 0.95, par_dir_model)
  opt_orth <- CalcNewOptimumAlongVector(startTraits, sigma, 0.95, orth_dir_model)

  RugRes_rand <- CalculateRuggednessParallel(parsMasked, model, "Randomised", opt_rand, sigma,
                                        nCores = future::availableCores(),
                                        seed = seed,
                                        path = DATA_PATH)

  RugRes_par <- CalculateRuggednessParallel(parsMasked, model, "Parallel", opt_par, sigma,
                                        nCores = future::availableCores(),
                                        seed = seed,
                                        path = DATA_PATH)
                                      
  RugRes_orth <- CalculateRuggednessParallel(parsMasked, model, "Orthogonal", opt_orth, sigma,
                                        nCores = future::availableCores(),
                                        seed = seed,
                                        path = DATA_PATH)

  RugRes <- rbind(RugRes_rand, RugRes_par, RugRes_orth)
  
  # Set identifiers
  RugRes$molComp <- rep(comps[(par_idx_range - 1) %% nComps + 1], each = NUM_STEPS+1)
  RugRes$bkg <- rep(c(rep(rep(1:NUM_BACKGROUNDS, each = nComps), times = REPS_PER_RUN)), each = NUM_STEPS+1)
  
  output_index <- match(model, models)
  d_ruggedness[[output_index]] <- RugRes
}

d_ruggedness <- data.table::rbindlist(d_ruggedness)

write_csv(d_ruggedness, paste0(SAVE_PATH, "d_ruggedness_", par_idx, ".csv"), col_names = F)
