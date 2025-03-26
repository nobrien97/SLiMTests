library(tidyverse)
library(deSolve)
library(mvtnorm)
library(broom)

DATA_PATH <- "/home/564/nb9894/tests/newMotifs/fitnessLandscape/ruggedness/R/" 
#DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/" 
SAVE_PATH <- "/scratch/ht96/nb9894/newMotifs/fitnessLandscape/ruggedness/"

setwd(DATA_PATH)

# Load functions
source("./fitnesslandscapefunctions.R")

CalculateRuggednessParallel <- function(g, model, optima, sigma, n = 10,
                                        width = 0.004,
                                        seed,
                                        path) {
  # g = genotypes (molecular components). Replicate starting points for the walk
  # w = fitnesses of the starting points
  # n = number of steps in the walk
  # seed = replicate seed for the run
  
  result <- data.frame(model = character(nrow(g)),
                         startW = numeric(nrow(g)),
                         endW = numeric(nrow(g)),
                         netChangeW = numeric(nrow(g)),
                         sumChangeW = numeric(nrow(g)),
                         numFitnessHoles = integer(nrow(g)))

  # Minimum fitness to be considered == 0
  fitness_epsilon <- 1e-10
  
  for(row_index in seq_len(nrow(g))) {
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
    changeFitnesses <- rollingFitnesses[rollingFitnesses > fitness_epsilon]
    netChangeW <- changeFitnesses[length(changeFitnesses)] - changeFitnesses[1]
    
    # 
    if (length(netChangeW) == 0) {
      netChangeW <- 0
    }
    
    result[row_index, ] <- c(model = model,
                         startW = rollingFitnesses[1],
                         endW = rollingFitnesses[n+1],
                         netChangeW = newChangeW,
                         sumChangeW = sum(abs(diff(changeFitnesses))),
                         numFitnessHoles = sum(rollingFitnesses < fitness_epsilon))
  }
  
  return(result)
}


# input arguments: index to load parameters
args <- commandArgs(trailingOnly = T)
par_idx <- as.numeric(args[1])


# Nosil method
# Generate Latin hypercube starting points
models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")

NUM_RUNS <- 10000
NUM_BACKGROUNDS <- 10
REPS_PER_RUN <- 10
MAX_COMP_SIZE <- log(3)
nComps <- length(comps)

# 10 replicates per run, each run will return a dataframe with 1200 * 5 = 6000 (5 per model) rows in it
# for 1000 total files to combine
ROWS_PER_RUN <- nComps * NUM_BACKGROUNDS * REPS_PER_RUN

# range of input rows to evaluate this run
par_idx_range <- (ROWS_PER_RUN * (par_idx - 1) + 1):(ROWS_PER_RUN * par_idx)

# Read in parameters:
# Data frame in blocks of 120 (nComps * NUM_BACKGROUNDS)
# each block is one replicate mutation applied in 10 backgrounds in 12 different molecular components
# 10000 total replicates for 1200000 applications of that replicate in the backgrounds and mol comps
# 
pars <- readRDS("pars.RDS")
pars <- pars[par_idx_range,]

seeds <- readRDS(paste0(DATA_PATH, "seeds.RDS"))
seed <- seeds[par_idx_range]

d_ruggedness <- data.frame(
  model = character(ROWS_PER_RUN * length(models)),
  startW = numeric(ROWS_PER_RUN * length(models)),
  endW = numeric(ROWS_PER_RUN * length(models)),
  netChangeW = numeric(ROWS_PER_RUN * length(models)),
  sumChangeW = numeric(ROWS_PER_RUN * length(models)),
  numFitnessHoles = integer(ROWS_PER_RUN * length(models)),
  molComp <- character(ROWS_PER_RUN * length(models)),
  bkg <- integer(ROWS_PER_RUN * length(models))
)

# opt_seed <- sample(1:.Machine$integer.max, 1)
# 423812551
opt_seed <- 423812551L

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
  opt <- CalcOptima(startTraits, sigma, 0.9)
  
  RugRes <- CalculateRuggednessParallel(parsMasked, model, opt, sigma,
                                        seed = seed,
                                        path = DATA_PATH)
  
  # Set identifiers
  RugRes$molComp <- comps[(par_idx_range - 1) %% nComps + 1]
  RugRes$bkg <- c(rep(rep(1:NUM_BACKGROUNDS, each = nComps), times = REPS_PER_RUN))
  
  output_index <- match(model, models)
  d_ruggedness[(ROWS_PER_RUN * (output_index - 1) + 1):(ROWS_PER_RUN * output_index),] <- RugRes
}

write_csv(d_ruggedness, paste0(SAVE_PATH, "d_ruggedness_", par_idx, ".csv"), col_names = F)
