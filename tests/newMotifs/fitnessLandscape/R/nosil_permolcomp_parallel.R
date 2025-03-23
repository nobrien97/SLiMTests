library(tidyverse)
#library(DoE.wrapper)
library(deSolve)
library(mvtnorm)
library(broom)

DATA_PATH <- "/home/564/nb9894/tests/newMotifs/fitnessLandscape/ruggedness/R/" 
#DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/" 
SAVE_PATH <- "/scratch/ht96/nb9894/newMotifs/fitnessLandscape/ruggedness/"

setwd(DATA_PATH)

# Load functions
source("./fitnesslandscapefunctions.R")

CalculateRuggednessSingle <- function(g, model, optima, sigma, n = 10,
                                        width = 0.004,
                                        path) {
  # g = genotypes (molecular components). Replicate starting points for the walk
  # w = fitnesses of the starting points
  # n = number of steps in the walk
  # seed = replicate seed for the run
  
  df_result <- data.frame(model = character(nrow(g)),
                          startW = numeric(nrow(g)),
                          endW = numeric(nrow(g)),
                          netChangeW = numeric(nrow(g)),
                          sumChangeW = numeric(nrow(g)),
                          numFitnessHoles = integer(nrow(g)))
  
  for(row_index in seq_len(nrow(g))) {
    nComps <- ncol(g)
    rollingGenotypes <- g[1:(n+1), ]
    rollingFitnesses <- numeric(n+1)
    
    # Set the seed for each walk
    #set.seed(seed)
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
    
    result <- data.frame(model = model,
                         startW = rollingFitnesses[1],
                         endW = rollingFitnesses[n+1],
                         netChangeW = changeFitnesses[length(changeFitnesses)] - changeFitnesses[1],
                         sumChangeW = sum(abs(diff(changeFitnesses))),
                         numFitnessHoles = sum(rollingFitnesses <= 0.0))
    df_result[row_index,] <- result
  }
  
  return(df_result)
}


# input arguments: seed and par
args <- commandArgs(trailingOnly = T)
seed_idx <- as.numeric(args[1])
par_idx <- as.numeric(args[2])

# Nosil method
# Generate Latin hypercube starting points
models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")

NUM_RUNS <- 10000
NUM_BACKGROUNDS <- 10
MAX_COMP_SIZE <- log(3)
nComps <- length(comps)

pars <- readRDS("pars.RDS")
pars <- pars[par_idx,]


seeds <- readRDS(paste0(DATA_PATH, "seeds.RDS"))
seed <- seeds[seed_idx]

d_ruggedness <- data.frame(
  model = character(length(models)),
  startW = numeric(length(models)),
  endW = numeric(length(models)),
  netChangeW = numeric(length(models)),
  sumChangeW = numeric(length(models)),
  numFitnessHoles = integer(length(models)),
  molComp <- character(length(models)),
  bkg <- integer(length(models))
)


# Iterate over models
for (model in models) {
  # randomly sample an optimum
  set.seed(seed)
  parsMasked <- ParsMask(pars, model)
  optMolComps <- as.data.frame(t(runif(ncol(parsMasked), 0, MAX_COMP_SIZE)))
  colnames(optMolComps) <- colnames(parsMasked)
  startSolution <- SolveModel(exp(optMolComps), model)
  startTraits <- GetTraitValues(startSolution, model, exp(optMolComps))
  sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
  opt <- CalcOptima(startTraits, sigma, 0.9)
  
  RugRes <- CalculateRuggednessSingle(parsMasked, model, opt, sigma, path = DATA_PATH)
  
  # Set identifiers
  RugRes$molComp <- comps[(par_idx - 1) %% nComps + 1]
  RugRes$bkg <- c(rep(rep(1:NUM_BACKGROUNDS, each = nComps), times = NUM_RUNS))[par_idx]
  
  output_index <- match(model, models)
  d_ruggedness[output_index,] <- RugRes
}

write_csv(d_ruggedness, paste0(SAVE_PATH, "d_ruggedness_", seed_idx, "_", par_idx, ".csv"), col_names = F)
