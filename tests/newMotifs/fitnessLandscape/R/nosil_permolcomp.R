# Runs latin hypercube samples using the Nosil method per model on each molecular component
# to see how each individually contributes to landscape ruggedness
# each molecular component is varied n times across m genetic backgrounds to measure interactions
library(tidyverse)
library(DoE.wrapper)
library(deSolve)
library(mvtnorm)
library(broom)
library(latex2exp)
library(paletteer)
library(ggbeeswarm)

library(future)
library(doParallel)
library(foreach)
source("./fitnesslandscapefunctions.R")

CalculateRuggednessParallel <- function(g, model, optima, sigma, n = 10, nCores = availableCores(),
                                        width = 0.004,
                                        seed = sample(1:.Machine$integer.max, nrow(g))) {
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
    
    setwd("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/")
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
    
    result <- data.frame(model = model,
                         startW = rollingFitnesses[1],
                         endW = rollingFitnesses[n+1],
                         netChangeW = changeFitnesses[length(changeFitnesses)] - changeFitnesses[1],
                         sumChangeW = sum(abs(diff(changeFitnesses))),
                         numFitnessHoles = sum(rollingFitnesses <= 0.0))
    return(result)
  }
  
  stopCluster(cl)
  return(df_result)
}

# Nosil method
# Generate Latin hypercube starting points
models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")

NUM_RUNS <- 10000
NUM_BACKGROUNDS <- 10
MAX_COMP_SIZE <- log(3)
nComps <- length(comps)

# Hypercube is NUM_RUNs per molecular component, per genetic background
pars <- lhs.design(NUM_RUNS, 1, type = "random")
pars <- pars * MAX_COMP_SIZE
pars <- pars[,1]

# Generate backgrounds
parBackgrounds <- matrix(runif((length(comps)) * NUM_BACKGROUNDS), ncol = (length(comps)))
parBackgrounds <- parBackgrounds * MAX_COMP_SIZE

# Replicate rows for the number of molecular components
parBackgrounds <-  parBackgrounds %x% rep(1, nComps)

# Replicate rows for the number of replicate values
parBackgrounds <-  rep(1, NUM_RUNS) %x% parBackgrounds

# replace diagonal components for each component and background
repeated_pars <- rep(pars, each = nComps)
for (i in seq_len(NUM_RUNS)) {
  # Fill NUM_BACKGROUNDS diagonals
  for (j in seq_len(NUM_BACKGROUNDS)) {
    k <- (i - 1) * nComps + j
    offset_start <- ((i - 1) * nComps * NUM_BACKGROUNDS) + ((j - 1) * nComps) + 1 # number of previous runs
    offset_end <-  offset_start + (nComps - 1)
    diag(parBackgrounds[offset_start:offset_end,]) <- pars[i]
  }
}

colnames(parBackgrounds) <- comps

pars <- as.data.frame(parBackgrounds)

# Run in parallel
#seeds <- sample(1:.Machine$integer.max, NUM_RUNS)
#saveRDS(seeds, "seeds.RDS")
seeds <- readRDS("seeds.RDS")

d_ruggedness <- data.frame(
  model = character(NUM_RUNS * length(models)),
  startW = numeric(NUM_RUNS * length(models)),
  endW = numeric(NUM_RUNS * length(models)),
  netChangeW = numeric(NUM_RUNS * length(models)),
  sumChangeW = numeric(NUM_RUNS * length(models)),
  numFitnessHoles = integer(NUM_RUNS * length(models))
)

for (model in models) {
  # randomly sample an optimum
  parsMasked <- ParsMask(pars, model)
  optMolComps <- as.data.frame(t(runif(ncol(parsMasked), 0, MAX_COMP_SIZE)))
  colnames(optMolComps) <- colnames(parsMasked)
  startSolution <- SolveModel(exp(optMolComps), model)
  startTraits <- GetTraitValues(startSolution, model, exp(optMolComps))
  sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
  opt <- CalcOptima(startTraits, sigma, 0.9)
  
  RugRes <- CalculateRuggednessParallel(parsMasked, model, opt, sigma, seed = seeds)
  
  output_index <- match(model, models)
  d_ruggedness[(NUM_RUNS * (output_index - 1) + 1):(NUM_RUNS * output_index),] <- RugRes
}