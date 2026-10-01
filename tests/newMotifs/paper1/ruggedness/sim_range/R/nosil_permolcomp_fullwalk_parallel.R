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


molComp_names <- list("NAR" = c(
  # NAR and PAR
  "aZ",
  "bZ",
  "KZ",
  "KXZ",
  "base", # baseline expression
  "Hilln", # hill coefficient
  "XMult" # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "base", # baseline expression
  "Hilln", # hill coefficient
  "XMult" # X multiplier
),
"FFBH" = c(
  # FFBH
  "aX",
  "KZX",
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "base", # baseline expression
  "Hilln", # hill coefficient
  "XMult" # X multiplier
)
)

molComp_names[["PAR"]] <- molComp_names[["NAR"]]
molComp_names[["FFLI1"]] <- molComp_names[["FFLC1"]]




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
    
    comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")

    nComps <- ncol(g)
    rollingGenotypes <- g[1:(n+1), ]
    rollingFitnesses <- numeric(n+1)
    
    # Set the seed for each walk
    set.seed(seed[row_index])
    # Sample n steps per genotype per a normal distribution with a given width
    # Assume width is split evenly across the components
    mutations <- rmvnorm(n, sigma = diag(nComps_total) * ( width / nComps ))
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
NUM_STEPS <- 2
REPS_PER_RUN <- 10

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
# 
pars <- readRDS(paste0(DATA_PATH, "pars.RDS"))

# Filter parameters to the current model
pars <- pars %>% filter(model == current_model) %>% select(-model) %>%
  select(all_of(active_comps))

pars <- pars[par_idx_range,]

# Read in seeds and choose the correct model
seeds <- readRDS(paste0(DATA_PATH, "seeds.RDS"))
seeds <- seeds[[current_model]]
seed <- as.integer(seeds[par_idx_range])

# Read in parallel/orthogonal/randomised directions
parallel_opt_dir <- read_csv(paste0(DATA_PATH, "parallel_traitdir.csv"), col_names = F)
orth_opt_dir <- read_csv(paste0(DATA_PATH, "orth_traitdir.csv"), col_names = F)


d_ruggedness <- list()

# Sample seed for optimum (different between replicates, same between backgrounds)
opt_seed <- sample(1:.Machine$integer.max, 1)

  # randomly sample an optimum
  set.seed(opt_seed)
  #parsMasked <- ParsMask(pars, model) 
  #model_comps <- CompsForModel(comps, model)
  
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

  RugRes_rand <- CalculateRuggednessParallel(pars, current_model, "Randomised", opt_rand, sigma,
                                        n = NUM_STEPS,
                                        nCores = future::availableCores(),
                                        seed = seed,
                                        path = DATA_PATH)

  RugRes_par <- CalculateRuggednessParallel(pars, current_model, "Parallel", opt_par, sigma,
                                        n = NUM_STEPS,
                                        nCores = future::availableCores(),
                                        seed = seed,
                                        path = DATA_PATH)
                                      
  RugRes_orth <- CalculateRuggednessParallel(pars, current_model, "Orthogonal", opt_orth, sigma,
                                        n = NUM_STEPS,
                                        nCores = future::availableCores(),
                                        seed = seed,
                                        path = DATA_PATH)

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
