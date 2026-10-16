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

REPS_PER_RUN <- 100
N_BACKGROUNDS <- 10 # Number of different fitness optima
ROWS_PER_RUN <- REPS_PER_RUN * N_BACKGROUNDS

# Load parameter combinations and seeds
pars_all <- readRDS(paste0(DATA_PATH, "pars.RDS"))
seeds_all <- readRDS(paste0(DATA_PATH, "seeds.RDS"))

# Read in parallel/orthogonal/randomised directions
parallel_opt_dir <- read_csv(paste0(DATA_PATH, "parallel_traitdir.csv"), col_names = F)
orth_opt_dir <- read_csv(paste0(DATA_PATH, "orth_traitdir.csv"), col_names = F)

# Load in phenotypic optima from simulations to use
SIM_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/"
rand_opt <- read_csv(paste0(SIM_PATH, "randomisedStartsM/slim_opt.csv"), 
                     col_types = paste0(c(rep("f", times = 2),
                                   rep("n", times = 12)), collapse = ""),
                     col_names = c("seed", "modelindex", paste0("opt.", 1:4), 
                                   paste0("sigma.", 1:4), paste0("traitShiftDir.", 1:4)))
par_opt <- read_csv(paste0(SIM_PATH, "parallelSel/slim_opt.csv"), 
                    col_types = paste0(c(rep("f", times = 2),
                                         rep("n", times = 12)), collapse = ""),
                    col_names = c("seed", "modelindex", paste0("opt.", 1:4), 
                                  paste0("sigma.", 1:4), paste0("traitShiftDir.", 1:4)))
orth_opt <- read_csv(paste0(SIM_PATH, "orthSel/slim_opt.csv"),
                     col_types = paste0(c(rep("f", times = 2),
                                          rep("n", times = 12)), collapse = ""),
                     col_names = c("seed", "modelindex", paste0("opt.", 1:4), 
                                   paste0("sigma.", 1:4), paste0("traitShiftDir.", 1:4)))

# Add combos
rand_opt <- AddCombosToDF(rand_opt)
# Only use recombination rate = 0.1 to match the parallel/orth sel case (modelindex 11:15)
rand_opt <- rand_opt %>%
  mutate(model = factor(model, levels = model_names, 
                        labels = model_names_noquote),
         dataset = "Randomised") %>%
  filter(r == 0.1)

par_opt <- AddCombosToDF(par_opt)


# output
d_ruggedness <- vector(mode = "list", length = length(models))

# Iterate over models
for (current_model in models) {
  output_index <- match(current_model, models)
  active_comps <- molComp_names[[current_model]] 
  nComps <- length(active_comps)

  # range of input rows to evaluate this run
  par_idx_range <- (ROWS_PER_RUN * (par_idx - 1) + 1):(ROWS_PER_RUN * par_idx)

  # Filter parameters to the current model
  pars <- pars_all %>% filter(model == current_model) %>% select(-model) %>%
    select(all_of(active_comps))
  
  pars <- pars[par_idx_range,]
  
  # Read in seeds and choose the correct model
  seeds <- as.integer(seeds_all[[current_model]])

  # get optimum at random from the PC space
  for (j in seq_len(N_BACKGROUNDS)) {
    set.seed(seeds[j])
    optMolComps <- as.data.frame(t(runif(nComps, min(pars), max(pars))))
    colnames(optMolComps) <- colnames(pars)
      
    # Write molcomps/starting values
    startInput <- tempfile(tmpdir = PBS_JOBFS, fileext = ".csv")
    startOpt <- tempfile(tmpdir = PBS_JOBFS, fileext = ".csv") 
    startOut <- tempfile(tmpdir = PBS_JOBFS, fileext = ".csv") 
    
    write.table(optMolComps, startInput, sep = ",", col.names = F, row.names = F)
    # Dummy optimum, fitness info doesn't matter here, but we want to feed the correct data
    # optimum data has per trait optimum values and sigmas combined, hence the * 2
    write.table(rbind(rep(1, times = NumTraits(current_model) * 2)), 
                startOpt, sep = ",", col.names = F, row.names = F)
    
    startTraits <- runLandscaper(startInput, startOut, startOpt, current_model, 1, useID = FALSE)
    startTraits <- unlist(startTraits[1,2:(NumTraits(current_model) + 1)])
    names(startTraits) <- NULL
    sigma <- CalcSelectionSigmas(startTraits, 0.05, 0.1, 0.1)
    par_dir_model <- unlist(parallel_opt_dir[which(models == current_model),c(1:length(startTraits), ncol(parallel_opt_dir))])
    orth_dir_model <- unlist(orth_opt_dir[which(models == current_model),c(1:length(startTraits), ncol(orth_opt_dir))])
  
    opt_rand <- CalcOptima(startTraits, sigma, 0.95)
    opt_par <- CalcNewOptimumAlongVector(startTraits, sigma, 0.95, par_dir_model)
    opt_orth <- CalcNewOptimumAlongVector(startTraits, sigma, 0.95, orth_dir_model)
    
    RugRes_rand <- CalculateRuggednessLandscaper(pars, current_model, "Randomised", opt_rand, sigma,
                                          n = 0,
                                          nCores = 2)
  
    RugRes_par <- CalculateRuggednessLandscaper(pars, current_model, "Parallel", opt_par, sigma,
                                          n = 0,
                                          nCores = 2)
                                        
    RugRes_orth <- CalculateRuggednessLandscaper(pars, current_model, "Orthogonal", opt_orth, sigma,
                                          n = 0,
                                          nCores = 2)
  
    RugRes <- rbind(RugRes_rand, RugRes_par, RugRes_orth)
    
    # Set identifiers
    RugRes$bkg <- j
    
    d_ruggedness[[output_index]][[j]] <- distinct(RugRes)
  }

  # Combine within models
  d_ruggedness[[output_index]] <- data.table::rbindlist(d_ruggedness[[output_index]])
}

# Combine across models
d_ruggedness <- data.table::rbindlist(d_ruggedness, fill = T)

# Rearrange columns
d_ruggedness <- d_ruggedness %>%
  select(step, model, dataset, fitness,
  aX, KZX, aY, bY, KY, KZ, KXZ,
  aZ, bZ, zZ, h, gX, bkg)

write_csv(d_ruggedness, paste0(SAVE_PATH, "d_ruggedness_", par_idx, ".csv"), col_names = F)