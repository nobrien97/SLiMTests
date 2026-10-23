library(tidyverse)
library(deSolve)
library(mvtnorm)
library(broom)

DATA_PATH <- "/home/564/nb9894/tests/newMotifs/paper1/ruggedness/sim_range/R/" 
# DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/" 
SAVE_PATH <- "/scratch/ht96/nb9894/newMotifs/paper1/ruggedness/sim_range/"

SAVE_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/"

setwd(DATA_PATH)

# Load functions
source("./fitnesslandscapefunctions.R")

# input arguments: index to load parameters
# args <- commandArgs(trailingOnly = T)
# par_idx <- as.numeric(args[1])

# REPS_PER_RUN <- 100
# N_BACKGROUNDS <- 10 # Number of different fitness optima
# ROWS_PER_RUN <- REPS_PER_RUN * N_BACKGROUNDS

# Load parameter combinations and seeds
pars_all <- readRDS(paste0(DATA_PATH, "d_pc_grid_molcomp.RDS"))

# Load in phenotypic optima from simulations to use
d_opt <- readRDS(paste0(DATA_PATH, "d_opt_grid.RDS"))


# output
d_ruggedness <- vector(mode = "list", length = length(models))

# Iterate over models
for (current_model in models) {
  output_index <- match(current_model, models)
  active_comps <- molComp_names[[current_model]] 
  nComps <- length(active_comps)

  # # range of input rows to evaluate this run
  # par_idx_range <- (ROWS_PER_RUN * (par_idx - 1) + 1):(ROWS_PER_RUN * par_idx)

  # Filter parameters to the current model
  # Split by dataset and isAdapted
  pars <- pars_all[[output_index]] %>% 
    mutate(rowID = row_number()) %>%
    relocate(rowID, .before = 1)

  pars_rand <- pars %>%
    filter(dataset == "Randomised")
  
  pars_orth <- pars %>%
    filter(dataset == "Orthogonal")
  
  pars_par <- pars %>%
    filter(dataset == "Parallel")
  
  opt_rand <- full_join(
    d_opt,
    pars_rand,
    by = c("seed", "model", "dataset")
  ) %>%
    drop_na() %>%
    select(rowID, starts_with("opt"), starts_with("sigma")) %>%
    select(rowID, where(~sum(.) != 0)) # filter out unused columns
  
  opt_par <- full_join(
    d_opt,
    pars_par,
    by = c("seed", "model", "dataset")
  ) %>%
    drop_na() %>%
    select(rowID, starts_with("opt"), starts_with("sigma")) %>%
    select(rowID, where(~sum(.) != 0)) # filter out unused columns
  
  
  opt_orth <- full_join(
    d_opt,
    pars_orth,
    by = c("seed", "model", "dataset")
  ) %>%
    drop_na() %>%
    select(rowID, starts_with("opt"), starts_with("sigma")) %>%
    select(rowID, where(~sum(.) != 0)) # filter out unused columns
  
  # Filter pars to correct rows
  pars_rand <- pars_rand %>%
    select(rowID, all_of(active_comps), -dataset)
  pars_orth <- pars_orth %>%
    select(rowID, all_of(active_comps), -dataset)
  pars_par <- pars_par %>%
    select(rowID, all_of(active_comps), -dataset)
  
  
  # get optimum at random from the PC space
    RugRes_rand <- CalculateRuggednessLandscaper(pars_rand, current_model, "Randomised", opt_rand,
                                          n = 0,
                                          nCores = 2)
  
    RugRes_par <- CalculateRuggednessLandscaper(pars_par, current_model, "Parallel", opt_par,
                                          n = 0,
                                          nCores = 2)
                                        
    RugRes_orth <- CalculateRuggednessLandscaper(pars_orth, current_model, "Orthogonal", opt_orth,
                                          n = 0,
                                          nCores = 2)
  
    RugRes <- rbind(RugRes_rand, RugRes_par, RugRes_orth)
    
    RugRes[match(RugRes$rowID, pars$rowID), 
           c("seed", "isAdapted")] <- pars[, c("seed", "isAdapted")]
    
    d_ruggedness[[output_index]] <- distinct(RugRes)
}

# Combine across models
d_ruggedness <- data.table::rbindlist(d_ruggedness, fill = T)

# Add additional descriptor columns

# Rearrange columns
d_ruggedness <- d_ruggedness %>%
  select(seed, model, dataset, isAdapted, fitness,
  aX, KZX, aY, bY, KY, KZ, KXZ,
  aZ, bZ, zZ, h, gX)

write_csv(d_ruggedness, paste0(SAVE_PATH, "d_ruggedness_sim.csv"), col_names = F)
