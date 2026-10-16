library(DoE.wrapper)
library(tidyverse)
library(factoextra)
library(cowplot)

# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 1649063102
seed <- 1649063102

set.seed(seed)

models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "zZ", "h", "gX")

molComp_names <- list("NAR" = c(
  # NAR and PAR
  "aZ",
  "bZ",
  "KZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
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
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
)
)

molComp_names[["PAR"]] <- molComp_names[["NAR"]]
molComp_names[["FFLI1"]] <- molComp_names[["FFLC1"]]

nComps <- length(comps)
NUM_RUNS <- 100000
N_BACKGROUNDS <- 10

d_pc_qg <- readRDS("/mnt/d/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/d_pc_qg.RDS")

# Sample N_BACKGROUNDS random seeds to generate from per model
# > sample(1:.Machine$integer.max, 1)
# [1] 1629634382
set.seed(1629634382)
landscape_candidates <- list()
for (model in models) {
  for (dataset in c("Parallel", "Orthogonal", "Randomised")) {
    seed_names <- names(d_pc_qg[[model]][[dataset]])
    n_samples <- min(length(seed_names), N_BACKGROUNDS)
    seeds <- sample(seed_names, n_samples)
    
    landscape_candidates[[model]][[dataset]] <- d_pc_qg[[model]][[dataset]][seeds]
  }
}


# For each model, get the minimum and maximum PC1 and PC2 loadings, generate a grid along those
# two directions
d_grid <- list()

# Yucky triple nested loop
for (cur_model in models) {
  for (cur_dataset in c("Parallel", "Orthogonal", "Randomised")) {
    seed_names <- names(landscape_candidates[[cur_model]][[cur_dataset]])
    for (cur_seed in seed_names) {
      # Get the right landscape
      prc <- landscape_candidates[[cur_model]][[cur_dataset]][[cur_seed]]
      
      mins <- unlist(apply(prc$x[,1:2], 2, min))
      maxs <- unlist(apply(prc$x[,1:2], 2, max))
      
      # Generate grid
      pc1_seq <- seq(mins[1], maxs[1], length.out = sqrt(NUM_RUNS))
      pc2_seq <- seq(mins[2], maxs[2], length.out = sqrt(NUM_RUNS))
      generated <- expand.grid(PC1 = pc1_seq,
                               PC2 = pc2_seq)
      
    # Translate the grid into component space
      d_traitspace <- (as.matrix(generated[,1:2]) %*% t(prc$rotation[,1:2])) + prc$center
      # Set any negative values to 1 (default value)
      d_traitspace[d_traitspace < 0] <- 1
      
      d_traitspace <- as.data.frame(d_traitspace)
      d_traitspace$model <- cur_model
      d_traitspace$seed <- cur_seed
      d_traitspace$dataset <- cur_dataset
    
      d_grid[[cur_model]][[cur_dataset]][[cur_seed]] <- d_traitspace
    }
    # Collapse the list after we're done with this model/dataset combo
    d_grid <- data.table::rbindlist(d_grid[[cur_model]][[cur_dataset]])
  }
  # Collapse again after we're done with this model
  d_grid <- data.table::rbindlist(d_grid[[cur_model]])
}

# Final collapse to combine across models
d_grid_test <- data.table::rbindlist(d_grid, fill = T)

d_grid <- d_grid %>% dplyr::select(model, seed, dataset, all_of(comps))

saveRDS(d_grid, "pars.RDS")


seeds <- sample(1:.Machine$integer.max, N_BACKGROUNDS)

seeds_list <- vector(mode = "list")
for (i in seq_along(models)) {
  current_model <- models[i]
  seeds_list[[current_model]] <- seeds
}
saveRDS(seeds_list, "seeds.RDS")
