library(DoE.wrapper)
library(tidyverse)
library(factoextra)

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
NUM_RUNS <- 1000000
N_BACKGROUNDS <- 10

d_pc_qg <- readRDS("/mnt/d/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/d_pc_qg.RDS")

# All motifs explain >=80% of variance in first two PCs
fviz_eig(d_pc_qg[[1]])
fviz_pca_biplot(d_pc_qg[[1]])

# For each model, get the minimum and maximum PC1 and PC2 loadings, generate a grid along those
# two directions
d_grid <- vector(mode = "list", length = length(d_pc_qg))

for (i in seq_along(d_pc_qg)) {
  prc <- d_pc_qg[[i]]
  
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
  d_traitspace$model <- models[i]

  d_grid[[i]] <- d_traitspace
}

d_grid <- data.table::rbindlist(d_grid, fill = T)

d_grid <- d_grid %>% dplyr::select(model, all_of(comps))

saveRDS(d_grid, "pars.RDS")


seeds <- sample(1:.Machine$integer.max, N_BACKGROUNDS)

seeds_list <- vector(mode = "list")
for (i in seq_along(models)) {
  current_model <- models[i]
  seeds_list[[current_model]] <- seeds
}
saveRDS(seeds_list, "seeds.RDS")
