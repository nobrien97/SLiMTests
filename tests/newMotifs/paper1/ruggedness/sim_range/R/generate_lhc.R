library(DoE.wrapper)


# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 1649063102
seed <- 1649063102

set.seed(seed)

models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")

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


nComps <- length(comps)
NUM_RUNS <- 10000
NUM_BACKGROUNDS <- 10

# Max comp size is set per component and per model via find_landscape_range.R
d_molcomp_maxvals <- readRDS("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/R/d_molcomp_maxvals.RDS")


output <- vector(mode = "list", length = length(models))

for (m in seq_along(models)) {
  current_model <- models[m]
  # Get the list of relevant parameters
  active_comps <- molComp_names[[current_model]] 
  n_active_comps <- length(active_comps)
  
  
  molcomp_min_values <- d_molcomp_maxvals %>% ungroup() %>%
    filter(model == current_model) %>% select(min_value) %>% unlist()
  names(molcomp_min_values) <- active_comps 
  
  molcomp_max_values <- d_molcomp_maxvals %>% ungroup() %>%
    filter(model == current_model) %>% select(max_value) %>% unlist()
  names(molcomp_max_values) <- active_comps 
  
  # Generate hypercube of parameters
  # Hypercube is NUM_RUNs per molecular component, per genetic background
  pars <- lhs.design(NUM_RUNS, 1, type = "random")
  pars <- pars * MAX_COMP_SIZE
  
  # Repeat for each component
  test_mat <- matrix(rep(pars$X1, times = n_active_comps), ncol = n_active_comps)
  
  # rescale each column by min/max vals
  pars <- sweep(sweep(test_mat, 2, molcomp_max_values - molcomp_min_values, "*"), 
                      2, molcomp_min_values, "+")
  
  # Generate backgrounds
  parBackgrounds <- matrix(runif(n_active_comps * NUM_BACKGROUNDS), ncol = (n_active_comps))
  
  scaled_backgrounds <- sweep(sweep(parBackgrounds, 2, molcomp_max_values - molcomp_min_values, "*"), 
                              2, molcomp_min_values, "+")
  
  # Replicate rows for the number of molecular components
  parBackgrounds <-  scaled_backgrounds %x% rep(1, n_active_comps)
  
  # Replicate rows for the number of replicate values
  parBackgrounds <-  rep(1, NUM_RUNS) %x% parBackgrounds
  
  # replace diagonal components for each component and background
  for (i in seq_len(NUM_RUNS)) {
    # Fill NUM_BACKGROUNDS diagonals
    for (j in seq_len(NUM_BACKGROUNDS)) {
      k <- (i - 1) * n_active_comps + j
      offset_start <- ((i - 1) * n_active_comps * NUM_BACKGROUNDS) + ((j - 1) * n_active_comps) + 1 # number of previous runs
      offset_end <-  offset_start + (n_active_comps - 1)
      diag(parBackgrounds[offset_start:offset_end,]) <- pars[i]
    }
  }
  
  colnames(parBackgrounds) <- active_comps
  pars_df <- as.data.frame(parBackgrounds)
  pars_df$model <- current_model
  output[[m]] <- pars_df
}

output <- data.table::rbindlist(output, fill = T)

saveRDS(output, "pars.RDS")


seeds <- sample(1:.Machine$integer.max, NUM_RUNS)

seeds_list <- vector(mode = "list")
for (i in seq_along(models)) {
  current_model <- models[i]
  seeds_list[[current_model]] <- rep(seeds, each = NUM_BACKGROUNDS * length(molComp_names[[current_model]]))
}
saveRDS(seeds_list, "seeds.RDS")
