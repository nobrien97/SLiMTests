library(DoE.wrapper)
library(tidyverse)

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

# Max comp size is set per component and per model via find_landscape_range.R
d_molcomp_maxvals <- readRDS("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/R/d_molcomp_maxvals.RDS")

d_molcomp_maxvals <- d_molcomp_maxvals %>%
  drop_na()


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
  
  # Generate uniform grid
  pars <- lhs.design(NUM_RUNS, 1, type = "random")

  pars_df <- as.data.frame(pars)
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
