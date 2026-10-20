library(tidyverse)
library(latex2exp)

setwd("/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/sim_range")

COMBO_PATH <- "~/tests/newMotifs/R/"
d_combos <- read.table(paste0(COMBO_PATH, "combos.csv"), header = F,
                       col.names = c("model", "r"))


# Adds the parameter combination to a dataframe
AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]])
}

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

model_names_noquote <- c("NAR", "PAR", "FFLC1", 
                         "FFLI1", "FFBH")

all_molcomp_features <- c(
  "aX" = TeX("$\\alpha_X$", output = "character"),
  "KZX" = TeX("$K_{ZX}$", output = "character"),
  "aY" = TeX("$\\alpha_Y$", output = "character"),
  "bY" = TeX("$\\beta_Y$", output = "character"),
  "KY" = TeX("$K_Y$", output = "character"),
  "aZ" = TeX("$\\alpha_Z$", output = "character"),
  "bZ" = TeX("$\\beta_Z$", output = "character"),
  "KZ" = TeX("$K_{Z}$", output = "character"),
  "KXZ" = TeX("$K_{XZ}$", output = "character"),
  "zZ" = TeX("$\\zeta_Z$", output = "character"), # baseline expression
  "h" = TeX("$h$", output = "character"), # hill coefficient
  "gX" = TeX("$\\gamma_X$", output = "character") # X multiplier
)

GetMolCompTable <- function(x, model) {
  # x is a dataframe, each row is a timepoint
  # we need each row to be a 1000xN dataframe
  colnames(x) <- NULL
  if (model == "NAR" | model == "PAR") {
    # Get molecular component dataframe
    result <- lapply(seq_len(nrow(x)), function(i) { 
      x_s <- unlist(x[i,])
      data.frame(row_id = i,
                 ind_id = 1:(length(x_s) / 7),
                               aZ    = log(x_s[seq(1, length(x_s), by = 7)]),
                               bZ    = log(x_s[seq(2, length(x_s), by = 7)]),
                               KZ    = log(x_s[seq(3, length(x_s), by = 7)]),
                               KXZ   = log(x_s[seq(4, length(x_s), by = 7)]),
                               zZ  = log(x_s[seq(5, length(x_s), by = 7)]),
                               h     = log(x_s[seq(6, length(x_s), by = 7)]),
                               gX = log(x_s[seq(7, length(x_s), by = 7)]))
      })
  } else if (model == "FFLC1" | model == "FFLI1") {
    # Get molecular component dataframe
    result <- lapply(seq_len(nrow(x)), function(i) { 
      x_s <- unlist(x[i,])
      data.frame(row_id = i,
                 ind_id = 1:(length(x_s) / 9),
                               aY    = log(x_s[seq(1, length(x_s), by = 9)]),
                               bY    = log(x_s[seq(2, length(x_s), by = 9)]),
                               KY    = log(x_s[seq(3, length(x_s), by = 9)]),
                               aZ    = log(x_s[seq(4, length(x_s), by = 9)]),
                               bZ    = log(x_s[seq(5, length(x_s), by = 9)]),
                               KXZ   = log(x_s[seq(6, length(x_s), by = 9)]),
                               zZ  = log(x_s[seq(7, length(x_s), by = 9)]),
                               h     = log(x_s[seq(8, length(x_s), by = 9)]),
                               gX = log(x_s[seq(9, length(x_s), by = 9)]))
    })
  } else if (model == "FFBH") {
    result <- lapply(seq_len(nrow(x)), function(i) { 
      x_s <- unlist(x[i,])
      data.frame(row_id = i,
                 ind_id = 1:(length(x_s) / 11),
                               aX    = log(x_s[seq(1, length(x_s), by = 11)]),
                               KZX   = log(x_s[seq(2, length(x_s), by = 11)]),
                               aY    = log(x_s[seq(3, length(x_s), by = 11)]),
                               bY    = log(x_s[seq(4, length(x_s), by = 11)]),
                               KY    = log(x_s[seq(5, length(x_s), by = 11)]),
                               aZ    = log(x_s[seq(6, length(x_s), by = 11)]),
                               bZ    = log(x_s[seq(7, length(x_s), by = 11)]),
                               KXZ   = log(x_s[seq(8, length(x_s), by = 11)]),
                               zZ  = log(x_s[seq(9, length(x_s), by = 11)]),
                               h     = log(x_s[seq(10, length(x_s), by = 11)]),
                               gX = log(x_s[seq(11, length(x_s), by = 11)]))
    })
  }
  
  # 
  return(data.table::rbindlist(result))
}

# Number of combinations per grid
NUM_RUNS <- 10000

DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/"

# Load data
# max number of cols is 11 molcomps * 1000 individuals + gen, seed, modelindex = 11003
max_cols <- 11003
d_molcomp_rand <- as_tibble(read.table(paste0(DATA_PATH, "randomisedStartsM/slim_sampled_moltrait.csv"),
                           colClasses = c("integer", rep("factor", times = 2),
                                                rep("numeric", times = 11000)),
                           sep = ",",
                           col.names = c("gen", "seed", "modelindex", paste0("X", 1:11000)),
                           fill = T, header = F
                           ))

d_molcomp_orth <- as_tibble(read.table(paste0(DATA_PATH, "orthSel/slim_sampled_moltrait.csv"),
                                       colClasses = c("integer", rep("factor", times = 2),
                                                      rep("numeric", times = 11000)),
                                       sep = ",",
                                       col.names = c("gen", "seed", "modelindex", paste0("X", 1:11000)),
                                       fill = T, header = F
))

d_molcomp_par <- as_tibble(read.table(paste0(DATA_PATH, "parallelSel/slim_sampled_moltrait.csv"),
                                       colClasses = c("integer", rep("factor", times = 2),
                                                      rep("numeric", times = 11000)),
                                       sep = ",",
                                       col.names = c("gen", "seed", "modelindex", paste0("X", 1:11000)),
                                       fill = T, header = F
))

# Match the indices of par and orth
d_molcomp_rand <- d_molcomp_rand %>%
  filter(as.numeric(as.character(modelindex)) > 10)

d_molcomp_rand <- AddCombosToDF(d_molcomp_rand) %>%
  mutate(dataset = "Randomised")

d_molcomp_orth <- AddCombosToDF(d_molcomp_orth) %>%
  mutate(dataset = "Orthogonal")

d_molcomp_par <- AddCombosToDF(d_molcomp_par) %>%
  mutate(dataset = "Parallel")

d_molcomp <- rbind(d_molcomp_par, d_molcomp_orth, d_molcomp_rand)

d_molcomp <- d_molcomp %>%
  mutate(dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
         model = factor(model, levels = model_names_noquote)) %>% select(-r)

# For each model/dataset, sample an adapted replicate and a maladapted replicate
# Load data
d_qg_tot <- readRDS("/g/data/ht96/nb9894/newMotifs/paper1/d_qg_tot.RDS")

# > seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 2081765552
set.seed(2081765552)

samples <- d_qg_tot %>%
  filter(gen == 60000, log10(r) == -1) %>%
  group_by(model, dataset, isAdapted) %>%
  slice_sample(n = 1) %>%
  select(model, dataset, isAdapted, seed) %>%
  mutate(model = factor(model, levels = c("'NAR'", "'PAR'",
                                          "'FFLC1'", "'FFLI1'",
                                          "'FFBH'"),
                        labels = model_names_noquote))

# filter to the sampled choices
d_molcomp_filtered <- inner_join(
  d_molcomp %>% filter(gen == 60000) %>% select(-c(gen, modelindex)),
  samples,
  by = c("seed", "model", "dataset")
)

# Convert from many columns to few
# Will have to split per model
d_molcomp_split <- d_molcomp_filtered %>%
  group_by(model) %>%
  group_split()

# Each row in d_molcomp_split is a grid
d_molcomp_table_list <- purrr::map(seq_along(model_names_noquote), function(i) {
  x <- d_molcomp_split[[i]]
  dataset <- rep(x$dataset, each = 1000) # Repeat each value by the number of individuals
  seed <- rep(x$seed, each = 1000)
  isAdapted <- rep(x$isAdapted, each = 1000)
  
  # Drop columns with NA values
  x <- x[,colSums(is.na(x)) == 0]
  
  # Remove label columns
  x <- x %>% select(-c(seed, model, dataset, isAdapted))
  
  result <- GetMolCompTable(x, model_names_noquote[i])
  result$model <- model_names_noquote[i]
  result$dataset <- dataset
  result$isAdapted <- isAdapted
  result$seed <- seed
  
  # Pivot longer so we have a molcomp table
  # result <- result %>%
  #   pivot_longer(cols = molComp_names[[model_names_noquote[[i]]]],
  #                names_to = "molComp_label",
  #                values_to = "molComp_value")
  
  return(as_tibble(result))
})

# We have a list of tibbles (a tibble per model), each with the individual data 
# for the example points. We need to do a PCA on this to transform to 2D, 
# then we need to build an even grid in PCA space and transform it back

d_pc_inds <- purrr::map(seq_along(model_names_noquote), function(i) {
  x <- d_molcomp_table_list[[i]]
  model_name <- model_names_noquote[i]
  molcomps <- molComp_names[[model_name]]
  
  # Split by dataset and isAdapted
  x_d_keys <- x %>%
    group_by(dataset, isAdapted, seed) %>%
    group_keys() %>%
    mutate(model = model_name)
  
  x_d <- x %>%
    group_by(dataset, isAdapted, seed) %>%
    group_split()
  
  # For each group in x_d, create a PCA from columns in the molComp_names vector
  grouped_pc <- lapply(seq_along(x_d), function(j) {
    pc_result <- list()
    x_d_j <- x_d[[j]] %>% select(all_of(molcomps))
    pc_result[["pca"]] <- prcomp(x_d_j)
    
    pc_result[["id"]] <- x_d_keys[j,]
    return(pc_result)
  })
  
  return(grouped_pc)
})

pc_grid <- purrr::map(seq_along(d_pc_inds), function(i) {
  model_pcs <- d_pc_inds[[i]]
  
  result <- purrr::map(seq_along(model_pcs), function(j) {
    cur_pc_element <- model_pcs[[j]]
    prc <- cur_pc_element$pca
    keys <- cur_pc_element$id
    
    mins <- unlist(apply(prc$x[,1:2], 2, min))
    maxs <- unlist(apply(prc$x[,1:2], 2, max))
    
    pc1_seq <- seq(mins[1], maxs[1], length.out = sqrt(NUM_RUNS))
    pc2_seq <- seq(mins[2], maxs[2], length.out = sqrt(NUM_RUNS))
    
    generated <- expand.grid(PC1 = pc1_seq,
                             PC2 = pc2_seq)
    
      # Translate the grid into component space
      d_traitspace <- (as.matrix(generated[,1:2]) %*% t(prc$rotation[,1:2])) + prc$center
      # Set any negative values to 1 (default value)
      d_traitspace[d_traitspace < 0] <- 1
      
      d_traitspace <- as_tibble(d_traitspace)
      d_traitspace$model <- model_names_noquote[[i]]
      d_traitspace$seed <- keys$seed
      d_traitspace$dataset <- keys$dataset
      d_traitspace$isAdapted <- keys$isAdapted
      
      return(d_traitspace)
  })
  
  data.table::rbindlist(result)
})

saveRDS(pc_grid, "d_pc_grid_molcomp.RDS")
saveRDS(d_pc_inds, "d_pc_molcomp.RDS")


# Exploring the whole space is very hard with 11 components, so we will construct our 
# fitness landscape along the PC1 and PC2 axes of the real data
# We can save the principal components and randomly generate points within the range of 
# PC1 and PC2 values, translate back to trait space, compute the fitness, then convert
# back to PC space
# This is done on a per simulation basis so we can align our fitness optima with the landscape
# and plot the individual walks

# Data is based on individual samples at end of walk, we can grab the mean values from d_qg_tot to 
# plot the entire walk on that surface
pc_grid <- readRDS("/mnt/d/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/d_pc_grid_molcomp.RDS")
d_pc_inds <- readRDS("/mnt/d/SLiMTests/tests/newMotifs/paper1/ruggedness/sim_range/d_pc_molcomp.RDS")


