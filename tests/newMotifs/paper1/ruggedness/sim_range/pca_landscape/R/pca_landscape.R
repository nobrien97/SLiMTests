# Use UMAP to plot fitness landscapes
# Helper functions, libraries etc.
HELPER_PATH <- "~/tests/newMotifs/paper1/ruggedness/umap/R/"

source(paste0(HELPER_PATH, "helperFns.R"))

# Read in model for this job
args <- commandArgs(trailingOnly = T)
ip <- args[1] # IP address of the node for connecting to the h2o java server

DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"
setwd(paste0(DATA_PATH, "sim_range"))
d_ruggedness <- data.table::fread(paste0(DATA_PATH, "sim_range/d_ruggedness_permolcomp.csv"), 
                                  header = F,
                                  colClasses = c("integer", "character", "character",
                                                 rep("numeric", times = 5),
                                                 "integer", "integer",
                                                 rep("numeric", times = 12),
                                                 "character", "integer"),
                                  col.names = c("step", "model", "dataset", "fitness", "startW", 
                                                "endW", "netChangeW", "sumChangeW", "numFitnessHoles", 
                                                "nSteps", "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                                                "aZ", "bZ", "h", "gX", "zZ",
                                                "molComp", "bkg"))
# Sample only the first step 
d_ruggedness <- d_ruggedness %>%
  filter(step == 1)

# Filter to only the columns we care about
d_ruggedness <- d_ruggedness %>%
  select(2:4, 11:22)

d_ruggedness_all <- d_ruggedness

# Now test PCA and autoencoder
h2o.init(ip = ip,
        port = 12345, 
        startH2O = FALSE)

landscape_list <- vector(mode = "list", length = 5)

for (model_name in model_names_noquote) {
  # Filter by model input
  d_ruggedness <- d_ruggedness_all %>%
    ungroup() %>%
    filter(model == model_name) %>%
    select(fitness, molComp_names[[model_name]])

  # > sample(1:.Machine$integer.max, 1)
  # [1] 1997523188
  seed <- 1997523188

  features <- as.h2o(d_ruggedness)
  n_features <- ncol(d_ruggedness) - 1

  pca_result <- h2o.prcomp(
    x = 2:ncol(d_ruggedness),
    training_frame = features,
    pca_method = "GramSVD",
    transform = "NONE", # No transformation, all traits on same scale already
    impute_missing = T,
    k = 1 - ncol(d_ruggedness),
    seed = seed
  )
  print("PCA done")
  print(pca_result)
  print("Eigenvectors:")
  print(pca_result@model$eigenvectors)

  # Save model so we can project adaptive walks onto it
  h2o.saveModel(pca_result, paste0(DATA_PATH, "sim_range/pca_", model_name, ".h2o"))

  pca_codings <- h2o.predict(pca_result, features)
  d_pca_codings <- pca_codings %>%
    as.data.frame() %>%
    select(PC1, PC2) %>%
    rename(LV1 = PC1,
          LV2 = PC2)

  saveRDS(d_pca_codings, paste0(DATA_PATH, "sim_range/d_pca_codings", model_name, ".RDS"))



  #d_pca_codings <- readRDS(paste0(DATA_PATH, "sim_range/d_pca_codings", model_name, ".RDS"))

  # Plot landscape
  d_landscape <- d_pca_codings %>%
    mutate(z = d_ruggedness$fitness) %>%
    rename(x = LV1,
          y = LV2) 

  # Fit surface with multilevel B-splines
  mba_landscape <- mba.surf(d_landscape, 
                              no.X = 1000, no.Y = 1000)

  d_landscape <- expand.grid(x = mba_landscape$xyz.est$x,
                            y = mba_landscape$xyz.est$y)
  d_landscape$z <- as.vector(mba_landscape$xyz.est$z)

  d_landscape <- d_landscape %>%
    mutate(z = if_else(is.na(z) | z < 0, 0, z),
    model = model_name)

  landscape_list[[match(model_name, model_names_noquote)]] <- d_landscape
}

out <- data.table::rbindlist(landscape_list)
saveRDS(out, paste0(DATA_PATH, "sim_range/pca_landscape/d_landscape_pca.RDS"))

# Adaptive walks
# Load in phenotype data
d_qg_tot <- readRDS("/g/data/ht96/nb9894/newMotifs/paper1/d_qg_tot.RDS")

# seed <- sample(1:.Machine$integer.max, 1)
# [1] 214952207
seed <- 214952207

# Sample some walks, split by model
d_qg_tot %>%
  mutate(model = factor(model, levels = levels(d_qg_tot$model),
                        labels = model_names_noquote)) %>%
  group_by(model, isAdapted) %>%
  filter(seed %in% sample(unique(seed), 3)) -> d_qg_sample





d_walks <- vector(mode = "list", length = length(model_names_noquote))
# Transform those walks into autoencoder space
for (i in seq_along(model_names_noquote)) {
  model_name <- model_names_noquote[i]
  pca_codings <- readRDS(paste0(DATA_PATH, "sim_range/d_pca_codings", model_name, ".RDS"))
  
  n_comps <- length(molComp_names[[model_name]])
  
  # Select molecular components
  d_qg_pca <- d_qg_sample %>% filter(model == model_name) %>%
    select(gen, seed, model, isAdapted, paste0("meanMC", 1:n_comps)) %>%
    mutate(across(starts_with("meanMC"), log))
  
  # Rename to match molcomp names
  colnames(d_qg_pca)[5:ncol(d_qg_pca)] <- molComp_names[[model_name]]
  
  # Load PCA model
  model_path <- list.files(paste0(DATA_PATH, "sim_range/pca_", model_name, ".h2o"), full.names = T)[1]
  
  pca <- h2o.loadModel(model_path)
  
  newpoints <- as.h2o(d_qg_pca %>% ungroup() %>% select(5:ncol(d_qg_pca)))
  #walk_predictions <- h2o.deepfeatures(ae, newpoints, layer = 2)
  walk_predictions <- h2o.predict(pca, newpoints)
  
  d_qg_pca[,c("LV1", "LV2")] <- as_tibble(walk_predictions)
  
  # Add missing molcomps and rearrange columns to match others
  extra_cols <- names(all_molcomp_features)[!(names(all_molcomp_features) %in% colnames(d_qg_pca))]
  d_qg_pca[,extra_cols] <- 0.0
  d_qg_pca <- d_qg_pca %>%
    relocate(all_of(names(all_molcomp_features)))
  
  d_walks[[i]] <- d_qg_pca
}

d_walks <- data.table::rbindlist(d_walks, fill = T)
saveRDS(d_walks, paste0(DATA_PATH, "sim_range/pca_landscape/d_walks.RDS"))

