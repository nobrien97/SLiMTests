# Use UMAP to plot fitness landscapes
# Helper functions, libraries etc.
HELPER_PATH <- "~/tests/newMotifs/paper1/ruggedness/umap/R/"

source(paste0(HELPER_PATH, "helperFns.R"))

# Read in model for this job
args <- commandArgs(trailingOnly = T)
model_name <- args[1]

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

# Filter by model input
d_ruggedness <- d_ruggedness %>%
  ungroup() %>%
  filter(model == model_name) %>%
  select(fitness, molComp_names[[model_name]])

# For this model, do a grid search to find the best n_neighbours to maximise
# correlation between distances between umap points and distances between 
# trait combinations
neighbours <- c(5, 15, 50, 100)
min_dists <- c(0.1, 0.25, 0.5, 0.99)
metrics <- c("euclidean", "manhattan")
combos <- expand.grid(neighbours, min_dists, metrics)
combos <- combos %>%
  rename(neighbour = Var1,
         min_dist = Var2,
         metric = Var3) %>%
  mutate(metric = as.character(metric))

# > sample(1:.Machine$integer.max, 1)
# [1] 1997523188
seed <- 1997523188

# Fit landscape on small sample so we can get a better idea of parameter differences
d_ruggedness_sample <- d_ruggedness %>%
  slice_sample(n = 1000)

umap_result <- vector(mode = "list", length = nrow(combos))
startTime <- as.numeric(Sys.time())
for (i in seq_len(nrow(combos))) {
  combo <- combos[i,]
  umap_result[[i]]$umap_data <- as_tibble(umap2(d_ruggedness_sample %>% select(-fitness),
                                                n_neighbors = combo$neighbour,
                                                min_dist = combo$min_dist,
                                                metric = combo$metric)) %>%
    rename(LV1 = V1, LV2 = V2)

  umap_result[[i]]$distances <- CalcDistancesUMAP(d_ruggedness_sample %>% select(-fitness),
                                                umap_result[[i]]$umap_data,
                                                n = 100)
}
endTime <- as.numeric(Sys.time())
print(paste("Seconds to finish UMAP grid search:", round(endTime - startTime, digits = 3)))

# Check which one is best
best_i <- which.max(unlist(lapply(umap_result, function(x) {
  return(x$distances$r2)
})))

print(paste0("best umap grid result = row ", best_i))

unlist(lapply(umap_result, function(x) {
  return(x$distances$r2)
}))

# Run with best parameter combination on the whole dataset
umap_big <- vector(mode = "list", length = 1)
set.seed(seed)

nThreads <- parallel::detectCores()

start_time <- as.numeric(Sys.time())
umap_model <- umap2(d_ruggedness %>% select(-fitness), 
      n_neighbors = combos[best_i,]$neighbour,
      min_dist = combos[best_i,]$min_dist,
      metric = combos[best_i,]$metric, ret_model = T,
      n_threads = nThreads)

umap_big[[1]]$umap_data <- as_tibble(umap_model$embedding) %>% 
  rename(LV1 = V1,LV2 = V2)

end_time <- as.numeric(Sys.time())
print(paste("Seconds to run big UMAP =", (end_time - start_time)))

# Save model so we can project adaptive walks onto it
save_uwot(umap_model, paste0(DATA_PATH, "sim_range/umap_", model_name))