# Use UMAP to plot fitness landscapes
# Helper functions, libraries etc.
source("./helperFns.R")

# Read in model for this job
args <- commandArgs(trailingOnly = T)
model_name <- args[1]

DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"
setwd(paste0(DATA_PATH, "log3"))
d_ruggedness <- data.table::fread(paste0(DATA_PATH, "log3/d_ruggedness_permolcomp.csv"), 
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
save_uwot(umap_model, paste0(DATA_PATH, "umap_", model_name))


# Now test PCA and autoencoder
h2o.init()

features <- as.h2o(d_ruggedness)
n_features <- ncol(d_ruggedness) - 1
ae <- h2o.deeplearning(x = 2:ncol(d_ruggedness),
                       training_frame = features,
                       autoencoder = T,
                       seed = seed,
                       hidden = c(n_features, 2, n_features),
                       epochs = 100,
                       activation = "Tanh",
                       sparse = F)
print("Finished training autoencoder")
print(ae)
ae_codings <- h2o.deepfeatures(ae, features, layer = 2)

# Save model so we can project adaptive walks onto it
saveRDS(ae, paste0(DATA_PATH, "ae_", model_name, ".RDS"))

d_ae <- as.data.frame(ae_codings) 
# LV for latent variable
colnames(d_ae) <- c("LV1", "LV2")

saveRDS(d_ae, paste0(DATA_PATH, "d_ae_codings", model_name, ".RDS"))

######
# PCA
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
saveRDS(pca_result, paste0(DATA_PATH, "pca_", model_name, ".RDS"))


pca_codings <- h2o.predict(pca_result, features)
d_pca_codings <- pca_codings %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rename(LV1 = PC1,
         LV2 = PC2)

saveRDS(d_pca_codings, paste0(DATA_PATH, "d_pca_codings", model_name, ".RDS"))

# Calculate distances to find which best fits data
d_pca_dist <- CalcDistancesUMAP(d_ruggedness %>% select(-fitness),
                                d_pca_codings,
                                n = 100000, seed = seed)
d_ae_dist <- CalcDistancesUMAP(d_ruggedness %>% select(-fitness),
                                d_ae,
                                n = 100000, seed = seed)
d_umap_dist <- CalcDistancesUMAP(d_ruggedness %>% select(-fitness),
                                 umap_big[[1]]$umap_data,
                                 n = 100000, seed = seed)
d_dr.dist <- rbind(d_pca_dist$dist.frame %>% mutate(id = row_number(), dr.method = "PCA",
                                                    r2 = d_pca_dist$r2),
                   d_ae_dist$dist.frame %>% mutate(id = row_number(), dr.method = "Autoencoder",
                                                   r2 = d_ae_dist$r2),
                   d_umap_dist$dist.frame %>% mutate(id = row_number(), dr.method = "UMAP",
                                                     r2 = d_umap_dist$r2))

d_dr.dist <- d_dr.dist %>%
  rename(lv_dist = umap_dist)

# Output combined distance data
saveRDS(d_dr.dist, paste0(DATA_PATH, "d_dr.dist_", model_name, ".RDS"))

d_dr.r2 <- d_dr.dist %>%
  group_by(dr.method) %>%
  filter(row_number() == 1)

# Plot distance between points in both trait and umap coordinates
ggplot(d_dr.dist,
       aes(x = trait_dist, y = lv_dist)) +
  facet_nested(. ~ "Method" + dr.method) +
  geom_point(shape = 21, alpha = 0.2) +
  theme_bw() +
  geom_text(data = d_dr.r2,
            mapping = aes(x = 3.5, y = 5.5,
                          label = TeX(paste0("$r^2 = ", 
                                             round(r2, digits = 3)), 
                                      output = "character")),
            parse = T) +
  labs(x = "Trait-space distance",
       y = "Latent-space distance") +
  theme(text = element_text(size = 14)) -> plt_dist
ggsave(paste0("plt_dist_", model_name, ".png"), 
       plt_dist, device = png, width = 12, height = 5,
       dpi = 600)

# Plot landscape using the method producing the greatest r2
d_landscape <- rbind(d_pca_codings %>% mutate(dr.method = "PCA"),
                     d_ae %>% mutate(dr.method = "Autoencoder"),
                     umap_big[[1]]$umap_data %>% mutate(dr.method = "UMAP"))

best_method <- (d_dr.r2 %>% ungroup() %>% filter(r2 == max(r2)))$dr.method[1]

# Plot landscape
d_landscape <- d_landscape %>%
  filter(dr.method == best_method) %>%
  select(-dr.method) %>%
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
  mutate(z = if_else(is.na(z) | z < 0, 0, z))

saveRDS(d_landscape, paste0(DATA_PATH, "d_landscape_", model_name, ".RDS"))

ggplot(d_landscape,
       aes(x = x, y = y, z = z)) +
  geom_raster(aes(fill = z)) +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  ggtitle(paste0("Fitness landscape for ", model_name)) +
  coord_equal() +
  labs(x = "LV1", y = "LV2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line')) -> plt_landscape
plt_landscape
ggsave(paste0("plt_landscape_", model_name, ".png"), 
       plt_landscape, device = png, width = 7, height = 7,
       dpi = 600)




# Plot without any interpolation by downsampling
downsample <- 1/20 
d_landscape_ds <- umap_big[[1]]$umap_data %>% mutate(fitness = d_ruggedness$fitness) %>%
  group_by(x = downsample * round(LV1 / downsample),
           y = downsample * round(LV2 / downsample)) %>%
  summarise(z = mean(fitness))

ggplot(d_landscape_ds,
       aes(x = x, y = y, fill = z, z = z, group = z)) +
  geom_raster(interpolate = F) +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "LV1", y = "LV2", 
       fill = "Fitness") +
  theme_bw() +
  coord_equal() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line')) -> plot_landscape_ds
ggsave(paste0("plt_landscape_ds_", model_name, ".png"), 
       plot_landscape_ds, device = png, width = 7, height = 7,
       dpi = 600)

# 30 million samples per group, how do we do that?
# fit UMAP on subset?
## Parametric UMAP is another option, fit on 