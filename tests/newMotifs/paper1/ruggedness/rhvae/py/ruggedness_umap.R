# Use UMAP to plot fitness landscapes
# Helper functions, libraries etc.
source("./helperFns.R")


# Read in model for this job
args <- commandArgs(trailingOnly = T)
model_name <- args[1]

DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"
d_ruggedness <- data.table::fread(paste0(DATA_PATH, "log3/d_ruggedness_permolcomp.csv"), 
                                  header = F)

setwd(paste0(DATA_PATH, "log3"))

colnames(d_ruggedness) <- c("step", "model", "dataset", "fitness", "startW", 
                            "endW", "netChangeW", "sumChangeW", "numFitnessHoles", 
                            "nSteps", "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                            "aZ", "bZ", "h", "gX", "zZ",
                            "molComp", "bkg")

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
    rename(UMAP1 = V1,UMAP2 = V2)
  
  umap_result[[i]]$distances <- CalcDistancesUMAP(d_ruggedness_sample %>% select(-fitness),
                                                umap_result[[i]]$umap_data,
                                                n = 100)
}
endTime <- as.numeric(Sys.time())
(endTime - startTime) / 60


# Check which one is best
best_i <- which.max(unlist(lapply(umap_result, function(x) {
  return(x$distances$r2)
})))
best_i

unlist(lapply(umap_result, function(x) {
  return(x$distances$r2)
}))

# With enough samples, there is hardly any difference between parameters - so choose sensible ones
best_i <- 1
umap_result <- vector(mode = "list", length = 1)
set.seed(42)
d_ruggedness_sample <- d_ruggedness %>%
  slice_sample(n = 100000)

umap_result[[best_i]]$umap_data <- as_tibble(umap2(d_ruggedness_sample %>% select(-fitness), 
                                              n_neighbors = 5,
                                              min_dist = 0.10,
                                              metric = "euclidean")) %>% 
  rename(UMAP1 = V1,UMAP2 = V2)
umap_result[[best_i]]$distances <- CalcDistancesUMAP(d_ruggedness_sample %>% select(-fitness),
                                                umap_result[[best_i]]$umap_data,
                                                n = 10000)

plot(umap_result[[best_i]]$distances$dist.frame, 
     main = paste("FFBH UMAP: r^2 =", round(umap_result[[best_i]]$distances$r2, digits = 3)),
     xlab = "Distance between trait coordinates", 
     ylab = "Distance between UMAP coordinates")

# Plot landscape
d_ffbh_umap <- umap_result[[best_i]]$umap_data %>% 
  mutate(z = d_ruggedness_sample$fitness) %>%
  rename(x = UMAP1,
         y = UMAP2)

# Create model to predict grid of values
# ffbh_mod <- loess(z ~ x + y, data = d_ffbh_umap)
# d_ffbh_umap <- with(d_ffbh_umap,
#                     expand.grid(x = seq(min(x), max(x), len = 100),
#                                 y = seq(min(y), max(y), len = 100)))
# d_ffbh_umap$z <- c(predict(ffbh_mod, newdata = d_ffbh_umap, se = F))


d_ffbh_umap <- akima::interp(d_ffbh_umap$x, 
                             d_ffbh_umap$y, 
                             d_ffbh_umap$z, duplicate = "strip",
                             nx = 100, ny = 100)
d_ffbh_umap <- akima::interp2xyz(d_ffbh_umap, data.frame = T) %>%
  drop_na() %>%
  mutate(z = if_else(z < 0, 0, z))

ggplot(d_ffbh_umap,
       aes(x = x, y = y, z = z)) +
  geom_raster(aes(fill = z)) +
  #scale_fill_viridis_c() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "UMAP1", y = "UMAP2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))


downsample <- 1/10
d_ffbh_umap <- umap_result[[best_i]]$umap_data %>% mutate(fitness = d_ruggedness_sample$fitness) %>%
  group_by(x = downsample * round(UMAP1 / downsample),
           y = downsample * round(UMAP2 / downsample)) %>%
  summarise(z = mean(fitness))


ggplot(d_ffbh_umap,
       aes(x = x, y = y, fill = z, z = z, group = z)) +
  geom_raster(interpolate = F) +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "UMAP1", y = "UMAP2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))

# Try with FFLC1
