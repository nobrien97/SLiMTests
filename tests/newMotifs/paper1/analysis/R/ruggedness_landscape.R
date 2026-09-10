# Plot dimensionality-reduced fitness landscapes 
library(tidyverse)
library(h2o)
library(paletteer)
library(latex2exp)
library(lattice)
library(latticeExtra)
library(cowplot)
library(scattermore)
library(ggh4x)
library(FNN)

source("helperFn.R")

# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 1271371565
seed <- 1271371565

DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"

d_ruggedness <- data.table::fread(paste0(DATA_PATH, "d_ruggedness_permolcomp.csv"), 
                                  header = F)

colnames(d_ruggedness) <- c("step", "model", "dataset", "fitness", "startW", 
                            "endW", "netChangeW", "sumChangeW", "numFitnessHoles", 
                            "nSteps", "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                            "aZ", "bZ", "Hilln", "XMult", "base",
                            "molComp", "bkg")

# Filter to only the columns we care about
d_ruggedness <- d_ruggedness %>%
  select(2:4, 11:22)

# Try to export a subset of ruggedness for local analysis
d_ruggedness_sbst <- d_ruggedness %>%
  filter(fitness >= 0) %>% # Remove invalid estimates
  group_by(model, dataset) %>% # 5 * 3 combinations
  slice_sample(n = 100000)

saveRDS(d_ruggedness_sbst, paste0(DATA_PATH, "d_ruggedness_sbst.RDS"))

d_ruggedness_sbst <- readRDS("/mnt/d/SLiMTests/tests/newMotifs/paper1/ruggedness/d_ruggedness_sbst.RDS")

d_ruggedness_sbst <- d_ruggedness_sbst %>%
  rename(gX = XMult,
         h = Hilln,
         zZ = base)

h2o.init()

d_ruggedness_nar <- d_ruggedness_sbst %>% ungroup() %>% filter(model == "NAR") %>%
  select(3, molComp_names[["NAR"]], 1:2)

features_nar <- as.h2o(d_ruggedness_nar)
ae_nar <- h2o.deeplearning(x = 2:8,
                 training_frame = features_nar,
                 autoencoder = T,
                 seed = seed,
                 hidden = c(7, 2, 7),
                 epochs = 100,
                 activation = "Tanh",
                 sparse = F)
ae_nar
ae_nar_codings <- h2o.deepfeatures(ae_nar, features_nar, layer = 2)

d_ruggedness_par <- d_ruggedness_sbst %>% ungroup() %>% filter(model == "PAR") %>%
  select(3, molComp_names[["PAR"]], 1:2)

features_par <- as.h2o(d_ruggedness_par)
ae_par <- h2o.deeplearning(x = 2:8,
                           training_frame = features_par,
                           autoencoder = T,
                           seed = seed,
                           hidden = c(7, 2, 7),
                           epochs = 100,
                           activation = "Tanh",
                           sparse = F)
ae_par
ae_par_codings <- h2o.deepfeatures(ae_par, features_par, layer = 2)

d_ruggedness_fflc1 <- d_ruggedness_sbst %>% ungroup() %>% filter(model == "FFLC1") %>%
  select(3, molComp_names[["FFLC1"]], 1:2)

features_fflc1 <- as.h2o(d_ruggedness_fflc1)
ae_fflc1 <- h2o.deeplearning(x = 2:10,
                           training_frame = features_fflc1,
                           autoencoder = T,
                           seed = seed,
                           hidden = c(9, 6, 2, 6, 9),
                           epochs = 100,
                           activation = "Tanh",
                           sparse = F)
ae_fflc1
ae_fflc1_codings <- h2o.deepfeatures(ae_fflc1, features_fflc1, layer = 3)

d_ruggedness_ffli1 <- d_ruggedness_sbst %>% ungroup() %>% filter(model == "FFLI1") %>%
  select(3, molComp_names[["FFLI1"]], 1:2)

features_ffli1 <- as.h2o(d_ruggedness_ffli1)
ae_ffli1 <- h2o.deeplearning(x = 2:10,
                             training_frame = features_ffli1,
                             autoencoder = T,
                             seed = seed,
                             hidden = c(9, 6, 2, 6, 9),
                             epochs = 100,
                             activation = "Tanh",
                             sparse = F)
ae_ffli1
ae_ffli1_codings <- h2o.deepfeatures(ae_ffli1, features_ffli1, layer = 3)

d_ruggedness_ffbh <- d_ruggedness_sbst %>% ungroup() %>% filter(model == "FFBH") %>%
  select(3, molComp_names[["FFBH"]], 1:2)

features_ffbh <- as.h2o(d_ruggedness_ffbh)
ae_ffbh <- h2o.deeplearning(x = 2:12,
                             training_frame = features_ffbh,
                             autoencoder = T,
                             seed = seed,
                             hidden = c(11, 6, 2, 6, 11),
                             epochs = 100,
                             activation = "Tanh",
                             sparse = F)
ae_ffbh
ae_ffbh_codings <- h2o.deepfeatures(ae_ffbh, features_ffbh, layer = 3)

ae_codings <- list(ae_nar_codings,
                   ae_par_codings,
                   ae_fflc1_codings,
                   ae_ffli1_codings,
                   ae_ffbh_codings)

list_ruggedness <- list(d_ruggedness_nar,
                        d_ruggedness_par,
                        d_ruggedness_fflc1,
                        d_ruggedness_ffli1,
                        d_ruggedness_ffbh)

d_codings <- purrr::map(seq_along(list_ruggedness), function(i) {
  codings <- as.data.frame(ae_codings[[i]])
  list_ruggedness[[i]] %>%
    mutate(DF1 = codings[,1],
           DF2 = codings[,2],
           dataset = factor(list_ruggedness[[i]]$dataset, 
                            levels = c("Parallel",
                                       "Orthogonal",
                                       "Randomised")),
           model = factor(list_ruggedness[[i]]$model,
                          levels = rev(model_names_noquote)))
  }, .progress = T)

d_codings <- data.table::rbindlist(d_codings, fill = T)

saveRDS(d_codings, "/mnt/i/SLiMTests/tests/newMotifs/paper1/ruggedness/d_codings.RDS")

d_codings <- readRDS("/mnt/i/SLiMTests/tests/newMotifs/paper1/ruggedness/d_codings.RDS")
d_codings <- d_codings %>% mutate(model = factor(model, levels = model_names_noquote))


design <- "
AABBCC
#DDEE#
"

# ggplot scattermore
plt_codings <- ggplot(d_codings,
       aes(x = DF1, y = DF2, colour = fitness)) +
  facet_manual(model ~ ., design) +
  geom_scattermore(pointsize = 3.2, 
                   interpolate = T, alpha = 1, pixels = c(1024, 1024)) +
  scale_colour_viridis_c(breaks = seq(0, 1, by = 0.25),
                       labels = seq(0, 1, by = 0.25),
                       limits = c(0.0, 1.0)) + 
  labs(x = TeX("$Z_1$"), y = TeX("$Z_2$"), colour = "Fitness") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black"),#paletteer_d("ggprism::viridis"), 1),
        legend.key.width = unit(3.5, "lines"))
plt_codings
ggsave("plt_landscape_autoencoder.png", width = 9, height = 7, device = png)

# Lattice contour plot
d_codings <- d_codings %>%
  group_by(model, dataset) %>%
  slice_sample(n = 1000)

# Plot
grid_codings <- levelplot(fitness ~ DF1*DF2 | dataset*model, d_codings, 
                            col.regions = paletteer_c("viridis::viridis", 100, 1),
                            #panel = panel.levelplot.points, 
                            xlab = list(label = "DF1", cex = 1.2), 
                            ylab = list(label = "DF2", cex = 1.2), 
                          as.table = F,
                            colorkey = list(space = "bottom",
                                            title = "Fitness",
                                            labels = list(cex = 1.2)),
                            #par.settings = list(layout.heights = list(xlab.key.padding = 4)), 
                            pretty = T) + 
  layer_(panel.2dsmoother(..., n = 200))

contours_codings <- 
  contourplot(
    fitness ~ DF1*DF2 | dataset*model, d_codings, 
    panel=panel.2dsmoother, col = "white",
    as.table = F,
    labels = list(col = "white",
                  cex = 1.2))

plt_codings <- grid_codings + contours_codings
plot_grid(plt_codings)
ggsave("plt_codings.png", device = png, width = 7, height = 7.4, bg = "white")



# PCA
pca_nar <- h2o.prcomp(
  x = 2:8,
  training_frame = features_nar,
  pca_method = "GramSVD",
  transform = "NONE",
  impute_missing = T,
  k = 7,
  seed = seed
)
pca_nar

pca_par <- h2o.prcomp(
  x = 2:8,
  training_frame = features_par,
  pca_method = "GramSVD",
  transform = "NONE",
  impute_missing = T,
  k = 7,
  seed = seed
)
pca_par

pca_fflc1 <- h2o.prcomp(
  x = 2:10,
  training_frame = features_fflc1,
  pca_method = "GramSVD",
  transform = "NONE",
  impute_missing = T,
  k = 9,
  seed = seed
)
pca_fflc1

pca_ffli1 <- h2o.prcomp(
  x = 2:10,
  training_frame = features_ffli1,
  pca_method = "GramSVD",
  transform = "NONE",
  impute_missing = T,
  k = 9,
  seed = seed
)
pca_ffli1

pca_ffbh <- h2o.prcomp(
  x = 2:12,
  training_frame = features_ffbh,
  pca_method = "GramSVD",
  transform = "NONE",
  impute_missing = T,
  k = 11,
  seed = seed
)
pca_ffbh

pca_all <- list(pca_nar,
                pca_par,
                pca_fflc1,
                pca_ffli1,
                pca_ffbh)

d_pca_codings <- purrr::map(seq_along(pca_all), function(i){
  pca_result <- pca_all[[i]]
  ruggedness_model <- as.h2o(list_ruggedness[[i]])
  
  codings <- h2o.predict(pca_result, ruggedness_model)
  
  codings <- codings %>%
    as.data.frame() %>%
    select(PC1, PC2)
    
  return(cbind(codings, list_ruggedness[[i]]))
}, .progress = T)

d_pca_codings <- data.table::rbindlist(d_pca_codings, fill = T)

d_pca_codings <- d_pca_codings %>%
  mutate(model = factor(model, levels = model_names_noquote),
         dataset = factor(dataset, levels = c("Parallel",
                                              "Orthogonal",
                                              "Randomised")))

saveRDS(d_pca_codings, "/mnt/i/SLiMTests/tests/newMotifs/paper1/ruggedness/d_pca_codings.RDS")
d_pca_codings <- readRDS("/mnt/d/SLiMTests/tests/newMotifs/paper1/ruggedness/d_pca_codings.RDS")

plt_pca <- ggplot(d_pca_codings,
                      aes(x = PC1, y = PC2, colour = fitness)) +
  #facet_nested("Model" + model ~ "Trait/selection alignment" + dataset) +
  facet_manual(model ~ ., design) +
  geom_scattermore(pointsize = 3.2, 
                   interpolate = T, alpha = 1, pixels = c(1024, 1024)) +
  scale_colour_viridis_c(breaks = seq(0, 1, by = 0.25),
                         labels = seq(0, 1, by = 0.25),
                         limits = c(0.0, 1.0)) + 
  labs(x = TeX("$PC1$"), y = TeX("$PC2$"), colour = "Fitness") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black"),#paletteer_d("ggprism::viridis"), 1),
        legend.key.width = unit(3.5, "lines"))
plt_pca
ggsave("plt_landscape_pca.png", width = 9, height = 7, device = png)

# Plot as a contour map - downsample
downsample <- 1/20
d_pca_grid <- d_pca_codings %>% select(model, dataset, PC1, PC2, fitness) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample), 
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(fitness))


contour_pal <- c("#220022", paletteer_d("ggprism::viridis", 6)[-c(1)])

plt_pca_tile <- ggplot(d_pca_grid %>% drop_na(),
                     aes(x = x, y = y, fill = tanh(z), z = z, group = z)) +
  facet_manual(model ~ ., design = design) +
  geom_tile() +
  # scale_fill_gradientn(colours = contour_pal,
  #                   breaks = c(0, seq(0.1, 1.0, by = 0.3)),
  #                   #labels = seq(0, 1, by = 0.25),
  #                   limits = c(0, 1)) + 
  scale_fill_gradientn(colours = contour_pal) +
  #scale_fill_viridis_c() +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
plt_pca_tile
ggsave("plt_pca_tile.png", plt_pca_tile, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")


# Interpolated landscape
# Fast interpolation via k nearest neighbours
get.nn <- function(data, labels, query) {
  nns <- get.knnx(data, query, k=1)
  labels[nns$nn.index]
  cbind(query, z = labels[nns$nn.index])
}

GRID_DIM <- 1000
d_pca_grid <- d_pca_codings %>% select(model, dataset, PC1, PC2, fitness) %>%
  group_by(model) %>%
  mutate(maxPC1 = max(PC1),
         minPC1 = min(PC1),
         maxPC2 = max(PC1),
         minPC2 = min(PC2)) %>%
  group_map(~ get.nn(.x %>% select(x = PC1, y = PC2), 
                   .x$fitness, expand.grid(data.frame(x = seq(.x$minPC1[1], .x$maxPC1[1], length.out = GRID_DIM),
                                       y = seq(.x$minPC2[1], .x$maxPC2[1], length.out = GRID_DIM)))))
  
d_pca_grid <- purrr::map(seq_along(d_pca_grid), function(i) {
    x <- d_pca_grid[[i]]
    x$model <- model_names_noquote[i]
    return(x)
    })

d_pca_grid <- data.table::rbindlist(d_pca_grid)  


plt_pca_contour <- ggplot(d_pca_grid,
                          aes(x = x, y = y, z = z)) +
  facet_manual(model ~ ., design = design) +
  geom_contour_filled() +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
plt_pca_contour


##########################################
### HPC
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



# Try PCA on whole dataset
PCPerModel <- function(data, z, columns, should.scale = T) {
  data <- data %>% select(all_of(columns))
  pc <- prcomp(data, scale = should.scale)
  result <- data.frame(PC1 = pc$x[,1],
                       PC2 = pc$x[,2],
                       z = z)
  return(result)
}

pc_ruggedness <- purrr::map(seq_along(model_names_noquote), function(i) {
  x <- d_ruggedness %>%
    filter(fitness >= 0.0) %>%
    rename(h = Hilln,
           gX = XMult,
           zZ = base) %>%
    mutate(model = factor(model, levels = model_names_noquote)) %>%
    filter(model == model_names_noquote[i])
  
  z <- x$fitness
  column_names <- (molComp_names[model_names_noquote[i]])[[1]]
  result <- PCPerModel(x, z, column_names)
  result$model <- model_names_noquote[i]
  return(result)
  
}, .progress = T)

d_pc_ruggedness <- data.table::rbindlist(pc_ruggedness, fill = T)

downsample <- 1/15
d_pc_ruggedness_ds <- d_pc_ruggedness %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

nrow(d_pc_ruggedness_ds)

# Plot
contour_pal <- c("#220022", paletteer_d("ggprism::viridis", 6)[-c(1)])
design <- "
AABBCC
#DDEE#
"

plt_pca_tile <- ggplot(d_pc_ruggedness_ds %>% drop_na() %>%
                         mutate(model = factor(model, 
                                               levels = model_names_noquote)),
                       aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design) +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                    breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                    limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_full_log10.png", plt_pca_tile, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")


# Do a zoom-in PCA on a smaller part of the space
# focus on area of maximum fitness

# ds_zoom should be lower resolution than the final downsample, 
# we want to identify the area of largest average fitness
ds_zoom <- 1/4
zoom_coords <- d_pc_ruggedness %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = ds_zoom * round(PC1 / ds_zoom),
           y = ds_zoom * round(PC2 / ds_zoom)) %>%
  summarise(z = mean(z)) %>%
  ungroup() %>%
  group_by(model) %>%
  # Identify which areas have greatest fitness and filter around there
  summarise(bigZ = max(z),
        bigZCoord = which(z == bigZ),
         bigPC1 = x[bigZCoord],
         bigPC2 = y[bigZCoord])

zoom_coords
  
# Now filter to area with greatest fitness
d_pc_ruggedness_zoom_ds <- left_join(d_pc_ruggedness,
                                     zoom_coords,
                                     by = "model") %>% 
  select(model, PC1, PC2, bigPC1, bigPC2, z) %>%
  group_by(model) %>%
  filter(between(PC1, bigPC1 - 1, bigPC1 + 1),
         between(PC2, bigPC2 - 1, bigPC2 + 1)) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))
  
plt_pca_tile_zoom <- ggplot(d_pc_ruggedness_zoom_ds %>% drop_na() %>%
                         mutate(model = factor(model, 
                                               levels = model_names_noquote)),
                       aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design, scales = "free") +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_zoom_log10.png", plt_pca_tile_zoom, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")

# This was for the sampling of trait combinations from 0 to log(10)
# Maybe not realistic for most cases, what about around the starting conditions?
# Those are sampled log(2), we look at 0 to log(3)
d_ruggedness_log3 <- data.table::fread(paste0(DATA_PATH, "log3/d_ruggedness_permolcomp.csv"), 
                                  header = F)

colnames(d_ruggedness_log3) <- c("step", "model", "dataset", "fitness", "startW", 
                            "endW", "netChangeW", "sumChangeW", "numFitnessHoles", 
                            "nSteps", "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                            "aZ", "bZ", "Hilln", "XMult", "base",
                            "molComp", "bkg")

# Filter to only the columns we care about
d_ruggedness_log3 <- d_ruggedness_log3 %>%
  select(2:4, 11:22)

pc_ruggedness_log3 <- purrr::map(seq_along(model_names_noquote), function(i) {
  x <- d_ruggedness_log3 %>%
    filter(fitness >= 0.0) %>%
    rename(h = Hilln,
           gX = XMult,
           zZ = base) %>%
    mutate(model = factor(model, levels = model_names_noquote)) %>%
    filter(model == model_names_noquote[i])
  
  z <- x$fitness
  column_names <- (molComp_names[model_names_noquote[i]])[[1]]
  result <- PCPerModel(x, z, column_names)
  result$model <- model_names_noquote[i]
  return(result)
  
}, .progress = T)

d_pc_ruggedness_log3 <- data.table::rbindlist(pc_ruggedness_log3, fill = T)

downsample <- 1/15
d_pc_ruggedness_log3_ds <- d_pc_ruggedness_log3 %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

nrow(d_pc_ruggedness_log3_ds)

# Plot
plt_pca_log3_tile <- ggplot(d_pc_ruggedness_log3_ds %>% drop_na() %>%
                         mutate(model = factor(model, 
                                               levels = model_names_noquote)),
                       aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design) +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_full_log3.png", plt_pca_log3_tile, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")

# ds_zoom should be lower resolution than the final downsample, 
# we want to identify the area of largest average fitness
ds_zoom <- 1/4
zoom_coords_log3 <- d_pc_ruggedness_log3 %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = ds_zoom * round(PC1 / ds_zoom),
           y = ds_zoom * round(PC2 / ds_zoom)) %>%
  summarise(z = mean(z)) %>%
  ungroup() %>%
  group_by(model) %>%
  # Identify which areas have greatest fitness and filter around there
  summarise(bigZ = max(z),
            bigZCoord = which(z == bigZ),
            bigPC1 = x[bigZCoord],
            bigPC2 = y[bigZCoord])

zoom_coords_log3

# Now filter to area with greatest fitness
d_pc_ruggedness_log3_zoom_ds <- left_join(d_pc_ruggedness_log3,
                                     zoom_coords_log3,
                                     by = "model") %>% 
  select(model, PC1, PC2, bigPC1, bigPC2, z) %>%
  group_by(model) %>%
  filter(between(PC1, bigPC1 - 1, bigPC1 + 1),
         between(PC2, bigPC2 - 1, bigPC2 + 1)) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

plt_pca_log3_tile_zoom <- ggplot(d_pc_ruggedness_log3_zoom_ds %>% drop_na() %>%
                              mutate(model = factor(model, 
                                                    levels = model_names_noquote)),
                            aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design, scales = "free") +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_zoom_log3.png", plt_pca_log3_tile_zoom, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")




#####
# Repeat but with no PC scaling (since the components are all on the same scale anyway)
pc_ruggedness <- purrr::map(seq_along(model_names_noquote), function(i) {
  x <- d_ruggedness %>%
    filter(fitness >= 0.0) %>%
    rename(h = Hilln,
           gX = XMult,
           zZ = base) %>%
    mutate(model = factor(model, levels = model_names_noquote)) %>%
    filter(model == model_names_noquote[i])
  
  z <- x$fitness
  column_names <- (molComp_names[model_names_noquote[i]])[[1]]
  result <- PCPerModel(x, z, column_names, should.scale = F)
  result$model <- model_names_noquote[i]
  return(result)
  
}, .progress = T)

d_pc_ruggedness <- data.table::rbindlist(pc_ruggedness, fill = T)

d_pc_ruggedness_ds <- d_pc_ruggedness %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

# Plot
plt_pca_tile <- ggplot(d_pc_ruggedness_ds %>% drop_na() %>%
                         mutate(model = factor(model, 
                                               levels = model_names_noquote)),
                       aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design) +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_full_log10_noscale.png", plt_pca_tile, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")


# Do a zoom-in PCA on a smaller part of the space
# focus on area of maximum fitness

# ds_zoom should be lower resolution than the final downsample, 
# we want to identify the area of largest average fitness
zoom_coords <- d_pc_ruggedness %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = ds_zoom * round(PC1 / ds_zoom),
           y = ds_zoom * round(PC2 / ds_zoom)) %>%
  summarise(z = mean(z)) %>%
  ungroup() %>%
  group_by(model) %>%
  # Identify which areas have greatest fitness and filter around there
  summarise(bigZ = max(z),
            bigZCoord = which(z == bigZ),
            bigPC1 = x[bigZCoord],
            bigPC2 = y[bigZCoord])

# Now filter to area with greatest fitness
d_pc_ruggedness_zoom_ds <- left_join(d_pc_ruggedness,
                                     zoom_coords,
                                     by = "model") %>% 
  select(model, PC1, PC2, bigPC1, bigPC2, z) %>%
  group_by(model) %>%
  filter(between(PC1, bigPC1 - 1, bigPC1 + 1),
         between(PC2, bigPC2 - 1, bigPC2 + 1)) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

plt_pca_tile_zoom <- ggplot(d_pc_ruggedness_zoom_ds %>% drop_na() %>%
                              mutate(model = factor(model, 
                                                    levels = model_names_noquote)),
                            aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design, scales = "free") +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_zoom_log10_noscale.png", plt_pca_tile_zoom, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")


pc_ruggedness_log3 <- purrr::map(seq_along(model_names_noquote), function(i) {
  x <- d_ruggedness_log3 %>%
    filter(fitness >= 0.0) %>%
    rename(h = Hilln,
           gX = XMult,
           zZ = base) %>%
    mutate(model = factor(model, levels = model_names_noquote)) %>%
    filter(model == model_names_noquote[i])
  
  z <- x$fitness
  column_names <- (molComp_names[model_names_noquote[i]])[[1]]
  result <- PCPerModel(x, z, column_names, should.scale = F)
  result$model <- model_names_noquote[i]
  return(result)
  
}, .progress = T)

d_pc_ruggedness_log3 <- data.table::rbindlist(pc_ruggedness_log3, fill = T)

downsample <- 1/15
d_pc_ruggedness_log3_ds <- d_pc_ruggedness_log3 %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

nrow(d_pc_ruggedness_log3_ds)

# Plot
plt_pca_log3_tile <- ggplot(d_pc_ruggedness_log3_ds %>% drop_na() %>%
                              mutate(model = factor(model, 
                                                    levels = model_names_noquote)),
                            aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design) +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_full_log3_noscale.png", plt_pca_log3_tile, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")

# ds_zoom should be lower resolution than the final downsample, 
# we want to identify the area of largest average fitness
ds_zoom <- 1/4
zoom_coords_log3 <- d_pc_ruggedness_log3 %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = ds_zoom * round(PC1 / ds_zoom),
           y = ds_zoom * round(PC2 / ds_zoom)) %>%
  summarise(z = mean(z)) %>%
  ungroup() %>%
  group_by(model) %>%
  # Identify which areas have greatest fitness and filter around there
  summarise(bigZ = max(z),
            bigZCoord = which(z == bigZ),
            bigPC1 = x[bigZCoord],
            bigPC2 = y[bigZCoord])

zoom_coords_log3

# Now filter to area with greatest fitness
d_pc_ruggedness_log3_zoom_ds <- left_join(d_pc_ruggedness_log3,
                                          zoom_coords_log3,
                                          by = "model") %>% 
  select(model, PC1, PC2, bigPC1, bigPC2, z) %>%
  group_by(model) %>%
  filter(between(PC1, bigPC1 - 1, bigPC1 + 1),
         between(PC2, bigPC2 - 1, bigPC2 + 1)) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

plt_pca_log3_tile_zoom <- ggplot(d_pc_ruggedness_log3_zoom_ds %>% drop_na() %>%
                                   mutate(model = factor(model, 
                                                         levels = model_names_noquote)),
                                 aes(x = x, y = y, fill = z, z = z, group = z)) +
  facet_manual(model ~ ., design = design, scales = "free") +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "PC1", y = "PC2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
ggsave("plt_pca_tile_zoom_log3_noscale.png", plt_pca_log3_tile_zoom, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")




##########################################