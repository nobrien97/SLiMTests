# Split data for pythae
library(tidyverse)

# Load data, split into groups
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


d_ruggedness <- readRDS("/mnt/i/SLiMTests/tests/newMotifs/paper1/ruggedness/d_ruggedness_sbst.RDS")
d_ruggedness <- d_ruggedness %>%
  rename(zZ = base,
         gX = XMult,
         h = Hilln)

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


# Use K-medoids to sample n centroids for fitting riemannian geo model
N_CENTROIDS <- 256
N_SAMPLES <- 1

d_ruggedness_nar <- d_ruggedness %>% filter(model == "NAR") %>% ungroup() %>%
  select(fitness, molComp_names[["NAR"]])
d_ruggedness_par <- d_ruggedness %>% filter(model == "PAR") %>% ungroup() %>%
  select(fitness, molComp_names[["PAR"]]) 
d_ruggedness_fflc1 <- d_ruggedness %>% filter(model == "FFLC1") %>% ungroup() %>%
  select(fitness, molComp_names[["FFLC1"]]) 
d_ruggedness_ffli1 <- d_ruggedness %>% filter(model == "FFLI1") %>% ungroup() %>%
  select(fitness, molComp_names[["FFLI1"]]) %>%
    cluster::clara(., N_CENTROIDS, samples = N_SAMPLES, pamLike = T))$medoids
d_ruggedness_ffbh <- d_ruggedness %>% filter(model == "FFBH") %>% ungroup() %>%
  select(fitness, molComp_names[["FFBH"]]) 


d_nar_centroids <- cluster::clara(d_ruggedness_nar, N_CENTROIDS, samples = N_SAMPLES, pamLike = T,
                 metric = "euclidean")$medoids
d_par_centroids <- cluster::clara(d_ruggedness_par, N_CENTROIDS, samples = N_SAMPLES, pamLike = T,
                                  metric = "euclidean")$medoids
d_fflc1_centroids <- cluster::clara(d_ruggedness_fflc1, N_CENTROIDS, samples = N_SAMPLES, pamLike = T,
                                  metric = "euclidean")$medoids
d_ffli1_centroids <- cluster::clara(d_ruggedness_ffli1, N_CENTROIDS, samples = N_SAMPLES, pamLike = T,
                                  metric = "euclidean")$medoids
d_ffbh_centroids <- cluster::clara(d_ruggedness_ffbh, N_CENTROIDS, samples = N_SAMPLES, pamLike = T,
                                  metric = "euclidean")$medoids

# Centroids
write_csv(as.data.frame(d_nar_centroids), "d_centroids_NAR.csv")
write_csv(as.data.frame(d_par_centroids), "d_centroids_PAR.csv")
write_csv(as.data.frame(d_fflc1_centroids), "d_centroids_FFLC1.csv")
write_csv(as.data.frame(d_ffli1_centroids), "d_centroids_FFLI1.csv")
write_csv(as.data.frame(d_ffbh_centroids), "d_centroids_FFBH.csv")



write_csv(d_ruggedness_nar, "d_ruggedness_NAR.csv")
write_csv(d_ruggedness_par, "d_ruggedness_PAR.csv")
write_csv(d_ruggedness_fflc1, "d_ruggedness_FFLC1.csv")
write_csv(d_ruggedness_ffli1, "d_ruggedness_FFLI1.csv")
write_csv(d_ruggedness_ffbh, "d_ruggedness_FFBH.csv")


##
# See if model properly keeps scale/distance between points relative to PCA and phenotype space
##
d_ruggedness_rh <- read_csv("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/rhvae/py/d_ruggedness_NAR_rh.csv")

# Attach dataset
d_ruggedness_rh$dataset <- d_ruggedness[d_ruggedness$model == "NAR",]$dataset

# Plot RH
downsample <- 1/60
d_ruggedness_nar_rh_ds <- d_ruggedness_rh %>% select(model, RH1, RH2, fitness) %>%
  group_by(model,
           x = downsample * round(RH1 / downsample),
           y = downsample * round(RH2 / downsample)) %>%
  summarise(z = mean(fitness))

library(ggh4x)
library(paletteer)
contour_pal <- c("#220022", paletteer_d("ggprism::viridis", 6)[-c(1)])
design <- "
AABBCC
#DDEE#
"


# Plot
plt_rhvae_test <- ggplot(d_ruggedness_nar_rh_ds %>% drop_na() %>%
                              mutate(model = factor(model, 
                                                    levels = model_names_noquote)),
                            aes(x = x, y = y, fill = z, z = z, group = z)) +
  #facet_manual(model ~ ., design = design) +
  geom_raster() +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  labs(x = "RH1", y = "RH2", 
       fill = "Fitness") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line'))
plt_rhvae_test

# distance matrix of samples
SAMPLE_SIZE <- 10000
d_distances <- data.frame(trait_dist = numeric(SAMPLE_SIZE),
                          latent_dist = numeric(SAMPLE_SIZE))

set.seed(42)
samples_i <- sample(sample_space, SAMPLE_SIZE * 2, replace = F)
samples_j <- samples_i[(SAMPLE_SIZE+1):(SAMPLE_SIZE * 2)]
samples_i <- samples_i[1:SAMPLE_SIZE]

# Transform to same space
d_ruggedness_trait_scale <- d_ruggedness_rh %>%
  select(2:8) %>%
  mutate(across(everything(), scale))

d_ruggedness_rh_scale <- d_ruggedness_rh %>%
  select(9:10) %>%
  mutate(across(everything(), scale))


d_distances$latent_dist <- sqrt(rowSums((d_ruggedness_rh_scale[samples_i, ] - d_ruggedness_rh_scale[samples_j, ])^2))
d_distances$trait_dist <- sqrt(rowSums((d_ruggedness_trait_scale[samples_i, ] - d_ruggedness_trait_scale[samples_j, ])^2))
  
# Calculate correlation
dist_r2 <- cor(d_distances$trait_dist, d_distances$latent_dist, method = "pearson")^2

ggplot(d_distances,
       aes(x = trait_dist, y = latent_dist)) +
  geom_point(shape = 21, alpha = 0.3) +
  ggtitle(paste("RHVAE | R^2 =", dist_r2)) +
  theme_bw() +
  labs(x = "Phenotype space distance", y = "Latent space distance")



# compare to PCA
z <- d_ruggedness_nar$fitness
column_names <- (molComp_names[["NAR"]])
pc_ruggedness_log3 <- PCPerModel(d_ruggedness_nar, z, column_names, should.scale = F)

# Scale so that trait space comparable
pc_ruggedness_log3 <- pc_ruggedness_log3 %>%
  mutate(across(1:2, scale))

d_distances$pc_dist <- sqrt(rowSums((pc_ruggedness_log3[samples_i, ] - pc_ruggedness_log3[samples_j, ])^2))
dist_pc_r2 <- cor(d_distances$trait_dist, d_distances$pc_dist, method = "pearson")^2

ggplot(d_distances,
       aes(x = trait_dist, y = pc_dist)) +
  geom_point(shape = 21, alpha = 0.3) +
  ggtitle(paste("PCA | R^2 =", dist_pc_r2)) +
  theme_bw() +
  labs(x = "Phenotype space distance", y = "Latent space distance")

# Quite a good fit, NAR landscape quite simple


# What about for the rugged FFLC1 or the holey FFBH? Does PCA still hold up?
z <- d_ruggedness_fflc1$fitness
column_names <- (molComp_names[["FFLC1"]])
pc_ruggedness_log3 <- PCPerModel(d_ruggedness_fflc1, z, column_names, should.scale = F)

pc_ruggedness_log3 <- pc_ruggedness_log3 %>%
  mutate(across(1:2, scale))

d_distances$pc_dist <- sqrt(rowSums((pc_ruggedness_log3[samples_i, ] - pc_ruggedness_log3[samples_j, ])^2))
dist_pc_r2 <- cor(d_distances$trait_dist, d_distances$pc_dist, method = "pearson")^2

ggplot(d_distances,
       aes(x = trait_dist, y = pc_dist)) +
  geom_point(shape = 21, alpha = 0.3) +
  ggtitle(paste("PCA | R^2 =", dist_pc_r2)) +
  theme_bw() +
  labs(x = "Phenotype space distance", y = "Latent space distance")

# HAHA NO

# FFBH model
z <- d_ruggedness_ffbh$fitness
column_names <- (molComp_names[["FFBH"]])
pc_ruggedness_log3 <- PCPerModel(d_ruggedness_ffbh, z, column_names, should.scale = F)

pc_ruggedness_log3 <- pc_ruggedness_log3 %>%
  mutate(across(1:2, scale))

d_distances$pc_dist <- sqrt(rowSums((pc_ruggedness_log3[samples_i, ] - pc_ruggedness_log3[samples_j, ])^2))
dist_pc_r2 <- cor(d_distances$trait_dist, d_distances$pc_dist, method = "pearson")^2

ggplot(d_distances,
       aes(x = trait_dist, y = pc_dist)) +
  geom_point(shape = 21, alpha = 0.3) +
  ggtitle(paste("PCA | R^2 =", dist_pc_r2)) +
  theme_bw() +
  labs(x = "Phenotype space distance", y = "Latent space distance")

# ALSO BAD
