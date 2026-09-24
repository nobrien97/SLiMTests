library(scattermore)

# Load in relevant data
HELPER_PATH <- "~/tests/newMotifs/paper1/ruggedness/umap/R/"
source(paste0(HELPER_PATH, "helperFns.R"))

# Read landscape data
DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"
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


d_landscape <- vector(mode = "list", length = length(model_names_noquote))
d_landscape_ds <- vector(mode = "list", length = length(model_names_noquote))

d_dist <- vector(mode = "list", length = length(model_names_noquote))

for (i in seq_along(model_names_noquote)) {
  model <- model_names_noquote[i]
  # MBA interpolated landscape
  d_landscape[[i]] <- readRDS(paste0(DATA_PATH, "d_landscape_", model, ".RDS"))
  d_landscape[[i]]$model <- model
  
  # Construct downsampled landscape
  pc_codings <- readRDS(paste0(DATA_PATH, "d_pca_codings", model, ".RDS"))
  ae_codings <- readRDS(paste0(DATA_PATH, "d_ae_codings", model, ".RDS"))
  umap_codings <- as_tibble(load_uwot(paste0(DATA_PATH, "umap_", model))$embedding)
  d_dr <- readRDS(paste0(DATA_PATH, "d_dr.dist_", model, ".RDS"))
  
  d_landscape_all <- rbind(pc_codings %>% mutate(dr.method = "PCA"),
                           ae_codings %>% mutate(dr.method = "Autoencoder"),
                           umap_codings %>% rename(LV1 = V1, LV2 = V2) %>%
                             mutate(dr.method = "UMAP"))
  
  d_dist[[i]] <- d_dr %>% mutate(model = model)
  
  best_method <- (d_dr %>%
                    group_by(dr.method) %>%
                    filter(row_number() == 1) %>% 
                    ungroup() %>% filter(r2 == max(r2)))$dr.method[1]
  
  fitnessvals <- d_ruggedness %>%
    ungroup() %>%
    filter(model == model_names_noquote[i])
  fitnessvals <- fitnessvals$fitness
  
  downsample <- 1/50 
  d_landscape_ds[[i]] <- d_landscape_all %>% 
    filter(dr.method == best_method) %>%
    select(-dr.method) %>%
    mutate(fitness = fitnessvals) %>%
    group_by(x = downsample * round(LV1 / downsample),
             y = downsample * round(LV2 / downsample)) %>%
    summarise(z = mean(fitness))
  
  d_landscape_ds[[i]]$model <- model
}

d_landscape <- data.table::rbindlist(d_landscape)
d_landscape_ds <- data.table::rbindlist(d_landscape_ds)
d_dist <- data.table::rbindlist(d_dist)

saveRDS(d_landscape, paste0(DATA_PATH, "d_landscape_combined.RDS"))
saveRDS(d_landscape_ds, paste0(DATA_PATH, "d_landscape_ds_combined.RDS"))
saveRDS(d_dist, paste0(DATA_PATH, "d_dist_combined.RDS"))

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

h2o.init(ip = "10.6.22.10",
         port = 12345, 
         startH2O = FALSE)

d_walks <- vector(mode = "list", length = length(model_names_noquote))
# Transform those walks into autoencoder space
for (i in seq_along(model_names_noquote)) {
  model_name <- model_names_noquote[i]
  ae_codings <- readRDS(paste0(DATA_PATH, "d_ae_codings", model_name, ".RDS"))
  
  n_comps <- length(molComp_names[[model_name]])
  
  # Select molecular components
  d_qg_ae <- d_qg_sample %>% filter(model == model_name) %>%
    select(gen, seed, model, isAdapted, paste0("meanMC", 1:n_comps)) %>%
    mutate(across(starts_with("meanMC"), log))
  
  # Rename to match molcomp names
  colnames(d_qg_ae)[5:ncol(d_qg_ae)] <- molComp_names[[model_name]]
  
  # Load AE model
  model_path <- list.files(paste0(DATA_PATH, "ae_", model_name, ".h2o"), full.names = T)[1]
  
  ae <- h2o.loadModel(model_path)
  
  newpoints <- as.h2o(d_qg_ae %>% ungroup() %>% select(5:ncol(d_qg_ae)))
  walk_predictions <- h2o.deepfeatures(ae, newpoints, layer = 2)
  #walk_predictions <- h2o.predict(ae, newpoints)
  
  d_qg_ae[,c("LV1", "LV2")] <- as_tibble(walk_predictions)
  
  # Add missing molcomps and rearrange columns to match others
  extra_cols <- names(all_molcomp_features)[!(names(all_molcomp_features) %in% colnames(d_qg_ae))]
  d_qg_ae[,extra_cols] <- 0.0
  d_qg_ae <- d_qg_ae %>%
    relocate(all_of(names(all_molcomp_features)))
  
  d_walks[[i]] <- d_qg_ae
}

d_walks <- data.table::rbindlist(d_walks)

saveRDS(d_walks, "d_walks.RDS")

# Plot as arrows on the landscape
d_walks <- readRDS("./d_walks.RDS")

# Plot

d_landscape <- readRDS("./d_landscape_combined.RDS")
d_landscape_ds <- readRDS("./d_landscape_ds_combined.RDS")
d_dist <- readRDS("./d_dist_combined.RDS")

d_r2 <- d_dist %>%
  group_by(model, dr.method) %>%
  summarise(r2 = r2[1])

design <- "
AABBCC
#DDEE#
"

d_walks <- d_walks %>%
  mutate(model = factor(model, levels = model_names_noquote))

ggplot(d_landscape %>%
         mutate(model = factor(model, levels = model_names_noquote)),
       aes(x = x, y = y, z = z)) +
  geom_raster(aes(fill = z)) +
  facet_manual(model ~ ., design = design) +
  scale_fill_gradientn(colours = contour_pal,
                       breaks = c(0, seq(0.1, 1.0, by = 0.3)),
                       limits = c(0, 1)) +
  scale_colour_manual(values = c("#990000", "#009900")) +
  # Plot some adaptive walks
  geom_segment(data = d_walks %>% filter(gen > 49500 & gen < 51000), inherit.aes = F,
            mapping = aes(x = LV1, xend = lead(LV1),
                          y = LV2, yend = lead(LV2),
                          colour = isAdapted, group = seed), 
            arrow = arrow(length = unit(0.3, "cm"), type = "closed"), linewidth = 0.5) +
  coord_equal() +
  labs(x = "LV1", y = "LV2", 
       fill = "Fitness", colour = "Population adapted?") +
  theme_bw() +
  theme(text = element_text(size=12), 
        legend.position = "bottom",
        legend.key.width = unit(3.5, 'line')) -> plt_landscape
plt_landscape
ggsave(paste0("plt_landscape.png"), 
       plt_landscape, device = png, width = 8.5, height = 7,
       dpi = 600)

ggplot(d_landscape_ds %>%
         mutate(model = factor(model, levels = model_names_noquote)),
       aes(x = x, y = y, fill = z, z = z, group = z)) +
  geom_raster(interpolate = F) +
  facet_manual(model ~ ., design = design) +
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
plot_landscape_ds
ggsave(paste0("plt_landscape_ds.png"), 
       plot_landscape_ds, device = png, width = 7, height = 7,
       dpi = 600)


# Distances
d_dist <- d_dist %>%
  mutate(model = factor(model, levels = model_names_noquote))

d_r2 <- d_r2 %>%
  mutate(model = factor(model, levels = model_names_noquote))

ggplot(d_dist,
       aes(x = trait_dist, y = lv_dist, colour = model)) +
  facet_nested("Model" + model ~ "Method" + dr.method) +
  geom_scattermore(shape = 21, alpha = 0.2, pixels = c(512, 512)) +
  theme_bw() +
  geom_text(data = d_r2,
            mapping = aes(x = 4, y = 4.25,
                          label = TeX(paste0("$r^2 = ",
                                             round(r2, digits = 3)),
                                      output = "character")),
            parse = T, inherit.aes = F) +
  scale_colour_manual(values = model_pal,
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names_noquote) +
  labs(x = "Component-space distance",
       y = "Latent-space distance",
       colour = "Model") +
  theme(text = element_text(size = 14),
        legend.position = "none") -> plt_dist
plt_dist
ggsave(paste0("plt_dist.png"), 
       plt_dist, device = png, width = 8, height = 8,
       dpi = 600)

