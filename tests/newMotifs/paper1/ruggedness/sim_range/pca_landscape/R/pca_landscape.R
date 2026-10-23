# Use UMAP to plot fitness landscapes
# Helper functions, libraries etc.
HELPER_PATH <- "~/tests/newMotifs/paper1/ruggedness/umap/R/"
HELPER_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/umap/R/"

source(paste0(HELPER_PATH, "helperFns.R"))

clamp <- function(x, y = c(0,1), na.rm = F) {
  if (na.rm) {
    x[!is.na(x)]
  }
  
  pmax(pmin(x, y[2]), y[1])
}



DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/"

setwd(paste0(DATA_PATH, "sim_range"))
d_ruggedness <- data.table::fread(paste0(DATA_PATH, "sim_range/d_ruggedness_sim.csv"), 
                                  header = F,
                                  colClasses = c("character", "character", "character", "logical", 
                                                 rep("numeric", times = 13)),
                                  col.names = c("seed", "model", "dataset", "isAdapted", "fitness",
                                                "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                                                "aZ", "bZ", "gX", "h", "zZ"))

d_ruggedness <- d_ruggedness %>%
  mutate(seed = factor(seed),
         model = factor(model, levels = model_names_noquote),
         dataset = factor(dataset, levels = c("Parallel",
                                              "Orthogonal",
                                              "Randomised")),
         isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"),
                            labels = c("Adapted", "Maladapted")))

# Filter to only the columns we care about
d_ruggedness <- d_ruggedness %>%
  select(2:4, 5:17)

d_ruggedness_all <- d_ruggedness

d_ruggedness_all %>%
  group_by(model) %>%
  summarise(across(aX:zZ, min))

d_ruggedness_all %>%
  group_by(model) %>%
  summarise(across(aX:zZ, max))

landscape_list <- vector(mode = "list", length = 5)

# Construct landscape for each example
d_pc <- d_ruggedness %>%
  group_by(model, dataset, isAdapted) %>%
  group_map(~ prcomp(.x %>% select(all_of(molComp_names[[as.character(.y$model)]])),
            scale = F))

trait_values <- as.matrix(d_ruggedness %>%
                            select(-fitness))

centered <- scale(trait_values, center = prc$center, scale = F)
pc_values <- as.data.frame(centered %*% prc$rotation[,1:2])


d_rug_pc <- d_ruggedness %>%
  group_by(model, dataset, isAdapted) %>%
  group_map(~ {
    prc <- prcomp(.x %>% select(all_of(molComp_names[[as.character(.y$model)]])),
                 scale = F)
    trait_values <- .x %>% select(all_of(molComp_names[[as.character(.y$model)]]))
    centered <- scale(trait_values, center = prc$center, scale = F)
    pc_values <- as.data.frame(centered %*% prc$rotation[,1:2])
    colnames(pc_values) <- c("PC1", "PC2")
    return(cbind(.y, seed = .x$seed, fitness = .x$fitness, pc_values))
  })

d_rug_pc <- data.table::rbindlist(d_rug_pc)
saveRDS(d_rug_pc, paste0(DATA_PATH, "sim_range/d_rug_pc_sim.RDS"))
d_rug_pc <- readRDS(paste0(DATA_PATH, "sim_range/d_rug_pc_sim.RDS"))


d_landscape <- d_rug_pc %>%
  rename(x = PC1, y = PC2, z = fitness) %>%
  group_by(model, dataset) %>% # TODO: Could add isAdapted in here
  group_map(~ {
    landscape_in <- .x %>%
      select(x, y, z)
    mba_landscape <- mba.surf(landscape_in, 
                              no.X = 100, no.Y = 100)
    result <- expand.grid(x = mba_landscape$xyz.est$x,
                          y = mba_landscape$xyz.est$y)
    result$z <- as.vector(mba_landscape$xyz.est$z)
    result <- result %>%
      mutate(z = if_else(is.na(z), 0, clamp(z)), 
             model = .y$model,
             dataset = .y$dataset,
             #isAdapted = .y$isAdapted,
             seed = .x$seed[1])
  })

d_landscape <- data.table::rbindlist(d_landscape)
saveRDS(d_landscape, paste0(DATA_PATH, "sim_range/d_landscape_pca_sim.RDS"))



