library(tidyverse)
library(paletteer)
library(latex2exp)
library(ggh4x)

DATA_PATH <- "/path/to/ruggedness/"

d_ruggedness <- data.table::fread(paste0(DATA_PATH, "d_ruggedness_permolcomp.csv"), 
                                  header = F)

colnames(d_ruggedness) <- c("step", "model", "dataset", "fitness", "startW", 
                            "endW", "netChangeW", "sumChangeW", "numFitnessHoles", 
                            "nSteps", "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                            "aZ", "bZ", "Hilln", "XMult", "base",
                            "molComp", "bkg")

# Filter to only the component values and model info
d_ruggedness <- d_ruggedness %>%
  select(2:4, 11:22)


# Define molecular components for each motif
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

# Motif names
model_names_noquote <- c("NAR", "PAR", "FFLC1", 
                         "FFLI1", "FFBH")



# Try PCA on whole dataset
PCPerModel <- function(data, z, columns) {
  data <- data %>% select(all_of(columns))
  pc <- prcomp(data, scale = T)
  result <- data.frame(PC1 = pc$x[,1],
                       PC2 = pc$x[,2],
                       z = z)
  return(result)
}

# Do PCA per motif
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

# Downsample by averaging across neighbouring PC1/PC2
downsample <- 1/15
d_pc_ruggedness_ds <- d_pc_ruggedness %>% select(model, PC1, PC2, z) %>%
  group_by(model,
           x = downsample * round(PC1 / downsample),
           y = downsample * round(PC2 / downsample)) %>%
  summarise(z = mean(z))

# Plot setup
contour_pal <- c("#220022", paletteer_d("ggprism::viridis", 6)[-c(1)])
design <- "
AABBCC
#DDEE#
"

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
ggsave("plt_pca_tile_full.png", plt_pca_tile, device = png, 
       width = 11, height = 8, dpi = 600, bg = "white")
