library(tidyverse)
library(uwot)
library(ggh4x)
library(paletteer)
library(latex2exp)
library(h2o)
library(MBA)

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


contour_pal <- c("#220022", paletteer_d("ggprism::viridis", 6)[-c(1)])
design <- "
AABBCC
#DDEE#
"


CalcDistancesUMAP <- function(trait_data, umap_data, n = 10000, seed = 42, should.scale = T) {
  # Trait data and umap data should correspond to each other (i.e. same row)
  if (nrow(trait_data) != nrow(umap_data)) {
    stop("UMAP rows should correspond to trait rows!")
  }
  
  if (n > nrow(trait_data)) {
    warning("n is greater than number of rows, clamping to number of rows")
    n <- nrow(trait_data)
  }
  
  sample_space <- 1:nrow(trait_data)
  result <- data.frame(trait_dist = numeric(n),
                       umap_dist = numeric(n))
  
  set.seed(seed)
  # Randomly sample pairs of rows
  samples_i <- sample(sample_space, n * 2, replace = F)
  samples_j <- samples_i[(n+1):(n * 2)]
  samples_i <- samples_i[1:n]
  
  if (should.scale) {
    # Transform both the umap and traits to same space
    trait_data <- trait_data %>%
      mutate(across(everything(), scale))
    
    umap_data <- umap_data %>%
      mutate(across(everything(), scale))
  }
  
  # Euclidean distances between each of i and j
  result$umap_dist <- sqrt(rowSums((umap_data[samples_i, ] - umap_data[samples_j, ])^2))
  result$trait_dist <- sqrt(rowSums((trait_data[samples_i, ] - trait_data[samples_j, ])^2))
  
  # Calculate correlation
  dist_r2 <- cor(result$trait_dist, result$umap_dist, method = "pearson")^2
  
  return(list("dist.frame" = result,
              "r2" = dist_r2))
}
