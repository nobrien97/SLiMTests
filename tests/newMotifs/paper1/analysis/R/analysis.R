library(tidyverse)
library(paletteer)
library(latex2exp)
library(brms)

# Helper functions
source("helperFn.R")


# 1) Network topology shapes M

## Eigenstructure of M for each motif
# Load M
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"

d_m <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex", "meanH", 
                                          paste0("phenomean_", 1:4),
                                          paste0("phenovar_", 1:4),
                                          paste0("phenocor_", c(12, 13, 14, 23, 24, 34)),
                                          paste0("molTrait_", 1:11)))



## Neutral trait correlations

DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"
d_qg <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex", "meanH", 
                                          paste0("phenomean_", 1:4),
                                  paste0("phenovar_", 1:4),
                                  paste0("phenocor_", c(12, 13, 14, 23, 24, 34)),
                                  paste0("molTrait_", 1:11)))

d_qg_traitCor <- d_qg %>%
    mutate(model = ModelFromIndex(modelindex)) %>%
    pivot_longer(cols = starts_with("phenocor"),
                 names_to = "traitCombo",
                 values_to = "cor", names_prefix = "phenocor_") %>%
    filter(!is.infinite(cor), !is.nan(cor))

# Fisher transformation (https://en.wikipedia.org/wiki/Fisher_transformation)
d_qg_traitCor_mdl <- d_qg_traitCor %>%
  select(gen, seed, modelindex, model, traitCombo, cor) %>%
  mutate(z = atanh(pmin(pmax(cor, -1 + 1e-8), 1 - 1e-8)),
         se = 1 / sqrt(5000 - 3)) %>%
  filter(gen == 25000)

hist(d_qg_traitCor_mdl$z, breaks = 100)

fit <- brm(
  z ~ model * traitCombo + (1 | seed),
  data = d_qg_traitCor_mdl,
  family = student(),
  chains = 4, cores = 4, iter = 4000
)

summary(fit)

z_post <- posterior_epred(fit)
r_post <- tanh(z_post)

# Summarise
r_long <- as.data.frame(t(r_post)) %>%
  setNames(paste0("draw_", seq_len(ncol(.)))) %>%
  bind_cols(d_qg_traitCor_mdl) %>%
  pivot_longer(starts_with("draw_"), names_to = "draw", values_to = "r_post") %>%
  group_by(gen, modelindex, model, traitCombo) %>%
  summarise(r_post_mean = mean(r_post),
            r_post_ci_lower = r_post_mean - CI(r_post),
            r_post_ci_upper = r_post_mean + CI(r_post)) %>%
  filter(!(model == "NAR" & grepl("[3-4]", traitCombo)), # remove invalid trait combinations
         !(model == "PAR" & grepl("[3-4]", traitCombo)),
         !(model == "FFLC1" & grepl("[4]", traitCombo)),
         !(model == "FFLI1" & grepl("[4]", traitCombo)))
  

# Setup correlation matrices
make_matrix <- function(x) {
  # triangular number to get length
  t <- nrow(x) 
  n <- ((-1 + sqrt(1 + 8 * t)) / 2) + 1
  
  cor_mat <- array(1, dim = c(n, n, 3))
  
  # lower CI
  cor_mat[,,1][lower.tri(cor_mat[,,1])] <- x$r_post_ci_lower
  cor_mat[,,1][upper.tri(cor_mat[,,1])] <- t(cor_mat[,,1])[upper.tri(cor_mat[,,1])]
  # mean
  cor_mat[,,2][lower.tri(cor_mat[,,2])] <- x$r_post_mean
  cor_mat[,,2][upper.tri(cor_mat[,,2])] <- t(cor_mat[,,2])[upper.tri(cor_mat[,,2])]
  # upper CI
  cor_mat[,,3][lower.tri(cor_mat[,,3])] <- x$r_post_ci_upper
  cor_mat[,,3][upper.tri(cor_mat[,,3])] <- t(cor_mat[,,3])[upper.tri(cor_mat[,,3])]
  
  return(cor_mat)
}

cor_mats <- r_long %>%
  group_by(model) %>%
  group_map(~ make_matrix(.x))


C_FFBH_post <- array(1, dim = c(4, 4, 3))
# lower CI
C_FFBH_post[,,1][lower.tri(C_FFBH_post[,,1])] <- r_long[r_long$model == "FFBH",]$r_post_ci_lower
C_FFBH_post[,,1][upper.tri(C_FFBH_post[,,1])] <- t(C_FFBH_post[,,1])[upper.tri(C_FFBH_post[,,1])]
# mean
C_FFBH_post[,,2][lower.tri(C_FFBH_post[,,2])] <- r_long[r_long$model == "FFBH",]$r_post_mean
C_FFBH_post[,,2][upper.tri(C_FFBH_post[,,2])] <- t(C_FFBH_post[,,2])[upper.tri(C_FFBH_post[,,2])]
# upper CI
C_FFBH_post[,,3][lower.tri(C_FFBH_post[,,3])] <- r_long[r_long$model == "FFBH",]$r_post_ci_upper
C_FFBH_post[,,3][upper.tri(C_FFBH_post[,,3])] <- t(C_FFBH_post[,,3])[upper.tri(C_FFBH_post[,,3])]




# Convert to matrix form
make_cor_mat <- function(x) {
  n <- length(x)
  
  matrix(x, nrow = n)
}

cor_mats <- lapply()



d_qg_traitCor_sum <- d_qg_traitCor %>%
  group_by(gen, model, traitCombo) %>%
  summarise(meanCor = mean(cor),
            var = var(cor))

d_qg_traitCor_sum %>%
  filter(!(model == "NAR" & grepl("[3-4]", traitCombo)),
         !(model == "PAR" & grepl("[3-4]", traitCombo)),
         !(model == "FFLC1" & grepl("[4]", traitCombo)),
         !(model == "FFLI1" & grepl("[4]", traitCombo))) %>%
  select(-var) -> d_traitCor_tab
d_traitCor_tab

# Get 

## Eigenvalue dispersion - do more complex motifs produce more anisotropic M?

## Leading eigenvector directions


# 2) Alignment with M predicts evolvability

## Evolvability against alignment of M with direction of selection

### Calculate alignment of M and beta

### Get evolvability along beta bTGb


# 3) Positive and negative controls confirm developmental bias

## Confirm populations evolving toward optima parallel to M adapt faster than random-direction populations (i.e. above)

## Populations orthogonal to M adapt most slowly and contain the least evolvability

## Difference between parallel and orthogonal is the cost of developmental constraint, compare among motifs


# 4) G and M under drift and selection

## Angle between leading eigenvectors of M and G

## expect G = M under drift


# 5) Evolvability, autonomy and V_A through M

## FFBH had high V_A but low adaptation - was VA concentrated along M's leading eigenvectors and misaligned with selection?

## Conditional evolvability in the direction of selection predicts adaptive success



