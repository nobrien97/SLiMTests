library(tidyverse)
library(paletteer)
library(latex2exp)
library(brms)
library(betareg)

# Helper functions
source("helperFn.R")


# 1) Network topology shapes M

## Eigenstructure of M for each motif
# Load M
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"
#DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"

d_m <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                          paste0("mean_", 1:4),
                                          paste0("var_", 1:4),
                                          paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m <- d_m %>%
  mutate(model = ModelFromIndexWithR(modelindex))
  
# Convert a row to a matrix
row_to_m <- function(x) {
  # Get the number of traits and covariance terms
  n <- 2
  cov_terms <- 12
  
  if (x$model == "FFLC1" | x$model == "FFLI1") {
    n <- 3
    cov_terms <- c(12, 13, 23)
  }
  
  if (x$model == "FFBH") {
    n <- 4
    cov_terms <- c(12, 13, 14, 23, 24, 34)
  }
  
  # Triangular number for number of covariance terms
  n_cov <- ((n-1) * n) / 2
  
  m <- matrix(NA_real_, nrow = n, ncol = n)
  
  # Variances
  diag(m) <- unlist(x[1,paste0("var_", 1:n)])
  
  # Covariances
  m[lower.tri(m)] <- unlist(x[1,paste0("cov_", cov_terms)])
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  
  return(m)
}

# get matrices
m_matrices <- d_m %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))

# Get eigenvectors of each M
e_m <- lapply(m_matrices, eigen)
saveRDS(e_m, "eigen_randomised_m.RDS")

# e_m <- readRDS("eigen_randomised_m.RDS")

# Calculate relative eigenvalue dispersion
Vrel <- function(l) {
  p <- length(l)
  avg_l <- mean(l)
  
  sum((l - avg_l)^2) / (p * (p-1) * avg_l^2)
}

vrel_m <- unlist(lapply(e_m, function(x) { Vrel(x$values) }))

# Add to data
d_vrel <- d_m %>%
  mutate(r = RFromIndex(modelindex),
         vrel = vrel_m) %>%
  select(gen, seed, modelindex, model, r, vrel)

# Plot
d_vrel_sum <- d_vrel %>%
  group_by(gen, model) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))

ggplot(d_vrel_sum,
       aes(x = gen - 50000, y = vrel_mean, colour = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = vrel_mean - vrel_CI, ymax = vrel_mean + vrel_CI, fill = model),
              colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1)) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# Stable over time, can take average?
d_vrel_tab <- d_vrel %>%
  group_by(model) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))
d_vrel_tab

lm_vrel <- lm(vrel ~ model, data = d_vrel)
summary(lm_vrel)

# beta distributed
hist(d_vrel$vrel, breaks = 100)

br_vrel <- betareg(vrel ~ model, data = d_vrel %>% filter(gen == 60000))

summary(br_vrel)

# SLOW
# beta_vrel <- brm(
#   bf(vrel ~ model),
#   data = d_vrel,
#   family = Beta(),
#   chains = 4, iter = 2000, 
#   cores = 4,
#   file = "beta_vrel"
# )
# 
# summary(beta_vrel)
# 
# vrel_post <- posterior_epred(beta_vrel)


## Neutral trait correlations

DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"
#DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"

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

hist(d_qg_traitCor_mdl$cor, breaks = 100)
hist(d_qg_traitCor_mdl$z, breaks = 100)

fit <- brm(
  z ~ model * traitCombo + (1 | seed),
  data = d_qg_traitCor_mdl,
  family = gaussian(),
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

# Store as correlation matrices - alphabetical order
cor_mats <- r_long %>%
  mutate(model = factor(model, levels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"))) %>%
  group_by(model) %>%
  group_map(~ make_matrix(.x))

# label
names(cor_mats) <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")




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



