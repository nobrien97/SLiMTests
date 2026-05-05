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


# Measure cosine similarity between M matrices and beta

d_opt <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)
# o = optimum, s = sigma, d = direction (-1, 1)
colnames(d_opt) <- c("seed", "modelindex", "o_t1", "o_t2", "o_t3", "o_t4", 
                     "s_t1", "s_t2", "s_t3", "s_t4", "d_t1", "d_t2", "d_t3",
                     "d_t4")



# Attach quant gen data
d_qg <- data.table::fread(paste0("/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/", 
                                 "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 29)), 
                          col.names = c("gen", "seed", "modelindex", "meanH",
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var", 
                                        "trait3_var", "trait4_var", "dist", 
                                        "dist1", "dist2", "dist3", "dist4", "mean_w",
                                        "var_w", "deltaPheno", "deltaW", 
                                        "meanMC1", "meanMC2", "meanMC3", "meanMC4", 
                                        "meanMC5", "meanMC6", "meanMC7", "meanMC8", 
                                        "meanMC9", "meanMC10", "meanMC11"), 
                          fill = T)

d_selvec_m <- d_qg %>%
  filter(gen >= 50000) %>%
  select(gen, seed, modelindex, isAdapted, ends_with("mean"))

d_selvec_m <- left_join(d_selvec_m, d_opt %>% 
                        select(seed, modelindex, starts_with("o_")) %>%
                        mutate(seed = factor(seed),
                               modelindex = factor(modelindex)), 
                      by = c("seed", "modelindex"))

d_selvec_m <- AddCombosToDF(d_selvec_m)

d_selvec_m <- d_selvec_m %>%
  mutate(modelindex = as.factor(modelindex),
         seed = as.factor(seed),
         model = factor(model, levels = model_names)) %>%
  rename(timePoint = gen)


id_m <- d_m %>% mutate(timePoint = gen) %>% select(timePoint, seed, modelindex)
id_m$clus <- 1
id_m$modelindex <- as.factor(id_m$modelindex)
id_m$seed <- as.factor(id_m$seed)

id_m <- AddCombosToDF(id_m)
id_m$model <- factor(id_m$model, levels = model_names)

id_m <- inner_join(id_m, d_qg %>% mutate(timePoint = gen) %>% 
                     select(timePoint, seed, modelindex, isAdapted),
                   by = c("timePoint", "seed", "modelindex"))

d_selvec_m <- inner_join(id_m, d_selvec_m, 
                         by = c("timePoint", "seed", "modelindex", "isAdapted", "model", "r"))


d_selvec_m <- d_selvec_m %>%
  rowwise() %>%
  mutate(t1_dir = o_t1 - trait1_mean,
         t2_dir = o_t2 - trait2_mean,
         t3_dir = o_t3 - trait3_mean,
         t4_dir = o_t4 - trait4_mean,
         norm = sqrt(sum(c_across(ends_with("dir"))^2)), # normalise
         t1_dir = t1_dir / norm,
         t2_dir = t2_dir / norm,
         t3_dir = t3_dir / norm,
         t4_dir = t4_dir / norm) %>%
  select(timePoint, seed, modelindex, isAdapted, model, r, norm, ends_with("dir"))


d_cossim_m <- GetCosineSimilarity(m_matrices, d_selvec_m %>% select(ends_with("dir")), id_m)

saveRDS(d_cossim_m, "d_cossim_m.RDS")

d_cossim_m <- AddCombosToDF(d_cossim_m)

d_cossim_m_sum <- d_cossim_m %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_m_sum$model <- as.factor(d_cossim_m_sum$model)


ggplot(d_cossim_m_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $m_{max}$ and $\\beta"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim_m_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Time point", 
       y = TeX("Evolvability ($\\beta^T M \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


## Neutral trait correlations

DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"
d_qg_traitCor <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex", "meanH", 
                                          paste0("phenomean_", 1:4),
                                  paste0("phenovar_", 1:4),
                                  paste0("phenocor_", c(12, 13, 14, 23, 24, 34)),
                                  paste0("molTrait_", 1:11)))

d_qg_traitCor <- d_qg_traitCor %>%
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

# Measure alignment of trait correlation matrices against G and M correlations
# G
h2_cor <- lapply(h2_pd, cov2cor)

# Add correlation matrix to id: index in cor_mats
id <- id %>%
  mutate(corMatIndex = match(model, model_names))

# Match correlation matrices with G matrix list
# First eigenvector of trait correlation matrix
cor_eig <- lapply(cor_mats, function(x) return(eigen(x[,,2])$vectors[,1]))
cor_eig <- lapply(cor_eig, function(x) {
  result <- numeric(4)
  result[1:length(x)] <- x
  return(result)
})

d_trait_cor <- as.data.frame(t(as.data.frame(cor_eig[id$corMatIndex])))

# Measure cosine similarity
d_cossim_R <- GetCosineSimilarity(h2_cor, d_trait_cor, id)

d_cossim_R <- AddCombosToDF(d_cossim_R)

d_cossim_R_sum <- d_cossim_R %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_R_sum$model <- as.factor(d_cossim_R_sum$model)


ggplot(d_cossim_R, 
       aes(x = timePoint, y = sqrt(cosSim^2), colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_R_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanCosSim, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $g_{max}$ and $r_{max}$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim_R, 
       aes(x = timePoint, y = bTMb, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_R_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanbTGb, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Evolvability ($\\beta^T G \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# M 
m_cor <- vector(mode = "list", length(m_matrices))
for (i in seq_along(m_matrices)) {
  result <- cov2cor(m_matrices[[i]])
  result[is.nan(result)] <- 0
  m_cor[[i]] <- result
}

# Add correlation matrix to id: index in cor_mats
id_m <- id_m %>%
  mutate(corMatIndex = match(model, model_names))

# Match correlation matrices with G matrix list
# First eigenvector of trait correlation matrix
cor_m_eig <- lapply(m_cor, function(x) return(eigen(x)$vectors[,1]))
cor_m_eig <- lapply(cor_m_eig, function(x) {
  result <- numeric(4)
  result[1:length(x)] <- x
  return(result)
})

d_m_cor <- as.data.frame(t(as.data.frame(cor_m_eig[id_m$corMatIndex])))

# Measure cosine similarity
d_cossim_m_traitcor <- GetCosineSimilarity(m_cor, d_m_cor, id_m)

saveRDS(d_cossim_m_traitcor, "d_cossim_m_traitcor.RDS")

d_cossim_m_traitcor <- AddCombosToDF(d_cossim_m_traitcor)

d_cossim_m_traitcor_sum <- d_cossim_m_traitcor %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_m_traitcor_sum$model <- as.factor(d_cossim_m_traitcor_sum$model)

ggplot(d_cossim_m_traitcor_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $m_{max}$ and $r_{max}$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim_m_traitcor_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Time point", 
       y = TeX("Evolvability ($r_{max}^T M r_{max}$)"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


# Contingency table: population adapted vs cosine similarity > threshold
## model type is collinear with this



## Eigenvalue dispersion - do more complex motifs produce more anisotropic M?

## Leading eigenvector directions


# 2) Alignment with M predicts evolvability

## Evolvability against alignment of M with direction of selection

# G matrix evolvability analysis



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



