library(tidyverse)
library(paletteer)
library(latex2exp)
library(brms)
library(betareg)
library(ggh4x)
library(ggbeeswarm)
library(cowplot)
library(nlme)
library(emmeans)

# Helper functions
source("helperFn.R")

# Combos
COMBO_PATH <- '/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
COMBO_PATH <- '/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
d_combos <- read_delim(COMBO_PATH, 
                       delim = " ", col_names = F)
names(d_combos) <- c("model", "r")




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


# Measure cosine similarity between M matrices and beta

d_opt <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)
d_opt <- read_csv("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)

# o = optimum, s = sigma, d = direction (-1, 1)
colnames(d_opt) <- c("seed", "modelindex", "o_t1", "o_t2", "o_t3", "o_t4", 
                     "s_t1", "s_t2", "s_t3", "s_t4", "d_t1", "d_t2", "d_t3",
                     "d_t4")



# Attach quant gen data
PATH_QG <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"

d_qg <- data.table::fread(PATH_QG, header = F, 
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

# Summarise adapted/maladapted at end of sim
d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98)) %>%
  mutate(model = factor(model, levels = model_names)) %>%
  ungroup() -> d_qg



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
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
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
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


## Neutral trait correlations

DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"
#DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"

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


# Read G matrices
G_DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/R/"
G_DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/R/"

d_h2_mrr <- read_csv(paste0(G_DATA_PATH, "out_h2_mrr.csv"), col_names = F)
d_h2_mkr <- read_csv(paste0(G_DATA_PATH, "out_h2_mkr.csv"), col_names = F)

d_h2_trait_mkr <- read_csv(paste0(G_DATA_PATH, "out_h2_trait_mkr.csv"), col_names = F)
d_h2_trait_mrr <- read_csv(paste0(G_DATA_PATH, "out_h2_trait_mrr.csv"), col_names = F)

colnames(d_h2_trait_mkr) <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_t1",
                              "VA_t2", "VA_t3", "VA_t4", "CVA_t1_t2", "CVA_t1_t3",
                              "CVA_t1_t4", "CVA_t2_t3", "CVA_t2_t4", "CVA_t3_t4",
                              "h2_t1", "h2_t2", "h2_t3", "h2_t4")

colnames(d_h2_trait_mrr) <- colnames(d_h2_trait_mkr)

# join
d_h2_trait_mkr$calcMode <- "mkr"
d_h2_trait_mrr$calcMode <- "mrr"

d_h2_trait <- rbind(d_h2_trait_mkr, d_h2_trait_mrr)

d_h2_trait %>% mutate(model = d_combos$model[.$modelindex],
                      model = factor(model, levels = model_names),
                      r = d_combos$r[.$modelindex]) -> d_h2_trait

d_h2_trait <- d_h2_trait %>%
  distinct(gen, seed, modelindex, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()


d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, isAdapted) %>% filter(gen >= 49500)


# inner join optPerc
d_h2_trait <- left_join(d_h2_trait, d_qg_optPerc, by = c("gen", "seed", "modelindex"))

# Counts for each model type:
table(d_h2_trait$model, d_h2_trait$isAdapted)

# Discretise generation
d_h2_trait <- d_h2_trait %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))

# summarise
d_h2_trait_sum <- d_h2_trait %>% 
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanH2w = mean(h2_w, na.rm = T),
                   seH2w = se(h2_w, na.rm = T),
                   meanVAw = mean(VA_w, na.rm = T),
                   seVAw = se(VA_w, na.rm = T))
d_h2_trait_sum$model <- as.factor(d_h2_trait_sum$model)

# Heritability distribution
ggplot(d_h2_trait %>% 
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?"),
       aes(x = timePoint, y = h2_w, colour = model)) +
  facet_nested(r_title + log10(r) ~ adapted_title + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_trait_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanH2w, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Narrow-sense heritability $(h^2)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Additive variance
# Small effects as separate figure
ggplot(d_h2_trait %>%
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?"),
       aes(x = timePoint, y = VA_w, colour = model)) +
  facet_nested(r_title + log10(r) ~ adapted_title + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_trait_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanVAw, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Additive variance in fitness $(VA)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names,
                      guide = guide_legend(override.aes = list(shape = 16,
                                                               size = 3))) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Split h2 into G matrices
d_h2_trait %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(5:8, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  group_by(modelindex, timePoint, isAdapted) %>%
  group_split(.) -> split_h2


# Separate into model indices
# each sublist is replicates of a model index
sourceCpp("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")
sourceCpp("/mnt/e/Documents/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")

lapply(split_h2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices


# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat <- unlist(cov_matrices, recursive = F)

# get ids from the matrix
cov_matrix_modelindex <- GetMatrixIDs(split_h2)

id <- data.table::rbindlist(cov_matrix_modelindex, 
                fill = T)
id$label <- as.character(1:nrow(id))
id$modelindex <- as.factor(id$modelindex)
id <- AddCombosToDF(id)
id$model <- factor(id$model, levels = model_names)

# Evolvability metrics

# First convert to nearest positive definite matrix
h2_pd <- lapply(h2_mat, function(x) {
  if (!is.positive.definite(x)) {return (as.matrix(nearPD(x)$mat))}
  return(x)
})

# Now find cosine similarity between selection vector and leading eigenvector of G
# Filter selvec to h2_pd matrices
d_selvec <- d_qg %>%
  filter(gen == 50000 | gen == 60000) %>%
  select(gen, seed, modelindex, ends_with("mean"))

d_selvec <- left_join(d_selvec, d_opt %>% 
                        select(seed, modelindex, starts_with("o_")) %>%
                        mutate(seed = factor(seed),
                               modelindex = factor(modelindex)), 
                      by = c("seed", "modelindex"))

d_selvec <- AddCombosToDF(d_selvec)

d_selvec <- d_selvec %>%
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
  select(gen, seed, modelindex, model, r, norm, ends_with("dir"))


d_selvec <- d_selvec %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
  select(-gen)

d_selvec2 <- inner_join(id, d_selvec, 
                        by = c("timePoint", "seed", "modelindex", "model", "r"))

d_cossim <- GetCosineSimilarity(h2_pd, d_selvec2 %>% select(ends_with("dir")), id)

d_cossim <- AddCombosToDF(d_cossim)

d_cossim_sum <- d_cossim %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_sum$model <- as.factor(d_cossim_sum$model)


ggplot(d_cossim, 
       aes(x = timePoint, y = sqrt(cosSim^2), colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanCosSim, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $g_{max}$ and $\\beta"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim, 
       aes(x = timePoint, y = bTMb, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_sum %>% ungroup() %>%
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
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")



d_ecr <- CalcECRATrait(h2_pd, id)
d_ecr <- AddCombosToDF(d_ecr)

# Refactor model
d_ecr <- d_ecr %>%
  mutate(model = factor(model, levels = model_names))

# Need to calculate cev means separately for the different models
# K- shouldn't mean over cev_KZ and KXZ
d_ecr_sum <- d_ecr %>%
  group_by(timePoint, model, r) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = cev, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = timePoint, y = cev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Time point", y = "Mean conditional evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = res, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = timePoint, y = res_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Time point", y = "Mean respondability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = aut, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = timePoint, y = aut_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Time point", y = "Mean autonomy",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = ev, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = timePoint, y = ev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Time point", y = "Mean evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                      plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO", label_size = 12)

plt_evol <- plot_grid(plt_evol,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol.png", device = png, bg = "white",
       width = 10, height = 7)




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

# Measure cosine similarity of g_max vs r_max
d_cossim_R <- GetCosineSimilarity(h2_cor, d_trait_cor, id)

d_cossim_R <- AddCombosToDF(d_cossim_R)

d_cossim_R_sum <- d_cossim_R %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_R_sum$model <- as.factor(d_cossim_R_sum$model)


ggplot(d_cossim_R, 
       aes(x = timePoint, y = abs(cosSim), colour = model)) +
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
       y = TeX("Evolvability ($r_{max}^T G ~r_{max}$)"),
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

# Match correlation matrices with M matrix list
# First eigenvector of trait correlation matrix
d_m_cor <- as.data.frame(t(as.data.frame(cor_eig[id_m$corMatIndex])))

# Measure cosine similarity between M matrix correlations and neutral trait correlations
d_cossim_m_traitcor <- GetCosineSimilarity(m_cor, d_m_cor, id_m)

saveRDS(d_cossim_m_traitcor, "d_cossim_m_traitcor.RDS")

d_cossim_m_traitcor <- AddCombosToDF(d_cossim_m_traitcor)

d_cossim_m_traitcor_sum <- d_cossim_m_traitcor %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
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
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
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
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Adapted populations look like they have less variation along the trait correlation axis
## need to get around the constraint
### can we look at how this is affected by the alignment between selection and trait corr?
### if selection is opposing trait corr, then adaptation will require M and G to become 
### misaligned with R

corBeta_mats <- cor_mats[id_m$corMatIndex]
corBeta_mats <- lapply(corBeta_mats, function(x) return(x[,,2]))

d_cossim_corBeta <- GetCosineSimilarity(corBeta_mats, d_selvec_m %>% select(ends_with("dir")), id_m)

saveRDS(d_cossim_corBeta, "d_cossim_corBeta.RDS")

d_cossim_corBeta <- AddCombosToDF(d_cossim_corBeta)

d_cossim_corBeta_sum <- d_cossim_corBeta %>%
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_corBeta_sum$model <- as.factor(d_cossim_corBeta_sum$model)

ggplot(d_cossim_corBeta_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $r_{max}$ and $\\beta$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

#  Across all points
ggplot(d_cossim_corBeta, 
       aes(x = abs(cosSim), fill = model)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted) +
  geom_histogram(bins = 100, show.legend = F) +
  labs(
       x = TeX("Absolute cosine similarity between $r_{max}$ and $\\beta$"),
       fill = "") +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")



ggplot(d_cossim_corBeta_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Time point", 
       y = TeX("Evolvability ($\\beta^T R \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                    labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                    breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


# Is the probability that you are adapted predicted well by the cosine similarity between r_max and beta,
# or by model type alone

d_cossim_corBeta_abs <- d_cossim_corBeta %>%
  mutate(absCosSim = abs(cosSim)) %>%
  filter(timePoint == 60000)

lm_rmax_beta <- glm(isAdapted ~ model * absCosSim,
                    data = d_cossim_corBeta_abs,
                    family = "binomial")
plot(lm_rmax_beta)
summary(lm_rmax_beta)

em_rmax_beta <- emmeans(lm_rmax_beta, ~ model * absCosSim, 
                        at = list(absCosSim = c(0.1, 0.5, 0.9)),
                        type = "response")
summary(em_rmax_beta)
test(em_rmax_beta)
emmip(em_rmax_beta, model ~ absCosSim)

# Or flip the model: what is the cosine similarity of maladapted/adapted pops

gls_cossim_rmax_beta <- gls(absCosSim ~ model * isAdapted,
                           data = d_cossim_corBeta_abs,
                           weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_cossim_rmax_beta)
em_cossim_rmax_beta <- emmeans(gls_cossim_rmax_beta, ~ model * isAdapted,
                        type = "response")
summary(em_cossim_rmax_beta)
test(em_cossim_rmax_beta)
emmip(em_cossim_rmax_beta, model ~ isAdapted)


# Prediction says that cosine similarity between beta and r_max decreases in adapted pops
## except for FFLI1, FFLC1 not significant
## So the pops that adapt tend to have lower or similar alignment of selection with trait corr
## selection might be shaping the trait correlations so the neutral expectation shifts?
## In this case alignment of G and M should be more important

## What about bTRb?
gls_bTRb_rmax_beta <- gls(bTMb ~ model * isAdapted,
                            data = d_cossim_corBeta_abs,
                            weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_bTRb_rmax_beta)
em_bTRb_rmax_beta <- emmeans(gls_bTRb_rmax_beta, ~ model * isAdapted,
                               type = "response")
summary(em_bTRb_rmax_beta)
test(em_bTRb_rmax_beta)
emmip(em_bTRb_rmax_beta, model ~ isAdapted)


## GLS for Gmax
d_cossim_gmax_abs <- d_cossim %>%
  filter(timePoint == "End") %>%
  mutate(absCosSim = abs(cosSim))

gls_cossim_gmax_beta <- gls(absCosSim ~ model * isAdapted,
                            data = d_cossim_gmax_abs,
                            weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_cossim_gmax_beta)
em_cossim_gmax_beta <- emmeans(gls_cossim_gmax_beta, ~ model * isAdapted,
                               type = "response")
summary(em_cossim_gmax_beta)
test(em_cossim_gmax_beta)
emmip(em_cossim_gmax_beta, model ~ isAdapted, CIs = T)

## bTGb
gls_bTGb_gmax_beta <- gls(bTMb ~ model * isAdapted,
                            data = d_cossim_gmax_abs,
                            weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_bTGb_gmax_beta)
em_bTGb_gmax_beta <- emmeans(gls_bTGb_gmax_beta, ~ model * isAdapted,
                               type = "response")
summary(em_bTGb_gmax_beta)
test(em_bTGb_gmax_beta)
emmip(em_bTGb_gmax_beta, model ~ isAdapted, CIs = T)

# M matrix
d_cossim_mmax_abs <- d_cossim_m %>%
  filter(timePoint == "60000") %>%
  mutate(absCosSim = abs(cosSim))

gls_cossim_mmax_beta <- gls(absCosSim ~ model * isAdapted,
                            data = d_cossim_mmax_abs,
                            weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_cossim_mmax_beta)
em_cossim_mmax_beta <- emmeans(gls_cossim_mmax_beta, ~ model * isAdapted,
                               type = "response")
summary(em_cossim_mmax_beta)
test(em_cossim_mmax_beta)
emmip(em_cossim_mmax_beta, model ~ isAdapted, CIs = T)

## bTGb
gls_bTGb_mmax_beta <- gls(bTMb ~ model * isAdapted,
                          data = d_cossim_mmax_abs,
                          weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_bTGb_mmax_beta)
em_bTGb_mmax_beta <- emmeans(gls_bTGb_mmax_beta, ~ model * isAdapted,
                             type = "response")
summary(em_bTGb_mmax_beta)
test(em_bTGb_mmax_beta)
emmip(em_bTGb_mmax_beta, model ~ isAdapted, CIs = T)

# Adapted populations tended to have more mutational variance along beta


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



