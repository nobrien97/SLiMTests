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
library(matrixcalc)
library(Matrix)
library(tidymodels)
library(randomForest)
library(caret)
library(iml)
library(sobolMDA)
library(ggalt)


# Helper functions
source("helperFn.R")

# Combos
COMBO_PATH <- '/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
COMBO_PATH <- '/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
d_combos <- read_delim(COMBO_PATH, 
                       delim = " ", col_names = F)
names(d_combos) <- c("model", "r")

# Attach quant gen data
PATH_QG <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_qg.csv"

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




# 1) Network topology shapes M

## Eigenstructure of M for each motif
# Load M
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"
DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar.csv"

d_m <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                          paste0("mean_", 1:4),
                                          paste0("var_", 1:4),
                                          paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m <- d_m %>%
  mutate(model = ModelFromIndexWithR(modelindex))
  

# get matrices
# m_matrices <- d_m %>%
#   rowwise() %>%
#   group_map(~ row_to_m(.x))

#saveRDS(m_matrices, "m_matrices.RDS")
m_matrices <- readRDS("m_matrices.RDS")


# Get eigenvectors of each M
#e_m <- lapply(m_matrices, eigen)
#saveRDS(e_m, "eigen_randomised_m.RDS")

e_m <- readRDS("eigen_randomised_m.RDS")

# Calculate relative eigenvalue dispersion

vrel_m <- unlist(lapply(e_m, function(x) { Vrel(x$values) }))

# Add to data
d_vrel <- d_m %>%
  mutate(model = factor(model, levels = model_names_noquote),
         r = RFromIndex(modelindex),
         vrel = vrel_m) %>%
  select(gen, seed, modelindex, model, r, vrel)

# join with quant gen data
d_vrel <- left_join(d_m %>% mutate(seed = factor(seed),
                                   modelindex = factor(modelindex),
                                   r = RFromIndex(modelindex),
                                   vrel = vrel_m),
                     d_qg %>% select(gen, seed, modelindex, model, r, isAdapted) %>%
                       mutate(model = str_remove_all(model, "'")),
                     by = c("gen", "seed", "modelindex", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote))


# Plot
d_vrel_sum <- d_vrel %>%
  group_by(gen, model, isAdapted) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))

ggplot(d_vrel_sum,
       aes(x = gen - 50000, y = vrel_mean, colour = model)) +
  facet_nested(.~"Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = vrel_mean - vrel_CI, ymax = vrel_mean + vrel_CI, fill = model),
              colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")



# Stable over time (except FFLC1, mirrors M evolvability), can take average?
d_vrel_tab <- d_vrel %>%
  group_by(model, isAdapted) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))
d_vrel_tab

ggplot(d_vrel,
       aes(x = model, y = vrel, colour = model)) +
  facet_nested(.~"Population adapted" + isAdapted) +
  geom_boxplot() +
  geom_point(data = d_vrel_tab, aes(y = vrel_mean), size = 2, stroke = 1,
             shape = 21,
             colour = "black", fill = "white") +
  labs(x = "Model", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_colour_manual(values = pal,
                      breaks = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
ggsave("plt_vrel.png", device = png, bg = "white",
       width = 8, height = 6)

# beta distributed
hist(d_vrel$vrel, breaks = 100)

br_vrel <- betareg(vrel ~ model + isAdapted, data = d_vrel %>% filter(gen == 60000))
plot(br_vrel)
summary(br_vrel)

em.vrel <- emmeans(br_vrel, ~ model + isAdapted, type = "response")
test(em.vrel)
pwpp(em.vrel, by = "model")

emmip(em.vrel, ~ model + isAdapted, CIs = T) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y= TeX("Predicted $V_{rel}$")) +
  theme(text = element_text(size = 12))
ggsave("plt_pred_vrel.png", device = png, bg = "white",
       width = 8, height = 6)


## V_rel shows that all models are quite anisotropic, but interestingly the FFBH is
## most isotropic of them all.
## There is a consistent difference between adapted and maladapted models: adapted
## have lower V_rel (more isotropy) than maladapted.

# Does the correlation structure among traits reflect correlations among molecular 
# components?
## Can look at covariance of trait M with mol comps?
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar_percomp.csv"
DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar_percomp.csv"
DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar_percomp.csv"

t_mc_combos <- expand.grid(1:11, 1:4)

d_m_molcomp <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                         paste0("cov_", t_mc_combos$Var2, "_", t_mc_combos$Var1)))

d_m_molcomp <- d_m_molcomp %>%
  mutate(model = ModelFromIndexWithR(modelindex))

# Split cov_ wider
d_m_molcomp <- d_m_molcomp %>% filter(gen == 60000) %>%
  pivot_longer(cols = starts_with("cov"),
               names_to = c("misc", "trait", "component"),
               names_sep = "_",
               values_to = "cov") %>%
  select(-misc)

d_m_molcomp <- left_join(d_m_molcomp %>% mutate(seed = factor(seed),
                                   modelindex = factor(modelindex),
                                   r = RFromIndex(modelindex)),
                    d_qg %>% select(gen, seed, modelindex, model, r, isAdapted) %>%
                      mutate(model = str_remove_all(model, "'")),
                    by = c("gen", "seed", "modelindex", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote))



d_m_molcomp_sum <- d_m_molcomp %>%
  group_by(model, isAdapted, trait, component) %>%
  summarise(meanCov = mean(abs(cov)),
            SECov = se(abs(cov)))

# Plot covariances as heatmap
ggplot(d_m_molcomp_sum, 
       aes(x = trait, y = component, fill = meanCov)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted) +
  geom_tile() +
  theme_bw() +
  labs(x = "Trait", y = "Component", fill = "Covariance")




# Measure cosine similarity between M matrices and beta
d_opt <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)
d_opt <- read_csv("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)
d_opt <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_opt.csv", col_names = F)

# o = optimum, s = sigma, d = direction (-1, 1)
colnames(d_opt) <- c("seed", "modelindex", "o_t1", "o_t2", "o_t3", "o_t4", 
                     "s_t1", "s_t2", "s_t3", "s_t4", "d_t1", "d_t2", "d_t3",
                     "d_t4")


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
  mutate(t1_dir = o_t1 - trait1_mean,
         t2_dir = o_t2 - trait2_mean,
         t3_dir = o_t3 - trait3_mean,
         t4_dir = o_t4 - trait4_mean,
         norm = sqrt(rowSums(pick(ends_with("dir"))^2)), # normalise
         t1_dir = t1_dir / norm,
         t2_dir = t2_dir / norm,
         t3_dir = t3_dir / norm,
         t4_dir = t4_dir / norm) %>%
  select(timePoint, seed, modelindex, isAdapted, model, r, norm, ends_with("dir"))


d_cossim_m <- GetCosineSimilarity(m_matrices, d_selvec_m %>% select(ends_with("dir")), id_m)

saveRDS(d_cossim_m, "d_cossim_m.RDS")
d_cossim_m <- readRDS("d_cossim_m.RDS")

d_cossim_m <- AddCombosToDF(d_cossim_m)

d_cossim_m_sum <- d_cossim_m %>%
  mutate(timePoint = timePoint - 50000) %>%
  group_by(timePoint, model, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_m_sum$model <- factor(d_cossim_m_sum$model, levels = model_names)


ggplot(d_cossim_m_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Absolute cosine similarity between $m_{max}$ and $\\beta"),
       colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal, labels = model_names_noquote) +
  scale_fill_manual(values = pal, labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim_m_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Model" + model~ "Population adapted" + isAdapted,
               scales = "free",
               labeller = labeller(model = model_names_labeller)) +
  geom_line() +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Evolvability ($\\beta^T M \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = pal, labels = model_names_noquote) +
  scale_fill_manual(values = pal, labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


## Neutral trait correlations
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"
DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/slim_qg.csv"
DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/neutralCorr/slim_qg.csv"

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
  dplyr::summarise(r_post_mean = mean(r_post),
            r_post_ci_lower = r_post_mean - CI(r_post),
            r_post_ci_upper = r_post_mean + CI(r_post)) %>%
  filter(!(model == "NAR" & grepl("[3-4]", traitCombo)), # remove invalid trait combinations
         !(model == "PAR" & grepl("[3-4]", traitCombo)),
         !(model == "FFLC1" & grepl("[4]", traitCombo)),
         !(model == "FFLI1" & grepl("[4]", traitCombo)))

d_traitCor_sum <- d_qg_traitCor_mdl %>%
  group_by(gen, modelindex, model, traitCombo) %>%
  dplyr::summarise(r_post_mean = mean(cor),
            r_post_ci = CI(cor),
            r_post_ci_lower = r_post_mean - CI(r_post),
            r_post_ci_upper = r_post_mean + CI(r_post)) %>%
  filter(!(model == "NAR" & grepl("[3-4]", traitCombo)), # remove invalid trait combinations
         !(model == "PAR" & grepl("[3-4]", traitCombo)),
         !(model == "FFLC1" & grepl("[4]", traitCombo)),
         !(model == "FFLI1" & grepl("[4]", traitCombo)))

  
# Store as correlation matrices - alphabetical order
cor_mats <- d_traitCor_sum %>%
  dplyr::mutate(model = factor(model, levels = model_names_noquote)) %>%
  group_by(model) %>%
  group_map(~ make_matrix(.x))

# label
names(cor_mats) <- model_names_noquote


# Table for supplementary

cor_mats_ci <- lapply(cor_mats, function(x) {
  return(x[,,2] - x[,,1])
})

cor_mats_table <- vector("list", 5)

for (i in seq_along(cor_mats_table)) {
  meanMat <- round(cor_mats[[i]][,,2], 3)
  n <- nrow(meanMat)
  result <- meanMat
  diag(result) <- "1"
  cor_mats_table[[i]] <- result
}

knitr::kable(cor_mats_table, format = "latex")



# Read G matrices
G_DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/R/"
G_DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/R/"
G_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/"

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
  scale_colour_manual(values = pal, labels = model_names) +
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
  scale_colour_manual(values = pal,
                      labels = model_names) +
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
  mutate(t1_dir = o_t1 - trait1_mean,
         t2_dir = o_t2 - trait2_mean,
         t3_dir = o_t3 - trait3_mean,
         t4_dir = o_t4 - trait4_mean,
         norm = sqrt(rowSums(pick(ends_with("dir"))^2)), # normalise
         t1_dir = t1_dir / norm,
         t2_dir = t2_dir / norm,
         t3_dir = t3_dir / norm,
         t4_dir = t4_dir / norm) %>%
  select(timePoint, seed, modelindex, dataset, isAdapted, model, r, norm, ends_with("dir"))


d_selvec <- d_selvec %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
  select(-gen)

d_selvec2 <- inner_join(id, d_selvec, 
                        by = c("timePoint", "seed", "modelindex", "model", "r"))

d_cossim <- GetCosineSimilarity(h2_pd, d_selvec2 %>% select(ends_with("dir")), id)

d_cossim <- AddCombosToDF(d_cossim)

d_cossim_sum <- d_cossim %>%
  group_by(timePoint, model, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_sum$model <- as.factor(d_cossim_sum$model)


ggplot(d_cossim, 
       aes(x = timePoint, y = sqrt(cosSim^2), colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanCosSim, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $g_{max}$ and $\\beta"),
       colour = "Model") +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim, 
       aes(x = timePoint, y = bTMb, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanbTGb, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Evolvability ($\\beta^T G \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = pal,
                      labels = model_names) +
  #coord_cartesian(ylim = c(0, 1)) +
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
  group_by(model, isAdapted) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = cev, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = cev_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean conditional evolvability (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev
plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = res, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = res_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean respondability (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_res
plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = aut, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = aut_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean autonomy (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_aut
plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = ev, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = ev_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean evolvability (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_ev
plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                      plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO", label_size = 12)

plt_evol <- plot_grid(plt_evol,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol_g.png", device = png, bg = "white",
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
d_cossim_R$model <- factor(d_cossim_R$model, levels = model_names)
d_cossim_R$dataset <- factor(id$dataset, levels = c("Orthogonal", "Parallel", "Randomised"))


d_cossim_R_sum <- d_cossim_R %>%
  group_by(model, dataset, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))


ggplot(d_cossim_R, 
       aes(x = model, y = abs(cosSim), colour = model)) +
  facet_nested("Selection/trait alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_cossim_R_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = model, y = meanCosSim, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  labs(x = "Model", 
       y = TeX("Absolute cosine similarity between $g_{max}$ and $r_{max}$"),
       colour = "Model") +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_cossim_gmax_R

ggplot(d_cossim_R, 
       aes(x = model, y = bTMb, colour = model)) +
  facet_nested("Selection/trait alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_cossim_R_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = model, y = meanbTGb, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  labs(x = "Model", 
       y = TeX("Evolvability ($r_{max}^T G ~r_{max}$)"),
       colour = "Model") +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_rtgr

plot_grid(plt_cossim_gmax_R,
          plt_rtgr,
          ncol = 2,
          labels = "AUTO")
ggsave("plt_cossim_g_rmax.png", device = png, width = 12, height = 5, bg = "white")


####################################


# M vs trait correlation
m_cor <- vector(mode = "list", length(m_matrices_tot))
for (i in seq_along(m_matrices_tot)) {
  result <- cov2cor(m_matrices_tot[[i]])
  result[is.nan(result)] <- 0
  m_cor[[i]] <- result
}

# Add correlation matrix to id: index in cor_mats
id_m_tot <- id_m_tot %>%
  mutate(corMatIndex = match(model, model_names))

# Match correlation matrices with M matrix list
# First eigenvector of trait correlation matrix
d_m_cor <- as.data.frame(t(as.data.frame(cor_eig[id_m_tot$corMatIndex])))

# Measure cosine similarity between M matrix correlations and neutral trait correlations
d_cossim_m_traitcor <- GetCosineSimilarity(m_cor, d_m_cor, id_m_tot)

#saveRDS(d_cossim_m_traitcor, "d_cossim_m_traitcor.RDS")

d_cossim_m_traitcor <- readRDS("d_cossim_m_traitcor.RDS")

d_cossim_m_traitcor <- AddCombosToDF(d_cossim_m_traitcor)
d_cossim_m_traitcor$model <- factor(d_cossim_m_traitcor$model, levels = model_names)
d_cossim_m_traitcor$timePoint <- d_cossim_m_traitcor$timePoint - 50000


d_cossim_m_traitcor_sum <- d_cossim_m_traitcor %>%
  group_by(timePoint, model, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))

ggplot(d_cossim_m_traitcor_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted,
               scales = "free",
               labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Absolute cosine similarity between $m_{max}$ and $r_{max}$"),
       colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  scale_fill_manual(values = pal,
                      labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_cossim_mmax_r
plt_cossim_mmax_r


ggplot(d_cossim_m_traitcor_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted,
               scales = "free",
               labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Evolvability ($r_{max}^T M r_{max}$)"),
       colour = "Model") +
  scale_y_continuous(breaks = equal_breaks(3, 0.05),
                     labels = scales::label_number(digits = 2)) +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  scale_fill_manual(values = pal,
                    labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_rtmr
plt_rtmr

plot_grid(plt_cossim_mmax_r,
          plt_rtmr,
          ncol = 2,
          labels = "AUTO")
ggsave("plt_cossim_m_rmax.png", device = png, width = 12, height = 5, bg = "white")



# Adapted populations look like they have less variation along the trait correlation axis
## need to get around the constraint
### can we look at how this is affected by the alignment between selection and trait corr?
### if selection is opposing trait corr, then adaptation will require M and G to become 
### misaligned with R

corBeta_mats <- cor_mats[id_m$corMatIndex]
corBeta_mats <- lapply(corBeta_mats, function(x) return(x[,,2]))

d_cossim_corBeta <- GetCosineSimilarity(corBeta_mats, d_selvec_m %>% select(ends_with("dir")), id_m)

saveRDS(d_cossim_corBeta, "d_cossim_corBeta.RDS")

d_cossim_corBeta <- readRDS("d_cossim_corBeta.RDS")

d_cossim_corBeta <- AddCombosToDF(d_cossim_corBeta)
d_cossim_corBeta$model <- factor(d_cossim_corBeta$model, levels = model_names)
d_cossim_corBeta$timePoint <- d_cossim_corBeta$timePoint - 50000

d_cossim_corBeta_sum <- d_cossim_corBeta %>%
  group_by(timePoint, model, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))

ggplot(d_cossim_corBeta_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted,
               scales = "free",
               labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Absolute cosine similarity between $r_{max}$ and $\\beta$"),
       colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  scale_fill_manual(values = pal,
                    labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_cossim_rmax_beta
plt_cossim_rmax_beta

ggsave("plt_cossim_rmax_beta.png", device = png, width = 10, height = 5, bg = "white")

## Increase in alignment of rmax with beta in FFLC1 adapted models corresponds to 
## increased conditional evolvability (M matrix) and reduced autonomy
## In other words, when trait correlations align with selection, mutational variance
## spikes along that direction

ggplot(d_cossim_corBeta_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted,
               scales = "free",
               labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Evolvability ($\\beta^T R \\beta$)"),
       colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  scale_fill_manual(values = pal,
                    labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_btrb
plt_btrb

# what is the cosine similarity of maladapted/adapted pops
gls_cossim_rmax_beta <- gls(absCosSim ~ model * isAdapted,
                           data = d_cossim_corBeta_abs,
                           weights = varIdent(form = ~ 1 | model * isAdapted))
summary(gls_cossim_rmax_beta)
plot(gls_cossim_rmax_beta)
em_cossim_rmax_beta <- emmeans(gls_cossim_rmax_beta, ~ model * isAdapted,
                        type = "response")
summary(em_cossim_rmax_beta)
test(em_cossim_rmax_beta)
emmip(em_cossim_rmax_beta, model ~ isAdapted)
pwpp(em_cossim_rmax_beta, by = "model")

emmip(em_cossim_rmax_beta, ~ model + isAdapted, CIs = T) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y= TeX("Predicted $abs(cos(\\theta)_\\beta^{rmax})$")) +
  theme(text = element_text(size = 12))
ggsave("plt_pred_cossim_beta_rmax.png", device = png, bg = "white",
       width = 8, height = 6)

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
pwpp(em_bTRb_rmax_beta, by = "model")


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
pwpp(em_cossim_gmax_beta, by = "model")


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
pwpp(em_bTGb_gmax_beta, by = "model")

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
## Evolvability = bTGb

# Join the evolvability estimates with M alignment
d_btgb_Malign <- left_join(d_cossim %>% 
                             select(timePoint, seed, modelindex, isAdapted,
                                    model, r, bTMb, cosSim) %>%
                             rename(bTGb = bTMb) %>%
                             mutate(absCS_Gb = abs(cosSim)) %>%
                             select(-cosSim),
                           d_cossim_m %>%
                             filter(timePoint == 60000 | timePoint == 50000) %>%
                             mutate(timePoint = if_else(timePoint == 50000, "Start", "End")) %>%
                             select(timePoint, seed, modelindex, isAdapted,
                                    model, r, bTMb, cosSim) %>%
                             mutate(absCS_Mb = abs(cosSim)) %>%
                             select(-cosSim),
                           by = c("timePoint", "seed", "modelindex", "isAdapted",
                                  "model", "r"))


ggplot(d_btgb_Malign,
       aes(x = absCS_Mb, y = bTGb, colour = model)) +
#  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Is adapted" + isAdapted) +
  geom_point(shape = 1) +
  labs(x = TeX("Absolute $cos(\\theta)_{M\\beta}$"), 
       y = TeX("Evolvability ($\\beta^T G \\beta$)"),
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

ggplot(d_btgb_Malign,
       aes(x = absCS_Mb, y = absCS_Gb, colour = model)) +
  #  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  facet_nested("Recombination rate (log10)" + log10(r)~ "Is adapted" + isAdapted) +
  geom_point(shape = 1) +
  labs(x = TeX("Alignment of M with $\\beta (abs(cos(\\theta)_\\beta^{M}))$"), 
       y = TeX("Alignment of G with $\\beta (abs(cos(\\theta)_\\beta^{G}))$"),
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



## Summary
d_btgb_Malign_sum <- d_btgb_Malign %>%
  group_by(model, isAdapted) %>%
  summarise(meanCosSim_Gb = mean(absCS_Gb),
            CICosSim_Gb = CI(absCS_Gb),
            meanCosSim_Mb = mean(absCS_Mb),
            CICosSim_Mb = CI(absCS_Mb),
            meanbTGb = mean(bTGb),
            CIbTGb = CI(bTGb))

ggplot(d_btgb_Malign_sum,
       aes(x = meanCosSim_Mb, y = meanCosSim_Gb, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_point(shape = 1) +
  geom_errorbar(aes(ymin = meanCosSim_Gb - CICosSim_Gb,
                    ymax = meanCosSim_Gb + CICosSim_Gb)) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("Mean $abs(cos(\\theta)_\\beta^{M})$"), 
       y = TeX("Mean $abs(cos(\\theta)_\\beta^{G})$"),
       colour = "Model") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

### cos(theta)M_beta vs bTGb
ggplot(d_btgb_Malign_sum,
       aes(x = meanCosSim_Mb, y = meanbTGb, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_point(shape = 1) +
  geom_errorbar(aes(ymin = meanbTGb - CIbTGb,
                    ymax = meanbTGb + CIbTGb)) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("M matrix/selection alignment ($abs(cos(\\theta)_\\beta^{M})$)"), 
       y = TeX("Evolvability ($\\beta^TG\\beta$)"),
       colour = "Model") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

### Alignment vs Selection probability
d_isAdapted <- d_h2_trait %>%
  filter(gen == 50000) %>%
  select(seed, modelindex, model, r, isAdapted) %>%
  distinct(.keep_all = T)
d_isAdapted <- as.data.frame(table(d_isAdapted$model, d_isAdapted$isAdapted))
d_isAdapted <- d_isAdapted %>%
  rename(model = Var1,
         isAdapted = Var2)

# Check we have 208 replicates * 3 r levels = 624 per model
d_isAdapted %>%
  group_by(model) %>%
  summarise(total = sum(Freq))

# join isAdapted counts
d_btgb_Malign_sum <- left_join(d_btgb_Malign_sum %>% 
                                 mutate(model = factor(model, levels = model_names),
                                        isAdapted = factor(isAdapted)),
                               d_isAdapted,
                               by = c("model", "isAdapted"))

ggplot(d_btgb_Malign_sum %>% 
         filter(isAdapted == T) %>% 
         mutate(Freq = Freq / 624),
       aes(x = meanCosSim_Mb, y = Freq, colour = model)) +
  #geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_point(shape = 1) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("M matrix/selection alignment ($abs(cos(\\theta)_\\beta^{M})$)"), 
       y = "Proportion of populations adapted",
       colour = "Model") +
  coord_flip(xlim = c(0, 1)) + 
  scale_colour_manual(values = pal,
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
ggsave("plt_Mcossim_probAdapt.png", device = png, width = 6, height = 6, bg = "white")



d_btgb_Malign_adapted <- d_btgb_Malign %>% filter(isAdapted)

gls_btmb_adapt <- gls(bTMb ~ model, data = d_btgb_Malign_adapted,
                      weights = varIdent(form = ~ 1 | model))
summary(gls_btmb_adapt)
em_btgb_adapt <- emmeans(gls_btmb_adapt, ~ model,
                          type = "response")
summary(em_btgb_adapt)
test(em_btgb_adapt)
emmip(em_btgb_adapt, ~ model, CIs = T) +
  theme_bw() +
  labs(x = TeX("Model"), y= TeX("Predicted $\\beta^TM\\beta$")) +
  theme(text = element_text(size = 12))
emmeans::contrast(em_btgb_adapt)


# Does M/beta alignment predict population adaptedness?
## use random forest
d_btgb_Malign_rf <- d_btgb_Malign %>%
  mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE")),
         model = factor(model, levels = model_names),
         timePoint = factor(timePoint, levels = c("Start", "End")),
         r = factor(log10(r), levels = c(-10, -5, -1))) %>%
  select(isAdapted, model, r, absCS_Mb, absCS_Gb, bTGb, bTMb)

# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215
set.seed(seed)
idx <- sample(2, nrow(d_btgb_Malign_rf), replace = T, prob = c(0.7, 0.3))
train_mbeta_adapted <- d_btgb_Malign_rf[idx == 1,]
test_mbeta_adapted <- d_btgb_Malign_rf[idx == 2,]

rf_mbeta_adapted <- randomForest(formula = isAdapted ~ model * r * absCS_Gb * absCS_Mb * bTGb * bTMb,
                                 data = train_mbeta_adapted,
                                 ntree = 500,
                                 proximity = T,
                                 importance = T,
                                 type = "classification")

print(rf_mbeta_adapted)

# Training data
p_train_mbeta_adapted <- predict(rf_mbeta_adapted, train_mbeta_adapted)
caret::confusionMatrix(p_train_mbeta_adapted, train_mbeta_adapted$isAdapted)

# Test data
p_test_mbeta_adapted <- predict(rf_mbeta_adapted, test_mbeta_adapted)
caret::confusionMatrix(p_test_mbeta_adapted, test_mbeta_adapted$isAdapted)


# Plot errors (black line = OOB, red = false positive, green = false negative)
plot(rf_mbeta_adapted)


# Number of nodes per tree
hist(treesize(rf_mbeta_adapted),
     main = "# Nodes for the RF trees",
     col = "forestgreen")

# variable importance - model, btmb, and cossim Mb important
varImpPlot(rf_mbeta_adapted,
           sort = T, type = 2, scale = F)
randomForest::importance(rf_mbeta_adapted, type = 1, scale = F)

# conditional
iscores_crf <- varimp(crf_mbeta_adapted, conditional = T)
iscores_crf

saveRDS(iscores_crf, "iscores_crf.RDS")
# iscores_crf <- readRDS("iscores_crf.RDS")

iscores <- varImp(rf_mbeta_adapted, conditional = T)
iscores <- iscores %>% tibble::rownames_to_column("feature") 
iscores$feature <- iscores$feature %>% as.factor()
iscores <- iscores[,-3]
colnames(iscores) <- c("feature", "score")

bor_mbeta <- Boruta::Boruta(isAdapted ~ ., data = d_btgb_Malign_rf)
bor_mbeta
plot(bor_mbeta)

d_bor_mbeta <- process_the_Boruta_data(bor_mbeta)

feature_names <- c(TeX("Shadow Variable (min)", output = "character"),
                   TeX("Shadow Variable (mean)", output = "character"),
                   TeX("Shadow Variable (max)", output = "character"),
                   TeX("Recombination rate", output = "character"),
                   TeX("G-selection alignment ($|cos(\\theta)|^G_\\beta$)",
                       output = "character"),
                   TeX("Evolvability ($\\beta^TG\\beta$)", output = "character"),
                   TeX("M-selection alignment ($|cos(\\theta)|^M_\\beta$)",
                       output = "character"),
                   TeX("M Evolvability ($\\beta^TM\\beta$)", output = "character"),
                   TeX("Motif", output = "character"))


ggplot(iscores %>%
         mutate(feature = factor(feature, levels = c("r", "absCS_Gb", "bTGb", 
                                                     "absCS_Mb", "bTMb", "model"))),
       aes(x = feature, y = score)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = parse(text = feature_names[4:9]),
                   guide = guide_axis(n.dodge = 2)) + 
  theme_bw() +
  labs(x = "Feature", y = "Importance") +
  theme(text = element_text(size = 12)) -> plt_rf_perm_imp


pal_boruta <- c(rep("#00C0EA", times = 3),
                rep("#00A000", times = 6))

ggplot(d_bor_mbeta %>% pivot_longer(everything()) %>%
         mutate(x = fct_reorder(name, value, median)),
       aes(x = x, y = value, fill = x)) +
  geom_boxplot(show.legend = F) +
  theme_bw() +
  scale_fill_manual(values = pal_boruta) +
  scale_x_discrete(labels = parse(text = feature_names),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Importance") +
  theme(text = element_text(size = 12)) -> plt_boruta_imp
plt_boruta_imp
ggsave("plt_boruta_import.png", device = png, bg = "white",
       width = 12, height = 8)

predictor <- iml::Predictor$new(rf_mbeta_adapted, 
                                data = test_mbeta_adapted[, 2:7], 
                                y = test_mbeta_adapted$isAdapted)

imp <- iml::FeatureImp$new(predictor,
                           loss = "ce",
                           n.repetitions = 100)

ggplot(imp$results %>% 
         mutate(feature = factor(feature, levels = c("r", "absCS_Gb", "bTGb", 
                                                     "absCS_Mb", "bTMb", "model"))),
       aes(x = feature, y = importance)) +
  geom_point() +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                width = 0.2) +
  scale_x_discrete(labels = parse(text = feature_names[4:9]),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Permutation Importance (Loss: CE)") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_perm_imp
plt_perm_imp
ggsave("plt_perm_feat_imp.png", device = png, width = 9, height = 5, bg = "white")

plot_grid(plt_boruta_imp,
          plt_perm_imp,
          nrow = 2,
          labels = "AUTO",
          align = "v")
ggsave("plt_feat_imp.png", device = png, width = 12, height = 9, bg = "white")


ale <- FeatureEffects$new(predictor, grid.size = 10)
ale$plot()

ale_plots <- vector(mode = "list", length = 7)

ale_labels <- feature_names[c(9, 4, 7, 5, 6, 8, 10)]
names(ale_labels) <- names(ale$results)

for (i in seq_along(ale_plots)) {
  d_ale <- ale$results[[i]]
  d_ale <- d_ale %>% filter(.class == "Adapted")
  x_label <- ale_labels[d_ale$.feature[1]] 
  
  if (d_ale$.feature[1] == "model") {
    d_ale$.borders <- factor(d_ale$.borders, 
                             levels = model_names,
                             labels = model_names_noquote)
  }
  
  if (!is.numeric(d_ale$.borders[1])) {
    geom_fn <- geom_lollipop
    scale_fn <- scale_x_discrete
  } else {
    geom_fn <- geom_line
    scale_fn <- scale_x_continuous
  }
  ale_plots[[i]] <- ggplot(d_ale,
         aes(x = .borders, y = .value)) +
    geom_fn() +
    scale_fn() +
    theme_bw() +
    labs(x = parse(text = x_label), y = "ALE of adaptation probability")
}

plot_grid(plotlist = ale_plots,
          labels= "AUTO")
ggsave("plt_ale.png", device = png, bg = "white",
       width = 9, height = 6)



# partial dependence
partialPlot(x = rf_mbeta_adapted,
            pred.data = train_mbeta_adapted,
            x.var = model)
partialPlot(x = rf_mbeta_adapted,
            pred.data = train_mbeta_adapted,
            x.var = absCS_Mb)


MDSplot(rf_mbeta_adapted, train_mbeta_adapted$model,
        pal = pal)
legend("topright", legend = model_names_noquote,
       fill = pal)
ggsave("plt_rf_mbeta_adapted_mds.png", device = png, 
       width = 5, height = 5, bg = "white")



#########################

gls_btgb_malign <- gls(bTGb ~ absCS_Mb * model,
                          data = d_btgb_Malign,
                          weights = varIdent(form = ~ 1 | model))
summary(gls_btgb_malign)
em_btgb_malign <- emmeans(gls_btgb_malign, ~ absCS_Mb * model,
                             type = "response", 
                          at = list(absCS_Mb = c(0, 0.5, 1)))
summary(em_btgb_malign)
test(em_btgb_malign)
emmip(em_btgb_malign, model ~ absCS_Mb, CIs = T) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  labs(x = TeX("Predicted $abs(cos(\\theta)_M^\\beta)$"), y= TeX("Predicted $abs(cos(\\theta)_G^\\beta)$")) +
  theme(text = element_text(size = 12))
emmeans::contrast(em_btgb_malign)
ggsave("plt_pred_cossimGbeta_cossimMbeta.png", device = png, bg = "white",
       width = 8, height = 6)

# Both G and M tend to be correlated with the direction of selection, not for NAR or PAR though
## But not a 1:1 relationship, Gmax is often misaligned with beta



## Positive relationship between bTGb and cosine similarity between M and b
## evolvability along the selection gradient increases when mutational variance is biased
## towards that direction - not really surprising


# 3) Positive and negative controls confirm developmental bias

# Load in data
PATH_QG_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"

d_qg_orth <- data.table::fread(PATH_QG_ORTH, header = F, 
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
d_qg_orth <- AddCombosToDF(d_qg_orth) 

d_qg_orth %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98)) %>%
  mutate(model = factor(model, levels = model_names),
         modelindex = factor(modelindex, levels = 1:15)) %>%
  ungroup() -> d_qg_orth


PATH_QG_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_qg.csv"

d_qg_par <- data.table::fread(PATH_QG_PAR, header = F, 
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
d_qg_par <- AddCombosToDF(d_qg_par) 

d_qg_par %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98)) %>%
  mutate(model = factor(model, levels = model_names),
         modelindex = factor(modelindex, levels = 1:15)) %>%
  ungroup() -> d_qg_par


DATA_PATH_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_mutvar.csv"
DATA_PATH_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_mutvar.csv"

d_m_orth <- read_csv(DATA_PATH_ORTH, col_names = c("gen", "seed", "modelindex",
                                         paste0("mean_", 1:4),
                                         paste0("var_", 1:4),
                                         paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m_par <- read_csv(DATA_PATH_PAR, col_names = c("gen", "seed", "modelindex",
                                                   paste0("mean_", 1:4),
                                                   paste0("var_", 1:4),
                                                   paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m_orth <- d_m_orth %>%
  mutate(model = ModelFromIndexWithR(modelindex))

d_m_par <- d_m_par %>%
  mutate(model = ModelFromIndexWithR(modelindex))

# get matrices
m_matrices_orth <- d_m_orth %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))
# 
m_matrices_par <- d_m_par %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))


#saveRDS(m_matrices_orth, "m_matrices_orth.RDS")
#saveRDS(m_matrices_par, "m_matrices_par.RDS")

m_matrices_orth <- readRDS("m_matrices_orth.RDS")
m_matrices_par <- readRDS("m_matrices_par.RDS")


# Get eigenvectors of each M
e_m_orth <- lapply(m_matrices_orth, eigen)
saveRDS(e_m_orth, "eigen_randomised_m_orth.RDS")

e_m_par <- lapply(m_matrices_par, eigen)
saveRDS(e_m_par, "eigen_randomised_m_par.RDS")

e_m_orth <- readRDS("eigen_randomised_m_orth.RDS")
e_m_par <- readRDS("eigen_randomised_m_par.RDS")

# Calculate relative eigenvalue dispersion
vrel_m_orth <- unlist(lapply(e_m_orth, function(x) { Vrel(x$values) }))
vrel_m_par <- unlist(lapply(e_m_par, function(x) { Vrel(x$values) }))

# join with quant gen data
d_vrel_orth <- left_join(d_m_orth %>% mutate(seed = factor(seed),
                                   modelindex = factor(modelindex, levels = 1:15),
                                   r = RFromIndex(modelindex),
                                   vrel = vrel_m_orth),
                    d_qg_orth %>% select(gen, seed, modelindex, model, r, isAdapted) %>%
                      mutate(model = str_remove_all(model, "'")),
                    by = c("gen", "seed", "modelindex", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote),
         dataset = "Orthogonal")

d_vrel_par <- left_join(d_m_par %>% mutate(seed = factor(seed),
                                             modelindex = factor(modelindex, levels = 1:15),
                                             r = RFromIndex(modelindex),
                                             vrel = vrel_m_par),
                         d_qg_par %>% select(gen, seed, modelindex, model, r, isAdapted) %>%
                           mutate(model = str_remove_all(model, "'")),
                         by = c("gen", "seed", "modelindex", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote),
         dataset = "Parallel")

d_vrel$dataset <- "Randomised"

d_vrel_tot <- rbind(d_vrel, d_vrel_orth, d_vrel_par)


# Plot
d_vrel_sum <- d_vrel_tot %>%
  group_by(gen, model, dataset, isAdapted) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))

ggplot(d_vrel_sum,
       aes(x = gen - 50000, y = vrel_mean, colour = model)) +
  facet_nested("Dataset" + dataset~"Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = vrel_mean - vrel_CI, ymax = vrel_mean + vrel_CI, fill = model),
              colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")



# Stable over time (except FFLC1, mirrors M evolvability), can take average?
d_vrel_tab <- d_vrel_tot %>%
  group_by(model, dataset, isAdapted) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))
d_vrel_tab

ggplot(d_vrel_tot,
       aes(x = model, y = vrel, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset~"Population adapted" + isAdapted) +
  geom_boxplot() +
  geom_point(data = d_vrel_tab, aes(y = vrel_mean), size = 2, stroke = 1,
             shape = 21,
             colour = "black", fill = "white") +
  labs(x = "Model", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_colour_manual(values = pal,
                      breaks = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_vrel_m
plt_vrel_m
ggsave("plt_vrel_m_alignment.png", device = png, bg = "white",
       width = 8, height = 6)

# No real difference in Vrel for M matrices, makes sense given selection not acting

#####################
# Measure cosine similarity between M matrices and beta
d_opt_par <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_opt.csv", col_names = F)
d_opt_orth <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_opt.csv", col_names = F)

# o = optimum, s = sigma, d = direction (-1, 1)
colnames(d_opt_par) <- c("seed", "modelindex", "o_t1", "o_t2", "o_t3", "o_t4", 
                     "s_t1", "s_t2", "s_t3", "s_t4", "d_t1", "d_t2", "d_t3",
                     "d_t4")
colnames(d_opt_orth) <- colnames(d_opt_par)

# Last value in d_opt_par and _orth is the eigenvalue of eigenvector r_i, we can ignore it
d_opt_par <- d_opt_par[,-15]
d_opt_orth <- d_opt_orth[,-15]

d_selvec_m_orth <- d_qg_orth %>%
  filter(gen >= 50000) %>%
  select(gen, seed, modelindex, isAdapted, dataset, ends_with("mean"))

d_selvec_m_par <- d_qg_par %>%
  filter(gen >= 50000) %>%
  select(gen, seed, modelindex, isAdapted, dataset, ends_with("mean"))


d_selvec_m_orth <- left_join(d_selvec_m_orth, d_opt_orth %>% 
                          select(seed, modelindex, dataset, starts_with("o_")) %>%
                          mutate(seed = factor(seed),
                                 modelindex = factor(modelindex)), 
                        by = c("seed", "modelindex", "dataset"))

d_selvec_m_par <- left_join(d_selvec_m_par, d_opt_par %>% 
                               select(seed, modelindex, dataset, starts_with("o_")) %>%
                               mutate(seed = factor(seed),
                                      modelindex = factor(modelindex)), 
                             by = c("seed", "modelindex", "dataset"))


d_selvec_m_orth <- AddCombosToDF(d_selvec_m_orth)
d_selvec_m_par <- AddCombosToDF(d_selvec_m_par)

d_selvec_m_orth <- d_selvec_m_orth %>%
  mutate(modelindex = factor(modelindex, levels = 1:15),
         seed = as.factor(seed),
         model = factor(model, levels = model_names)) %>%
  rename(timePoint = gen)

d_selvec_m_par <- d_selvec_m_par %>%
  mutate(modelindex = factor(modelindex, levels = 1:15),
         seed = as.factor(seed),
         model = factor(model, levels = model_names)) %>%
  rename(timePoint = gen)



id_m_orth <- d_m_orth %>% mutate(timePoint = gen) %>% select(timePoint, seed, modelindex)
id_m_orth$clus <- 1
id_m_orth$modelindex <- factor(id_m_orth$modelindex, levels = 1:15)
id_m_orth$seed <- as.factor(id_m_orth$seed)

id_m_orth <- AddCombosToDF(id_m_orth)
id_m_orth$model <- factor(id_m_orth$model, levels = model_names)

id_m_orth <- inner_join(id_m_orth, d_qg_orth %>% mutate(timePoint = gen) %>% 
                     select(timePoint, seed, modelindex, isAdapted),
                   by = c("timePoint", "seed", "modelindex"))


id_m_par <- d_m_par %>% mutate(timePoint = gen) %>% select(timePoint, seed, modelindex)
id_m_par$clus <- 1
id_m_par$modelindex <- factor(id_m_par$modelindex, levels = 1:15)
id_m_par$seed <- as.factor(id_m_par$seed)

id_m_par <- AddCombosToDF(id_m_par)
id_m_par$model <- factor(id_m_par$model, levels = model_names)

id_m_par <- inner_join(id_m_par, d_qg_par %>% mutate(timePoint = gen) %>% 
                          select(timePoint, seed, modelindex, isAdapted),
                        by = c("timePoint", "seed", "modelindex"))


d_selvec_m_orth <- inner_join(id_m_orth, d_selvec_m_orth, 
                         by = c("timePoint", "seed", "modelindex", "isAdapted", "model", "r"))

d_selvec_m_orth <- d_selvec_m_orth %>%
  mutate(t1_dir = o_t1 - trait1_mean,
         t2_dir = o_t2 - trait2_mean,
         t3_dir = o_t3 - trait3_mean,
         t4_dir = o_t4 - trait4_mean,
         norm = sqrt(rowSums(pick(ends_with("dir"))^2)), # normalise
         t1_dir = t1_dir / norm,
         t2_dir = t2_dir / norm,
         t3_dir = t3_dir / norm,
         t4_dir = t4_dir / norm) %>%
  select(timePoint, seed, modelindex, dataset, isAdapted, model, r, norm, ends_with("dir"))

d_selvec_m_par <- inner_join(id_m_par, d_selvec_m_par, 
                             by = c("timePoint", "seed", "modelindex", "isAdapted", "model", "r"))

d_selvec_m_par <- d_selvec_m_par %>%
  mutate(t1_dir = o_t1 - trait1_mean,
         t2_dir = o_t2 - trait2_mean,
         t3_dir = o_t3 - trait3_mean,
         t4_dir = o_t4 - trait4_mean,
         norm = sqrt(rowSums(pick(ends_with("dir"))^2)), # normalise
         t1_dir = t1_dir / norm,
         t2_dir = t2_dir / norm,
         t3_dir = t3_dir / norm,
         t4_dir = t4_dir / norm) %>%
  select(timePoint, seed, modelindex, isAdapted, model, r, norm, ends_with("dir"))


d_cossim_m_orth <- GetCosineSimilarity(m_matrices_orth, d_selvec_m_orth %>% select(ends_with("dir")), id_m_orth)
d_cossim_m_par <- GetCosineSimilarity(m_matrices_par, d_selvec_m_par %>% select(ends_with("dir")), id_m_par)

d_cossim_m_orth$dataset <- "Orthogonal"
d_cossim_m_par$dataset <- "Parallel"

d_cossim_m <- readRDS("d_cossim_m.RDS")
d_cossim_m$dataset <- "Randomised"

# Combine
d_cossim_m <- rbind(d_cossim_m, d_cossim_m_orth, d_cossim_m_par)
saveRDS(d_cossim_m, "d_cossim_m_datasets.RDS")
d_cossim_m <- readRDS("d_cossim_m_datasets.RDS")

d_cossim_m <- AddCombosToDF(d_cossim_m)

d_cossim_m_sum <- d_cossim_m %>%
  mutate(timePoint = timePoint - 50000) %>%
  group_by(timePoint, model, dataset, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(sqrt(cosSim^2), na.rm = T),
                   seCosSim = se(sqrt(cosSim^2), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_m_sum$model <- factor(d_cossim_m_sum$model, levels = model_names)


ggplot(d_cossim_m_sum, 
       aes(x = timePoint, y = meanCosSim, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset~"Population adapted" + isAdapted) +
  geom_line() +
  geom_ribbon(aes(ymin = meanCosSim - seCosSim, ymax = meanCosSim + seCosSim,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Absolute cosine similarity between $m_{max}$ and $\\beta"),
       colour = "Model") +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_manual(values = pal, labels = model_names_noquote) +
  scale_fill_manual(values = pal, labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim_m_sum, 
       aes(x = timePoint, y = meanbTGb, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset~"Population adapted" + isAdapted,
               scales = "free",
               labeller = labeller(model = model_names_labeller)) +
  geom_line() +
  geom_ribbon(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                  fill = model), colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generations post-optimum shift", 
       y = TeX("Evolvability ($\\beta^T M \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = pal, labels = model_names_noquote) +
  scale_fill_manual(values = pal, labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

## G matrix

G_ORTH_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"

d_h2_trait_mkr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_trait_mkr.csv"), col_names = F)
d_h2_trait_mrr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_trait_mrr.csv"), col_names = F)

d_h2_trait_mkr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_trait_mkr.csv"), col_names = F)
d_h2_trait_mrr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_trait_mrr.csv"), col_names = F)


colnames(d_h2_trait_mkr_orth) <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_t1",
                              "VA_t2", "VA_t3", "VA_t4", "CVA_t1_t2", "CVA_t1_t3",
                              "CVA_t1_t4", "CVA_t2_t3", "CVA_t2_t4", "CVA_t3_t4",
                              "h2_t1", "h2_t2", "h2_t3", "h2_t4")

colnames(d_h2_trait_mrr_orth) <- colnames(d_h2_trait_mkr_orth)
colnames(d_h2_trait_mkr_par) <- colnames(d_h2_trait_mkr_orth)
colnames(d_h2_trait_mrr_par) <- colnames(d_h2_trait_mkr_orth)

# join
d_h2_trait_mkr_orth$calcMode <- "mkr"
d_h2_trait_mrr_orth$calcMode <- "mrr"
d_h2_trait_mkr_par$calcMode <- "mkr"
d_h2_trait_mrr_par$calcMode <- "mrr"

d_h2_trait_mkr_orth$dataset <- "Orthogonal"
d_h2_trait_mrr_orth$dataset <- "Orthogonal"
d_h2_trait_mkr_par$dataset <- "Parallel"
d_h2_trait_mrr_par$dataset <- "Parallel"
d_h2_trait_mkr$dataset <- "Randomised"
d_h2_trait_mrr$dataset <- "Randomised"


d_h2_trait <- rbind(d_h2_trait_mkr, d_h2_trait_mrr, 
                    d_h2_trait_mkr_orth, d_h2_trait_mrr_orth,
                    d_h2_trait_mkr_par, d_h2_trait_mrr_par)

d_h2_trait %>% mutate(model = d_combos$model[.$modelindex],
                      model = factor(model, levels = model_names),
                      r = d_combos$r[.$modelindex]) -> d_h2_trait

d_h2_trait <- d_h2_trait %>%
  distinct(gen, seed, modelindex, dataset, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()

# Join qg together
d_qg$dataset <- "Randomised"
d_qg_orth$dataset <- "Orthogonal"
d_qg_par$dataset <- "Parallel"

d_qg_tot <- rbind(d_qg, d_qg_orth, d_qg_par)

d_qg_optPerc <- d_qg_tot %>% select(gen, seed, modelindex, dataset, isAdapted) %>% filter(gen >= 49500)

# inner join optPerc
d_h2_trait <- left_join(d_h2_trait, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))

# Counts for each model type:
table(d_h2_trait$model, d_h2_trait$isAdapted)
table(d_h2_trait$model, d_h2_trait$dataset, d_h2_trait$isAdapted)

# table of prop adapted

tab_propadapted_aligned <- d_qg_tot %>%
  filter(gen == 50000, log10(r) == -1) %>%
  mutate(model = factor(model, levels = model_names,
                        labels = model_names_noquote)) %>%
  group_by(model, dataset) %>%
  dplyr::summarise(propAdapted = sum(isAdapted) / n())

stargazer::stargazer(as.data.frame(tab_propadapted_aligned),
                     summary = F, rownames = F)


# Discretise generation
d_h2_trait <- d_h2_trait %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))

# summarise
d_h2_trait_sum <- d_h2_trait %>% 
  group_by(timePoint, model, dataset, r, isAdapted) %>%
  dplyr::summarise(meanH2w = mean(h2_w, na.rm = T),
                   seH2w = se(h2_w, na.rm = T),
                   meanVAw = mean(VA_w, na.rm = T),
                   seVAw = se(VA_w, na.rm = T))
d_h2_trait_sum$model <- as.factor(d_h2_trait_sum$model)

# Split h2 into G matrices
d_h2_trait %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(5:8, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(modelindex, timePoint, dataset, isAdapted) %>%
  group_split(.) -> split_h2


# Separate into model indices
# each sublist is replicates of a model index
sourceCpp("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")

lapply(split_h2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices


# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat <- unlist(cov_matrices, recursive = F)

# get ids from the matrix
cov_matrix_modelindex <- GetMatrixIDsWithDataset(split_h2)

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
d_selvec <- d_qg_tot %>%
  filter(gen == 50000 | gen == 60000) %>%
  select(gen, seed, modelindex, dataset, ends_with("mean"))

# join opt
d_opt$dataset <- "Randomised"
d_opt_orth$dataset <- "Orthogonal"
d_opt_par$dataset <- "Parallel"

d_opt_tot <- rbind(d_opt, d_opt_orth, d_opt_par)

d_selvec <- left_join(d_selvec, d_opt_tot %>% 
                        select(seed, modelindex, dataset, starts_with("o_")) %>%
                        mutate(seed = factor(seed),
                               modelindex = factor(modelindex)), 
                      by = c("seed", "modelindex", "dataset"))

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
  select(gen, seed, modelindex, dataset, model, r, norm, ends_with("dir"))


d_selvec <- d_selvec %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
  select(-gen)

d_selvec2 <- inner_join(id, d_selvec, 
                        by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim <- GetCosineSimilarity(h2_pd, d_selvec2 %>% select(ends_with("dir")), id)

d_cossim <- AddCombosToDF(d_cossim)

d_cossim$dataset <- d_selvec2$dataset

d_cossim_sum <- d_cossim %>%
  group_by(timePoint, model, dataset, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))
d_cossim_sum$model <- as.factor(d_cossim_sum$model)


ggplot(d_cossim, 
       aes(x = timePoint, y = sqrt(cosSim^2), colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanCosSim, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Absolute cosine similarity between $g_{max}$ and $\\beta"),
       colour = "Model") +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_cossim, 
       aes(x = timePoint, y = bTMb, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_cossim_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanbTGb, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Evolvability ($\\beta^T G \\beta$)"),
       colour = "Model") +
  scale_colour_manual(values = pal,
                      labels = model_names) +
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

d_ecr <- CalcECRATrait(h2_pd, id)
d_ecr <- AddCombosToDF(d_ecr)
d_ecr$dataset <- d_selvec2$dataset

# Refactor model
d_ecr <- d_ecr %>%
  mutate(model = factor(model, levels = model_names))

# Need to calculate cev means separately for the different models
# K- shouldn't mean over cev_KZ and KXZ
d_ecr_sum <- d_ecr %>%
  group_by(model, dataset, isAdapted) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = cev, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = cev_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean conditional evolvability (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev
plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = res, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = res_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean respondability (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_res
plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = aut, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = aut_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean autonomy (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_aut
plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = ev, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F, show.legend = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = ev_mean, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  labs(x = "Model", y = "Mean evolvability (G matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_ev
plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                      plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO", label_size = 12)

plt_evol <- plot_grid(plt_evol,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol_g_alignment.png", device = png, bg = "white",
       width = 12, height = 8)


## Evolvability against alignment of M with direction of selection
## Evolvability = bTGb

# Join the evolvability estimates with M alignment
d_btgb_Malign_tot <- left_join(d_cossim %>% 
                             select(timePoint, seed, modelindex, dataset, isAdapted,
                                    model, r, bTMb, cosSim) %>%
                             rename(bTGb = bTMb) %>%
                             mutate(absCS_Gb = abs(cosSim)) %>%
                             select(-cosSim),
                           d_cossim_m %>%
                             filter(timePoint == 60000 | timePoint == 50000) %>%
                             mutate(timePoint = if_else(timePoint == 50000, "Start", "End"),
                                    timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
                             select(timePoint, seed, modelindex, dataset, isAdapted,
                                    model, r, bTMb, cosSim) %>%
                             mutate(absCS_Mb = abs(cosSim)) %>%
                             select(-cosSim),
                           by = c("timePoint", "seed", "modelindex", "dataset", "isAdapted",
                                  "model", "r"))


ggplot(d_btgb_Malign_tot,
       aes(x = absCS_Mb, y = bTGb, colour = model)) +
  #  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_point(shape = 1) +
  labs(x = TeX("Absolute $cos(\\theta)_{M\\beta}$"), 
       y = TeX("Evolvability ($\\beta^T G \\beta$)"),
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

ggplot(d_btgb_Malign_tot,
       aes(x = absCS_Mb, y = absCS_Gb, colour = model)) +
  #  facet_nested("Recombination rate (log10)" + log10(r)~ "Population adapted" + isAdapted) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_point(shape = 1) +
  labs(x = TeX("Alignment of M with $\\beta (abs(cos(\\theta)_\\beta^{M}))$"), 
       y = TeX("Alignment of G with $\\beta (abs(cos(\\theta)_\\beta^{G}))$"),
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



## Summary
d_btgb_Malign_sum <- d_btgb_Malign_tot %>%
  group_by(model, dataset, isAdapted) %>%
  summarise(meanCosSim_Gb = mean(absCS_Gb),
            CICosSim_Gb = CI(absCS_Gb),
            meanCosSim_Mb = mean(absCS_Mb),
            CICosSim_Mb = CI(absCS_Mb),
            meanbTGb = mean(bTGb),
            CIbTGb = CI(bTGb))

ggplot(d_btgb_Malign_sum,
       aes(x = meanCosSim_Mb, y = meanCosSim_Gb, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_point(shape = 1) +
  geom_errorbar(aes(ymin = meanCosSim_Gb - CICosSim_Gb,
                    ymax = meanCosSim_Gb + CICosSim_Gb)) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("Mean $abs(cos(\\theta)_\\beta^{M})$"), 
       y = TeX("Mean $abs(cos(\\theta)_\\beta^{G})$"),
       colour = "Model") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

### cos(theta)M_beta vs bTGb
ggplot(d_btgb_Malign_sum,
       aes(x = meanCosSim_Mb, y = meanbTGb, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset ~ "Population adapted" + isAdapted) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_point(shape = 1) +
  geom_errorbar(aes(ymin = meanbTGb - CIbTGb,
                    ymax = meanbTGb + CIbTGb)) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("M matrix/selection alignment ($|(cos(\\theta)|_\\beta^{M}$)"), 
       y = TeX("Evolvability ($\\beta^TG\\beta$)"),
       colour = "Model") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

### Alignment vs Selection probability
d_isAdapted <- d_h2_trait %>%
  filter(gen == 50000) %>%
  select(seed, modelindex, model, dataset, r, isAdapted) %>%
  distinct(.keep_all = T)
d_isAdapted <- as.data.frame(table(d_isAdapted$model, d_isAdapted$isAdapted))
d_isAdapted <- d_isAdapted %>%
  rename(model = Var1,
         isAdapted = Var2)

# Check we have 208 replicates * 3 r levels = 624 per model
# plus the orth and parallel sets of another 208 replicates = 1040 per model
d_isAdapted %>%
  group_by(model) %>%
  summarise(total = sum(Freq))

# join isAdapted counts
d_btgb_Malign_sum <- left_join(d_btgb_Malign_sum %>% 
                                 mutate(model = factor(model, levels = model_names),
                                        isAdapted = factor(isAdapted)),
                               d_isAdapted,
                               by = c("model", "isAdapted"))

ggplot(d_btgb_Malign_sum %>% 
         filter(isAdapted == T) %>% 
         mutate(Freq = Freq / 1040),
       aes(x = meanCosSim_Mb, y = Freq, colour = model)) +
  #geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  facet_nested(.~"Selection/trait correlation alignment" + dataset) +
  geom_point(shape = 1) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("M matrix/selection alignment ($abs(cos(\\theta)_\\beta^{M})$)"), 
       y = "Proportion of populations adapted",
       colour = "Model") +
  coord_flip(xlim = c(0, 1)) + 
  scale_colour_manual(values = pal,
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
ggsave("plt_Mcossim_probAdapt_align.png", device = png, width = 6, height = 6, bg = "white")



# Does M/beta alignment predict population adaptedness?
## use random forest
d_btgb_Malign_rf <- d_btgb_Malign_tot %>%
  mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                            labels = c("Adapted", "Maladapted")),
         model = factor(model, levels = model_names),
         timePoint = factor(timePoint, levels = c("Start", "End")),
         dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
         r = factor(log10(r), levels = c(-10, -5, -1))) %>%
  select(isAdapted, model, dataset, r, absCS_Mb, absCS_Gb, bTGb, bTMb)

# TODO: ADD VREL

# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215
set.seed(seed)
# Sample per group to avoid unbalanced groups
adapted_counts <- table(d_btgb_Malign_rf$isAdapted)
total_counts <- sum(adapted_counts)
num_responses <- length(adapted_counts)
adapted_weights <- total_counts / (num_responses * adapted_counts)
names(adapted_weights) <- levels(d_btgb_Malign_rf$isAdapted)


idx <- sample(2, nrow(d_btgb_Malign_rf), replace = T, prob = c(0.7, 0.3))
train_mbeta_adapted <- d_btgb_Malign_rf[idx == 1,]
test_mbeta_adapted <- d_btgb_Malign_rf[idx == 2,]

# no balancing
rf_mbeta_adapted_nobal <- randomForest(formula = isAdapted ~ model * dataset * r * absCS_Gb * absCS_Mb * bTGb * bTMb,
                                 data = train_mbeta_adapted,
                                 ntree = 500,
                                 proximity = T,
                                 importance = T,
                                 type = "classification")

print(rf_mbeta_adapted_nobal)

# With balancing (class weights)
rf_mbeta_adapted <- randomForest(formula = isAdapted ~ model * dataset * r * absCS_Gb * absCS_Mb * bTGb * bTMb,
                                 data = train_mbeta_adapted,
                                 strata = train_mbeta_adapted$isAdapted,
                                classwt = adapted_weights,
                                 ntree = 500,
                                 proximity = T,
                                 importance = T,
                                 type = "classification")

print(rf_mbeta_adapted)

# Training data
p_train_mbeta_adapted <- predict(rf_mbeta_adapted, train_mbeta_adapted)
caret::confusionMatrix(p_train_mbeta_adapted, train_mbeta_adapted$isAdapted)

p_train_mbeta_adapted_nobal <- predict(rf_mbeta_adapted_nobal, train_mbeta_adapted)
caret::confusionMatrix(p_train_mbeta_adapted_nobal, train_mbeta_adapted$isAdapted)


# Test data
p_test_mbeta_adapted <- predict(rf_mbeta_adapted, test_mbeta_adapted)
p_test_mbeta_adapted_probs <- predict(rf_mbeta_adapted, test_mbeta_adapted,
                                      type = "prob")[,1]

caret::confusionMatrix(p_test_mbeta_adapted, test_mbeta_adapted$isAdapted)

p_test_mbeta_adapted_nobal <- predict(rf_mbeta_adapted_nobal, test_mbeta_adapted)
p_test_mbeta_adapted_nobal_probs <- predict(rf_mbeta_adapted_nobal, test_mbeta_adapted,
                                      type = "prob")[,1]

caret::confusionMatrix(p_test_mbeta_adapted_nobal, test_mbeta_adapted$isAdapted)



# roc
d_roc <- roc(response = test_mbeta_adapted$isAdapted,
    predictor = p_test_mbeta_adapted_probs)

d_roc_nobal <- roc(response = test_mbeta_adapted$isAdapted,
             predictor = p_test_mbeta_adapted_nobal_probs,
             levels = rev(levels(test_mbeta_adapted$isAdapted)))

d_rocs <- data.frame(model = c(rep("Weighted", times = length(d_roc$sensitivities)),
                               rep("Unbalanced", times = length(d_roc_nobal$sensitivities))),
                     sens = c(d_roc$sensitivities, d_roc_nobal$sensitivities),
                     spec = c(d_roc$specificities, d_roc_nobal$specificities))


roc_aucs <- c(pROC::auc(d_roc_nobal), pROC::auc(d_roc))


ggplot(d_rocs,
       aes(x = 1 - spec, y = sens, colour = model)) +
  geom_line() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  annotate("text", x = c(0.75, 0.75), y = c(0.375, 0.25), 
           label = paste("AUC:", round(roc_aucs, digits = 3)),
           colour = pal[1:2]) +
  theme_bw() +
  scale_colour_manual(values = pal) +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = "RF Model") +
  theme(legend.position = "bottom",
        text = element_text(size = 12))
ggsave("plt_RF_ROC.png", device = png, width = 4, height = 4, bg = "white")

## Balancing by class weights (according to freuqency of adapted vs maladapted) 
## is really a trade off: increases the specificity but decreases sensitivity
### seems to be reducing overfitting, increases ability to classify maladapted
### cases at the cost of decreasing classification of adapted pops.
### ROC decreases though, so doesn't seem worth it?


# Plot errors (black line = OOB, red = false positive, green = false negative)
plot(rf_mbeta_adapted)


# Number of nodes per tree
hist(treesize(rf_mbeta_adapted),
     main = "# Nodes for the RF trees",
     col = "forestgreen")

# variable importance - model, btmb, and cossim Mb important
varImpPlot(rf_mbeta_adapted,
           sort = T, type = 1, scale = F)
randomForest::importance(rf_mbeta_adapted, type = 1, scale = F)

# conditional
iscores_crf <- varimp(crf_mbeta_adapted, conditional = T)
iscores_crf

saveRDS(iscores_crf, "iscores_crf.RDS")
# iscores_crf <- readRDS("iscores_crf.RDS")

iscores <- varImp(rf_mbeta_adapted, conditional = T)
iscores <- iscores %>% tibble::rownames_to_column("feature") 
iscores$feature <- iscores$feature %>% as.factor()
iscores <- iscores[,-3]
colnames(iscores) <- c("feature", "score")

bor_mbeta <- Boruta::Boruta(isAdapted ~ ., data = d_btgb_Malign_rf)
bor_mbeta
plot(bor_mbeta)

d_bor_mbeta <- process_the_Boruta_data(bor_mbeta)

feature_names <- c(TeX("Shadow Variable (min)", output = "character"),
                   TeX("Shadow Variable (mean)", output = "character"),
                   TeX("Shadow Variable (max)", output = "character"),
                   TeX("Recombination rate", output = "character"),
                   TeX("$G/\\beta$ alignment",
                       output = "character"),
                   TeX("$R_{max} / \\beta$ alignment", output = "character"),
                   TeX("G Evolvability", output = "character"),
                   TeX("$M/\\beta$ alignment",
                       output = "character"),
                   TeX("M Evolvability", output = "character"),
                   TeX("Motif", output = "character"))


ggplot(iscores %>%
         mutate(feature = factor(feature, levels = c("r", "absCS_Gb", "dataset", "bTGb", 
                                                     "absCS_Mb", "bTMb", "model"))),
       aes(x = feature, y = score)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = parse(text = feature_names[4:10]),
                   guide = guide_axis(n.dodge = 2)) + 
  theme_bw() +
  labs(x = "Feature", y = "Importance") +
  theme(text = element_text(size = 12)) -> plt_rf_perm_imp


pal_boruta <- c(rep("#00C0EA", times = 3),
                rep("#00A000", times = 7))

ggplot(d_bor_mbeta %>% pivot_longer(everything()) %>%
         mutate(x = fct_reorder(name, value, median)),
       aes(x = x, y = value, fill = x)) +
  geom_boxplot(show.legend = F, linewidth = 0.25) +
  theme_bw() +
  scale_fill_manual(values = pal_boruta) +
  scale_x_discrete(labels = parse(text = feature_names),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Boruta Importance") +
  theme(text = element_text(size = 12)) -> plt_boruta_imp
plt_boruta_imp
ggsave("plt_boruta_import_align.png", device = png, bg = "white",
       width = 12, height = 8)

predictor <- iml::Predictor$new(rf_mbeta_adapted_nobal, 
                                data = test_mbeta_adapted[, 2:8], 
                                y = test_mbeta_adapted$isAdapted)

# Need to set the option future globals maxsize
options(future.globals.maxSize = 3221225472)
imp <- iml::FeatureImp$new(predictor,
                           loss = "ce",
                           n.repetitions = 100)

ggplot(imp$results %>% 
         mutate(feature = factor(feature, levels = c("r", "absCS_Gb", "dataset", "bTGb", 
                                                     "absCS_Mb", "bTMb", "model"))),
       aes(x = feature, y = importance)) +
  geom_point() +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                width = 0.2) +
  scale_x_discrete(labels = parse(text = feature_names[4:10]),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Permutation Importance") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_perm_imp
plt_perm_imp
ggsave("plt_perm_feat_imp_align.png", device = png, width = 9, height = 5, bg = "white")

# Sobol MDA
rf_sob_mbeta_adapted <- sobolMDA::ranger(isAdapted ~ .,
                                         data = train_mbeta_adapted, num.trees = 500, 
                                         importance = "sobolMDA")

sob_mbeta_adapted <- rf_sob_mbeta_adapted$variable.importance
d_sob_mbeta_adapted <- data.frame(feature = names(sob_mbeta_adapted),
                                  sobelMDA = sob_mbeta_adapted)

d_sob_mbeta_adapted$feature <- factor(d_sob_mbeta_adapted$feature,
                                      levels = c("r",
                                                 "absCS_Gb",
                                                 "dataset",
                                                 "bTGb",
                                                 "absCS_Mb",
                                                 "bTMb",
                                                 "model"))

ggplot(d_sob_mbeta_adapted,
       aes(x = feature, y = sobelMDA)) +
  geom_point() +
  geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
               linewidth = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = parse(text = feature_names[4:10]),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Sobel MDA") +
  theme(text = element_text(size = 12)) -> plt_sob_mbeta_adapted
plt_sob_mbeta_adapted


layout <- "
AAAA
AAAA
AAAA
BBCC
BBCC
"
plt_boruta_imp +
  plt_perm_imp +
  plt_sob_mbeta_adapted +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

ggsave("plt_feat_imp_align.png", device = png, width = 12, height = 10, bg = "white")


ale <- FeatureEffects$new(predictor, grid.size = 10)
ale$plot()

ale_plots <- vector(mode = "list", length = 7)

ale_labels <- feature_names[c(10, 6, 4, 8, 5, 7, 9)]
names(ale_labels) <- names(ale$results)

for (i in seq_along(ale_plots)) {
  d_ale <- ale$results[[i]]
  d_ale <- d_ale %>% filter(.class == "Adapted")
  x_label <- ale_labels[d_ale$.feature[1]] 
  
  if (d_ale$.feature[1] == "model") {
    d_ale$.borders <- factor(d_ale$.borders, 
                             levels = model_names,
                             labels = model_names_noquote)
  }
  
  if (!is.numeric(d_ale$.borders[1])) {
    geom_fn <- geom_lollipop
    scale_fn <- scale_x_discrete
  } else {
    geom_fn <- geom_line
    scale_fn <- scale_x_continuous
  }
  ale_plots[[i]] <- ggplot(d_ale,
                           aes(x = .borders, y = .value)) +
    geom_fn() +
    scale_fn() +
    theme_bw() +
    labs(x = parse(text = x_label), y = "ALE of adaptation probability")
}

plot_grid(plotlist = ale_plots,
          labels= "AUTO")
ggsave("plt_ale_align.png", device = png, bg = "white",
       width = 12, height = 9)



# partial dependence
partialPlot(x = rf_mbeta_adapted,
            pred.data = test_mbeta_adapted,
            x.var = model)
partialPlot(x = rf_mbeta_adapted,
            pred.data = test_mbeta_adapted,
            x.var = absCS_Mb)
partialPlot(x = rf_mbeta_adapted,
            pred.data = test_mbeta_adapted,
            x.var = dataset)

MDSplot(rf_mbeta_adapted, test_mbeta_adapted$model,
        pal = pal)
legend("topright", legend = model_names_noquote,
       fill = pal)
ggsave("plt_rf_mbeta_adapted_mds_align.png", device = png, 
       width = 5, height = 5, bg = "white")



## Confirm populations evolving toward optima parallel to M adapt faster than random-direction populations (i.e. above)

## Populations orthogonal to M adapt most slowly and contain the least evolvability

## Difference between parallel and orthogonal is the cost of developmental constraint, compare among motifs


# 4) G and M under drift and selection

# Measure cosine similarity between M matrix correlations and G matrix correlations
## Match G matrix rows to M matrix rows
id_m$dataset <- "Randomised"
id_m_orth$dataset <- "Orthogonal"
id_m_par$dataset <- "Parallel"


id_m_for_g_match <- id_m_tot %>%
  mutate(m_idx = row_number()) %>%
  filter(timePoint == 50000 | timePoint == 60000) %>%
  mutate(timePoint = if_else(timePoint == 50000, "Start", "End"))

id_m_g_matched <- left_join(id %>%
                              select(timePoint, seed, modelindex,
                                     dataset, isAdapted), 
                            id_m_for_g_match %>%
                              select(timePoint, seed, modelindex,
                                     dataset, isAdapted, m_idx),
                            by = c("timePoint", "seed", "modelindex", "dataset",
                                   "isAdapted"))

id_m_g_matched <- AddCombosToDF(id_m_g_matched)

# covariance instead of correlation matrices
m_cov_for_g_match <- m_matrices_tot[id_m_g_matched$m_idx]

m_cov_for_g_match <- lapply(m_cov_for_g_match, function(x) {
  result <- matrix(0, 4, 4)
  result[1:nrow(x), 1:nrow(x)] <- x
  return(result)
})



## Angle between leading eigenvectors of M and G
d_cossim_m_vs_g <- GetCosineSimilarityTwoMats(h2_pd, m_cov_for_g_match, id_m_g_matched)
d_cossim_m_vs_g <- AddCombosToDF(d_cossim_m_vs_g)
d_cossim_m_vs_g$model <- factor(d_cossim_m_vs_g$model, levels = model_names)
d_cossim_m_vs_g$dataset <- factor(id_m_g_matched$dataset, levels = c("Orthogonal",
                                                                     "Parallel",
                                                                     "Randomised"))
                                    
d_cossim_m_vs_g_sum <- d_cossim_m_vs_g %>%
  group_by(model, isAdapted) %>%
  dplyr::summarise(meanCosSim = mean(abs(cosSim), na.rm = T),
                   seCosSim = se(abs(cosSim), na.rm = T),
                   meanbTGb = mean(bTMb, na.rm = T),
                   sebTGb = se(bTMb, na.rm = T))

# similarity between g_max and m_max
ggplot(d_cossim_m_vs_g, 
       aes(x = model, y = abs(cosSim), colour = model)) +
  facet_nested(. ~ "Population adapted" + isAdapted) +
  geom_quasirandom(show.legend = F) +
  geom_point(data = d_cossim_m_vs_g_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = model, y = meanCosSim, group = model), colour = "black",
             fill = "white", shape = 21, size = 2, stroke = 1, position = position_dodge(0.9)) +
  labs(x = "Model", y = TeX("$G_{max}/M_{max}$ alignment ($|cos(\\theta)|^M_G$)"),
       colour = "Model") +
  scale_x_discrete(labels = model_names_noquote) +
  scale_colour_manual(values = pal,
                      labels = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggsave("plt_cossim_gmax_mmax_tot.png", device = png, bg = "white",
       width = 8, height = 4)


# Amount of genetic variance shared along the m_max
ggplot(d_cossim_m_vs_g_sum, 
       aes(x = model, y = meanbTGb, colour = model)) +
  facet_nested(.~ "Population adapted" + isAdapted) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanbTGb - sebTGb, ymax = meanbTGb + sebTGb,
                    colour = model), show.legend = F,
                position = position_dodge(0.9)) +
  labs(x = "Model", 
       y = TeX("Evolvability along $m_{max} (m_{max}^T~G~m_{max})$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


## expect G = M under drift

m_matrices_tot <- c(m_matrices, m_matrices_orth, m_matrices_par)
id_m_tot <- rbind(id_m %>% select(-corMatIndex), id_m_orth, id_m_par)

# 5) Evolvability, autonomy and V_A through M
d_ecr_m <- CalcECRATrait(m_matrices_tot, id_m_tot)
d_ecr_m <- AddCombosToDF(d_ecr_m)


saveRDS(d_ecr_m, "d_ecr_m.RDS")
d_ecr_m <- readRDS("d_ecr_m.RDS")

# Refactor model
d_ecr_m <- d_ecr_m %>%
  mutate(model = factor(model, levels = model_names))

d_ecr_m$dataset <- id_m_tot$dataset

# d_ecr_m_plt <- d_ecr_m %>%
#   filter(timePoint == 50000 | timePoint == 60000) %>%
#   mutate(timePoint = if_else(timePoint == 50000, "Start", "End"),
#          timePoint = factor(timePoint, levels = c("Start", "End")))

# Need to calculate cev means separately for the different models
# K- shouldn't mean over cev_KZ and KXZ
d_ecr_m_sum <- d_ecr_m %>%
  filter(!is.nan(cev)) %>%
  mutate(timePoint = timePoint - 50000) %>%
  group_by(timePoint, model, isAdapted) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr_m_sum %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = cev_mean, colour = model)) +
  facet_nested("Model" + model~"Population adapted" + isAdapted,
               scales = "free", labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = cev_mean - cev_se, ymax = cev_mean + cev_se, 
                  fill = model), alpha = 0.2, colour = NA, show.legend = F) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(breaks = equal_breaks(3, 0.05),
                     labels = scales::label_scientific(digits = 2)) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post-optimum shift", y = "Mean conditional evolvability (M matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        panel.spacing.y = unit(1, "lines"),
        text = element_text(size = 10)) -> plt_cev_m
plt_cev_m

ggplot(d_ecr_m_sum %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = res_mean, colour = model)) +
  facet_nested("Model" + model~"Population adapted" + isAdapted,
               scales = "free", labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = res_mean - res_se, ymax = res_mean + res_se, 
                  fill = model), alpha = 0.2, colour = NA, show.legend = F) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(breaks = equal_breaks(3, 0.05),
                     labels = scales::label_scientific(digits = 2)) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post-optimum shift", y = "Mean respondability (M matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        panel.spacing.y = unit(1, "lines"),
        text = element_text(size = 10)) -> plt_res_m
plt_res_m

ggplot(d_ecr_m_sum %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = aut_mean, colour = model)) +
  facet_nested("Model" + model~"Population adapted" + isAdapted,
               scales = "free", labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = aut_mean - aut_se, ymax = aut_mean + aut_se, 
                  fill = model), alpha = 0.2, colour = NA, show.legend = F) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(breaks = equal_breaks(3, 0.05),
                     labels = scales::label_scientific(digits = 2)) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post-optimum shift", y = "Mean autonomy (M matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        panel.spacing.y = unit(1, "lines"),
        text = element_text(size = 10)) -> plt_aut_m
plt_aut_m


ggplot(d_ecr_m_sum %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = timePoint, y = ev_mean, colour = model)) +
  facet_nested("Model" + model~"Population adapted" + isAdapted,
               scales = "free", labeller = labeller(model = model_names_labeller)) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin = ev_mean - ev_se, ymax = ev_mean + ev_se, 
                  fill = model), alpha = 0.2, colour = NA, show.legend = F) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(breaks = equal_breaks(3, 0.05),
                     labels = scales::label_scientific(digits = 2)) +
  labs(x = "Generations post-optimum shift", y = "Mean evolvability (M matrix)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(-5, 0, 0, 0),
        panel.spacing.y = unit(1, "lines"),
        text = element_text(size = 10)) -> plt_ev_m

plt_ev_m

leg <- get_legend(plt_ev_m)

plt_evol_m <- plot_grid(plt_ev_m + theme(legend.position = "none"),
                      plt_cev_m + theme(legend.position = "none"),
                      plt_res_m + theme(legend.position = "none"),
                      plt_aut_m + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO", label_size = 12)

plt_evol_m <- plot_grid(plt_evol_m,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol_m
ggsave("plt_evol_m.png", device = png, bg = "white",
       width = 10, height = 8) 

## Variation in the M matrix (across all selection gradients) differs between motifs
## All models see increases in conditional evolvability in adapted vs maladapted pops
## NAR, PAR, and FFLI1 models behave similarly in the other three: increases to evolvability
## and respondability in adapted vs maladapted, decrease in autonomy (leveraging mutational
## correlations, interesting). FFLC1 and FFBH behave similarly in the opposite direction:
## adapted pops have reduced evolvability, increased autonomy and reduced respondability

## increase in autonomy suggests that trait correlations severely limit adaptation in those
## two, there is a requirement to reduce developmental constraints to increase the conditional
## evolvability. A side effect is a decrease in net evolvability and respondability however.
## how does this occur? Adapted populations have an M matrix which is uncorrelated, but can only make
## small changes as a result? Whereas for other populations, correlations are less problematic,
## so larger pleiotropic effects are tolerable, increasing overall evolvability and conditional
## ev without strongly hampering fitness?

## This is average across all directions



## FFBH had high V_A but low adaptation - was VA concentrated along M's leading 
# eigenvectors and misaligned with selection?

## Conditional evolvability in the direction of selection predicts adaptive success

# G matrix Vrel
e_g <- lapply(h2_pd, eigen)
saveRDS(e_g, "eigen_randomised_g_tot.RDS")

vrel_g <- unlist(lapply(e_g, function(x) { Vrel(x$values) }))

# join with quant gen data
d_vrel_g <- left_join(id %>% mutate(seed = factor(seed),
                                   modelindex = factor(modelindex),
                                   r = RFromIndex(modelindex),
                                   vrel = vrel_g) %>%
                        select(-label),
                    d_qg_tot %>% select(gen, seed, modelindex, dataset, isAdapted, model, r) %>%
                      filter(gen == 50000 | gen == 60000) %>%
                      mutate(timePoint = if_else(gen == 50000, "Start", "End"),
                             timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
                      select(-gen),
                    by = c("timePoint", "seed", "modelindex", "dataset", "isAdapted", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names,
                        labels = model_names_noquote))

# Plot
d_vrel_sum <- d_vrel_g %>%
  group_by(timePoint, model, dataset, isAdapted) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))

ggplot(d_vrel_sum,
       aes(x = timePoint, y = vrel_mean, colour = model)) +
  facet_nested("Dataset" + dataset~"Population adapted" + isAdapted) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = vrel_mean - vrel_CI, ymax = vrel_mean + vrel_CI, colour = model),
              show.legend = F, position = position_dodge(0.9)) +
  labs(x = "Generations post-optimum shift", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom") 



# Stable over time (except FFLC1, mirrors M evolvability), can take average?
d_vrel_g_tab <- d_vrel_g %>%
  group_by(model, dataset, isAdapted) %>%
  summarise(vrel_mean = mean(vrel),
            vrel_CI = CI(vrel))
d_vrel_g_tab

ggplot(d_vrel_g,
       aes(x = model, y = vrel, colour = model)) +
  facet_nested("Selection/trait correlation alignment" + dataset~"Population adapted" + isAdapted) +
  geom_boxplot() +
  geom_point(data = d_vrel_g_tab, aes(y = vrel_mean), size = 2, stroke = 1,
             shape = 21,
             colour = "black", fill = "white") +
  labs(x = "Model", y = TeX("$V_{rel}$"), colour = "Model") +
  scale_colour_manual(values = pal,
                      breaks = model_names_noquote) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_vrel_g
plt_vrel_g

ggsave("plt_vrel_g_alignment.png", device = png, bg = "white",
       width = 8, height = 6)

# Combine with M matrix version

plot_grid(plt_vrel_g + labs(y = TeX("$V_{rel}$ (G)")),
          plt_vrel_m + labs(y = TeX("$V_{rel}$ (M)")),
          labels = "AUTO",
          nrow = 2)
ggsave("plt_vrel_alignment.png", device = png, bg = "white",
       width = 8, height = 10)
