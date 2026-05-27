# Random forest model of adaptedness
library(tidyverse)
library(paletteer)
library(Rcpp)
library(ggh4x)
library(latex2exp)
library(randomForest)
library(pROC)
library(patchwork)
library(ggalt)
library(cowplot)

source("helperFn.R")

# Load in data from the selection experiments

# Combos
COMBO_PATH <- '/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
COMBO_PATH <- '/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
d_combos <- read_delim(COMBO_PATH, 
                       delim = " ", col_names = F)
names(d_combos) <- c("model", "r")


PATH_QG <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_qg.csv"
PATH_QG <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_qg.csv"

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




# Load in data
PATH_QG_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"
PATH_QG_ORTH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"

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
PATH_QG_PAR <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_qg.csv"

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


#########
# Mutation data


DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"
DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar.csv"
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar.csv"
DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar.csv"

d_m <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                         paste0("mean_", 1:4),
                                         paste0("var_", 1:4),
                                         paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m <- d_m %>%
  mutate(model = ModelFromIndexWithR(modelindex))


# get matrices
m_matrices <- d_m %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))

#saveRDS(m_matrices, "m_matrices.RDS")
m_matrices <- readRDS("m_matrices.RDS")


# Get eigenvectors of each M
#e_m <- lapply(m_matrices, eigen)
#saveRDS(e_m, "eigen_randomised_m.RDS")

e_m <- readRDS("eigen_randomised_m.RDS")


DATA_PATH_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_mutvar.csv"
DATA_PATH_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_mutvar.csv"

DATA_PATH_ORTH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/orthSel/slim_mutvar.csv"
DATA_PATH_PAR <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_mutvar.csv"


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

d_m$dataset <- "Randomised"
d_m_orth$dataset <- "Orthogonal"
d_m_par$dataset <- "Parallel"

d_m_tot <- rbind(d_m, d_m_orth, d_m_par)


# get matrices
m_matrices_orth <- d_m_orth %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))
# 
m_matrices_par <- d_m_par %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))


m_matrices_tot <- c(m_matrices, m_matrices_orth, m_matrices_par)

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
vrel_m <- unlist(lapply(e_m, function(x) { Vrel(x$values) }))
vrel_m_orth <- unlist(lapply(e_m_orth, function(x) { Vrel(x$values) }))
vrel_m_par <- unlist(lapply(e_m_par, function(x) { Vrel(x$values) }))

# join with quant gen data
d_vrel <- left_join(d_m %>% mutate(seed = factor(seed),
                                   modelindex = factor(modelindex),
                                   r = RFromIndex(modelindex),
                                   vrel = vrel_m),
                    d_qg %>% select(gen, seed, modelindex, model, r, isAdapted) %>%
                      mutate(model = str_remove_all(model, "'")),
                    by = c("gen", "seed", "modelindex", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote))


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

# Combine
d_vrel_tot <- rbind(d_vrel, d_vrel_orth, d_vrel_par)


# Now the same for G matrices
# G matrix Vrel
e_g <- lapply(h2_pd, eigen)
saveRDS(e_g, "eigen_randomised_g_tot.RDS")

vrel_g <- unlist(lapply(e_g, function(x) { Vrel(x$values) }))

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



######################################################################
# Load optima
d_opt <- read_csv("/mnt/j/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_opt.csv", col_names = F)
d_opt_par <- read_csv("/mnt/j/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_opt.csv", col_names = F)
d_opt_orth <- read_csv("/mnt/j/SLiMTests/tests/newMotifs/paper1/orthSel/slim_opt.csv", col_names = F)

# o = optimum, s = sigma, d = direction (-1, 1)
colnames(d_opt) <- c("seed", "modelindex", "o_t1", "o_t2", "o_t3", "o_t4", 
                         "s_t1", "s_t2", "s_t3", "s_t4", "d_t1", "d_t2", "d_t3",
                         "d_t4")
colnames(d_opt_par) <- colnames(d_opt)
colnames(d_opt_orth) <- colnames(d_opt_par)

# Last value in d_opt_par and _orth is the eigenvalue of eigenvector r_i, we can ignore it
d_opt_par <- d_opt_par[,-15]
d_opt_orth <- d_opt_orth[,-15]





# Now G matrix
G_DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/R/"
G_DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/R/"
G_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/"
G_DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/"


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


G_ORTH_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"
G_ORTH_DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"


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
sourceCpp("/mnt/e/Documents/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")

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

# First convert to nearest positive definite matrix
h2_pd <- lapply(h2_mat, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
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
  mutate(t1_dir = o_t1 - trait1_mean,
         t2_dir = o_t2 - trait2_mean,
         t3_dir = o_t3 - trait3_mean,
         t4_dir = o_t4 - trait4_mean,
         norm = sqrt(rowSums(pick(ends_with("dir"))^2)), # normalise
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


###########################################################################
# M matrix
d_selvec_m_tot <- d_qg_tot %>%
  filter(gen >= 50000) %>%
  select(gen, seed, modelindex, dataset, isAdapted, ends_with("mean"))

d_selvec_m_tot <- left_join(d_selvec_m_tot, d_opt_tot %>% 
                          select(seed, modelindex, dataset, starts_with("o_")) %>%
                          mutate(seed = factor(seed),
                                 modelindex = factor(modelindex)), 
                        by = c("seed", "modelindex", "dataset"))


d_selvec_m_tot <- AddCombosToDF(d_selvec_m_tot)

d_selvec_m_tot <- d_selvec_m_tot %>%
  mutate(modelindex = as.factor(modelindex),
         seed = as.factor(seed),
         model = factor(model, levels = model_names),
         dataset = factor(dataset, levels = c("Orthogonal", "Parallel", "Randomised"))) %>%
  rename(timePoint = gen)

id_m_tot <- d_m_tot %>% mutate(timePoint = gen) %>% select(timePoint, seed, modelindex, dataset)
id_m_tot$clus <- 1
id_m_tot$modelindex <- as.factor(id_m_tot$modelindex)
id_m_tot$seed <- as.factor(id_m_tot$seed)

id_m_tot <- AddCombosToDF(id_m_tot)
id_m_tot$model <- factor(id_m_tot$model, levels = model_names)

id_m_tot <- inner_join(id_m_tot, d_qg_tot %>% mutate(timePoint = gen) %>% 
                     select(timePoint, seed, modelindex, dataset, isAdapted),
                   by = c("timePoint", "seed", "modelindex", "dataset"))

d_selvec_m_tot <- inner_join(id_m_tot, d_selvec_m_tot, 
                         by = c("timePoint", "seed", "modelindex", "dataset", "isAdapted", "model", "r"))


d_selvec_m_tot <- d_selvec_m_tot %>%
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


d_cossim_m_tot <- GetCosineSimilarity(m_matrices_tot, d_selvec_m_tot %>% select(ends_with("dir")), id_m_tot)
d_cossim_m_tot$dataset <- id_m_tot$dataset

saveRDS(d_cossim_m_tot, "d_cossim_m_datasets.RDS")
d_cossim_m_tot <- readRDS("d_cossim_m_datasets.RDS")

d_cossim_m_tot <- AddCombosToDF(d_cossim_m_tot)


## Evolvability against alignment of M with direction of selection
## Evolvability = bTGb

# Join the evolvability estimates with M alignment
d_btgb_Malign_tot <- left_join(d_cossim %>% 
                                 select(timePoint, seed, modelindex, dataset, isAdapted,
                                        model, r, bTMb, cosSim) %>%
                                 rename(bTGb = bTMb) %>%
                                 mutate(absCS_Gb = abs(cosSim)) %>%
                                 select(-cosSim),
                               d_cossim_m_tot %>%
                                 filter(timePoint == 60000 | timePoint == 50000) %>%
                                 mutate(timePoint = if_else(timePoint == 50000, "Start", "End"),
                                        timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
                                 select(timePoint, seed, modelindex, dataset, isAdapted,
                                        model, r, bTMb, cosSim) %>%
                                 mutate(absCS_Mb = abs(cosSim)) %>%
                                 select(-cosSim),
                               by = c("timePoint", "seed", "modelindex", "dataset", "isAdapted",
                                      "model", "r"))

## Summary
d_btgb_Malign_sum <- d_btgb_Malign_tot %>%
  group_by(model, dataset, isAdapted) %>%
  summarise(meanCosSim_Gb = mean(absCS_Gb),
            CICosSim_Gb = CI(absCS_Gb),
            meanCosSim_Mb = mean(absCS_Mb),
            CICosSim_Mb = CI(absCS_Mb),
            meanbTGb = mean(bTGb),
            CIbTGb = CI(bTGb))

### Alignment vs Selection probability
d_isAdapted <- d_h2_trait %>%
  filter(gen == 50000) %>%
  select(seed, modelindex, model, dataset, r, isAdapted) %>%
  distinct(.keep_all = T)

# Dataset does change the adaptation outcome
fisher.test(table(d_isAdapted$dataset, d_isAdapted$isAdapted))

d_isAdapted <- as.data.frame(table(d_isAdapted$model, d_isAdapted$dataset, d_isAdapted$isAdapted))
d_isAdapted <- d_isAdapted %>%
  rename(model = Var1,
         dataset = Var2,
         isAdapted = Var3,
         count = Freq) %>%
  group_by(model, dataset) %>%
  mutate(groupTotal = sum(count)) %>%
  ungroup() %>%
  mutate(Freq = count / groupTotal)

# Check we have 208 replicates * 3 r levels = 624 per model
# plus the orth and parallel sets of another 208 replicates = 1040 per model

# join isAdapted counts
d_btgb_Malign_sum <- left_join(d_btgb_Malign_sum %>% 
                                 mutate(model = factor(model, levels = model_names),
                                        isAdapted = factor(isAdapted)),
                               d_isAdapted,
                               by = c("model", "isAdapted", "dataset"))

ggplot(d_btgb_Malign_sum %>% 
         filter(isAdapted == T),
       aes(x = meanCosSim_Mb, y = Freq, colour = model)) +
  #geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  facet_nested(.~"Selection/trait correlation alignment" + dataset) +
  geom_point(shape = 1) +
  geom_errorbar(aes(xmin = meanCosSim_Mb - CICosSim_Mb, 
                    xmax = meanCosSim_Mb + CICosSim_Mb)) +
  labs(x = TeX("$M_{max} / \\beta$ alignment ($abs(cos(\\theta)_\\beta^{M})$)"), 
       y = "Proportion of populations adapted",
       colour = "Model") +
  coord_flip(ylim = c(0, 1), xlim = c(0,1)) + 
  scale_colour_manual(values = pal,
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
ggsave("plt_Mcossim_probAdapt_align.png", device = png, width = 6, height = 6, bg = "white")


# Does M/beta alignment predict population adaptedness?

# Add Vrel
d_btgb_Malign_tot_vrel <- left_join(d_btgb_Malign_tot %>%
                                      mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                                                labels = c("Adapted", "Maladapted")),
                                             model = factor(model, levels = model_names, labels = model_names_noquote),
                                             timePoint = factor(timePoint, levels = c("Start", "End")),
                                             dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                                             r = factor(log10(r), levels = c(-10, -5, -1))),
                                    d_vrel_tot %>%
                                      filter(gen == 50000 | gen == 60000) %>%
                                      mutate(timePoint = if_else(gen == 50000, "Start", "End"),
                                             isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                                                labels = c("Adapted", "Maladapted")),
                                             model = factor(model, levels = model_names_noquote),
                                             timePoint = factor(timePoint, levels = c("Start", "End")),
                                             dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                                             r = factor(log10(r), levels = c(-10, -5, -1))) %>%
                                      select(timePoint, seed, modelindex, model, dataset, isAdapted, r,
                                             vrel),
                                    by = c("timePoint", "seed", "modelindex",
                                           "model", "dataset", "isAdapted", "r"))

# Add vrel for G matrix
d_btgb_Malign_tot_vrel <- left_join(d_btgb_Malign_tot_vrel %>%
                                      rename(vrel_m = vrel),
                                    d_vrel_g %>%
                                      rename(vrel_g = vrel) %>%
                                      mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                                                labels = c("Adapted", "Maladapted")),
                                             dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                                             r = factor(log10(r), levels = c(-10, -5, -1))) %>%
                                      select(timePoint, seed, modelindex, model, dataset, isAdapted, r,
                                             vrel_g),
                                    by = c("timePoint", "seed", "modelindex",
                                           "model", "dataset", "isAdapted", "r"))



## use random forest
d_btgb_Malign_rf <- d_btgb_Malign_tot_vrel %>%
  select(isAdapted, model, dataset, r, absCS_Mb, absCS_Gb, bTGb, bTMb, vrel_g, vrel_m)

saveRDS(d_btgb_Malign_rf, "d_btgb_Malign_rf.RDS")



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
rf_mbeta_adapted_nobal <- randomForest(formula = isAdapted ~ model * dataset * r * absCS_Gb * absCS_Mb * 
                                         bTGb * bTMb * vrel_g * vrel_m,
                                       data = train_mbeta_adapted,
                                       ntree = 500,
                                       proximity = T,
                                       importance = T,
                                       type = "classification")

print(rf_mbeta_adapted_nobal)

# With balancing (class weights)
rf_mbeta_adapted_bal <- randomForest(formula = isAdapted ~ model * dataset * r * absCS_Gb * absCS_Mb * 
                                       bTGb * bTMb * vrel_g * vrel_m,
                                 data = train_mbeta_adapted,
                                 strata = train_mbeta_adapted$isAdapted,
                                 classwt = adapted_weights,
                                 ntree = 500,
                                 proximity = T,
                                 importance = T,
                                 type = "classification")

print(rf_mbeta_adapted_bal)

# Training data
p_train_mbeta_adapted_bal <- predict(rf_mbeta_adapted_bal, train_mbeta_adapted)
caret::confusionMatrix(p_train_mbeta_adapted_bal, train_mbeta_adapted$isAdapted)

p_train_mbeta_adapted_nobal <- predict(rf_mbeta_adapted_nobal, train_mbeta_adapted)
caret::confusionMatrix(p_train_mbeta_adapted_nobal, train_mbeta_adapted$isAdapted)


# Test data
p_test_mbeta_adapted_bal <- predict(rf_mbeta_adapted_bal, test_mbeta_adapted)
p_test_mbeta_adapted_bal_probs <- predict(rf_mbeta_adapted_bal, test_mbeta_adapted,
                                      type = "prob")[,1]

caret::confusionMatrix(p_test_mbeta_adapted_bal, test_mbeta_adapted$isAdapted)

p_test_mbeta_adapted_nobal <- predict(rf_mbeta_adapted_nobal, test_mbeta_adapted)
p_test_mbeta_adapted_nobal_probs <- predict(rf_mbeta_adapted_nobal, test_mbeta_adapted,
                                            type = "prob")[,1]

caret::confusionMatrix(p_test_mbeta_adapted_nobal, test_mbeta_adapted$isAdapted)



# roc
d_roc_bal <- roc(response = test_mbeta_adapted$isAdapted,
             predictor = p_test_mbeta_adapted_bal_probs)

d_roc_nobal <- roc(response = test_mbeta_adapted$isAdapted,
                   predictor = p_test_mbeta_adapted_nobal_probs,
                   levels = rev(levels(test_mbeta_adapted$isAdapted)))

d_rocs <- data.frame(model = c(rep("Weighted", times = length(d_roc_bal$sensitivities)),
                               rep("Unbalanced", times = length(d_roc_nobal$sensitivities))),
                     sens = c(d_roc_bal$sensitivities, d_roc_nobal$sensitivities),
                     spec = c(d_roc_bal$specificities, d_roc_nobal$specificities))


roc_aucs <- c(pROC::auc(d_roc_nobal), pROC::auc(d_roc_bal))


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
plot(rf_mbeta_adapted_bal)
plot(rf_mbeta_adapted_nobal)


# Number of nodes per tree
hist(treesize(rf_mbeta_adapted_nobal),
     main = "# Nodes for the RF trees",
     col = "forestgreen")

#############################
# Importance measures
## Boruta, permutation importance, sobol MDA

bor_mbeta <- Boruta::Boruta(isAdapted ~ ., data = d_btgb_Malign_rf)
bor_mbeta
plot(bor_mbeta)

d_bor_mbeta <- process_the_Boruta_data(bor_mbeta)

feature_names <- c(TeX("Shadow Variable (min)", output = "character"),
                   TeX("Shadow Variable (mean)", output = "character"),
                   TeX("Shadow Variable (max)", output = "character"),
                   TeX("Recombination rate", output = "character"),
                   TeX("$V_{rel}$ (G)", output = "character"),
                   TeX("$G/\\beta$ alignment",
                       output = "character"),
                   TeX("$R_{max} / \\beta$ alignment", output = "character"),
                   TeX("G Evolvability", output = "character"),
                   TeX("$M/\\beta$ alignment",
                       output = "character"),
                   TeX("$V_{rel}$ (M)", output = "character"),
                   TeX("M Evolvability", output = "character"),
                   TeX("Motif", output = "character"))


pal_boruta <- c(rep("#00C0EA", times = 3),
                rep("#00A000", times = 9))

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


# Permutation
predictor <- iml::Predictor$new(rf_mbeta_adapted_nobal, 
                                data = test_mbeta_adapted[, 2:10], 
                                y = test_mbeta_adapted$isAdapted)

# Need to set the option future globals maxsize
options(future.globals.maxSize = 3221225472)
imp <- iml::FeatureImp$new(predictor,
                           loss = "ce",
                           n.repetitions = 100)

ggplot(imp$results %>% 
         mutate(feature = factor(feature, levels = c("r", "vrel_g", "absCS_Gb", "dataset", "bTGb", 
                                                     "absCS_Mb", "vrel_m", "bTMb", "model"))),
       aes(x = feature, y = importance)) +
  geom_point() +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                width = 0.2) +
  scale_x_discrete(labels = parse(text = feature_names[4:12]),
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
                                                 "vrel_g",
                                                 "absCS_Gb",
                                                 "dataset",
                                                 "bTGb",
                                                 "absCS_Mb",
                                                 "vrel_m",
                                                 "bTMb",
                                                 "model"))

ggplot(d_sob_mbeta_adapted,
       aes(x = feature, y = sobelMDA)) +
  geom_point() +
  geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
               linewidth = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = parse(text = feature_names[4:12]),
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


# Accumulated local effects
ale <- FeatureEffects$new(predictor, grid.size = 10)
ale$plot()

ale_plots <- vector(mode = "list", length = 9)

ale_labels <- feature_names[c(12, 7, 4, 9, 6, 8, 11, 5, 10)]
names(ale_labels) <- names(ale$results)

for (i in seq_along(ale_plots)) {
  d_ale <- ale$results[[i]]
  d_ale <- d_ale %>% filter(.class == "Adapted")
  x_label <- ale_labels[d_ale$.feature[1]] 
  
  if (d_ale$.feature[1] == "model") {
    d_ale$.borders <- factor(d_ale$.borders, 
                             levels = model_names_noquote,
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

