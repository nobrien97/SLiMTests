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
library(ggbeeswarm)
library(xgboost)
library(caret)
library(nlme)
library(emmeans)

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
PATH_QG <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_qg.csv"

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
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98),
         timeToAdapt = first(gen[gen > 50000 & mean_w > 0.98]) - 50000) %>%
  mutate(model = factor(model, levels = model_names)) %>%
  ungroup() -> d_qg




# Load in data
PATH_QG_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"
PATH_QG_ORTH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"

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
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98),
         timeToAdapt = first(gen[gen > 50000 & mean_w > 0.98]) - 50000) %>%
  mutate(model = factor(model, levels = model_names),
         modelindex = factor(modelindex, levels = 1:15)) %>%
  ungroup() -> d_qg_orth


PATH_QG_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_qg.csv"
PATH_QG_PAR <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_qg.csv"

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
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98),
         timeToAdapt = first(gen[gen > 50000 & mean_w > 0.98]) - 50000) %>%
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

DATA_PATH_ORTH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/orthSel/slim_mutvar.csv"
DATA_PATH_PAR <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_mutvar.csv"


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



######################################################################
# Load optima
d_opt <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_opt.csv", col_names = F)
d_opt_par <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_opt.csv", col_names = F)
d_opt_orth <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_opt.csv", col_names = F)

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
G_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/"


d_h2_mrr <- read_csv(paste0(G_DATA_PATH, "out_h2_mrr.csv"), col_names = F)
d_h2_mkr <- read_csv(paste0(G_DATA_PATH, "out_h2_mkr.csv"), col_names = F)

colnames(d_h2_mrr) <- h2_colnames
colnames(d_h2_mkr) <- h2_colnames

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


d_h2_mkr$calcMode <- "mkr"
d_h2_mrr$calcMode <- "mrr"


G_ORTH_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"
G_ORTH_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"


d_h2_mkr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_mkr.csv"), col_names = F)
d_h2_mrr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_mrr.csv"), col_names = F)

d_h2_mkr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_mkr.csv"), col_names = F)
d_h2_mrr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_mrr.csv"), col_names = F)


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

colnames(d_h2_mrr_par) <- h2_colnames
colnames(d_h2_mkr_par) <- h2_colnames
colnames(d_h2_mrr_orth) <- h2_colnames
colnames(d_h2_mkr_orth) <- h2_colnames


# join
d_h2_mkr_orth$calcMode <- "mkr"
d_h2_mrr_orth$calcMode <- "mrr"
d_h2_mkr_par$calcMode <- "mkr"
d_h2_mrr_par$calcMode <- "mrr"

d_h2_trait_mkr_orth$calcMode <- "mkr"
d_h2_trait_mrr_orth$calcMode <- "mrr"
d_h2_trait_mkr_par$calcMode <- "mkr"
d_h2_trait_mrr_par$calcMode <- "mrr"

d_h2_mkr_orth$dataset <- "Orthogonal"
d_h2_mrr_orth$dataset <- "Orthogonal"
d_h2_mkr_par$dataset <- "Parallel"
d_h2_mrr_par$dataset <- "Parallel"
d_h2_mkr$dataset <- "Randomised"
d_h2_mrr$dataset <- "Randomised"

d_h2_trait_mkr_orth$dataset <- "Orthogonal"
d_h2_trait_mrr_orth$dataset <- "Orthogonal"
d_h2_trait_mkr_par$dataset <- "Parallel"
d_h2_trait_mrr_par$dataset <- "Parallel"
d_h2_trait_mkr$dataset <- "Randomised"
d_h2_trait_mrr$dataset <- "Randomised"

d_h2 <- rbind(d_h2_mkr, d_h2_mrr,
              d_h2_mkr_orth, d_h2_mrr_orth,
              d_h2_mkr_par, d_h2_mrr_par)

d_h2_trait <- rbind(d_h2_trait_mkr, d_h2_trait_mrr, 
                    d_h2_trait_mkr_orth, d_h2_trait_mrr_orth,
                    d_h2_trait_mkr_par, d_h2_trait_mrr_par)

d_h2_trait %>% mutate(model = d_combos$model[.$modelindex],
                      model = factor(model, levels = model_names),
                      r = d_combos$r[.$modelindex]) -> d_h2_trait

d_h2 %>% mutate(model = d_combos$model[.$modelindex],
                model = factor(model, levels = model_names),
                r = d_combos$r[.$modelindex]) -> d_h2

d_h2_trait <- d_h2_trait %>%
  distinct(gen, seed, modelindex, dataset, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()

d_h2 <- d_h2 %>%
  distinct(gen, seed, modelindex, dataset, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()


# Join qg together
d_qg$dataset <- "Randomised"
d_qg_orth$dataset <- "Orthogonal"
d_qg_par$dataset <- "Parallel"

d_qg_tot <- rbind(d_qg, d_qg_orth, d_qg_par)

saveRDS(d_qg_tot, "d_qg_tot.RDS")

d_qg_optPerc <- d_qg_tot %>% select(gen, seed, modelindex, dataset, isAdapted) %>% filter(gen >= 49500)

# inner join optPerc
d_h2_trait <- left_join(d_h2_trait, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))

d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))


# Counts for each model type:
table(d_h2_trait$model, d_h2_trait$isAdapted)
table(d_h2_trait$model, d_h2_trait$dataset, d_h2_trait$isAdapted)

table(d_h2$model, d_h2$isAdapted)
table(d_h2$model, d_h2$dataset, d_h2$isAdapted)

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
#saveRDS(d_opt_tot, "/mnt/i/SLiMTests/tests/newMotifs/paper1/d_opt_tot.RDS")

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

#saveRDS(d_selvec, "/mnt/i/SLiMTests/tests/newMotifs/paper1/d_selvec.RDS")

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

# Vrel for G matrices
# Now the same for G matrices
# G matrix Vrel
e_g <- lapply(h2_pd, eigen)
saveRDS(e_g, "eigen_randomised_g_tot.RDS")

absCS_Mb <- unlist(lapply(e_g, function(x) { Vrel(x$values) }))

d_absCS_Mb <- left_join(id %>% mutate(seed = factor(seed),
                                    modelindex = factor(modelindex),
                                    r = RFromIndex(modelindex),
                                    vrel = absCS_Mb) %>%
                        select(-label),
                      d_qg_tot %>% select(gen, seed, modelindex, dataset, isAdapted, model, r) %>%
                        filter(gen == 50000 | gen == 60000) %>%
                        mutate(timePoint = if_else(gen == 50000, "Start", "End"),
                               timePoint = factor(timePoint, levels = c("Start", "End"))) %>%
                        select(-gen),
                      by = c("timePoint", "seed", "modelindex", "dataset", "isAdapted", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names,
                        labels = model_names_noquote))

d_ecr <- CalcECRATrait(h2_pd, id)
d_ecr <- AddCombosToDF(d_ecr)
d_ecr$dataset <- d_selvec2$dataset

# Refactor model
d_ecr <- d_ecr %>%
  mutate(model = factor(model, levels = model_names))



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

# 5) Evolvability, autonomy and V_A through M
d_ecr_m <- CalcECRATrait(m_matrices_tot, id_m_tot)
d_ecr_m <- AddCombosToDF(d_ecr_m)


saveRDS(d_ecr_m, "d_ecr_m.RDS")
d_ecr_m <- readRDS("d_ecr_m.RDS")

# Refactor model
d_ecr_m <- d_ecr_m %>%
  mutate(model = factor(model, levels = model_names))

d_ecr_m$dataset <- id_m_tot$dataset



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
                                    d_absCS_Mb %>%
                                      rename(absCS_Mb = vrel) %>%
                                      mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                                                labels = c("Adapted", "Maladapted")),
                                             dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                                             r = factor(log10(r), levels = c(-10, -5, -1))) %>%
                                      select(timePoint, seed, modelindex, model, dataset, isAdapted, r,
                                             absCS_Mb),
                                    by = c("timePoint", "seed", "modelindex",
                                           "model", "dataset", "isAdapted", "r"))

# Add on conditional evolvability
d_cev <- left_join(d_ecr %>%
                     select(timePoint, seed, modelindex, isAdapted, model, r, dataset, cev) %>%
                     rename(cev_g = cev) %>%
                     mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                               labels = c("Adapted", "Maladapted")),
                            dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                            r = factor(log10(r), levels = c(-10, -5, -1)),
                            model = factor(model, levels = model_names,
                                           labels = model_names_noquote)),
                   d_ecr_m %>% filter(timePoint == 50000 | timePoint == 60000) %>%
                     mutate(timePoint = if_else(timePoint == 50000, "Start", "End"),
                            isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                               labels = c("Adapted", "Maladapted")),
                            model = factor(model, levels = model_names,
                                           labels = model_names_noquote),
                            timePoint = factor(timePoint, levels = c("Start", "End")),
                            dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                            r = factor(log10(r), levels = c(-10, -5, -1))
                            ) %>%
                     select(timePoint, seed, modelindex, isAdapted, model, r, dataset, cev) %>%
                     rename(cev_m = cev),
                   by = c("timePoint", "seed", "modelindex", "isAdapted", "model", "r",
                          "dataset")
                   )



d_btgb_Malign_tot_vrel <- left_join(d_btgb_Malign_tot_vrel,
                                    d_cev,
                                    by = c("timePoint", "seed", "modelindex", "isAdapted", "model", "r",
                                           "dataset"))

# Add time to adaptation
d_btgb_Malign_tot_vrel <- left_join(d_btgb_Malign_tot_vrel,
                                    d_qg_tot %>% filter(gen == 50000 | gen == 60000) %>%
                                      mutate(timePoint = if_else(gen == 50000, "Start", "End"),
                                             isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                                                labels = c("Adapted", "Maladapted")),
                                             model = factor(model, levels = model_names,
                                                            labels = model_names_noquote),
                                             timePoint = factor(timePoint, levels = c("Start", "End")),
                                             dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")),
                                             r = factor(log10(r), levels = c(-10, -5, -1))
                                      ) %>% select(timePoint, seed, modelindex, dataset, timeToAdapt),
                                    by = c("timePoint", "seed", "modelindex", "dataset"))

# Save
saveRDS(d_btgb_Malign_tot_vrel, "d_btgb_Malign_tot_vrel.RDS")
d_btgb_Malign_tot_vrel <- readRDS("d_btgb_Malign_tot_vrel.RDS")
d_btgb_Malign_tot_vrel <- readRDS("/mnt/i/SLiMTests/tests/newMotifs/paper1/d_btgb_Malign_tot_vrel.RDS")

## use random forest
d_btgb_Malign_rf <- d_btgb_Malign_tot_vrel %>% filter(timePoint == "End") %>%
  select(isAdapted, model, dataset, timeToAdapt, r, absCS_Mb, 
         absCS_Gb, bTGb, bTMb, absCS_Mb, vrel_m,
         cev_g, cev_m)

saveRDS(d_btgb_Malign_rf, "d_btgb_Malign_rf.RDS")
d_btgb_Malign_rf <- readRDS("d_btgb_Malign_rf.RDS")


# Filter out r < 0.1, makes analysis simpler
d_btgb_Malign_rf_nor <- d_btgb_Malign_rf %>% filter(r == -1) %>%
  select(-r)

# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215
set.seed(seed)
# Sample per group to avoid unbalanced groups
adapted_counts <- table(d_btgb_Malign_rf_nor$isAdapted)
total_counts <- sum(adapted_counts)
num_responses <- length(adapted_counts)
adapted_weights <- total_counts / (num_responses * adapted_counts)
names(adapted_weights) <- levels(d_btgb_Malign_rf_nor$isAdapted)


idx <- sample(2, nrow(d_btgb_Malign_rf_nor), replace = T, prob = c(0.7, 0.3))
train_mbeta_adapted <- d_btgb_Malign_rf_nor[idx == 1,]
test_mbeta_adapted <- d_btgb_Malign_rf_nor[idx == 2,]

# no balancing
rf_mbeta_adapted_nobal <- randomForest(formula = isAdapted ~ model * dataset * absCS_Gb * absCS_Mb * 
                                         bTGb * bTMb * absCS_Mb * vrel_m * cev_g * cev_m,
                                       data = train_mbeta_adapted,
                                       ntree = 500,
                                       proximity = T,
                                       importance = T,
                                       type = "classification")

print(rf_mbeta_adapted_nobal)

# With balancing (class weights)
rf_mbeta_adapted_bal <- randomForest(formula = isAdapted ~ model * dataset * absCS_Gb * absCS_Mb * 
                                       bTGb * bTMb * absCS_Mb * vrel_m * cev_g * cev_m,
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

# Tune the threshold
best_threshold <- pROC::coords(d_roc_bal, "best", best.method = "youden")


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

bor_mbeta <- Boruta::Boruta(isAdapted ~ ., data = d_btgb_Malign_rf_nor)
bor_mbeta
plot(bor_mbeta)

d_bor_mbeta <- process_the_Boruta_data(bor_mbeta)

feature_names <- c(TeX("Shadow Variable (min)", output = "character"),
                   TeX("Shadow Variable (mean)", output = "character"),
                   TeX("Shadow Variable (max)", output = "character"),
                   TeX("$G/\\beta$ alignment",
                       output = "character"),
                   TeX("$e_\\beta^G (\\beta^TG\\beta)$", output = "character"),
                   TeX("$V_{rel}$ (G)", output = "character"),
                   TeX("$R_{max} / \\beta$ alignment", output = "character"),
                   TeX("$e_{c}^G$", output = "character"),
                   TeX("$M/\\beta$ alignment",
                       output = "character"),
                   TeX("$e_\\beta^M (\\beta^TM\\beta)$", output = "character"),
                   TeX("$e_{c}^M$", output = "character"),
                   TeX("$V_{rel}$ (M)", output = "character"),
                   TeX("Motif", output = "character"))


pal_boruta <- c(rep("#00C0EA", times = 3),
                rep("#00A000", times = 10))

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
                                data = test_mbeta_adapted[, 2:11], 
                                y = test_mbeta_adapted$isAdapted,
                                type = "prob")

# Need to set the option future globals maxsize
options(future.globals.maxSize = 3221225472)
imp <- iml::FeatureImp$new(predictor,
                           loss = "ce",
                           n.repetitions = 100)

ggplot(imp$results %>% 
         mutate(feature = factor(feature, levels = c("absCS_Gb", "bTGb", "absCS_Mb", "dataset",  
                                                     "cev_g", "absCS_Mb", "bTMb", "cev_m", 
                                                     "vrel_m",  "model"))),
       aes(x = feature, y = importance)) +
  geom_point() +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                width = 0.2) +
  scale_x_discrete(labels = parse(text = feature_names[4:13]),
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
                                      levels = c("absCS_Gb",
                                                 "bTGb",
                                                 "absCS_Mb",
                                                 "dataset",
                                                 "cev_g",
                                                 "absCS_Mb",
                                                 "bTMb",
                                                 "cev_m",
                                                 "vrel_m",
                                                 "model"))

ggplot(d_sob_mbeta_adapted,
       aes(x = feature, y = sobelMDA)) +
  geom_point() +
  geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
               linewidth = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = parse(text = feature_names[4:13]),
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

ale_plots <- vector(mode = "list", length = 10)

ale_labels <- feature_names[c(13, 7, 9, 4, 5, 10, 6, 12, 8, 11)]
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
    # Remove outliers
    d_ale <- d_ale %>% filter(.borders < 1000)
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

#########################################################
# Per model
#seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215

rf_result <- RunRandomForestPerMotif(d_btgb_Malign_rf_nor, seed)
saveRDS(rf_result, "rf_result.RDS")
rf_result <- readRDS("rf_result.RDS")

rf_result[["NAR"]]

rf_result[["NAR"]]$cMat_bal
rf_result[["PAR"]]$cMat_bal
rf_result[["FFLC1"]]$cMat_bal
rf_result[["FFLI1"]]$cMat_bal
rf_result[["FFBH"]]$cMat_bal









#############################################################
# Per model comparisons for effects of molecular components
# Add on mean molecular component values
d_rf_molcomp <-  d_qg_tot %>% rename(timePoint = gen) %>%
                              filter(timePoint == 50000 | timePoint == 60000,
                                     log10(r) == -1) %>%
                              mutate(timePoint = if_else(timePoint == 50000, "Start", "End"),
                                     isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                                                        labels = c("Adapted", "Maladapted")),
                                     model = factor(model, levels = model_names,
                                                    labels = model_names_noquote),
                                     timePoint = factor(timePoint, levels = c("Start", "End")),
                                     dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised"))
                              ) %>%
                              select(isAdapted, model, dataset,
                                     starts_with("meanMC"))


# Per motif RF model
#seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215

rf_molcomps_result <- RunRandomForestMolCompPerMotif(d_rf_molcomp, seed)
saveRDS(rf_molcomps_result, "rf_molcomps_result.RDS")
rf_molcomps_result <- readRDS("rf_molcomps_result.RDS")


rf_molcomps_result[["NAR"]]$plt_featimp
rf_molcomps_result[["PAR"]]$plt_featimp
rf_molcomps_result[["FFLC1"]]$plt_featimp
rf_molcomps_result[["FFLI1"]]$plt_featimp
rf_molcomps_result[["FFBH"]]$plt_featimp

# Plot the most important features for each motif (Boruta importance)
d_bor_tot <- plyr::rbind.fill(rf_molcomps_result[["NAR"]]$bor %>% mutate(model = "NAR"),
                   rf_molcomps_result[["PAR"]]$bor %>% mutate(model = "PAR"),
                   rf_molcomps_result[["FFLC1"]]$bor %>% mutate(model = "FFLC1"),
                   rf_molcomps_result[["FFLI1"]]$bor %>% mutate(model = "FFLI1"),
                   rf_molcomps_result[["FFBH"]]$bor %>% mutate(model = "FFBH")
)

d_bor_top <- d_bor_tot %>%
  mutate(rowNum = row_number()) %>%
  pivot_longer(cols = -c(model, rowNum), names_to = "feature", values_to = "boruta") %>%
  filter(feature != "dataset", !grepl("shadow", feature)) %>% # Exclude trait/selection alignment and shadow vars
  mutate(boruta = if_else(is.na(boruta), 0, boruta)) %>% # Fill na values
  group_by(rowNum) %>%
  mutate(boruta = boruta / sum(boruta)) %>% # Normalise so max importance is 1
  ungroup() %>%
  select(-rowNum) %>%
  group_by(model, feature) %>%
  summarise(medBoruta = median(boruta, na.rm = T)) %>%
  filter(!is.na(medBoruta)) %>%
  arrange(desc(medBoruta), .by_group = T) %>%
  slice_head(n = 3)
d_bor_top

d_bor_tot %>%
  mutate(rowNum = row_number()) %>%
  pivot_longer(cols = -c(model, rowNum), names_to = "feature", values_to = "boruta") %>%
  filter(feature != "dataset", !grepl("shadow", feature)) %>% # Exclude trait/selection alignment and shadow vars
  #filter(interaction(model, feature) %in% d_bor_top$mod.feat) %>%
  #mutate(boruta = if_else(is.na(boruta), 0, boruta)) %>% # Fill na values
  group_by(rowNum) %>%
  mutate(boruta = (boruta / sum(boruta, na.rm = T))) %>% # Normalise so it sums to one
  ungroup() %>%
  select(-rowNum) %>%
  mutate(model = factor(model, levels = model_names_noquote),
         math_feat = c(all_molcomp_features,
                       "dataset" = TeX("Trait/selection alignment", output = "character"))[feature]) -> d_bor_plt
d_bor_plt %>%
  group_by(model, math_feat) %>%
  summarise(meanBoruta = mean(boruta),
            CIBoruta = CI(boruta)) -> d_bor_plt_sum

# Meadow plot
ggplot(d_bor_plt_sum,
       aes(x = model, y = math_feat, fill = meanBoruta)) +
  geom_tile() +
  geom_beeswarm(data = d_bor_plt, 
              mapping = aes(fill = boruta),
              shape = 21, size = 1.5) +
  theme_bw() +
  scale_y_discrete(labels = function(l) parse(text = l), limits = rev) +
  scale_fill_paletteer_c("viridis::viridis", na.value = "#444", limits = c(0.0, 0.27)) +
  guides(fill = guide_colourbar(theme = theme(legend.key.width = unit(dev.size()[1] / 2, "inches")),
                                title.vjust = 0.8)) +
  labs(x = "Model", y = "Molecular component", fill = "Normalised Boruta Importance") +
  theme(legend.position = "bottom",
        legend.title.align = 1,
        text = element_text(size = 12)) -> plt_model_molcomp_imp
plt_model_molcomp_imp
ggsave("plt_model_molcomp_imp.png", plt_model_molcomp_imp, device = png, 
       bg = "white", width = 7, height = 9)

# Reduce fig for presentation
ggplot(d_bor_plt_sum %>% filter(math_feat %in% c("alpha[Z]", "beta[Z]")),
       aes(x = model, y = math_feat, fill = meanBoruta)) +
  geom_tile() +
  theme_bw() +
  scale_y_discrete(labels = function(l) parse(text = l), limits = rev) +
  scale_fill_paletteer_c("viridis::viridis",
                         breaks = seq(0, 0.25, length.out = 7),
                         labels = round(seq(0, 0.25, length.out = 7), digits = 2)) +
  guides(fill = guide_colourbar(theme = theme(legend.key.width = unit(dev.size()[1] / 4, "inches")),
                                title.vjust = 0.8)) +
  labs(x = "Model", y = "Molecular component", fill = "Importance") +
  theme(legend.position = "bottom",
        legend.title.align = 1,
        text = element_text(size = 12)) -> plt_model_molcomp_imp_pres
plt_model_molcomp_imp_pres
ggsave("plt_model_molcomp_imp_pres.png", plt_model_molcomp_imp_pres, device = png, 
       bg = "white", width = 8, height = 5)


###############
# ALE for molcomp models - presentation figure
rf_molcomps_result[["NAR"]]$alePlot
rf_molcomps_result[["PAR"]]$alePlot
rf_molcomps_result[["FFLC1"]]$alePlot
rf_molcomps_result[["FFLI1"]]$alePlot
rf_molcomps_result[["FFBH"]]$alePlot

# Combine ALEs across motifs
d_ale <- lapply(rf_molcomps_result, function(x) {
  result <- x$ale$results[-1]
  bind_rows(result, .id = "column_label")
})

d_ale <- bind_rows(d_ale, .id = "model")
d_ale$.feature <- c(all_molcomp_features,
                    "dataset" = TeX("Trait/selection alignment", output = "character"))[d_ale$.feature]

# Remove outliers
d_ale <- d_ale %>% filter(.borders < 1000, .class == "Adapted") %>%
  mutate(model = factor(model, levels = model_names_noquote))


plt_ale_pres <- ggplot(d_ale %>% filter(column_label == "aZ" | column_label == "bZ"),
  aes(x = .borders, y = .value, colour = model, group = model)) +
  facet_nested("Molecular~component" + .feature~model,
               labeller = label_parsed, scales = "free", independent = "x") +
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_colour_manual(values = pal) +
  theme_bw() +
  labs(x = "Molecular component value", y = "Effect on adaptedness (ALE)",
       colour = "Model") +
  theme(legend.position = "bottom", text = element_text(size = 12))
plt_ale_pres
ggsave("plt_ale_molcomp_pres.png", plt_ale_pres, device = png, bg = "white",
       width = 10, height = 5)


# Full ALE
plt_ale_molcomp <- ggplot(d_ale,
                       aes(x = .borders, y = .value, colour = model, group = model)) +
  facet_nested("Molecular~component" + .feature~model,
               labeller = label_parsed, scales = "free", independent = "x") +
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_colour_manual(values = pal) +
  theme_bw() +
  labs(x = "Molecular component value", y = "Effect on adaptedness (ALE)",
       colour = "Model") +
  theme(legend.position = "bottom", text = element_text(size = 12))
plt_ale_molcomp
ggsave("plt_ale_molcomp.png", plt_ale_molcomp, device = png, bg = "white",
       width = 7.5, height = 12)

# Average ALE across all molcomps
d_ale %>% filter(.class == "Adapted") %>%
  group_by(model, .feature) %>%
  summarise(meanALE = mean(.value),
            CIALE = CI(.value)) -> d_ale_molcomp_sum

ggplot(d_ale_molcomp_sum,
       aes(x = model, y = meanALE, colour = model, group = model)) +
  facet_nested("Molecular~component" + .feature~.,
               labeller = label_parsed) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = meanALE - CIALE, ymax = meanALE + CIALE),
                width = 0.2) +
  scale_colour_manual(values = pal) +
  theme_bw() +
  labs(x = "Model", y = "Effect on adaptedness (ALE)",
       colour = "Model") +
  theme(legend.position = "bottom", text = element_text(size = 12)) -> plt_ale_molcomp_sum
plt_ale_molcomp_sum

kableExtra::kable(d_ale_molcomp_sum, "latex")




# Interactions
rf_molcomps_result[["NAR"]]$ia$plot()
rf_molcomps_result[["PAR"]]$ia$plot()
rf_molcomps_result[["FFLC1"]]$ia$plot()
rf_molcomps_result[["FFLI1"]]$ia$plot()
rf_molcomps_result[["FFBH"]]$ia$plot()

# Combine ALEs across motifs
d_int <- lapply(rf_molcomps_result, function(x) {
  x$ia$results %>% filter(.feature != "dataset", .class == "Adapted")
})

d_int <- bind_rows(d_int, .id = "model")
d_int$feat_math <- c(all_molcomp_features,
                    "dataset" = TeX("Trait/selection alignment", output = "character"))[d_int$.feature]

# Setup model
d_int <- d_int %>%
  mutate(model = factor(model, levels = model_names_noquote))

plt_int_pres <- ggplot(d_int %>% filter(.feature == "aZ" | .feature == "bZ" | .feature == "h"),
                       aes(x = .interaction, y = model, colour = model)) +
  facet_nested("Molecular component"+feat_math~., labeller = labeller(feat_math = label_parsed)) +
  geom_lollipop(horizontal = T, show.legend = F) +
  scale_colour_manual(values = pal) +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  labs(x = "Overall interaction strength", y = "Molecular component",
       colour = "Model") +
  theme(legend.position = "bottom", text = element_text(size = 12))
plt_int_pres
ggsave("plt_int_molcomp_pres.png", plt_int_pres, device = png, bg = "white",
       width = 4, height = 8)


# Full interactions plot
plt_int <- ggplot(d_int,
                       aes(x = .interaction, y = model, colour = model)) +
  facet_nested("Molecular component"+feat_math~., labeller = labeller(feat_math = label_parsed)) +
  geom_lollipop(horizontal = T, show.legend = F) +
  scale_colour_manual(values = pal) +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  labs(x = "Overall interaction strength", y = "Molecular component",
       colour = "Model") +
  theme(legend.position = "bottom", text = element_text(size = 12))
plt_int
ggsave("plt_int_molcomp.png", plt_int, device = png, bg = "white",
       width = 4, height = 10)

d_int %>%
  group_by(model) %>%
  summarise(mean(.interaction, na.rm = T))

# Plot interactions for alpha/beta/h with everything else
i_alpha <- lapply(rf_molcomps_result, function(x) {
  Interaction$new(x[["pred"]], feature = "aZ")
})
plot(i_alpha[["NAR"]])
plot(i_alpha[["PAR"]])
plot(i_alpha[["FFLC1"]])
plot(i_alpha[["FFLI1"]])
plot(i_alpha[["FFBH"]])

i_beta <- lapply(rf_molcomps_result, function(x) {
  Interaction$new(x[["pred"]], feature = "bZ")
})
plot(i_beta[["NAR"]])
plot(i_beta[["PAR"]])
plot(i_beta[["FFLC1"]])
plot(i_beta[["FFLI1"]])
plot(i_beta[["FFBH"]])

i_h <- lapply(rf_molcomps_result, function(x) {
  Interaction$new(x[["pred"]], feature = "h")
})
plot(i_h[["NAR"]])
plot(i_h[["PAR"]])
plot(i_h[["FFLC1"]])
plot(i_h[["FFLI1"]])
plot(i_h[["FFBH"]])






















# Sample per group to avoid unbalanced groups
adapted_counts <- table(d_rf_molcomp$isAdapted)
total_counts <- sum(adapted_counts)
num_responses <- length(adapted_counts)
adapted_weights <- total_counts / (num_responses * adapted_counts)
names(adapted_weights) <- levels(d_rf_molcomp$isAdapted)


idx <- sample(2, nrow(d_rf_molcomp), replace = T, prob = c(0.7, 0.3))
train_molcomp <- d_rf_molcomp[idx == 1,]
test_molcomp <- d_rf_molcomp[idx == 2,]

# no balancing
rf_molcomp_nobal <- randomForest(formula = isAdapted ~ .,
                                       data = train_molcomp,
                                       ntree = 500,
                                       proximity = T,
                                       importance = T,
                                       type = "classification")

print(rf_molcomp_nobal)

# With balancing (class weights)
rf_molcomp_bal <- randomForest(formula = isAdapted ~ .,
                                     data = train_molcomp,
                                     strata = train_molcomp$isAdapted,
                                     classwt = adapted_weights,
                                     ntree = 500,
                                     proximity = T,
                                     importance = T,
                                     type = "classification")

print(rf_molcomp_bal)

# Training data
p_train_molcomp_bal <- predict(rf_molcomp_bal, train_molcomp)
caret::confusionMatrix(p_train_molcomp_bal, train_molcomp$isAdapted)

p_train_molcomp_nobal <- predict(rf_molcomp_nobal, train_molcomp)
caret::confusionMatrix(p_train_molcomp_nobal, train_molcomp$isAdapted)


# Test data
p_test_molcomp_bal <- predict(rf_molcomp_bal, test_molcomp)
p_test_molcomp_bal_probs <- predict(rf_molcomp_bal, test_molcomp,
                                          type = "prob")[,1]

caret::confusionMatrix(p_test_molcomp_bal, test_molcomp$isAdapted)

p_test_molcomp_nobal <- predict(rf_molcomp_nobal, test_molcomp)
p_test_molcomp_nobal_probs <- predict(rf_molcomp_nobal, test_molcomp,
                                            type = "prob")[,1]

caret::confusionMatrix(p_test_molcomp_nobal, test_molcomp$isAdapted)

# roc
d_roc_bal <- roc(response = test_molcomp$isAdapted,
                 predictor = p_test_molcomp_bal_probs)

# Tune the threshold
best_threshold <- pROC::coords(d_roc_bal, "best", best.method = "youden")


d_roc_nobal <- roc(response = test_molcomp$isAdapted,
                   predictor = p_test_molcomp_nobal_probs,
                   levels = rev(levels(test_molcomp$isAdapted)))

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
ggsave("plt_RF_ROC_molcomp.png", device = png, width = 4, height = 4, bg = "white")



# Plot errors (black line = OOB, red = false positive, green = false negative)
plot(rf_molcomp_bal)
plot(rf_molcomp_nobal)


# Number of nodes per tree
hist(treesize(rf_molcomp_bal),
     main = "# Nodes for the RF trees",
     col = "forestgreen")

#############################
# Importance measures
## Boruta, permutation importance, sobol MDA

bor_mbeta <- Boruta::Boruta(isAdapted ~ ., data = d_rf_molcomp)
bor_mbeta
plot(bor_mbeta)

d_bor_mbeta <- process_the_Boruta_data(bor_mbeta)

feature_names <- c(TeX("Shadow Variable (min)", output = "character"),
                   TeX("Shadow Variable (mean)", output = "character"),
                   TeX("Shadow Variable (max)", output = "character"),
                   TeX("$G/\\beta$ alignment",
                       output = "character"),
                   TeX("$e_\\beta^G (\\beta^TG\\beta)$", output = "character"),
                   TeX("$V_{rel}$ (G)", output = "character"),
                   TeX("$R_{max} / \\beta$ alignment", output = "character"),
                   TeX("$e_{c}^G$", output = "character"),
                   TeX("$M/\\beta$ alignment",
                       output = "character"),
                   TeX("$e_\\beta^M (\\beta^TM\\beta)$", output = "character"),
                   TeX("$e_{c}^M$", output = "character"),
                   TeX("$V_{rel}$ (M)", output = "character"),
                   TeX("Motif", output = "character"))


pal_boruta <- c(rep("#00C0EA", times = 3),
                rep("#00A000", times = 10))

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
predictor <- iml::Predictor$new(rf_molcomp_bal, 
                                data = test_molcomp[, 2:14], 
                                y = test_molcomp$isAdapted,
                                type = "prob")

# Need to set the option future globals maxsize
options(future.globals.maxSize = 3221225472)
imp <- iml::FeatureImp$new(predictor,
                           loss = "ce",
                           n.repetitions = 100)

ggplot(imp$results,
       aes(x = feature, y = importance)) +
  geom_point() +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                width = 0.2) +
  # scale_x_discrete(labels = parse(text = feature_names[4:13]),
  #                  guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Permutation Importance") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_perm_imp
plt_perm_imp
ggsave("plt_perm_feat_imp_molcomp.png", device = png, width = 9, height = 5, bg = "white")

# Sobol MDA
rf_sob_mbeta_adapted <- sobolMDA::ranger(isAdapted ~ .,
                                         data = train_molcomp, num.trees = 500, 
                                         importance = "sobolMDA")
sob_mbeta_adapted <- rf_sob_mbeta_adapted$variable.importance
d_sob_mbeta_adapted <- data.frame(feature = names(sob_mbeta_adapted),
                                  sobelMDA = sob_mbeta_adapted)

d_sob_mbeta_adapted$feature <- factor(d_sob_mbeta_adapted$feature,
                                      levels = c("absCS_Gb",
                                                 "bTGb",
                                                 "absCS_Mb",
                                                 "dataset",
                                                 "cev_g",
                                                 "absCS_Mb",
                                                 "bTMb",
                                                 "cev_m",
                                                 "vrel_m",
                                                 "model"))

ggplot(d_sob_mbeta_adapted,
       aes(x = feature, y = sobelMDA)) +
  geom_point() +
  geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
               linewidth = 0.5) +
  theme_bw() +
  # scale_x_discrete(labels = parse(text = feature_names[4:13]),
  #                  guide = guide_axis(n.dodge = 2)) +
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

ale_plots <- vector(mode = "list", length = 10)

ale_labels <- feature_names[c(13, 7, 9, 4, 5, 10, 6, 12, 8, 11)]
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
    # Remove outliers
    d_ale <- d_ale %>% filter(.borders < 1000)
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





# xgboost
## What is it about the models that predicts adaptedness?
## Why are NAR vs PAR more adapted than others
## Two parts: 
### 1) a model that across all datasets shows relationship between
### variation features and adaptedness
### 2) a model which shows differences in those features between models
### 3) a model that shows how variability in the molecular components contributes
###### to those features on a per-model basis

## 1) Use xgboost without dataset, model, only the features and outcome
d_xgb <- d_btgb_Malign_rf_nor %>% select(-c(dataset, model, timeToAdapt))
#seed <- sample(1:.Machine$integer.max, 1)
#793952135
seed <- 793952135
set.seed(seed)
train.test = c(0.7, 0.3)
idx_xgb <- sample(2, nrow(d_xgb), replace = T, prob = train.test)
train_xgb <- d_xgb[idx_xgb == 1,]
test_xgb <- d_xgb[idx_xgb == 2,]

# convert to proper format
train_matrix <- as.data.frame(train_xgb %>% mutate(across(where(~is.factor(.)),
                                                      ~as.numeric(.))))
train_matrix <- xgb.DMatrix(data = as.matrix(train_matrix %>% select(-isAdapted)),
                            label = train_matrix$isAdapted) 

test_matrix <- as.data.frame(test_xgb %>% mutate(across(where(~is.factor(.)),
                                                    ~as.numeric(.))))

# Tune XGB  hyperparameters
opt_weight <- sum(train_xgb$isAdapted == "Maladapted") / sum(train_xgb$isAdapted == "Adapted")

hyper_grid <- expand.grid(
  eta = 0.05,
  max_depth = 3,
  min_child_weight = 3,
  scale_pos_weight = c(opt_weight, 0.0001, 0.000001),
  subsample = 0.5,
  colsample_bytree = 0.5,
  gamma = 0.0,
  lambda = c(0.01, 0.1, 1),
  alpha = c(0.01, 0.1, 1),
  RMSE = NA,
  trees = NA
)

pb <- progress::progress_bar$new(
  format = "  Running [:bar] :percent in :elapsedfull",
  total = nrow(hyper_grid), clear = FALSE, width = 60)

for (i in seq_len(nrow(hyper_grid))) {
  set.seed(seed)
  xgb_model <- xgb.cv(data = train_matrix,
                      nrounds = 1000,
                      early_stopping_rounds = 50,
                      nfold = 10,
                      verbose = 0,
                      params = list(
                        objective = "reg:squarederror",
                        scale_pos_weight = hyper_grid$scale_pos_weight[i],
                        eta = hyper_grid$eta[i],
                        max_depth = hyper_grid$max_depth[i],
                        min_child_weight = hyper_grid$min_child_weight[i],
                        subsample = hyper_grid$subsample[i],
                        colsample_bytree = hyper_grid$colsample_bytree[i],
                        gamma = hyper_grid$gamma[i],
                        lambda = hyper_grid$lambda[i],
                        alpha = hyper_grid$alpha[i]
                      )
  )
  
  
  hyper_grid$RMSE[i] <- min(xgb_model$evaluation_log$test_rmse_mean)
  hyper_grid$trees[i] <- xgb_model$early_stop$best_iteration
  
  pb$tick()
}

# Pick the best learning rate to continue with
hyper_grid %>%
  filter(RMSE > 0) %>%
  arrange(RMSE) %>%
  glimpse()

best_params <- as.list(dplyr::arrange(hyper_grid, RMSE)[1,] %>% 
                         select(-c(trees, RMSE)))
best_params[["objective"]] <- "reg:squarederror"
best_nrounds <- dplyr::arrange(hyper_grid, RMSE)[1, "trees"]

xgb_model_final <- xgb.train(
  params = best_params,
  data = train_matrix,
  nrounds = best_nrounds,
  verbose = 0
)

# Confusion matrix w/ test set
pred <- predict(xgb_model_final, as.matrix(test_matrix %>% select(-isAdapted)))
pred <- as.numeric(pred > 0.5)

xgb_confusion <- caret::confusionMatrix(factor(pred, labels = c("Adapted", "Maladapted")), 
                                         factor(test_matrix$isAdapted,
                                                labels = c("Adapted", "Maladapted")))
xgb_confusion

# Feature importance - SHAP values
shap_values <- predict(xgb_model_final, as.matrix(test_matrix %>% select(-isAdapted)),
                       predcontrib = T)
d_shap_xgb <- as_tibble(shap_values) %>% select(-BIAS) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "shap") %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(shap))) %>%
  ungroup() %>%
  mutate(scaled_imp = mean_abs_shap / max(mean_abs_shap))

ggplot(d_shap_xgb,
       aes(x = reorder(feature, scaled_imp), y = scaled_imp)) +
  geom_lollipop()  +
  coord_flip() +
  labs(x = "Feature", y = "Normalised Shapley importance") +
  theme_bw() +
  theme(text = element_text(size = 12))


## XGB having a lot of trouble with false positives, try random forest
table(train_xgb$isAdapted)

ctrl <- trainControl(method = "cv",
                    classProbs = T,
                    savePredictions = "final",
                    summaryFunction = twoClassSummary)
nmin <- sum(train_xgb$isAdapted == "Maladapted")

set.seed(seed)
rf_model <- train(isAdapted ~ .,
                  data = train_xgb,
                  method = "rf",
                  ntree = 1500,
                  tuneLength = 5,
                  metric = "ROC",
                  trControl = ctrl,
                  strata = train_xgb$isAdapted,
                  sampsize = rep(nmin, 2))

rf_prediction <- predict(rf_model, test_xgb)

confusionMatrix(rf_prediction, reference = test_xgb$isAdapted)

set.seed(seed)
rf_model_unbal <- train(isAdapted ~ .,
                        data = train_xgb,
                        method = "rf",
                        ntree = 1500,
                        tuneLength = 5,
                        metric = "ROC",
                        trControl = ctrl)
confusionMatrix(rf_model_unbal)

rf_unbal_prediction <- predict(rf_model_unbal, test_xgb)

confusionMatrix(rf_unbal_prediction, reference = test_xgb$isAdapted)


rf_probs <- predict(rf_model, test_xgb, type = "prob")[,1]
rf_unbal_probs <- predict(rf_model_unbal, test_xgb, type = "prob")[,1]
rf_roc <- roc(response = test_xgb$isAdapted,
              predictor = rf_probs,
              levels = rev(levels(test_xgb$isAdapted)))
rf_unbal_roc <- roc(response = test_xgb$isAdapted,
              predictor = rf_unbal_probs,
              levels = rev(levels(test_xgb$isAdapted)))

plot(rf_roc, col = rgb(1, 0, 0, 0.5), lwd = 2)
plot(rf_unbal_roc, col = rgb(0, 0, 1, 0.5), lwd = 2, add = T)

ggplot(d_xgb,
       aes(x = isAdapted)) +
  geom_bar(stat = "count") +
  theme_bw()


# Boruta for most important predictors
bor_adapted <- Boruta::Boruta(isAdapted ~ ., d_btgb_Malign_rf %>% select(-timeToAdapt))
plot(bor_adapted)

d_bor <- process_the_Boruta_data(bor_adapted)

shadow_names <- c("shadowMin" = TeX("Shadow Variable (min)", output = "character"),
                  "shadowMean" = TeX("Shadow Variable (mean)", output = "character"),
                  "shadowMax" = TeX("Shadow Variable (max)", output = "character"))

# Sort variables by Boruta median for all other importance plots
bor_order <- names(sort(unlist(d_bor %>%
                                 summarise_all(median))))

# Axis labels
motif_labels <- c("r" = TeX("Recombination rate", output = "character"),
                  "cev_g" = TeX("$e_c^G$", output = "character"),             
                  "dataset" = TeX("Trait/selection alignment", output = "character"),
                  "bTGb" = TeX("$\\beta^TG\\beta$", output = "character"),
                  "absCS_Mb" = TeX("$V_{rel}(G)$"), output = "character",
                  "absCS_Gb" = TeX("$|cos(\\theta)_\\beta^G|$", output = "character"),
                  "absCS_Mb" = TeX("$|cos(\\theta)_\\beta^M|$", output = "character"),
                  "vrel_m" = TeX("$V_{rel}(M)$"), output = "character",
                  "bTMb" = TeX("$\\beta^TM\\beta$", output = "character"),
                  "cev_m" = TeX("$e_c^M$", output = "character"),             
                  "model" = TeX("Motif", output = "character"),             
                  
                  shadow_names)[bor_order]

pal_boruta <- generateCol(bor_adapted, colCode=c("#00A000","#EECC00","#DD0000","#00C0EA"),
                          col = NULL)

ggplot(d_bor %>% pivot_longer(everything()) %>%
         mutate(x = fct_reorder(name, value, median)),
       aes(x = x, y = value, fill = x)) +
  geom_boxplot(show.legend = F, linewidth = 0.25) +
  theme_bw() +
  scale_fill_manual(values = pal_boruta) +
  scale_x_discrete(labels = parse(text = motif_labels),
    guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Boruta Importance") +
  theme(text = element_text(size = 12)) -> plt_boruta_imp
plt_boruta_imp
ggsave("plt_boruta_adapt_pred.png", plt_boruta_imp, 
       device = png, bg = "white",
       width = 12, height = 6)
# Most important predictors are the M parameters


# logistic regression model
log_mod_all <- lme4::glmer(isAdapted ~ bTMb + bTGb + absCS_Mb + absCS_Gb +
                         absCS_Mb + vrel_m + cev_g + cev_m +
                         (1 | model) + (1 | dataset) + (1 | r),
                       data = d_btgb_Malign_rf, family = "binomial")
log_mod <- lme4::glmer(isAdapted ~ bTMb + absCS_Mb + 
                         vrel_m + cev_m +
                         (1 | model) + (1 | dataset) + (1 | r), 
                       data = d_btgb_Malign_rf, family = "binomial") 

performance::compare_performance(log_mod_all, log_mod, rank = T)
performance::check_collinearity(log_mod_all)

summary(log_mod_all)
plot(log_mod_all)
# odds-ratio scale coefficients
exp(log_mod_all@beta)
ci.wald <- confint(log_mod_all, method = "Wald")
exp(ci.wald)
report::report(log_mod_all)

# Export GLM table to latex
stargazer::stargazer(log_mod_all)
texreg::texreg(log_mod_all)

# Significant effects on adaptation - cossim M vs beta, vrel(G) and vrel(M)

# Confusion matrix
pred_logmod_all <- factor(as.numeric(predict(log_mod_all, type = "response") > 0.5) + 1,
                      labels = c("Adapted", "Maladapted"))
cm_logmod <- as.data.frame(table(pred_logmod_all, d_btgb_Malign_rf$isAdapted)) %>%
  rename(pred = pred_logmod_all,
         obs = Var2) %>%
  group_by(obs) %>%
  mutate(Freq = Freq / sum(Freq))

ggplot(cm_logmod,
       aes(x = obs, y = pred, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, digits = 3), colour = Freq), size = 6) +
  labs(x = "Observed outcome", y = "Predicted outcome", fill = "Frequency") +
  scale_fill_viridis_c() +
  scale_colour_gradientn(colours = c("#FDE725FF", "#5DC863FF", "#21908CFF",
                                     "#3B528BFF", "#440154FF"),
                         values = scales::rescale(c(0, 0.49, 0.5, 0.51 ,1)),
                         guide = "none") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  guides(fill = guide_colourbar(barwidth=15)) +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
ggsave("plt_gls_pred_adaptedness.png", width = 6, height = 5, bg = "white",
       device = png)

caret::confusionMatrix(pred_logmod_all, reference = d_btgb_Malign_rf$isAdapted)


# 2) Most important features for predicting adaptation across all datasets/models
# were absCS_Mb, Vrel(G), and Vrel(M). Now how do the model * dataset combos vary in these elements?
# in adapted populations only?
d_btgb_adapted <- d_btgb_Malign_rf %>%
  filter(isAdapted == "Adapted") %>%
  select(model, dataset, absCS_Mb, absCS_Mb, vrel_m)

gls.CS.model <- gls(absCS_Mb ~ model * dataset, data = d_btgb_adapted,
                    weights = varComb(varIdent(form=~1|model),
                                      varIdent(form=~1|dataset)))
summary(gls.CS.model)
plot(gls.CS.model)
# predict
gls.cs.pred <- predict(gls.CS.model)
plot(d_btgb_adapted$absCS_Mb, gls.cs.pred)
performance::check_distribution(gls.CS.model)

# Beta regression: ml type is best
beta.cs.model_ml <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "ML")
beta.cs.model_bc <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "BC")
beta.cs.model_br <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "BR")
performance::compare_performance(beta.cs.model_ml, beta.cs.model_bc, beta.cs.model_br)

beta.cs.model <- beta.cs.model_ml
summary(beta.cs.model)
plot(beta.cs.model)
# predict
beta.cs.pred <- predict(beta.cs.model, type = "response")
plot(d_btgb_adapted$absCS_Mb, beta.cs.pred)
performance::model_performance(beta.cs.model)

# Beta mixed effects model - control for the other variables
glmmbeta.cs.model <- GLMMadaptive::mixed_model(fixed = absCS_Mb ~ model * dataset,
                                               random = ~1|vrel_m + 1|absCS_Mb, 
                                               family = GLMMadaptive::beta.fam(),
                                      data = d_btgb_adapted)
summary(glmmbeta.cs.model)
performance::compare_performance(beta.cs.model, glmmbeta.cs.model, rank = T)

# Use beta model - again ML is best
beta.vrelg.model_ml <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "ML")
beta.vrelg.model_bc <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "BC")
beta.vrelg.model_br <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "BR")
performance::compare_performance(beta.vrelg.model_ml, beta.vrelg.model_ml, beta.vrelg.model_br)

beta.vrelg.model <- beta.vrelg.model_ml
summary(beta.vrelg.model)
plot(beta.vrelg.model)
# predict
beta.vrelg.pred <- predict(beta.vrelg.model, type = "response")
plot(d_btgb_adapted$absCS_Mb, beta.vrelg.pred)
performance::model_performance(beta.vrelg.model)

# ML is best again
beta.vrelm.model_ml <- betareg::betareg(vrel_m ~ model * dataset, data = d_btgb_adapted, type = "ML")
beta.vrelm.model_bc <- betareg::betareg(vrel_m ~ model * dataset, data = d_btgb_adapted, type = "BC")
beta.vrelm.model_br <- betareg::betareg(vrel_m ~ model * dataset, data = d_btgb_adapted, type = "BR")
performance::compare_performance(beta.vrelm.model_ml, beta.vrelm.model_ml, beta.vrelm.model_br)

beta.vrelm.model <- beta.vrelm.model_ml
summary(beta.vrelm.model)
plot(beta.vrelm.model)

# predict
beta.vrelm.pred <- predict(beta.vrelm.model, type = "response")
plot(d_btgb_adapted$vrel_m, beta.vrelm.pred)
performance::model_performance(beta.vrelm.model)


texreg::texreg(list(beta.cs.model, beta.vrelg.model, beta.vrelm.model))


# Emmeans for per-model differences
em.cs.model <- emmeans(beta.cs.model, spec = ~ model * dataset, type = "response")
em.vrelg.model <- emmeans(beta.vrelg.model, spec = ~ model * dataset, type = "response")
em.vrelm.model <- emmeans(beta.vrelm.model, spec = ~ model * dataset, type = "response")

test(em.cs.model)

pwpp(em.cs.model, by = "model")

d_em.cs.model <- emmip(em.cs.model, ~ model | dataset, CIs = T, plotit = F)


plt_em_cs <- ggplot(  d_em.cs.model,
  aes(x = model, y = yvar, colour = dataset, group = dataset)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = dataset), width = 0.2,
                alpha = 0.5, show.legend = F) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_paletteer_d("nationalparkcolors::Denali", direction = -1) +
  labs(x = "Model", y= TeX("Predicted $\\cos(\\theta)^M_\\beta$"),
       colour = "Trait/selection alignment") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
plt_em_cs
ggsave("plt_pred_cs.png", plt_em_cs, device = png, bg = "white",
       width = 8, height = 6)

test(em.vrelg.model)
pwpp(em.vrelg.model, by = "model")

d_em.vrelg.model <- emmip(em.vrelg.model, ~ model | dataset, CIs = T, plotit = F)


plt_em_vrelg <- ggplot(d_em.vrelg.model,
  aes(x = model, y = yvar, colour = dataset, group = dataset)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = dataset), width = 0.2,
                alpha = 0.5, show.legend = F) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_paletteer_d("nationalparkcolors::Denali", direction = -1) +
  theme_bw() +
  labs(x = "Model", y= TeX("Predicted $V_{rel}^G$"), colour = "Trait/selection alignment") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
plt_em_vrelg
ggsave("plt_pred_vrelg.png", plt_em_vrelg, device = png, bg = "white",
       width = 8, height = 6)

test(em.vrelm.model)
pwpp(em.vrelm.model, by = "model")

d_em.vrelm.model <- emmip(em.vrelm.model, ~ model | dataset, CIs = T, plotit = F)

plt_em_vrelm <- ggplot(d_em.vrelm.model,
                       aes(x = model, y = yvar, colour = dataset, group = dataset)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = dataset), width = 0.2,
                alpha = 0.5, show.legend = F) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_paletteer_d("nationalparkcolors::Denali", direction = -1) +
  labs(x = "Model", y= TeX("Predicted $V_{rel}^M$"), 
       colour = "Trait/selection alignment") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
plt_em_vrelm
ggsave("plt_pred_vrelm.png", plt_em_vrelm, device = png, bg = "white",
       width = 8, height = 6)

leg <- get_legend(plt_em_vrelm)

plt_em <- plot_grid(plt_em_cs + theme(legend.position = "none"),
          plt_em_vrelg + theme(legend.position = "none"),
          plt_em_vrelm + theme(legend.position = "none"),
          labels = "AUTO",
          ncol = 3)

plot_grid(plt_em, leg, nrow = 2,
          rel_heights = c(1, 0.1))

ggsave("plt_pred_vrel_cs.png", device = png, bg = "white",
       width = 12, height = 4)

d_em_cs <- as.data.frame(summary(em.cs.model))[,-c(4,5)]
d_em_vrelg <- as.data.frame(summary(em.vrelg.model))[,-c(4,5)]
d_em_vrelm <- as.data.frame(summary(em.vrelm.model))[,-c(4,5)]

# Make tables
kableExtra::kable(rbind(d_em_cs, d_em_vrelg, d_em_vrelm), format = "latex")


pairs(em.cs.model, by = "model")
pairs(em.cs.model, by = "dataset")
plot(em.cs.model, comparisons = F)

pairs(em.vrelg.model, by = "model")
pairs(em.vrelg.model, by = "dataset")
plot(em.vrelg.model, comparisons = F)

pairs(em.vrelm.model, by = "model")
pairs(em.vrelm.model, by = "dataset")
plot(em.vrelm.model, comparisons = F)

## Average CS over dataset for model average effects
summary(em.cs.model.mod_only <- emmeans(beta.cs.model, spec = ~ model, type = "response"))

## Average absCS_Mb over models for dataset average effects
summary(em.vrelg.model.dat_only <- emmeans(em.vrelg.model, spec = ~ dataset, type = "response"))
pairs(em.vrelg.model.dat_only)

## Average absCS_Mb over dataset for model average effects
summary(em.vrelg.model.mod_only <- emmeans(em.vrelg.model, spec = ~ model, type = "response"))
pairs(em.vrelg.model.mod_only)
contrast(em.vrelg.model.mod_only)

## Average over everything for absCS_Mb mean estimate
summary(em.vrelg.model.everything <- emmeans(em.vrelg.model, spec = ~ 1, type = "response"))

## Average over everything for vrel_m mean estimate
summary(em.vrelm.model.everything <- emmeans(em.vrelm.model, spec = ~ 1, type = "response"))

## Average vrel_m over dataset for model average effects
summary(em.vrelm.model.mod_only <- emmeans(em.vrelg.model, spec = ~ model, type = "response"))
pairs(em.vrelm.model.mod_only)
contrast(em.vrelm.model.mod_only)



# How do molecular components produce this genetic/mutational variance?
beta.cs.model_ml <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "ML")
beta.cs.model_bc <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "BC")
beta.cs.model_br <- betareg::betareg(absCS_Mb ~ model * dataset, data = d_btgb_adapted, type = "BR")
performance::compare_performance(beta.cs.model_ml, beta.cs.model_bc, beta.cs.model_br)

beta.cs.model <- beta.cs.model_ml
summary(beta.cs.model)
plot(beta.cs.model)
# predict
beta.cs.pred <- predict(beta.cs.model, type = "response")
plot(d_btgb_adapted$absCS_Mb, beta.cs.pred)
performance::model_performance(beta.cs.model)

PMmat_test <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_PMmat_.csv",
                       col_names = F)

colnames(PMmat_test)[1:3] <- c("gen", "seed", "modelindex")
PMmat_test$gen <- as.integer(PMmat_test$gen)
PMmat_test$seed <- factor(PMmat_test$seed)
PMmat_test$modelindex <- factor(PMmat_test$modelindex)

PMmat_test <- AddCombosToDF(PMmat_test)

PMmat_test$model <- factor(PMmat_test$model,
                           levels = model_names,
                           labels = model_names_noquote)


test_var <- diag(readMCMatrix(PMmat_test[1,] %>% select(-c(gen, seed, modelindex, r, model)), 
             as.character(unlist(PMmat_test[1,"model"]))))

d_var_test <- getMCVar(PMmat_test)

# Load up actual data
PMmat_par <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_PMmat.csv",
                       col_names = F)
PMmat_orth <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_PMmat.csv",
                      col_names = F)
PMmat_rand <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_PMmat.csv",
                      col_names = F)

PMmat_par$dataset <- "Parallel"
PMmat_orth$dataset <- "Orthogonal"
PMmat_rand$dataset <- "Randomised"

PMmat <- rbind(PMmat_par, PMmat_orth, PMmat_rand)
colnames(PMmat)[1:3] <- c("gen", "seed", "modelindex")
PMmat$gen <- as.integer(PMmat$gen)
PMmat$seed <- factor(PMmat$seed)
PMmat$modelindex <- factor(PMmat$modelindex)
PMmat <- AddCombosToDF(PMmat)
PMmat$model <- factor(PMmat$model,
                           levels = model_names,
                           labels = model_names_noquote)

# Get variances
## Look at molecular component mutational variance as the difference between
## gen 60000 and 59950 / 50

## Many very large values
## Some molecular components reach a saturation point, where any more variation
## is effectively neutral
## because of multiplicative scaling, this results in very large parameter values
## How do we deal with this?
## Can we scale the molecular components according to their effect on the phenotype?
## Could get the inds.ODEPars values sampled in  slim_sampled_moltrait.csv
## Then for each parameter transform according to their maximum value, which should be
## 1 / K^h

# PMmat <- PMmat %>%
#   filter(gen == 60000 | gen == 59950)
# 
# PMmat %>%
#   arrange(gen) %>% # 59950 first
#   group_by(seed, modelindex, dataset) %>% 
#   mutate(across(starts_with("X"), ~ .x - lag(.x))) %>% View(.)
  
PMmat <- PMmat %>%
  filter(gen == 60000)


# Get mean component values
d_mc_means <- d_qg_tot %>%
  filter(gen == 60000) %>%
  select(seed, model, r, dataset, starts_with("meanMC"))

# Align PMmat with d_mc_means
d_mc_means <- d_mc_means %>%
  arrange(seed, model, r, dataset)

PMmat <- PMmat %>%
  arrange(seed, model, r, dataset)


d_mc_means$model <- factor(d_mc_means$model,
                      levels = model_names,
                      labels = model_names_noquote)


d_mc_var <- getMCVar(PMmat, d_mc_means)

# Select rows matching the d_btgb_Malign dataset
d_mc_var_btgb <- left_join(d_btgb_Malign_tot_vrel %>% filter(timePoint == "End"), 
                           d_mc_var %>% select(-c(gen)) %>%
                             mutate(r = factor(log10(r))),
                           by = c("seed", "modelindex", "dataset", "model", "r"))

# Remove unnecessary columns
d_mc_var_btgb <- d_mc_var_btgb %>%
  select(-c(timePoint, modelindex, timeToAdapt))

d_mc_var_btgb$dataset <- factor(d_mc_var_btgb$dataset,
                                levels = c("Orthogonal", "Parallel", "Randomised")) 

# Very large values -> set an upper bound
## TODO: redo on data from slim_sampled_moltrait.csv, transformed to correspond
## to phenotypic change
## map ODEPar value to sigmoidal curve
MAX_VAL <- 1000
d_mc_var_btgb_filtered <- d_mc_var_btgb %>%
  mutate(across(starts_with("var"),
                ~if_else(.x > MAX_VAL, MAX_VAL, .x))) %>%
  filter(if_any(starts_with("var"), ~ .x >= MAX_VAL))


## How well does variance in each molecular component contribute to the important adaptation
## traits?
## more beta regressions
## This time per motif
d_mc_var_nar <- d_mc_var_btgb_filtered %>% filter(model == "NAR") %>%
  select(absCS_Mb, starts_with("var"), dataset) %>%
  select(where(~!any(is.na(.))))


beta.cs.mc.nar_ml <- betareg::betareg(absCS_Mb ~ ., 
                                  data = d_mc_var_nar, type = "ML")

beta.cs.mc.nar_bc <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_nar, type = "BC")
beta.cs.mc.nar_br <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_nar, type = "BR")
performance::compare_performance(beta.cs.mc.nar_ml, beta.cs.mc.nar_bc, beta.cs.mc.nar_br)

beta.cs.mc.nar <- beta.cs.mc.nar_ml
summary(beta.cs.mc.nar)
plot(beta.cs.mc.nar)
# predict
beta.cs.mc.nar.pred <- predict(beta.cs.mc.nar, type = "response")
plot(d_mc_var_nar$absCS_Mb, beta.cs.mc.nar.pred)
performance::model_performance(beta.cs.mc.nar)


# beta regression sucks, try a Boruta?
bor_nar_mc <- Boruta::Boruta(absCS_Mb ~ ., d_mc_var_nar)
plot(bor_nar_mc)
## It identifies beta, the only significant component in the regression, but still unsure

# Now try randomForest
#seed <- sample(1:.Machine$integer.max, 1)
#793952135
seed <- 793952135
set.seed(seed)
train.test = c(0.7, 0.3)
idx_nar <- sample(2, nrow(d_mc_var_nar), replace = T, prob = train.test)
train_nar <- d_mc_var_nar[idx_nar == 1,]
test_nar <- d_mc_var_nar[idx_nar == 2,]

#nmin <- sum(train_nar$isAdapted == "Maladapted")
ctrl <- trainControl(method = "cv",
                     savePredictions = "final")


set.seed(seed)
rf_mc_nar <- train(absCS_Mb ~ .,
                  data = train_nar,
                  method = "rf",
                  ntree = 500,
                  tuneLength = 5,
                  metric = "RMSE",
                  trControl = ctrl)

rf_nar_prediction <- predict(rf_mc_nar, test_nar)
plot(test_nar$absCS_Mb, rf_nar_prediction, ylim = c(0, 1))

rf_mc_nar$finalModel
plot(rf_mc_nar$finalModel)

# the random forest is even worse somehow
## the beta regression model is still significant, there is just a lot of 
## unexplained variation within models

# keep on going with betareg and boruta I guess?
d_mc_var_par <- d_mc_var_btgb_filtered %>% filter(model == "PAR") %>%
  select(absCS_Mb, starts_with("var"), dataset) %>%
  select(where(~!any(is.na(.))))

beta.cs.mc.par_ml <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_par, type = "ML")
beta.cs.mc.par_bc <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_par, type = "BC")
beta.cs.mc.par_br <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_par, type = "BR")

## These models don't even run, not enough variance?

performance::compare_performance(beta.cs.mc.par_ml, beta.cs.mc.par_bc, beta.cs.mc.par_br)

beta.cs.mc.par <- beta.cs.mc.par_ml
summary(beta.cs.mc.par)
plot(beta.cs.model)


bor_par_mc <- Boruta::Boruta(absCS_Mb ~ ., d_mc_var_par)
plot(bor_par_mc)


# FFLC1
d_mc_var_fflc1 <- d_mc_var_btgb_filtered %>% filter(model == "FFLC1") %>%
  select(absCS_Mb, starts_with("var"), dataset) %>%
  select(where(~!any(is.na(.))))

beta.cs.mc.fflc1_ml <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_fflc1, type = "ML")
beta.cs.mc.fflc1_bc <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_fflc1, type = "BC")
beta.cs.mc.fflc1_br <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_fflc1, type = "BR")

performance::compare_performance(beta.cs.mc.fflc1_ml, beta.cs.mc.fflc1_bc, beta.cs.mc.fflc1_br,
                                 rank = T)

beta.cs.mc.fflc1 <- beta.cs.mc.fflc1_ml
summary(beta.cs.mc.fflc1)
plot(beta.cs.mc.fflc1)


bor_fflc1_mc <- Boruta::Boruta(absCS_Mb ~ ., d_mc_var_fflc1)
plot(bor_fflc1_mc)

# FFLI1
d_mc_var_ffli1 <- d_mc_var_btgb_filtered %>% filter(model == "FFLI1") %>%
  select(absCS_Mb, starts_with("var"), dataset) %>%
  select(where(~!any(is.na(.))))

beta.cs.mc.ffli1_ml <- betareg::betareg(absCS_Mb ~ ., 
                                        data = d_mc_var_ffli1, type = "ML")
beta.cs.mc.ffli1_bc <- betareg::betareg(absCS_Mb ~ ., 
                                        data = d_mc_var_ffli1, type = "BC")
beta.cs.mc.ffli1_br <- betareg::betareg(absCS_Mb ~ ., 
                                        data = d_mc_var_ffli1, type = "BR")

performance::compare_performance(beta.cs.mc.ffli1_ml, beta.cs.mc.ffli1_bc, beta.cs.mc.ffli1_br,
                                 rank = T)

beta.cs.mc.ffli1 <- beta.cs.mc.ffli1_ml
summary(beta.cs.mc.ffli1)
plot(beta.cs.mc.ffli1)


bor_ffli1_mc <- Boruta::Boruta(absCS_Mb ~ ., d_mc_var_ffli1)
plot(bor_ffli1_mc)

# FFLI1
d_mc_var_ffbh <- d_mc_var_btgb_filtered %>% filter(model == "FFBH") %>%
  select(absCS_Mb, starts_with("var"), dataset) #%>%
  #select(where(~!any(is.na(.))))

beta.cs.mc.ffbh_ml <- betareg::betareg(absCS_Mb ~ ., 
                                        data = d_mc_var_ffbh, type = "ML")
beta.cs.mc.ffbh_bc <- betareg::betareg(absCS_Mb ~ ., 
                                        data = d_mc_var_ffbh, type = "BC")
beta.cs.mc.ffbh_br <- betareg::betareg(absCS_Mb ~ ., 
                                        data = d_mc_var_ffbh, type = "BR")

performance::compare_performance(beta.cs.mc.ffbh_ml, beta.cs.mc.ffbh_bc, beta.cs.mc.ffbh_br,
                                 rank = T)

beta.cs.mc.ffbh <- beta.cs.mc.ffbh_ml
summary(beta.cs.mc.ffbh)
plot(beta.cs.mc.ffbh)


bor_ffbh_mc <- Boruta::Boruta(absCS_Mb ~ ., d_mc_var_ffbh)
plot(bor_ffbh_mc)

# next is absCS_Mb
d_mc_var_vrelg_nar <- d_mc_var_btgb_filtered %>% filter(model == "NAR") %>%
  select(absCS_Mb, starts_with("var"), dataset) %>%
  select(where(~!any(is.na(.))))


beta.vrelg.mc.nar_ml <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_vrelg_nar, type = "ML")

beta.vrelg.mc.nar_bc <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_vrelg_nar, type = "BC")
beta.vrelg.mc.nar_br <- betareg::betareg(absCS_Mb ~ ., 
                                      data = d_mc_var_vrelg_nar, type = "BR")
performance::compare_performance(beta.vrelg.mc.nar_ml, beta.vrelg.mc.nar_bc, beta.vrelg.mc.nar_br)

beta.vrelg.mc.nar <- beta.vrelg.mc.nar_ml
summary(beta.vrelg.mc.nar)
plot(beta.vrelg.mc.nar)
# predict
beta.vrelg.mc.nar.pred <- predict(beta.vrelg.mc.nar, type = "response")
plot(d_mc_var_vrelg_nar$absCS_Mb, beta.vrelg.mc.nar.pred)
performance::model_performance(beta.vrelg.mc.nar)

# vrel_m
d_mc_var_vrelm_nar <- d_mc_var_btgb_filtered %>% filter(model == "NAR") %>%
  select(vrel_m, starts_with("var"), dataset) %>%
  select(where(~!any(is.na(.))))


beta.vrelm.mc.nar_ml <- betareg::betareg(vrel_m ~ ., 
                                         data = d_mc_var_vrelm_nar, type = "ML")

beta.vrelm.mc.nar_bc <- betareg::betareg(vrel_m ~ ., 
                                         data = d_mc_var_vrelm_nar, type = "BC")
beta.vrelm.mc.nar_br <- betareg::betareg(vrel_m ~ ., 
                                         data = d_mc_var_vrelm_nar, type = "BR")
performance::compare_performance(beta.vrelm.mc.nar_ml, beta.vrelm.mc.nar_bc, beta.vrelm.mc.nar_br)

beta.vrelm.mc.nar <- beta.vrelm.mc.nar_ml
summary(beta.vrelm.mc.nar)
plot(beta.vrelm.mc.nar)
# predict
beta.vrelm.mc.nar.pred <- predict(beta.vrelm.mc.nar, type = "response")
plot(d_mc_var_vrelm_nar$vrel_m, beta.vrelm.mc.nar.pred)
performance::model_performance(beta.vrelm.mc.nar)



#####################################################################################
# Read in M_C G matrix
H2_DATA_PATH_RAND <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/"
H2_DATA_PATH_ORTH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
H2_DATA_PATH_PAR <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"

max_molComps <- length(all_molcomp_features)
CV_names <- combn(names(all_molcomp_features), 2)
CV_names <- paste("CVA", CV_names[1,], CV_names[2,], sep = "_") 


d_h2_rand_mrr <- read_csv(paste0(H2_DATA_PATH_RAND, "out_h2_mrr.csv"), col_names = F)
d_h2_rand_mkr <- read_csv(paste0(H2_DATA_PATH_RAND, "out_h2_mkr.csv"), col_names = F)


colnames(d_h2_rand_mkr) <- c("gen", "seed", "modelindex", "VA_W", "h2_W",
                             paste0("VA_", names(all_molcomp_features)), CV_names,
                             paste0("h2_", names(all_molcomp_features)))
colnames(d_h2_rand_mrr) <- colnames(d_h2_rand_mkr)

d_h2_rand_mkr$calcMode <- "mkr"
d_h2_rand_mrr$calcMode <- "mrr"

d_h2_mc_rand <- rbind(d_h2_rand_mkr, d_h2_rand_mrr)

d_h2_mc_rand$modelindex <- factor(d_h2_mc_rand$modelindex)
d_h2_mc_rand$seed <- factor(d_h2_mc_rand$seed)

d_h2_mc_rand <- AddCombosToDF(d_h2_mc_rand)
d_h2_mc_rand$dataset <- "Randomised" 


d_h2_orth_mrr <- read_csv(paste0(H2_DATA_PATH_ORTH, "out_h2_mrr.csv"), col_names = F)
d_h2_orth_mkr <- read_csv(paste0(H2_DATA_PATH_ORTH, "out_h2_mkr.csv"), col_names = F)


colnames(d_h2_orth_mkr) <- c("gen", "seed", "modelindex", "VA_W", "h2_W",
                             paste0("VA_", names(all_molcomp_features)), CV_names,
                             paste0("h2_", names(all_molcomp_features)))
colnames(d_h2_orth_mrr) <- colnames(d_h2_orth_mkr)

d_h2_orth_mkr$calcMode <- "mkr"
d_h2_orth_mrr$calcMode <- "mrr"

d_h2_mc_orth <- rbind(d_h2_orth_mkr, d_h2_orth_mrr)

d_h2_mc_orth$modelindex <- factor(d_h2_mc_orth$modelindex)
d_h2_mc_orth$seed <- factor(d_h2_mc_orth$seed)

d_h2_mc_orth <- AddCombosToDF(d_h2_mc_orth)
d_h2_mc_orth$dataset <- "Orthogonal" 


d_h2_par_mrr <- read_csv(paste0(H2_DATA_PATH_PAR, "out_h2_mrr.csv"), col_names = F)
d_h2_par_mkr <- read_csv(paste0(H2_DATA_PATH_PAR, "out_h2_mkr.csv"), col_names = F)


colnames(d_h2_par_mkr) <- c("gen", "seed", "modelindex", "VA_W", "h2_W",
                            paste0("VA_", names(all_molcomp_features)), CV_names,
                            paste0("h2_", names(all_molcomp_features)))
colnames(d_h2_par_mrr) <- colnames(d_h2_par_mkr)

d_h2_par_mkr$calcMode <- "mkr"
d_h2_par_mrr$calcMode <- "mrr"

d_h2_mc_par <- rbind(d_h2_par_mkr, d_h2_par_mrr)

d_h2_mc_par$modelindex <- factor(d_h2_mc_par$modelindex)
d_h2_mc_par$seed <- factor(d_h2_mc_par$seed)

d_h2_mc_par <- AddCombosToDF(d_h2_mc_par)
d_h2_mc_par$dataset <- "Parallel" 

d_h2_mc <- rbind(d_h2_mc_rand, d_h2_mc_orth, d_h2_mc_par)

d_h2_mc$model <- factor(d_h2_mc$model, 
                        levels = model_names,
                        labels = model_names_noquote)

# Discretise generation
d_h2_mc <- d_h2_mc %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))


# Construct G matrices from rows
# inner join optPerc
d_h2_mc_adapt <- left_join(d_btgb_Malign_tot_vrel %>% select(timePoint, seed, modelindex,
                                                             dataset, isAdapted), 
                           d_h2_mc %>% select(-gen), 
                           by = c("timePoint", "seed", "modelindex", "dataset"))

# Counts for each model type:
table(d_h2_mc_adapt$model, d_h2_mc_adapt$isAdapted)

# Split h2 into G matrices
d_h2_mc_adapt %>%
  select(!VA_W) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(8:19, is.na), timePoint == "End", isAdapted == "Adapted") %>%  # Drop rows with no variance
  distinct(timePoint, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(modelindex, dataset, timePoint, isAdapted) %>%
  group_split(.) -> split_h2_mc


# Separate into model indices
# each sublist is replicates of a model index
sourceCpp("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")
sourceCpp("/mnt/e/Documents/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")

lapply(split_h2_mc, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_mc

# Select rows and columns by motif type
h2_mat <- unlist(cov_matrices_mc, recursive = F)

# get ids from the matrix
cov_matrix_modelindex <- GetMatrixIDsWithDataset(split_h2_mc)

h2_mat_motif <- h2_mat

for (i in seq_len(length(h2_mat))) {
  x <- h2_mat[[i]]
  id <- cov_matrix_modelindex[[i]]
  model <- as.character(d_combos$model[as.numeric(levels(id$modelindex))[id$modelindex]])
  model <- str_replace_all(model, "'", "")
  mcs <- unlist(molComp_names[model])

  h2_mat_motif[[i]] <- x[mcs,mcs]  
}

id <- data.table::rbindlist(cov_matrix_modelindex, 
                            fill = T)
id$label <- 1:nrow(id)
id$modelindex <- as.factor(id$modelindex)
id <- AddCombosToDF(id)
id$model <- factor(id$model, levels = model_names, labels = model_names_noquote)

# First convert to nearest positive definite matrix
h2_pd <- lapply(h2_mat_motif, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

# Split per model
id %>%
  group_by(model) %>%
  group_split() -> id_permotif

# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

g_mc_nar <- h2_pd[as.integer(id_permotif[[1]]$label)]
g_mc_par <- h2_pd[as.integer(id_permotif[[2]]$label)]
g_mc_fflc1 <- h2_pd[as.integer(id_permotif[[3]]$label)]
g_mc_ffli1 <- h2_pd[as.integer(id_permotif[[4]]$label)]
g_mc_ffbh <- h2_pd[as.integer(id_permotif[[5]]$label)]

# Eigentensors
# TODO: Bootstrap estimate
# Sample matrix pairs from within/between treatments
# Look at similarity/shared components and eigenvalues
# Return list of traits and percent inclusion per model
# null distribution across all treatments

# per motif labels
id_nar <- id_permotif[[1]]
id_nar$label <- 1:nrow(id_nar)

etd_nar <- EigenTensorExperiment(g_mc_nar, id_nar, n = 100)


## Another option is to do PCASim to look at variability in usage of molecular components
## Or to run eigentensor decomposition to look at most important components amongst all
## simulations -> along with PCASim as a measure of how frequently these are used?
### This might be the best option
### Gather list of all motif matrices
### find the average scaled variance per trait across all matrices -> 
### if this is high, relative to others it is an important source of variance
### decompose the tensor of matrices and find the smallest eigentensor: this 
### describes the trait-space direction in which variation is most stable across
### simulations
### If a component has high mean relative variance and a large contribution to the
### smallest PC, it is a repeated contributor to the adapted populations
##
## would need to bootstrap this method and compare to a null distribution taken from
## maladapted populations

## Eigentensors take forever, too much data :(
## Look at model of conditional evolvability of each matrix
### Problem is that the molecular components aren't the direct targets of selection, 
### would need to estimate the strength of selection on each. Unless we do an average
### estimate across all selection gradients, mean conditional evolvability?

# Plot mean G matrices for each group
d_h2_mc_adapt %>%
  select(!VA_W) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(8:19, is.na), timePoint == "End", isAdapted == "Adapted") %>%  # Drop rows with no variance
  distinct(timePoint, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(model, dataset) %>%
  group_map(~ .x %>% 
              summarise(across(starts_with("VA") | starts_with("CVA"), 
                               list(mean), .names = "{.col}")) %>%
              mutate(model = .y$model,
                     dataset = .y$dataset),
            .keep = T) -> split_h2_mc_mean
View(split_h2_mc_mean[[1]])

lapply(split_h2_mc_mean, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_mean_matrices_mc

# Select rows and columns by motif type
h2_mean_mat <- unlist(cov_mean_matrices_mc, recursive = F)

# get ids from the matrix
cov_mean_matrix_modelindex <- GetMeanMatrixIDs(split_h2_mc_mean)

h2_mean_mat_motif <- h2_mean_mat

for (i in seq_len(length(h2_mean_mat))) {
  x <- h2_mean_mat[[i]]
  id <- cov_mean_matrix_modelindex[[i]]
  mcs <- unlist(molComp_names[as.character(id$model)])
  
  h2_mean_mat_motif[[i]] <- x[mcs,mcs]  
}

# First convert to nearest positive definite matrix
h2_mean_pd <- lapply(h2_mean_mat_motif, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

id_mean <- data.table::rbindlist(cov_mean_matrix_modelindex, 
                            fill = T)
id_mean$label <- 1:nrow(id_mean)

# From each covariance matrix, construct an adjacency matrix
## convert to correlation matrix,
cor_mean <- lapply(h2_mean_pd, function(x) {
  cov2cor(x)
})


THRESHOLD <- 0.01
graph_g_mean <- list()

for (i in seq_along(cor_mean)) {
  weights <- cor_mean[[i]]
  # unweighted
  adj <- AdjMatFromCorMat(weights, THRESHOLD)
  diag(adj) <- 0
  rownames(adj) <- rownames(weights)
  colnames(adj) <- rownames(adj)
  
  variances <- diag(h2_mean_pd[[i]])
  
  names(variances) <- all_molcomp_features[names(variances)]
  
  # Fix floating point issues
  weights[lower.tri(weights)] <- weights[upper.tri(weights)]
  
  
  # Scale weights from [-1, 1] to [0, 1]
  weights <- (weights + 1) / 2
  
  # set diag to 0
  diag(weights) <- 0.5
  
  
  graph <- (igraph::graph_from_adjacency_matrix(adj,
                                                mode = "undirected", 
                                                weighted = F))
  
  # Get nodes and edges
  d_graph <- igraph::as_data_frame(graph)

  # Add group name
  d_graph[,c("model", "dataset")] <- cov_mean_matrix_modelindex[[i]]
  
  # Reshape into tidy graph
  graph <- d_graph %>%
    mutate(model = factor(model, levels = model_names_noquote),
           dataset = factor(dataset, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised"))) %>%
    group_by(model, dataset) %>% 
    group_map(~tbl_graph(nodes = all_molcomp_features,
                         edges = tibble(
                           from = as.integer(
                             factor(.x$from, levels = names(all_molcomp_features))),
                           to = as.integer(
                             factor(.x$to, levels = names(all_molcomp_features)))
                           #from = .x$from,
                           #to = .x$to
                         ),
                         directed = F) %>%
                # activate("edges") %>%
                # mutate(weight = .x$weight) %>%
                activate("nodes") %>%
                mutate(degree = centrality_degree(#weights = .x$weight,
                                                  normalized = T),
                       betweenness = centrality_betweenness(directed = F,
                                                            normalized = T),
                       eigen = centrality_eigen(),
              variance = variances[nodes]))#weights = .x$weight)))
  
  graph_g_mean[[i]] <- (graph)
}

names(graph_g_mean) <- interaction(id_mean$model, 
                                   id_mean$dataset)


# Create shared layout for everything
set.seed(seed)
gr_layout <- create_layout(graph_g_mean[[15]][[1]], layout = "fr")
xlims <- c(min(gr_layout$x), max(gr_layout$x))
ylims <- c(min(gr_layout$y), max(gr_layout$y))

#Separate figure for each motif
mean_graph_figs <- lapply(graph_g_mean, function(x)
{
  x <- x[[1]] %>%
    activate("nodes") %>%
    filter(degree > 0)
  
  # Consistent layout
  layout_sbst <- gr_layout[gr_layout$nodes %in% igraph::V(x)$nodes,]
  layout_sbst$degree <- x %>% activate(nodes) %>% pull(degree)
  layout_sbst$eigen <- x %>% activate(nodes) %>% pull(eigen)
  layout_sbst$variance <- x %>% activate(nodes) %>% pull(variance)
  attr(layout_sbst, "graph") <- attr(layout_sbst, "graph") %>%
    activate(edges) %>%
    filter(F) %>% # Remove all edges
    bind_edges(x %>% activate(edges) %>% data.frame(), node_key = "nodes")
  
  ggraph(graph = layout_sbst) +
    geom_edge_link(colour = "black",  start_cap = circle(0.2),
                       end_cap = circle(0.2)) +
                       #sep = unit(1, "lines")) +
    geom_node_point(shape = 21, aes(fill = degree, size = variance)) +
    geom_node_text(aes(label = nodes), parse = T) +
    scale_fill_gradient(low = "#8E8FEE", high = "#CD2626",
                        limits = c(0, 0.8), n.breaks = 4) +
    expand_limits(x = xlims, y = ylims) +
    scale_size(range = c(6, 12), limits = c(0, 1)) +
    #facet_edges(dataset~model, nrow = 3) +
    labs(fill = "Degree centrality", size = "Genetic variance") +
    theme_graph() +
    theme(legend.position = "bottom")
}) 

leg <- get_legend(mean_graph_figs[[1]])

plt_mean_graphs <- plot_grid(plotlist = lapply(mean_graph_figs, function(x) {x + theme(legend.position = "none")}),
          nrow = 3,
          ncol = 5,
          byrow = F)

plot_grid(plt_mean_graphs,
          leg, 
          nrow = 2,
          rel_heights = c(1, 0.1))


# Plot each separate matrix, then combine into a mean plot rather than mean across matrices 
# first
cor_pd <- lapply(h2_pd, function(x){
  cov2cor(x)
})

graph_g <- list()

# Choose threshold 0.1
THRESHOLD <- 0.2

pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total (:percent eta: :eta)", total = length(cor_pd))

for (i in seq_along(cor_pd)) {
  pb$tick()
  weights <- cor_pd[[i]]
  # unweighted
  adj <- AdjMatFromCorMat(weights, THRESHOLD)
  diag(adj) <- 0
  rownames(adj) <- rownames(weights)
  colnames(adj) <- rownames(adj)
  
  variances <- diag(h2_pd[[i]])
  
  names(variances) <- all_molcomp_features[names(variances)]
  
  # Fix floating point issues
  weights[lower.tri(weights)] <- weights[upper.tri(weights)]
  
  
  # Scale weights from [-1, 1] to [0, 1]
  weights <- (weights + 1) / 2
  
  # set diag to 0
  diag(weights) <- 0.5
  
  
  graph <- (igraph::graph_from_adjacency_matrix(adj,
                                                mode = "undirected", 
                                                weighted = F))
  
  # Get nodes and edges
  d_graph <- igraph::as_data_frame(graph)
  
  if (nrow(d_graph > 0)) {
    # Add group name
    #d_graph[,c("timePoint", "model", "dataset", "isAdapted")] <- id[i, c("timePoint", "model", "dataset", "isAdapted")]
    
    # Reshape into tidy graph
    graph <- tbl_graph(nodes = all_molcomp_features,
                           edges = tibble(
                             from = as.integer(
                               factor(d_graph$from, levels = names(all_molcomp_features))),
                             to = as.integer(
                               factor(d_graph$to, levels = names(all_molcomp_features)))
                             #from = .x$from,
                             #to = .x$to
                           ),
                           directed = F) %>%
                  # activate("edges") %>%
                  # mutate(weight = .x$weight) %>%
                  activate("nodes") %>%
                  mutate(degree = centrality_degree(#weights = .x$weight,
                    normalized = T),
                    betweenness = centrality_betweenness(directed = F,
                                                         normalized = T),
                    eigen = centrality_eigen(),
                    variance = variances[nodes])#weights = .x$weight)))
  } 
  else
  {
    graph <- as_tbl_graph(graph) %>%
      activate("nodes") %>%
      mutate(degree = centrality_degree(#weights = .x$weight,
        normalized = T),
        betweenness = centrality_betweenness(directed = F,
                                             normalized = T),
        eigen = centrality_eigen(),
        nodes = all_molcomp_features[name],
        variance = variances[name])
  }
  graph_g[[i]] <- graph
}

# Calculate average graph per group
# First create table
all_edges <- purrr::map(seq_along(graph_g), function(i) {
  x <- graph_g[[i]]
  x %>%
    activate(edges) %>%
    as_tibble() %>%
    mutate(node1 = pmin(from, to),
           node2 = pmax(from, to),
           graph_id = i,
           model = id[i,]$model,
           dataset = id[i,]$dataset,
           isAdapted = id[i,]$isAdapted)
}, .progress = T) %>%
  data.table::rbindlist(.)

## Edge weight is the proportion of graphs with the edge (correlation > THRESHOLD)
## 
summary_edges <- all_edges %>%
  filter(isAdapted == "Adapted") %>%
  group_by(model, dataset) %>%
  mutate(num_graphs = n_distinct(graph_id)) %>%
  group_by(node1, node2, model, dataset) %>% # Proportion of edges above threshold vs total in that model
  summarise(proportion = n_distinct(graph_id) / num_graphs[1],
            num_graphs = num_graphs[1]) %>%
  ungroup() %>%
  rename(from = node1, to = node2) 

all_nodes <- purrr::map(seq_along(graph_g), function(i) {
  x <- graph_g[[i]]
  x %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(model = id[i,]$model,
           dataset = id[i,]$dataset,
           isAdapted = id[i,]$isAdapted,
           graph_id = i)
}, .progress = T) %>%
  data.table::rbindlist(., fill = T)

summary_nodes <- all_nodes %>%
  filter(isAdapted == "Adapted") %>%
  mutate(nodes = factor(nodes, levels = all_molcomp_features)) %>%
  group_by(model, dataset, nodes) %>%
  summarise(meanVariance = mean(variance, na.rm = T),
            CIVariance = CI(variance, na.rm = T),
            meanDegree = mean(degree),
            CIDegree = CI(degree)) %>%
  ungroup()


# Split edges/nodes by group
summary_edges <- summary_edges %>%
  group_by(model, dataset) %>%
  group_split()

summary_nodes <- summary_nodes %>%
  group_by(model, dataset) %>%
  group_split()

# Create graph
graph_summary <- purrr::map(seq_along(summary_nodes), function(i) {
  x_n <- summary_nodes[[i]]
  x_e <- summary_edges[[i]]
  tbl_graph(
    nodes = x_n,
    edges = x_e,
    directed = F,
    node_key = "nodes"
  )
  })

# Plot average graphs
# Create shared layout for everything
set.seed(seed)
gr_layout <- create_layout(graph_summary[[15]], layout = "fr")
xlims <- c(min(gr_layout$x), max(gr_layout$x))
ylims <- c(min(gr_layout$y), max(gr_layout$y))

#Separate figure for each motif
summary_graph_figs <- purrr::map(seq_along(graph_summary), function(i)
{
  x <- graph_summary[[i]] %>%
    activate("nodes") %>%
    filter(meanDegree > 0)
  
  # Consistent layout
  layout_sbst <- gr_layout[gr_layout$nodes %in% igraph::V(x)$nodes,]
  layout_sbst$degree <- x %>% activate(nodes) %>% pull(meanDegree)
  layout_sbst$variance <- x %>% activate(nodes) %>% pull(meanVariance)
  attr(layout_sbst, "graph") <- attr(layout_sbst, "graph") %>%
    activate(edges) %>%
    filter(F) %>% # Remove all edges, stick on appropriate edges
    bind_edges(x %>% activate(edges) %>% data.frame(), node_key = "nodes")
  
  ggraph(graph = layout_sbst) +
    geom_edge_link(aes(edge_colour = proportion),  start_cap = circle(0.2),
                   end_cap = circle(0.2)) +
    #sep = unit(1, "lines")) +
    geom_node_point(shape = 21, aes(fill = degree, size = variance)) +
    geom_node_text(aes(label = nodes), parse = T) +
    scale_fill_viridis_c(limits = c(0, 1), n.breaks = 5) +
    # scale_fill_gradient2(low = "#8E8FEE", high = "#CD2626",
    #                    limits = c(-3, 0), n.breaks = 5) +
    scale_edge_colour_gradient(low = "#999", high = "#000",
                               limits = c(0, 1), n.breaks = 5) +
    expand_limits(x = xlims, y = ylims) +
    scale_size(range = c(6, 12), limits = c(0, 1)) +
    #facet_edges(dataset~model, nrow = 3) +
    labs(fill = "Degree centrality", size = "Mean genetic variance",
         edge_colour = "Confidence") +
    theme_graph() +
    theme(legend.position = "bottom")
}, .progress = T) 

summary_graph_figs[[15]]
summary_graph_figs[[1]]

leg <- get_legend(summary_graph_figs[[1]])

plt_sum_graphs <- plot_grid(plotlist = lapply(summary_graph_figs, function(x) {x + theme(legend.position = "none")}),
                             nrow = 3,
                             ncol = 5,
                             byrow = F)

plot_grid(plt_sum_graphs,
          leg, 
          nrow = 2,
          rel_heights = c(1, 0.1))







# Run conditional evolvability
test_g_cevol <- cEvolPerTrait(g_mc_nar[[1]])
cevol_mc <- ConditionalEvolvabilityExperiment(h2_pd, id)

# Join with V_rel dataset
d_btgb_cev_mc <- inner_join(d_btgb_Malign_tot_vrel %>%
                             filter(timePoint == "End") %>%
                             select(seed, modelindex, dataset, isAdapted,
                                    model, r, absCS_Mb, absCS_Mb, vrel_m),
                           cevol_mc,
                           by = c("seed", "modelindex", "dataset"))
saveRDS(d_btgb_cev_mc, "/mnt/d/SLiMTests/tests/newMotifs/paper1/d_btgb_cev_mc.RDS")

# Summary statistics for each model
d_btgb_cev_mc %>%
  rename_with(~sub("^(cev|aut)", "values\\1", .)) %>%
  pivot_longer(cols = starts_with(("values")), 
               names_pattern = "(.*)(cev|aut)_(.*)",
               names_to = c(".value", "names", "molComp")) %>%
  pivot_wider(names_from = names,
              values_from = values) -> d_mc_cev

d_mc_cev_nar <- d_mc_cev %>% filter(model == "NAR") %>% drop_na() %>%
  mutate(molComp = factor(molComp))
# Does integration differ between molecular components in each model?
beta.int.mc.nar <- betareg::betareg(1 - aut ~ molComp * dataset,
                                    data = d_mc_cev_nar)
summary(beta.int.mc.nar)

rlm.cev.mc.nar <- MASS::rlm(cev ~ molComp * dataset,
                                    data = d_mc_cev_nar)
plot(rlm.cev.mc.nar)
summary(rlm.cev.mc.nar)

pred.rlm.cev.nar <- predict(rlm.cev.mc.nar, d_mc_cev_nar$cev)

em.rlm.cev.mc.nar <- emmeans(rlm.cev.mc.nar, ~ molComp + dataset + r)
summary(em.rlm.cev.mc.nar)
pairs(em.rlm.cev.mc.nar, by = "molComp")

# Recombination has a strong effect on conditional evolvability - larger values mean
# more conditional evolvability

report::report(rlm.cev.mc.nar)
  
d_mc_cev_sum <- d_mc_cev %>%
  group_by(model, dataset, isAdapted, molComp) %>%
  summarise(meanCEV = mean(cev),
            CICEV = CI(cev),
            meanAut = mean(aut),
            CIAut = CI(aut),
            meanInt = mean(1 - aut),
            CIInt = CI(1 - aut))

ggplot(d_mc_cev_sum %>% filter(isAdapted == "Adapted"),
       aes(x = molComp, y = meanCEV, colour = model)) +
  facet_nested("Selection/trait alignment" + dataset ~ "Motif" + model) +
  geom_point() +
  geom_errorbar(aes(ymin = meanCEV - CICEV, ymax = meanCEV + CICEV)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5)) +
  scale_x_discrete(labels = function(x) parse(text = all_molcomp_features[x])) +
  labs(x = "Molecular component", y = "Mean conditional evolvability") +
  guides(colour = guide_none()) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom") -> plt_mean_cev_molcomp
plt_mean_cev_molcomp


ggplot(d_mc_cev_sum %>% filter(isAdapted == "Adapted"),
       aes(x = molComp, y = meanInt, colour = model)) +
  facet_nested("Selection/trait alignment" + dataset ~ "Motif" + model) +
  geom_point() +
  geom_errorbar(aes(ymin = meanInt - CIInt, ymax = meanInt + CIInt)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5)) +
  scale_x_discrete(labels = function(x) parse(text = all_molcomp_features[x])) +
  labs(x = "Molecular component", y = "Mean integration") +
  guides(colour = guide_none()) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom") -> plt_mean_int_molcomp
plt_mean_int_molcomp


ggplot(d_mc_cev %>% filter(isAdapted == "Adapted", seed == "750019483", r == -1),
       aes(x = molComp, y = 1 - aut, colour = model)) +
  facet_nested("Selection/trait alignment" + dataset ~ "Motif" + model) +
  geom_point() +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5)) +
  scale_x_discrete(labels = function(x) parse(text = all_molcomp_features[x])) +
  labs(x = "Molecular component", y = "Integration") +
  guides(colour = guide_none()) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom")

# Plot distribution of contribution to total integration?
## If it is very high, that component contributes a lot to constraint
## Want a measure of which molComps covary most often and which ones don't
## Cooccurrence betwork? Measure how often integration values are similar
## split up range of integration values into a factor 0, 1, 2, 3, 4, with each
# capturing a fifth of the total space (0 = [0, 0.2), 1 = [0.2, 0.4), 2 = [0.4, 0.6),
# 3 = [0.6, 0.8), 4 = [0.8, 1.0])
# count matrix of each, build cooccurence matrix across all simulations
# Use coocure
library(cooccure)
# Requires one row per entry
d_coocure <- d_btgb_cev_mc %>%
  filter(isAdapted == "Adapted") %>%
  select(model, dataset, 21:32)
  
names(d_coocure)[3:14] <- sub('aut_', '', names(d_coocure)[3:14])

d_coocure <- d_coocure %>%
  mutate(across(3:14, ~ 1 - .x,
                .names = "int_{.col}")) %>%
  select(1:2, starts_with("int")) %>%
  mutate(across(starts_with("int"), ~ cut(.x, breaks = seq(from = 0.0, to = 1.0, by = 0.2))),
       id = interaction(model, dataset)) %>%
  select(-c(1:2))

d_coocure <- d_mc_cev %>%
  filter(isAdapted == "Adapted") %>%
  select(seed, model, dataset, molComp, aut) %>%
  mutate(int = (1 - aut),
         id = interaction(seed, model, dataset),
         model.dataset = interaction(model, dataset)) %>%
  select(-aut)
d_coocure$int[is.na(d_coocure$int)] <- 0.0

co_all <- cooccurrence(d_coocure, group = "model.dataset",
                       field = "molComp", by = "id",
                       weight_by = "int",
                       similarity = "association")

library(ggraph)
library(tidygraph)
co_all_graph <- co_all %>%
  separate_wider_delim(group, ".", names = c("model", "dataset")) %>%
  mutate(model = factor(model, levels = model_names_noquote),
         dataset = factor(dataset, levels = c("Parallel",
                                              "Orthogonal",
                                              "Randomised"))) %>%
  group_by(model, dataset) %>% 
  group_map(~tbl_graph(nodes = all_molcomp_features,
                        edges = tibble(
                          from = as.integer(
                          factor(.x$from, levels = names(all_molcomp_features))),
                          to = as.integer(
                            factor(.x$to, levels = names(all_molcomp_features)))
                          ),
                       directed = F) %>%
              activate("edges") %>%
              mutate(weight = .x$weight) %>%
              activate("nodes") %>%
              mutate(degree = centrality_degree(weights = .x$weight,
                                                normalized = T),
                     betweenness = centrality_betweenness(directed = F,
                                                          normalized = T),
                     eigen = centrality_eigen(weights = .x$weight)))

group_names <- co_all %>%
  separate_wider_delim(group, ".", names = c("model", "dataset")) %>%
  mutate(model = factor(model, levels = model_names_noquote),
         dataset = factor(dataset, levels = c("Parallel",
                                              "Orthogonal",
                                              "Randomised"))) %>%
  group_by(model, dataset) %>% 
  group_keys() %>%
  select(model, dataset)
  
names(co_all_graph) <- interaction(group_names$model, group_names$dataset)


#Separate figure for each motif
co_figs <- lapply(co_all_graph, function(x)
  {
  x <- x %>%
    activate("nodes") %>%
    filter(degree > 0)
  ggraph(x, layout = "fr") +
    geom_edge_parallel(aes(colour = weight),  start_cap = circle(0.2),
                       end_cap = circle(0.2),
                       
                       sep = unit(1, "lines")) +
    geom_node_point(shape = 21, aes(fill = degree, size = ), size = 6) +
    geom_node_text(aes(label = nodes), parse = T) +
    scale_fill_gradient(low = "#8E8FEE", high = "#CD2626") +
    #scale_size(range = c(6, 12)) +
    scale_edge_colour_gradient(low = "#0066DD", high = "#33BBAA") +
    #facet_edges(dataset~model, nrow = 3) +
    labs(fill = "Degree centrality", colour = "Association") +
    theme_graph() +
    theme(legend.position = "bottom")
}) 

plot_grid(plotlist = co_figs,
          nrow = 3,
          ncol = 5,
           byrow = F)



nodes <- co_all_graph %>% select(5:6)
edges <- co_all_graph %>% select(1:4) %>%
  mutate(from = as.numeric(factor(from, levels = names(all_molcomp_features))),
         to = as.numeric(factor(to, levels = names(all_molcomp_features))))

co_all_graph <- 
  tbl_graph(nodes = nodes,
            edges = edges) %>%
  activate("edges") %>%
  mutate(from_c = all_molcomp_features[from],
         to_c = all_molcomp_features[to]) 
  
co_all_graph <- co_all_graph  %>%
  activate("edges") %>%
  mutate(model = factor(model, levels = model_names_noquote),
         dataset = factor(dataset, levels = c("Parallel",
                                            "Orthogonal",
                                            "Randomised")))
co_all_graph <- co_all_graph %>%
  activate("nodes") %>%
  mutate(name = all_molcomp_features[name],
         degree = centrality_degree(),
         betweenness = centrality_betweenness())

ggraph(co_all_graph) +
  geom_edge_parallel(aes(colour = log(weight)),  start_cap = circle(0.2),
                     end_cap = circle(0.2),
                     arrow = arrow(length = unit(0.5, "lines")),
                                   sep = unit(1, "lines")) +
  geom_node_point(shape = 21, aes(size = degree, fill = betweenness)) +
  geom_node_text(aes(label = name), parse = T) +
  scale_fill_gradient(low = "#8E8FEE", high = "#CD2626") +
  scale_size(range = c(6, 12)) +
  scale_edge_colour_gradient(low = "#0066DD", high = "#33BBAA") +
  facet_edges(dataset~model, nrow = 3) +
  labs(fill = "Betweenness", size = "Degree centrality") +
  theme_graph() +
  theme(legend.position = "bottom")
  


d_int_cooccurrence <- d_mc_cev %>%
  group_by(seed, modelindex, dataset) %>%
  mutate(int = 1 - aut,
         int_fac = cut(int, breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       right = F)) %>%
  ungroup() %>%
  filter(isAdapted == "Adapted")

net_cmx <- crossprod(table(d_int_cooccurrence %>% select(dataset, model, molComp, int_fac)))


d_mc_int_sum <- d_mc_cev %>%
  group_by(seed, modelindex, dataset) %>%
  mutate(total = sum(1 - aut, na.rm = T),
         contribution = (1 - aut) / total) %>%
  ungroup()
  # group_by(model, dataset, isAdapted, molComp) %>%
  # summarise(meanInt = mean(1 - aut),
  #           CIInt = CI(1 - aut),
  #           meanContrib = mean(contribution),
  #           CIContrib = CI(contribution))



ggplot(d_mc_int_sum %>% filter(isAdapted == "Adapted"),
       aes(x = molComp, y = contribution, colour = model)) +
  facet_nested("Selection/trait alignment" + dataset ~ "Motif" + model) +
  geom_boxplot() +
  #geom_errorbar(aes(ymin = meanContrib - CIContrib, ymax = meanContrib + CIContrib)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5)) +
  scale_x_discrete(labels = function(x) parse(text = all_molcomp_features[x])) +
  labs(x = "Molecular component", y = "Contribution to integration") +
  guides(colour = guide_none()) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom")

# Run betareg of conditional evolvability vs V_rel etc. per model
beta.cs.mc.nar <- betareg::betareg(absCS_Mb ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                                   cev_zZ * cev_h * cev_gX, 
                                   train_cs.mc.nar)
summary(beta.cs.mc.nar)
plot(beta.cs.mc.nar)

#seed <- sample(1:.Machine$integer.max, 1)
# seed
# [1] 162176257
seed <- 162176257


### NAR
MCEffects_NAR <- CalculateMCEffects(d_btgb_cev_mc_split[[1]], 
                   list(as.formula(absCS_Mb ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                     cev_zZ * cev_h * cev_gX),
                     as.formula(vrel_g ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                                  cev_zZ * cev_h * cev_gX),
                     as.formula(vrel_m ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                                  cev_zZ * cev_h * cev_gX)
                     ),
                   seed)

## Cos similarity
# Beta model no different to a null model
summary(MCEffects_NAR$absCS_Mb$beta_model)
MCEffects_NAR$absCS_Mb$beta_model_lrtest

# Random forest isn't good either
MCEffects_NAR$absCS_Mb$rf
MCEffects_NAR$absCS_Mb$rf_rmse

plot(MCEffects_NAR$absCS_Mb$boruta)

## Vrel_G
# Significant beta model
summary(MCEffects_NAR$vrel_g$beta_model)
MCEffects_NAR$vrel_g$beta_model_lrtest
MCEffects_NAR$vrel_g$shapley
MCEffects_NAR$vrel_g$lmg

# Much worse random forest
(MCEffects_NAR$vrel_g$rf)
# Test prediction accuracy
MCEffects_NAR$vrel_g$rf_rmse

plot(MCEffects_NAR$vrel_g$boruta)

MCEffects_NAR$vrel_g$interact$plot()
MCEffects_NAR$vrel_g$imp$plot()
MCEffects_NAR$vrel_g$ale$plot()

## Vrel_M
# Bad beta model
summary(MCEffects_NAR$vrel_m$beta_model)
MCEffects_NAR$vrel_m$beta_model_lrtest
MCEffects_NAR$vrel_m$shapley
MCEffects_NAR$vrel_m$lmg

plot(MCEffects_NAR$vrel_m$boruta)


# Similarly bad random forest
(MCEffects_NAR$vrel_m$rf)
MCEffects_NAR$vrel_m$interact$plot()
MCEffects_NAR$vrel_m$imp$plot()
MCEffects_NAR$vrel_m$ale$plot()

### PAR
MCEffects_PAR <- CalculateMCEffects(d_btgb_cev_mc_split[[2]], 
                                    list(as.formula(absCS_Mb ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                                                      cev_zZ * cev_h * cev_gX),
                                         as.formula(vrel_g ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                                                      cev_zZ * cev_h * cev_gX),
                                         as.formula(vrel_m ~ cev_aZ * cev_bZ * cev_KZ * cev_KXZ *
                                                      cev_zZ * cev_h * cev_gX)
                                    ),
                                    seed)

# Bad beta model
summary(MCEffects_PAR$absCS_Mb$beta_model)
MCEffects_PAR$absCS_Mb$beta_model_lrtest
MCEffects_PAR$absCS_Mb$shapley
MCEffects_PAR$absCS_Mb$lmg

# Poor RF model as well
(MCEffects_PAR$absCS_Mb$rf)
MCEffects_PAR$absCS_Mb$rf_rmse
plot(MCEffects_PAR$absCS_Mb$boruta)

MCEffects_PAR$absCS_Mb$interact$plot()
MCEffects_PAR$absCS_Mb$imp$plot()
MCEffects_PAR$absCS_Mb$ale$plot()

# Significant beta model
summary(MCEffects_PAR$vrel_g$beta_model)
MCEffects_PAR$vrel_g$beta_model_lrtest
MCEffects_PAR$vrel_g$shapley
MCEffects_PAR$vrel_g$lmg

# Good RF model
(MCEffects_PAR$vrel_g$rf)
MCEffects_PAR$vrel_g$rf_rmse

plot(MCEffects_PAR$vrel_g$boruta)

MCEffects_PAR$vrel_g$interact$plot()
MCEffects_PAR$vrel_g$imp$plot()
MCEffects_PAR$vrel_g$ale$plot()

# Bad beta model
summary(MCEffects_PAR$vrel_m$beta_model)
MCEffects_PAR$vrel_m$beta_model_lrtest
MCEffects_PAR$vrel_m$shapley
MCEffects_PAR$vrel_m$lmg

# Bad RF model
(MCEffects_PAR$vrel_m$rf)
(MCEffects_PAR$vrel_m$rf_rmse)

plot(MCEffects_PAR$vrel_m$boruta)
MCEffects_PAR$vrel_m$interact$plot()
MCEffects_PAR$vrel_m$imp$plot()
MCEffects_PAR$vrel_m$ale$plot()


## FFLC1
MCEffects_FFLC1 <- CalculateMCEffects(d_btgb_cev_mc_split[[3]], 
                                    list(as.formula(absCS_Mb ~ cev_aY * cev_bY * cev_KY * 
                                                      cev_aZ * cev_bZ * cev_KXZ *
                                                      cev_zZ * cev_h * cev_gX),
                                         as.formula(vrel_g ~ cev_aY * cev_bY * cev_KY * 
                                                      cev_aZ * cev_bZ * cev_KXZ *
                                                      cev_zZ * cev_h * cev_gX),
                                         as.formula(vrel_m ~ cev_aY * cev_bY * cev_KY * 
                                                      cev_aZ * cev_bZ * cev_KXZ *
                                                      cev_zZ * cev_h * cev_gX)
                                    ),
                                    seed)

 # Beta model could not be built (too few adapted pops?)
summary(MCEffects_FFLC1$absCS_Mb$beta_model)
MCEffects_FFLC1$absCS_Mb$beta_model_lrtest
MCEffects_FFLC1$absCS_Mb$shapley
MCEffects_FFLC1$absCS_Mb$lmg

# Poor RF model
(MCEffects_FFLC1$absCS_Mb$rf)
MCEffects_FFLC1$absCS_Mb$rf_rmse
plot(MCEffects_FFLC1$absCS_Mb$boruta)



# Beta could not be fit
summary(MCEffects_FFLC1$vrel_g$beta_model)
MCEffects_FFLC1$vrel_g$beta_model_lrtest
MCEffects_FFLC1$vrel_g$shapley
MCEffects_FFLC1$vrel_g$lmg

# Bad RF model
(MCEffects_FFLC1$vrel_g$rf)
MCEffects_FFLC1$vrel_g$rf_rmse
plot(MCEffects_FFLC1$vrel_g$boruta)

MCEffects_FFLC1$vrel_g$interact$plot()
MCEffects_FFLC1$vrel_g$imp$plot()
MCEffects_FFLC1$vrel_g$ale$plot()

# Could not be fit
summary(MCEffects_FFLC1$vrel_m$beta_model)
MCEffects_FFLC1$vrel_m$beta_model_lrtest
MCEffects_FFLC1$vrel_m$shapley
MCEffects_FFLC1$vrel_m$lmg

# Bad RF
(MCEffects_FFLC1$vrel_m$rf)
MCEffects_FFLC1$vrel_m$rf_rmse
plot(MCEffects_FFLC1$vrel_m$boruta)
MCEffects_FFLC1$vrel_m$interact$plot()
MCEffects_FFLC1$vrel_m$imp$plot()
MCEffects_FFLC1$vrel_m$ale$plot()


MCEffects_FFLI1 <- CalculateMCEffects(d_btgb_cev_mc_split[[4]], 
                                      list(as.formula(absCS_Mb ~ cev_aY * cev_bY * cev_KY * 
                                                        cev_aZ * cev_bZ * cev_KXZ *
                                                        cev_zZ * cev_h * cev_gX),
                                           as.formula(vrel_g ~ cev_aY * cev_bY * cev_KY * 
                                                        cev_aZ * cev_bZ * cev_KXZ *
                                                        cev_zZ * cev_h * cev_gX),
                                           as.formula(vrel_m ~ cev_aY * cev_bY * cev_KY * 
                                                        cev_aZ * cev_bZ * cev_KXZ *
                                                        cev_zZ * cev_h * cev_gX)
                                      ),
                                      seed)

# Could not be fit
summary(MCEffects_FFLI1$absCS_Mb$beta_model)
MCEffects_FFLI1$absCS_Mb$beta_model_lrtest
MCEffects_FFLI1$absCS_Mb$shapley
MCEffects_FFLI1$absCS_Mb$lmg

# Poor RF fit
(MCEffects_FFLI1$absCS_Mb$rf)
MCEffects_FFLI1$absCS_Mb$rf_rmse
plot(MCEffects_FFLI1$absCS_Mb$boruta)

MCEffects_FFLI1$absCS_Mb$interact$plot()
MCEffects_FFLI1$absCS_Mb$imp$plot()
MCEffects_FFLI1$absCS_Mb$ale$plot()

# Could not be fit
summary(MCEffects_FFLI1$vrel_g$beta_model)
MCEffects_FFLI1$vrel_g$beta_model_lrtest
MCEffects_FFLI1$vrel_g$shapley
MCEffects_FFLI1$vrel_g$lmg

# Reasonable RF fit
(MCEffects_FFLI1$vrel_g$rf)
MCEffects_FFLI1$vrel_g$rf_rmse
plot(MCEffects_FFLI1$vrel_g$boruta)
MCEffects_FFLI1$vrel_g$interact$plot()
MCEffects_FFLI1$vrel_g$imp$plot()
MCEffects_FFLI1$vrel_g$ale$plot()

# Could not be fit
summary(MCEffects_FFLI1$vrel_m$beta_model)
MCEffects_FFLI1$vrel_m$beta_model_lrtest
MCEffects_FFLI1$vrel_m$shapley
MCEffects_FFLI1$vrel_m$lmg

# Poor RF fit
(MCEffects_FFLI1$vrel_m$rf)
MCEffects_FFLI1$vrel_m$rf_rmse
plot(MCEffects_FFLI1$vrel_m$boruta)
MCEffects_FFLI1$vrel_m$interact$plot()
MCEffects_FFLI1$vrel_m$imp$plot()
MCEffects_FFLI1$vrel_m$ale$plot()


MCEffects_FFBH <- CalculateMCEffects(d_btgb_cev_mc_split[[5]], 
                                      list(as.formula(absCS_Mb ~ cev_aX * cev_KZX * cev_aY * 
                                                        cev_bY * cev_KY * 
                                                        cev_aZ * cev_bZ * cev_KXZ *
                                                        cev_zZ * cev_h * cev_gX),
                                           as.formula(vrel_g ~ cev_aX * cev_KZX * cev_aY * 
                                                        cev_bY * cev_KY * 
                                                        cev_aZ * cev_bZ * cev_KXZ *
                                                        cev_zZ * cev_h * cev_gX),
                                           as.formula(vrel_m ~ cev_aX * cev_KZX * cev_aY * 
                                                        cev_bY * cev_KY * 
                                                        cev_aZ * cev_bZ * cev_KXZ *
                                                        cev_zZ * cev_h * cev_gX)
                                      ),
                                      seed)

# Could not be fit
summary(MCEffects_FFBH$absCS_Mb$beta_model)
MCEffects_FFBH$absCS_Mb$beta_model_lrtest
MCEffects_FFBH$absCS_Mb$shapley
MCEffects_FFBH$absCS_Mb$lmg

# Poor fit
(MCEffects_FFBH$absCS_Mb$rf)
MCEffects_FFBH$absCS_Mb$rf_rmse
plot(MCEffects_FFBH$absCS_Mb$boruta)
MCEffects_FFBH$absCS_Mb$interact$plot()
MCEffects_FFBH$absCS_Mb$imp$plot()
MCEffects_FFBH$absCS_Mb$ale$plot()

# Could not be fit
summary(MCEffects_FFBH$vrel_g$beta_model)
MCEffects_FFBH$vrel_g$beta_model_lrtest
MCEffects_FFBH$vrel_g$shapley
MCEffects_FFBH$vrel_g$lmg

# Good fit
(MCEffects_FFBH$vrel_g$rf)
MCEffects_FFBH$vrel_g$rf_rmse
plot(MCEffects_FFBH$vrel_g$boruta)

MCEffects_FFBH$vrel_g$interact$plot()
MCEffects_FFBH$vrel_g$imp$plot()
MCEffects_FFBH$vrel_g$ale$plot()

# Could not be fit
summary(MCEffects_FFBH$vrel_m$beta_model)
MCEffects_FFBH$vrel_m$beta_model_lrtest
MCEffects_FFBH$vrel_m$shapley
MCEffects_FFBH$vrel_m$lmg

# Poor fit
(MCEffects_FFBH$vrel_m$rf)
MCEffects_FFBH$vrel_m$rf_rmse
plot(MCEffects_FFBH$vrel_m$boruta)

MCEffects_FFBH$vrel_m$interact$plot()
MCEffects_FFBH$vrel_m$imp$plot()
MCEffects_FFBH$vrel_m$ale$plot()


# No clear correlation between conditional evolvabilities and the cosine similarity
# between M and beta - maybe to be expected?
# Same with Vrel_m

# Vrel_g seems to be most consistently predictable by the conditional evolvabilities
# 

# Report three most important components per boruta, use that to fit a beta model?
# Doesn't really improve things, the data is just too noisy - not good predictors

# Could also just report means of cev for each model, as we know that these evolvabilities
# are components of G anyway
# Could also measure cev over a trait




