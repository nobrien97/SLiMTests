# Random forest model of adaptedness
library(tidyverse)

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

# Combine
d_vrel_tot <- rbind(d_vrel, d_vrel_orth, d_vrel_par)

######################################################################
# Load optima
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

# Measure cosine similarity between M matrices and beta
d_opt <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)
d_opt <- read_csv("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_opt.csv", col_names = F)
d_opt <- read_csv("/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_opt.csv", col_names = F)

# o = optimum, s = sigma, d = direction (-1, 1)
colnames(d_opt) <- c("seed", "modelindex", "o_t1", "o_t2", "o_t3", "o_t4", 
                     "s_t1", "s_t2", "s_t3", "s_t4", "d_t1", "d_t2", "d_t3",
                     "d_t4")





# Now G matrix
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
d_selvec_m <- d_qg %>%
  filter(gen >= 50000) %>%
  select(gen, seed, modelindex, isAdapted, ends_with("mean"))

d_selvec_m <- left_join(d_selvec_m, d_opt %>% 
                          select(seed, modelindex, starts_with("o_")) %>%
                          mutate(seed = factor(seed),
                                 modelindex = factor(modelindex)), 
                        by = c("seed", "modelindex"))

d_selvec_m$dataset <- "Randomised"

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

d_selvec_m_tot <- rbind(d_selvec_m, d_selvec_orth, d_selvec_par)

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

d_selvec_m_tot <- inner_join(id_m_tot, d_selvec_m, 
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
saveRDS(d_cossim_m_tot, "d_cossim_m_datasets.RDS")
d_cossim_m_tot <- readRDS("d_cossim_m_datasets.RDS")

