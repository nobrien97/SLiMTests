source("helperFn.R")

COMBO_PATH <- '/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
d_combos <- read_delim(COMBO_PATH, 
                       delim = " ", col_names = F)
names(d_combos) <- c("model", "r")


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

saveRDS(d_h2, "d_h2.RDS")
saveRDS(d_h2_trait, "d_h2_trait.RDS")

d_h2 <- readRDS("d_h2.RDS")
d_h2_trait <- readRDS("d_h2_trait.RDS")


# Join qg together
d_qg$dataset <- "Randomised"
d_qg_orth$dataset <- "Orthogonal"
d_qg_par$dataset <- "Parallel"

d_qg_tot <- rbind(d_qg, d_qg_orth, d_qg_par)

saveRDS(d_qg_tot, "d_qg_tot.RDS")

d_qg_tot <- readRDS("/mnt/i/SLiMTests/tests/newMotifs/paper1/d_qg_tot.RDS")

d_qg_optPerc <- d_qg_tot %>% select(gen, seed, modelindex, dataset, isAdapted) %>% filter(gen >= 49500)

# inner join optPerc
d_h2_trait <- left_join(d_h2_trait, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))

d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))

# Discretise generation
d_h2_trait <- d_h2_trait %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))

d_h2 <- d_h2 %>%
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
id$model <- factor(id$model, levels = model_names, labels = model_names_noquote)

# Resize h2_pd according to n traits
h2_mat <- purrr::map(seq_along(h2_mat), function(i) {
  mat <- h2_mat[[i]]
  
  n <- GetMotifTraitRange(as.character(id$model[i]))
  
  mat[n, n]
})
  

# First convert to nearest positive definite matrix
h2_pd <- lapply(h2_mat, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

saveRDS(h2_pd, "h2_pd_trait.RDS")
saveRDS(id, "id_trait.RDS")

# Get mean G matrix
d_h2_trait %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(5:8, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(timePoint, model, dataset, isAdapted, .keep_all = T) %>%
  group_map(~ .x %>% 
              summarise(across(starts_with("VA") | starts_with("CVA"), 
                               list(mean), .names = "{.col}")) %>%
              mutate(timePoint = .y$timePoint,
                     model = .y$model,
                     dataset = .y$dataset,
                     isAdapted = .y$isAdapted),
            .keep = T) -> split_h2_mean

lapply(split_h2_mean, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_mean_matrices

h2_mean_mat <- unlist(cov_mean_matrices, recursive = F)

cov_mean_matrix_modelindex <- GetMeanMatrixIDs(split_h2_mean)
id_mean <- data.table::rbindlist(cov_mean_matrix_modelindex, 
                                 fill = T)
id_mean$label <- 1:nrow(id_mean)

# Resize h2_pd according to n traits
h2_mean_mat <- purrr::map(seq_along(h2_mean_mat), function(i) {
  mat <- h2_mean_mat[[i]]
  
  n <- GetMotifTraitRange(as.character(id_mean$model[i]))
  
  mat[n, n]
})

# First convert to nearest positive definite matrix
h2_mean_pd <- lapply(h2_mean_mat, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

eig_mean_pd <- purrr::map(seq_along(h2_mean_pd), function(i){
  eig <- eigen(h2_mean_pd[[i]])
  
  circle <- cbind(cos(seq(0, 2*pi, length.out = 100)),
                  sin(seq(0, 2*pi, length.out = 100)))
  
  ellipse_pc <- circle %*% diag(sqrt(eig$values[1:2])) %*% t(eig$vectors[,1:2])
  
  # Project onto ellipse for the major/minor axes
  V12 <- eig$vectors[1:2,1:2]
  
  Sigma_proj <- V12 %*%
    diag(eig$values[1:2]) %*%
    t(V12)
  
  eig_proj <- eigen(Sigma_proj)
  
  result <- id_mean[rep(i, times = 100),]
  result$pc1 <- ellipse_pc[,1]
  result$pc2 <- ellipse_pc[,2]
  result$pc1_exp <- eig$values[1] / sum(eig$values)
  result$pc2_exp <- eig$values[2] / sum(eig$values)
  result$a <- sqrt(eig$values[1])
  result$b <- sqrt(eig$values[2])
  major <- sqrt(eig_proj$values[1]) * eig_proj$vectors[,1]
  minor <- sqrt(eig_proj$values[2]) * eig_proj$vectors[,2]
  
  result$major_x <- major[1]
  result$major_y <- major[2]
  result$minor_x <- minor[1]
  result$minor_y <- minor[2]
  
  result
})

d_pc_mean <- data.table::rbindlist(eig_mean_pd, 
                                   fill = T)

d_pc_mean$model <- factor(d_pc_mean$model, 
                          levels = model_names,
                          labels = model_names_noquote)

d_pc_mean$dataset <- factor(d_pc_mean$dataset,
                            levels = c("Parallel",
                                       "Orthogonal",
                                       "Randomised"))

d_pc_mean_labels <- d_pc_mean %>%
  group_by(model, dataset, timePoint, isAdapted) %>%
  summarise(pc1_exp = pc1_exp[1],
            pc2_exp = pc2_exp[2])

d_pc_mean_lines <- d_pc_mean %>%
  group_by(model, dataset, timePoint, isAdapted) %>%
  slice_head(n = 1)

ggplot(d_pc_mean %>% filter(isAdapted == T),
       aes(x = pc1, y = pc2, colour = timePoint,
           group = label)) +
  facet_nested("Model" + model ~ 
                 "Trait/selection alignment / Measurement time" + dataset + timePoint) +
  geom_path() +
  geom_segment(data = d_pc_mean_lines %>% filter(isAdapted == T), inherit.aes = F,
              aes(x = -minor_x, xend = minor_x, y = -minor_y, yend = minor_y,
                  colour = timePoint, group = label)) +
  geom_segment(data = d_pc_mean_lines %>% filter(isAdapted == T), inherit.aes = F,
               aes(x = -major_x, xend = major_x, y = -major_y, yend = major_y,
                   colour = timePoint, group = label)) +
  geom_label(data = d_pc_mean_labels %>% filter(isAdapted == T), inherit.aes = F,
             aes(x = 0, y = 0.5, 
                 label = paste0("PC1: ", round(pc1_exp * 100, digits = 2), "%"))) +
  geom_label(data = d_pc_mean_labels %>% filter(isAdapted == T), inherit.aes = F,
             aes(x = 0, y = -0.5, 
                 label = paste0("PC2: ", round(pc2_exp * 100, digits = 2), "%"))) +
  scale_colour_manual(values = c("#001889FF", "#E98935FF")) +
  coord_fixed() +
  labs(x = "PC1", y = "PC2", colour = "Time") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))
ggsave("plt_gmat_trait.png", device = png, width = 17, height = 9)

#################################################
# Repeat with mc matrices
d_h2 %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(6:17, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(timePoint, modelindex, dataset, isAdapted) %>%
  group_split(.) -> split_h2_mc

lapply(split_h2_mc, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_mc


# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat_mc <- unlist(cov_matrices_mc, recursive = F)

# get ids from the matrix
cov_matrix_modelindex_mc <- GetMatrixIDsWithDataset(split_h2_mc)

id_mc <- data.table::rbindlist(cov_matrix_modelindex_mc, 
                            fill = T)
id_mc$label <- as.character(1:nrow(id_mc))
id_mc$modelindex <- as.factor(id_mc$modelindex)
id_mc <- AddCombosToDF(id_mc)
id_mc$model <- factor(id_mc$model, levels = model_names, labels = model_names_noquote)

# Resize h2_pd according to n traits
h2_mat_mc <- purrr::map(seq_along(h2_mat_mc), function(i) {
  mat <- h2_mat_mc[[i]]
  
  n <- molComp_names[[as.character(id_mc$model[i])]]
  mat[n, n]
})


# First convert to nearest positive definite matrix
h2_pd_mc <- lapply(h2_mat_mc, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

saveRDS(h2_pd_mc, "h2_pd_mc.RDS")
saveRDS(id_mc, "id_trait_mc.RDS")



# Mean matrix
d_h2 %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(6:17, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(timePoint, model, dataset, isAdapted, .keep_all = T) %>%
  group_map(~ .x %>% 
              summarise(across(starts_with("VA") | starts_with("CVA"), 
                               list(mean), .names = "{.col}")) %>%
              mutate(timePoint = .y$timePoint,
                     model = .y$model,
                     dataset = .y$dataset,
                     isAdapted = .y$isAdapted),
            .keep = T) -> split_h2_mc_mean

lapply(split_h2_mc_mean, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_mean_matrices_mc

h2_mean_mat_mc <- unlist(cov_mean_matrices_mc, recursive = F)

cov_mean_matrix_modelindex_mc <- GetMeanMatrixIDs(split_h2_mc_mean)
id_mean_mc <- data.table::rbindlist(cov_mean_matrix_modelindex_mc, 
                                 fill = T)
id_mean_mc$label <- 1:nrow(id_mean_mc)

# Resize h2_pd according to n traits
h2_mean_mat_mc <- purrr::map(seq_along(h2_mean_mat_mc), function(i) {
  mat <- h2_mean_mat_mc[[i]]
  
  n <- GetMotifTraitRange(as.character(id_mean_mc$model[i]))
  
  mat[n, n]
})

# First convert to nearest positive definite matrix
h2_mean_pd_mc <- lapply(h2_mean_mat_mc, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

eig_mean_pd_mc <- purrr::map(seq_along(h2_mean_pd_mc), function(i){
  eig <- eigen(h2_mean_pd_mc[[i]])
  
  circle <- cbind(cos(seq(0, 2*pi, length.out = 100)),
                  sin(seq(0, 2*pi, length.out = 100)))
  
  ellipse_pc <- circle %*% diag(sqrt(eig$values[1:2])) %*% t(eig$vectors[,1:2])
  
  # Project onto ellipse for the major/minor axes
  V12 <- eig$vectors[1:2,1:2]
  
  Sigma_proj <- V12 %*%
    diag(eig$values[1:2]) %*%
    t(V12)
  
  eig_proj <- eigen(Sigma_proj)
  
  result <- id_mean_mc[rep(i, times = 100),]
  result$pc1 <- ellipse_pc[,1]
  result$pc2 <- ellipse_pc[,2]
  result$pc1_exp <- eig$values[1] / sum(eig$values)
  result$pc2_exp <- eig$values[2] / sum(eig$values)
  result$a <- sqrt(eig$values[1])
  result$b <- sqrt(eig$values[2])
  major <- sqrt(eig_proj$values[1]) * eig_proj$vectors[,1]
  minor <- sqrt(eig_proj$values[2]) * eig_proj$vectors[,2]
  
  result$major_x <- major[1]
  result$major_y <- major[2]
  result$minor_x <- minor[1]
  result$minor_y <- minor[2]
  
  result
})

d_pc_mean_mc <- data.table::rbindlist(eig_mean_pd_mc, 
                                   fill = T)

d_pc_mean_mc$model <- factor(d_pc_mean_mc$model, 
                          levels = model_names,
                          labels = model_names_noquote)

d_pc_mean_mc$dataset <- factor(d_pc_mean_mc$dataset,
                            levels = c("Parallel",
                                       "Orthogonal",
                                       "Randomised"))

d_pc_mean_mc_labels <- d_pc_mean_mc %>%
  group_by(model, dataset, timePoint, isAdapted) %>%
  summarise(pc1_exp = pc1_exp[1],
            pc2_exp = pc2_exp[2])

d_pc_mean_mc_lines <- d_pc_mean_mc %>%
  group_by(model, dataset, timePoint, isAdapted) %>%
  slice_head(n = 1)

ggplot(d_pc_mean_mc %>% filter(isAdapted == T),
       aes(x = pc1, y = pc2, colour = timePoint,
           group = label)) +
  facet_nested("Model" + model ~ 
                 "Trait/selection alignment / Measurement time" + dataset + timePoint) +
  geom_path() +
  geom_segment(data = d_pc_mean_mc_lines %>% filter(isAdapted == T), inherit.aes = F,
               aes(x = -minor_x, xend = minor_x, y = -minor_y, yend = minor_y,
                   colour = timePoint, group = label)) +
  geom_segment(data = d_pc_mean_mc_lines %>% filter(isAdapted == T), inherit.aes = F,
               aes(x = -major_x, xend = major_x, y = -major_y, yend = major_y,
                   colour = timePoint, group = label)) +
  geom_label(data = d_pc_mean_mc_labels %>% filter(isAdapted == T), inherit.aes = F,
             aes(x = 0, y = 0.5, 
                 label = paste0("PC1: ", round(pc1_exp * 100, digits = 2), "%"))) +
  geom_label(data = d_pc_mean_mc_labels %>% filter(isAdapted == T), inherit.aes = F,
             aes(x = 0, y = -0.5, 
                 label = paste0("PC2: ", round(pc2_exp * 100, digits = 2), "%"))) +
  scale_colour_manual(values = c("#001889FF", "#E98935FF")) +
  coord_fixed() +
  labs(x = "PC1", y = "PC2", colour = "Time") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))
ggsave("plt_gmat_mc.png", device = png, width = 17, height = 9)


# Now for M matrices
m_matrices <- readRDS("m_matrices.RDS")
m_matrices_orth <- readRDS("m_matrices_orth.RDS")
m_matrices_par <- readRDS("m_matrices_par.RDS")

m_matrices_tot <- c(m_matrices, m_matrices_orth, m_matrices_par)


DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar.csv"

d_m <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                         paste0("mean_", 1:4),
                                         paste0("var_", 1:4),
                                         paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m <- d_m %>%
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex))




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
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex))

d_m_par <- d_m_par %>%
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex))

d_m$dataset <- "Randomised"
d_m_orth$dataset <- "Orthogonal"
d_m_par$dataset <- "Parallel"

d_m_tot <- rbind(d_m, d_m_orth, d_m_par)

# Average
d_m_valid <- inner_join(d_m_tot %>%
                         mutate(model = factor(model, levels = model_names_noquote),
                                seed = factor(seed),
                                modelindex = factor(modelindex)),
                       d_qg_tot %>%
                         select(gen, seed, modelindex, r, model, dataset, isAdapted, timeToAdapt) %>%
                         mutate(model = factor(model, levels = model_names,
                                               labels = model_names_noquote))) %>%
  filter(log10(r) == -1,
         isAdapted == T,
         gen == 50000 | gen == 60000)

d_m_valid %>%
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  mutate(timePoint = factor(if_else(gen == 50000, "Start", "End"))) %>%
  group_by(gen, timePoint, model, dataset, isAdapted) %>%
  summarise(across(starts_with("var") | starts_with("cov"), 
                             list(mean), .names = "{.col}")) -> d_m_mean

# Means
m_matrices_valid_mean <- d_m_mean %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))


id_m_valid_mean <- d_m_mean %>%
  select(1:5)

# Non-means
m_matrices_valid <- d_m_valid %>%
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  mutate(timePoint = factor(if_else(gen == 50000, "Start", "End"))) %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))

id_m_valid <- d_m_valid %>%
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  mutate(timePoint = factor(if_else(gen == 50000, "Start", "End"))) %>%
  select(timePoint, seed, modelindex, isAdapted, dataset, model, r) %>%
  mutate(label = 1:n())

m_matrices_valid_pd <- lapply(m_matrices_valid, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})



# First convert to nearest positive definite matrix
m_matrices_valid_mean_pd <- lapply(m_matrices_valid_mean, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})


eig_mean_m_pd <- purrr::map(seq_along(m_matrices_valid_mean_pd), function(i){
  eig <- eigen(m_matrices_valid_mean_pd[[i]])
  
  circle <- cbind(cos(seq(0, 2*pi, length.out = 100)),
                  sin(seq(0, 2*pi, length.out = 100)))
  
  ellipse_pc <- circle %*% diag(sqrt(eig$values[1:2])) %*% t(eig$vectors[,1:2])
  
  # Project onto ellipse for the major/minor axes
  V12 <- eig$vectors[1:2,1:2]
  
  Sigma_proj <- V12 %*%
    diag(eig$values[1:2]) %*%
    t(V12)
  
  eig_proj <- eigen(Sigma_proj)
  
  result <- id_m_valid_mean[rep(i, times = 100),]
  result$label <- i
  result$pc1 <- ellipse_pc[,1]
  result$pc2 <- ellipse_pc[,2]
  result$pc1_exp <- eig$values[1] / sum(eig$values)
  result$pc2_exp <- eig$values[2] / sum(eig$values)
  result$a <- sqrt(eig$values[1])
  result$b <- sqrt(eig$values[2])
  major <- sqrt(eig_proj$values[1]) * eig_proj$vectors[,1]
  minor <- sqrt(eig_proj$values[2]) * eig_proj$vectors[,2]
  
  result$major_x <- major[1]
  result$major_y <- major[2]
  result$minor_x <- minor[1]
  result$minor_y <- minor[2]
  
  result
}, .progress = T)

d_pc_m_mean <- data.table::rbindlist(eig_mean_m_pd, 
                                      fill = T)

d_pc_m_mean$dataset <- factor(d_pc_m_mean$dataset,
                               levels = c("Parallel",
                                          "Orthogonal",
                                          "Randomised"))

d_pc_m_mean$timePoint <- factor(d_pc_m_mean$timePoint,
                                levels = c("Start",
                                           "End"))


d_pc_mean_m_labels <- d_pc_m_mean %>%
  group_by(model, dataset, timePoint, isAdapted) %>%
  summarise(pc1_exp = pc1_exp[1],
            pc2_exp = pc2_exp[2])

d_pc_mean_m_lines <- d_pc_m_mean %>%
  group_by(model, dataset, timePoint, isAdapted) %>%
  slice_head(n = 1)

ggplot(d_pc_m_mean %>% filter(isAdapted == T),
       aes(x = pc1, y = pc2, colour = timePoint,
           group = label)) +
  facet_nested("Model" + model ~ 
                 "Trait/selection alignment / Measurement time" + dataset + timePoint,
               scales = "free") +
  geom_path() +
  geom_segment(data = d_pc_mean_m_lines %>% filter(isAdapted == T), inherit.aes = F,
               aes(x = -minor_x, xend = minor_x, y = -minor_y, yend = minor_y,
                   colour = timePoint, group = label)) +
  geom_segment(data = d_pc_mean_m_lines %>% filter(isAdapted == T), inherit.aes = F,
               aes(x = -major_x, xend = major_x, y = -major_y, yend = major_y,
                   colour = timePoint, group = label)) +
  # geom_label(data = d_pc_mean_m_labels %>% filter(isAdapted == T), inherit.aes = F,
  #            aes(x = 0, y = 0.5, 
  #                label = paste0("PC1: ", round(pc1_exp * 100, digits = 2), "%"))) +
  # geom_label(data = d_pc_mean_m_labels %>% filter(isAdapted == T), inherit.aes = F,
  #            aes(x = 0, y = -0.5, 
  #                label = paste0("PC2: ", round(pc2_exp * 100, digits = 2), "%"))) +
  scale_colour_manual(values = c("#001889FF", "#E98935FF")) +
  #coord_fixed() +
  labs(x = "PC1", y = "PC2", colour = "Time") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))
ggsave("plt_mmat.png", device = png, width = 17, height = 9)


# Vrel_GT, Vrel_Gc, Vrel_M
## Measure how isomorphically variation is distributed across traits/components
vrel_gt <- purrr::map(seq_along(h2_pd), function(i) {
  g <- h2_pd[[i]]
  
  vrel_i <- Vrel(eigen(g)$values)
  
  result <- id[i,]
  
  result$vrel <- vrel_i
  return(result)
})

d_vrel_gt <- data.table::rbindlist(vrel_gt)

vrel_gc <- purrr::map(seq_along(h2_pd_mc), function(i) {
  g <- h2_pd_mc[[i]]
  
  vrel_i <- Vrel(eigen(g)$values)
  
  result <- id_mc[i,]
  
  result$vrel <- vrel_i
  return(result)
})
d_vrel_gc <- data.table::rbindlist(vrel_gc)

vrel_m <- purrr::map(seq_along(m_matrices_valid_pd), function(i) {
  g <- m_matrices_valid_pd[[i]]
  
  vrel_i <- Vrel(eigen(g)$values)
  
  result <- id_m_valid[i,]
  
  result$vrel <- vrel_i
  return(result)
})
d_vrel_m <- data.table::rbindlist(vrel_m)

# Combine
d_vrel_gt$matrix <- TeX("$G_T$", output = "character")
d_vrel_gc$matrix <- TeX("$G_C$", output = "character")
d_vrel_m$matrix <- TeX("$M$", output = "character")

d_vrel <- rbind(d_vrel_gt, d_vrel_gc, d_vrel_m)
d_vrel$dataset <- factor(d_vrel$dataset, levels = c("Parallel",
                                                    "Orthogonal",
                                                    "Randomised"))

d_vrel_sum <- d_vrel %>%
  filter(isAdapted == T, log10(r) == -1) %>%
  group_by(model, matrix, dataset) %>%
  summarise(meanVrel = mean(vrel))

ggplot(d_vrel %>% 
         filter(isAdapted == T, log10(r) == -1),
       aes(x = model, y = vrel, colour = model)) +
  facet_nested("Matrix" + matrix ~ "Trait/selection alignment" + dataset, 
               labeller = labeller(matrix = label_parsed)) +
  geom_quasirandom(show.legend = F) +
  geom_point(data = d_vrel_sum, aes(y = meanVrel), colour = "black", 
             fill = "white", stroke = 1,
             shape = 21, inherit.aes = T, show.legend = F) +
  scale_colour_manual(values = pal) +
  labs(x = "Model", y = TeX("$V_{rel}")) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_vrel.png", device = png, bg = "white", width = 10, height = 8)


# bTGb, bTCb, bTMb
## Measure how much variance is along the selection gradient
d_selvec <- readRDS("/mnt/i/SLiMTests/tests/newMotifs/paper1/d_selvec.RDS")
d_selvec$model <- factor(d_selvec$model, levels = model_names, labels = model_names_noquote)

d_selvec2_gt <- inner_join(id, d_selvec, 
                        by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim_gt <- GetCosineSimilarity(h2_pd, d_selvec2_gt %>% select(ends_with("dir")), id)

d_selvec2_mc <- inner_join(id_mc, d_selvec, 
                           by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim_gc <- GetCosineSimilarity(h2_pd_mc, d_selvec2_mc %>% select(ends_with("dir")), id_mc)

d_selvec2_m <- inner_join(id_m_valid, d_selvec, 
                           by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim_m <- GetCosineSimilarity(m_matrices_valid_pd, d_selvec2_m %>% select(ends_with("dir")), id_m_valid)

# Combine
d_cossim_gt$matrix <- TeX("$G_T$", output = "character")
d_cossim_gc$matrix <- TeX("$G_C$", output = "character")
d_cossim_m$matrix <- TeX("$M$", output = "character")

d_cossim <- rbind(d_cossim_gt, d_cossim_gc, d_cossim_m)

d_cossim <- AddCombosToDF(d_cossim)

d_cossim$model <- factor(d_cossim$model, levels = model_names, labels = model_names_noquote)

d_cossim_sum <- d_cossim %>%
  filter(isAdapted == T, log10(r) == -1) %>%
  group_by(model, matrix, dataset) %>%
  summarise(meanLogbTMb = mean(log10(bTMb)))



ggplot(d_cossim %>% 
         filter(isAdapted == T, log10(r) == -1),
       aes(x = model, y = log10(bTMb), colour = model)) +
  facet_nested("Matrix" + matrix ~ "Trait/selection alignment" + dataset, 
               labeller = labeller(matrix = label_parsed)) +
  geom_quasirandom(show.legend = F) +
  geom_point(data = d_cossim_sum, aes(y = meanLogbTMb), colour = "black", 
             fill = "white", stroke = 1,
             shape = 21, inherit.aes = T, show.legend = F) +
  scale_colour_manual(values = pal) +
  labs(x = "Model", y = TeX("$log_{10}(\\beta^T X \\beta)$")) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_btxb.png", device = png, bg = "white", width = 10, height = 8)

# PCASim
## Measure how similar the covariance within/between components/traits is between replicates
krz_in_gt <- id %>%
  mutate(g = h2_pd,
         group = interaction(model, dataset, log10(r)))

krz_in_gc <- id_mc %>%
  mutate(g = h2_pd_mc,
         group = interaction(model, dataset, log10(r)))

krz_in_m <- id_m_valid %>%
  mutate(g = m_matrices_valid_pd,
         group = interaction(model, dataset, log10(r)))


# Save krz_in: run this part on HPC
saveRDS(krz_in_gt, "/mnt/d/SLiMTests/tests/newMotifs/paper1/pca_in_gt.RDS")
saveRDS(krz_in_gc, "/mnt/d/SLiMTests/tests/newMotifs/paper1/pca_in_gc.RDS")
saveRDS(krz_in_m, "/mnt/d/SLiMTests/tests/newMotifs/paper1/pca_in_m.RDS")

# Bootstrap in ten parts for RAM reasons
# This is slow: uncomment to run, otherwise read in precalculated data
# Generate seeds
# newseed <- sample(1:.Machine$integer.max, 10)
# [1]  314145285 2009911717  267335506  231424073 1190700402 1189454198  395819651  848071181 1762114410
# [10] 1739509036
newseed <- c(314145285L, 2009911717L, 267335506L, 231424073L, 1190700402L, 
             1189454198L, 395819651L, 848071181L, 1762114410L, 1739509036L)
bootPCASim <- vector(mode = "list", length = length(newseed))

# Per model inputs
krz_in_gc_NAR <- krz_in_gc %>% filter(model == "NAR", log10(r) == -1)
krz_in_gc_PAR <- krz_in_gc %>% filter(model == "PAR", log10(r) == -1)
krz_in_gc_FFLC1 <- krz_in_gc %>% filter(model == "FFLC1", log10(r) == -1)
krz_in_gc_FFLI1 <- krz_in_gc %>% filter(model == "FFLI1", log10(r) == -1)
krz_in_gc_FFBH <- krz_in_gc %>% filter(model == "FFBH", log10(r) == -1)


for (i in seq_along(newseed)) {
  # Set seed
  set.seed(newseed[i])
  # Run replicate but only within models
  res_NAR <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_NAR, "group", T))
  res_PAR <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_PAR, "group", T))
  res_FFLC1 <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_FFLC1, "group", T))
  res_FFLI1 <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_FFLI1, "group", T))
  res_FFBH <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_FFBH, "group", T))
  
  # To data.frame
  res_NAR <- unnest(as.data.frame(t(res_NAR)), cols = everything()) %>%
    mutate(model = "NAR")
  res_PAR <- unnest(as.data.frame(t(res_PAR)), cols = everything()) %>%
    mutate(model = "PAR")
  res_FFLC1 <- unnest(as.data.frame(t(res_FFLC1)), cols = everything()) %>%
    mutate(model = "FFLC1")
  res_FFLI1 <- unnest(as.data.frame(t(res_FFLI1)), cols = everything()) %>%
    mutate(model = "FFLI1")
  res_FFBH <- unnest(as.data.frame(t(res_FFBH)), cols = everything()) %>%
    mutate(model = "FFBH")
  
  # Combine to output
  bootPCASim[[i]] <- rbind(res_NAR, res_PAR, res_FFLC1, res_FFLI1, res_FFBH)
}

# Output list into combined df
bootPCASim2 <- bind_rows(bootPCASim)
bootPCASim <- bootPCASim2 %>%
  separate(group1, c("model1", "dataset1", "r1"), "\\.",
           extra = "merge") %>%
  separate(group2, c("model2", "dataset2", "r2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = as.numeric(r1),
         r2 = as.numeric(r2),
         dataset1 = factor(dataset1, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised")),
         dataset2 = factor(dataset2, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised"))) %>%
  rename(PCASim = krzCor)

saveRDS(bootPCASim, paste0("/mnt/d/SLiMTests/tests/newMotifs/paper1/d_bootPCASim.RDS"))
bootPCASim <- readRDS(paste0("/mnt/d/SLiMTests/tests/newMotifs/paper1/d_bootPCASim.RDS"))

# Plot

# Get model comparisons for labelling
bootPCASim <- bootPCASim %>%
  mutate(datasetCombo = GetModelComparison(dataset1, dataset2, c("Parallel",
                                                                 "Orthogonal",
                                                                 "Randomised")),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)),
         model = factor(model, levels = model_names_noquote),
         nMCs = unlist(lapply(molComp_names[as.character(model)], length)))

# recomb by modelCombo - we don't have all the recombination levels for the
# non-randomised datasets
bootPCASim_sum <- bootPCASim %>%
  group_by(model, dataset1, dataset2) %>%
  summarise(meanPCASim = mean(PCASim),
            ciPCASim = CI(PCASim),
            nMCs = length(molComp_names[as.character(model[1])][[1]]))

# Facet design
design <- c(
  "
  AABB
  CCDD
  #EE#
  "
)


ggplot(bootPCASim_sum, aes(
  x = (dataset1), y = (dataset2)
)) +
  #facet_nested("Model" + model ~ .) + 
  facet_manual(.~model, design = design, axes = T) + 
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = seq(0.25, 1.0, by = 0.25),
                       labels = seq(0.25, 1.0, by = 0.25),
                       limits = c(0.25, 1)) +
  # scale_fill_gradientn(colours = paletteer_c("ggthemes::Blue-Green Sequential", n = 6),
  #                      values = c(0.0, 0.7, 0.75, 0.8, 0.9, 1.0)) +
  labs(x = "Trait/selection alignment (Matrix 1)", y = "Trait selection alignment (Matrix 2)", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 12), 
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "bottom",
        legend.key.width = unit(3.5, "lines"))
ggsave("PCASim_dataset.png", device = png, width = 8, height = 10)


# Total variance explained
prop_var_mc <- purrr::map(seq_along(h2_pd_mc), function(i) {
  vars <- diag(h2_pd_mc[[i]])
  d_template <- data.frame(matrix(ncol = 13, nrow = 1))
  names(d_template) <- c(names(all_molcomp_features), "totalVar")
  d_template[1,names(vars)] <- vars / sum(vars)
  d_template[1,13] <- sum(vars)
  return(d_template)
})

# Attach ID
d_prop_vars <- data.table::rbindlist(prop_var_mc)
d_prop_vars <- cbind(d_prop_vars, id_mc)

# Pivot longer
d_prop_vars <- d_prop_vars %>%
  pivot_longer(cols = (1:12),
               names_to = "molComp",
               values_to = "varExpl")

d_prop_vars <- d_prop_vars %>%
  drop_na(varExpl) %>%
  filter(isAdapted == T) %>%
  mutate(dataset = factor(dataset, levels = c("Parallel",
                                              "Orthogonal",
                                              "Randomised")))

d_prop_vars_sum <- d_prop_vars %>%
  group_by(model, dataset, molComp) %>%
  summarise(meanVarExpl = mean(varExpl),
            CIVarExpl = CI(varExpl),
            meanVar = mean(totalVar),
            CIVar = CI(totalVar))

ggplot(d_prop_vars_sum,
       aes(x = molComp, y = meanVarExpl, fill = model)) +
  facet_nested("Model" + model ~ "Trait/selection alignment" + dataset) +
  geom_col() +
  geom_text(aes(x = 6.5, y = 0.2, label = 
                  paste("Total variance =", round(meanVar, digits = 3),
                        "±", round(CIVar, digits = 3))), size = 4) +
  scale_x_discrete(labels = function(x) {parse(text = all_molcomp_features[x])}) +
  scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = meanVarExpl - CIVarExpl, ymax = meanVarExpl + CIVarExpl), width = 0.2) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                    guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      guide = "none") +
  labs(x = "Molecular component", y = "Mean genetic variance explained (%)") +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_var_expl.png", device = png, width = 13, height = 10)

# {Print most importance features}
print(xtable::xtable(d_prop_vars_sum %>%
                       group_by(model, dataset) %>%
                       slice_max(meanVarExpl, n = 4) %>%
                       ungroup() %>%
                       select(-model),
                     
                     digits = 3), include.rownames = F)


# Repeat with unscaled variance estimates
G_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/getH2/"


d_h2_noscale_mrr <- read_csv(paste0(G_DATA_PATH, "out_h2_noscale_mrr.csv"), col_names = F)
d_h2_noscale_mkr <- read_csv(paste0(G_DATA_PATH, "out_h2_noscale_mkr.csv"), col_names = F)

colnames(d_h2_noscale_mrr) <- h2_colnames
colnames(d_h2_noscale_mkr) <- h2_colnames

d_h2_noscale_trait_mkr <- read_csv(paste0(G_DATA_PATH, "out_h2_noscale_trait_mkr.csv"), col_names = F)
d_h2_noscale_trait_mrr <- read_csv(paste0(G_DATA_PATH, "out_h2_noscale_trait_mrr.csv"), col_names = F)

colnames(d_h2_noscale_trait_mkr) <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_t1",
                              "VA_t2", "VA_t3", "VA_t4", "CVA_t1_t2", "CVA_t1_t3",
                              "CVA_t1_t4", "CVA_t2_t3", "CVA_t2_t4", "CVA_t3_t4",
                              "h2_t1", "h2_t2", "h2_t3", "h2_t4")

colnames(d_h2_noscale_trait_mrr) <- colnames(d_h2_noscale_trait_mkr)

# join
d_h2_noscale_trait_mkr$calcMode <- "mkr"
d_h2_noscale_trait_mrr$calcMode <- "mrr"


d_h2_noscale_mkr$calcMode <- "mkr"
d_h2_noscale_mrr$calcMode <- "mrr"


G_ORTH_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"
G_ORTH_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/orthSel/getH2/"
G_PAR_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/"


d_h2_noscale_mkr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_noscale_mkr.csv"), col_names = F)
d_h2_noscale_mrr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_noscale_mrr.csv"), col_names = F)

d_h2_noscale_mkr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_noscale_mkr.csv"), col_names = F)
d_h2_noscale_mrr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_noscale_mrr.csv"), col_names = F)


d_h2_noscale_trait_mkr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_noscale_trait_mkr.csv"), col_names = F)
d_h2_noscale_trait_mrr_orth <- read_csv(paste0(G_ORTH_DATA_PATH, "out_h2_noscale_trait_mrr.csv"), col_names = F)

d_h2_noscale_trait_mkr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_noscale_trait_mkr.csv"), col_names = F)
d_h2_noscale_trait_mrr_par <- read_csv(paste0(G_PAR_DATA_PATH, "out_h2_noscale_trait_mrr.csv"), col_names = F)


colnames(d_h2_noscale_trait_mkr_orth) <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_t1",
                                   "VA_t2", "VA_t3", "VA_t4", "CVA_t1_t2", "CVA_t1_t3",
                                   "CVA_t1_t4", "CVA_t2_t3", "CVA_t2_t4", "CVA_t3_t4",
                                   "h2_t1", "h2_t2", "h2_t3", "h2_t4")

colnames(d_h2_noscale_trait_mrr_orth) <- colnames(d_h2_noscale_trait_mkr_orth)
colnames(d_h2_noscale_trait_mkr_par) <- colnames(d_h2_noscale_trait_mkr_orth)
colnames(d_h2_noscale_trait_mrr_par) <- colnames(d_h2_noscale_trait_mkr_orth)

colnames(d_h2_noscale_mrr_par) <- h2_colnames
colnames(d_h2_noscale_mkr_par) <- h2_colnames
colnames(d_h2_noscale_mrr_orth) <- h2_colnames
colnames(d_h2_noscale_mkr_orth) <- h2_colnames


# join
d_h2_noscale_mkr_orth$calcMode <- "mkr"
d_h2_noscale_mrr_orth$calcMode <- "mrr"
d_h2_noscale_mkr_par$calcMode <- "mkr"
d_h2_noscale_mrr_par$calcMode <- "mrr"

d_h2_noscale_trait_mkr_orth$calcMode <- "mkr"
d_h2_noscale_trait_mrr_orth$calcMode <- "mrr"
d_h2_noscale_trait_mkr_par$calcMode <- "mkr"
d_h2_noscale_trait_mrr_par$calcMode <- "mrr"

d_h2_noscale_mkr_orth$dataset <- "Orthogonal"
d_h2_noscale_mrr_orth$dataset <- "Orthogonal"
d_h2_noscale_mkr_par$dataset <- "Parallel"
d_h2_noscale_mrr_par$dataset <- "Parallel"
d_h2_noscale_mkr$dataset <- "Randomised"
d_h2_noscale_mrr$dataset <- "Randomised"

d_h2_noscale_trait_mkr_orth$dataset <- "Orthogonal"
d_h2_noscale_trait_mrr_orth$dataset <- "Orthogonal"
d_h2_noscale_trait_mkr_par$dataset <- "Parallel"
d_h2_noscale_trait_mrr_par$dataset <- "Parallel"
d_h2_noscale_trait_mkr$dataset <- "Randomised"
d_h2_noscale_trait_mrr$dataset <- "Randomised"

d_h2_noscale <- rbind(d_h2_noscale_mkr, d_h2_noscale_mrr,
              d_h2_noscale_mkr_orth, d_h2_noscale_mrr_orth,
              d_h2_noscale_mkr_par, d_h2_noscale_mrr_par)

d_h2_noscale_trait <- rbind(d_h2_noscale_trait_mkr, d_h2_noscale_trait_mrr, 
                    d_h2_noscale_trait_mkr_orth, d_h2_noscale_trait_mrr_orth,
                    d_h2_noscale_trait_mkr_par, d_h2_noscale_trait_mrr_par)

d_h2_noscale_trait %>% mutate(model = d_combos$model[.$modelindex],
                      model = factor(model, levels = model_names),
                      r = d_combos$r[.$modelindex]) -> d_h2_noscale_trait

d_h2_noscale %>% mutate(model = d_combos$model[.$modelindex],
                model = factor(model, levels = model_names),
                r = d_combos$r[.$modelindex]) -> d_h2_noscale

d_h2_noscale_trait <- d_h2_noscale_trait %>%
  distinct(gen, seed, modelindex, dataset, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()

d_h2_noscale <- d_h2_noscale %>%
  distinct(gen, seed, modelindex, dataset, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()

saveRDS(d_h2_noscale, "/mnt/i/SLiMTests/tests/newMotifs/paper1/d_h2_noscale.RDS")
saveRDS(d_h2_noscale_trait, "/mnt/i/SLiMTests/tests/newMotifs/paper1/d_h2_noscale_trait.RDS")

d_h2_noscale <- left_join(d_h2_noscale, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))


# Discretise generation
d_h2_noscale <- d_h2_noscale %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))



d_h2_noscale %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(6:17, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(timePoint, modelindex, dataset, isAdapted) %>%
  group_split(.) -> split_h2_noscale_mc

lapply(split_h2_noscale_mc, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_noscale_mc


h2_mat_noscale_mc <- unlist(cov_matrices_noscale_mc, recursive = F)

# get ids from the matrix
cov_matrix_modelindex_noscale_mc <- GetMatrixIDsWithDataset(split_h2_noscale_mc)

id_noscale_mc <- data.table::rbindlist(cov_matrix_modelindex_noscale_mc, 
                               fill = T)
id_noscale_mc$label <- as.character(1:nrow(id_noscale_mc))
id_noscale_mc$modelindex <- as.factor(id_noscale_mc$modelindex)
id_noscale_mc <- AddCombosToDF(id_noscale_mc)
id_noscale_mc$model <- factor(id_noscale_mc$model, levels = model_names, labels = model_names_noquote)

# Resize h2_pd according to n traits
h2_mat_noscale_mc <- purrr::map(seq_along(h2_mat_noscale_mc), function(i) {
  mat <- h2_mat_noscale_mc[[i]]
  
  n <- molComp_names[[as.character(id_noscale_mc$model[i])]]
  mat[n, n]
})


# First convert to nearest positive definite matrix
h2_pd_noscale_mc <- lapply(h2_mat_noscale_mc, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

saveRDS(h2_pd_noscale_mc, "h2_pd_noscale_mc.RDS")
saveRDS(id_noscale_mc, "id_trait_noscale_mc.RDS")


prop_var_noscale_mc <- purrr::map(seq_along(h2_pd_noscale_mc), function(i) {
  vars <- diag(h2_pd_noscale_mc[[i]])
  d_template <- data.frame(matrix(ncol = 13, nrow = 1))
  names(d_template) <- c(names(all_molcomp_features), "totalVar")
  d_template[1,names(vars)] <- vars / sum(vars)
  d_template[1,13] <- sum(vars)
  return(d_template)
})

# Attach ID
d_prop_vars_noscale <- data.table::rbindlist(prop_var_noscale_mc)
d_prop_vars_noscale <- cbind(d_prop_vars_noscale, id_noscale_mc)

# Pivot longer
d_prop_vars_noscale <- d_prop_vars_noscale %>%
  pivot_longer(cols = (1:12),
               names_to = "molComp",
               values_to = "varExpl")

d_prop_vars_noscale <- d_prop_vars_noscale %>%
  drop_na(varExpl) %>%
  filter(isAdapted == T, log10(r) == -1, timePoint == "End") %>%
  mutate(dataset = factor(dataset, levels = c("Parallel",
                                              "Orthogonal",
                                              "Randomised")))

d_prop_vars_noscale_sum <- d_prop_vars_noscale %>%
  group_by(model, dataset, molComp) %>%
  summarise(meanVarExpl = mean(varExpl),
            CIVarExpl = CI(varExpl),
            meanVar = mean(totalVar),
            CIVar = CI(totalVar))

ggplot(d_prop_vars_noscale_sum,
       aes(x = molComp, y = meanVarExpl, fill = model)) +
  facet_nested("Model" + model ~ "Trait/selection alignment" + dataset) +
  geom_col() +
  geom_text(aes(x = 6.5, y = 0.85, label = 
                  paste("Total variance =", round(meanVar, digits = 3),
                        "±", round(CIVar, digits = 3))), size = 4) +
  scale_x_discrete(labels = function(x) {parse(text = all_molcomp_features[x])}) +
  scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = meanVarExpl - CIVarExpl, ymax = meanVarExpl + CIVarExpl), width = 0.2) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                    guide = "none") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      guide = "none") +
  labs(x = "Molecular component", y = "Mean genetic variance explained (%)") +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_var_expl_noscale.png", device = png, width = 13, height = 10)

# {Print most importance features}
print(xtable::xtable(d_prop_vars_noscale_sum %>%
                       group_by(model, dataset) %>%
                       slice_max(meanVarExpl, n = 4) %>%
                       ungroup() %>%
                       select(-model),
                     
                     digits = 3), include.rownames = F)


# Now for traits
# inner join optPerc
d_h2_noscale_trait <- left_join(d_h2_noscale_trait, d_qg_optPerc, by = c("gen", "seed", "modelindex", "dataset"))

# Discretise generation
d_h2_noscale_trait <- d_h2_noscale_trait %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))



# Split h2 into G matrices
d_h2_noscale_trait %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(!if_all(5:8, is.na)) %>%  # Drop rows with no variance
  distinct(gen, seed, modelindex, dataset, .keep_all = T) %>%
  group_by(modelindex, timePoint, dataset, isAdapted) %>%
  group_split(.) -> split_h2_noscale_trait


# Separate into model indices
# each sublist is replicates of a model index

lapply(split_h2_noscale_trait, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_noscale_trait


# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat_noscale_trait <- unlist(cov_matrices_noscale_trait, recursive = F)

# get ids from the matrix
cov_matrix_modelindex_noscale_trait <- GetMatrixIDsWithDataset(split_h2_noscale_trait)

id_noscale_trait <- data.table::rbindlist(cov_matrix_modelindex_noscale_trait, 
                            fill = T)
id_noscale_trait$label <- as.character(1:nrow(id_noscale_trait))
id_noscale_trait$modelindex <- as.factor(id_noscale_trait$modelindex)
id_noscale_trait <- AddCombosToDF(id_noscale_trait)
id_noscale_trait$model <- factor(id_noscale_trait$model, levels = model_names, labels = model_names_noquote)

# Resize h2_pd according to n traits
h2_mat_noscale_trait <- purrr::map(seq_along(h2_mat_noscale_trait), function(i) {
  mat <- h2_mat_noscale_trait[[i]]
  
  n <- GetMotifTraitRange(as.character(id_noscale_trait$model[i]))
  
  mat[n, n]
})


# First convert to nearest positive definite matrix
h2_pd_noscale_trait <- lapply(h2_mat_noscale_trait, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})

saveRDS(h2_pd_noscale_trait, "/mnt/i/SLiMTests/tests/newMotifs/paper1/h2_pd_noscale_trait.RDS")
saveRDS(id_noscale_trait, "/mnt/i/SLiMTests/tests/newMotifs/paper1/id_noscale_trait.RDS")





# Vrel_GT, Vrel_Gc, Vrel_M
## Measure how isomorphically variation is distributed across traits/components
vrel_noscale_gt <- purrr::map(seq_along(h2_pd_noscale_trait), function(i) {
  g <- h2_pd_noscale_trait[[i]]
  
  vrel_i <- Vrel(eigen(g)$values)
  
  result <- id_noscale_trait[i,]
  
  result$vrel <- vrel_i
  return(result)
})

d_vrel_noscale_gt <- data.table::rbindlist(vrel_noscale_gt)

vrel_noscale_gc <- purrr::map(seq_along(h2_pd_noscale_mc), function(i) {
  g <- h2_pd_noscale_mc[[i]]
  
  vrel_i <- Vrel(eigen(g)$values)
  
  result <- id_noscale_mc[i,]
  
  result$vrel <- vrel_i
  return(result)
})
d_vrel_noscale_gc <- data.table::rbindlist(vrel_noscale_gc)

vrel_m <- purrr::map(seq_along(m_matrices_valid_pd), function(i) {
  g <- m_matrices_valid_pd[[i]]
  
  vrel_i <- Vrel(eigen(g)$values)
  
  result <- id_m_valid[i,]
  
  result$vrel <- vrel_i
  return(result)
})
d_vrel_m <- data.table::rbindlist(vrel_m)

# Combine
d_vrel_noscale_gt$matrix <- TeX("$G_T$", output = "character")
d_vrel_noscale_gc$matrix <- TeX("$G_C$", output = "character")
d_vrel_m$matrix <- TeX("$M$", output = "character")

d_vrel_noscale <- rbind(d_vrel_noscale_gt, d_vrel_noscale_gc, d_vrel_m)
d_vrel_noscale$dataset <- factor(d_vrel_noscale$dataset, levels = c("Parallel",
                                                    "Orthogonal",
                                                    "Randomised"))

d_vrel_noscale_sum <- d_vrel_noscale %>%
  filter(isAdapted == T, log10(r) == -1) %>%
  group_by(model, matrix, dataset) %>%
  summarise(meanVrel = mean(vrel))

ggplot(d_vrel_noscale %>% 
         filter(isAdapted == T, log10(r) == -1),
       aes(x = model, y = vrel, colour = model)) +
  facet_nested("Matrix" + matrix ~ "Trait/selection alignment" + dataset, 
               labeller = labeller(matrix = label_parsed)) +
  geom_quasirandom(show.legend = F) +
  geom_point(data = d_vrel_noscale_sum, aes(y = meanVrel), colour = "black", 
             fill = "white", stroke = 1,
             shape = 21, inherit.aes = T, show.legend = F) +
  scale_colour_manual(values = pal) +
  labs(x = "Model", y = TeX("$V_{rel}")) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_vrel_noscale.png", device = png, bg = "white", width = 10, height = 8)


# bTGb, bTCb, bTMb
## Measure how much variance is along the selection gradient
d_selvec <- readRDS("/mnt/i/SLiMTests/tests/newMotifs/paper1/d_selvec.RDS")
d_selvec$model <- factor(d_selvec$model, levels = model_names, labels = model_names_noquote)

d_selvec2_noscale_gt <- inner_join(id_noscale_trait, d_selvec, 
                           by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim_noscale_gt <- GetCosineSimilarity(h2_pd_noscale_trait, d_selvec2_noscale_gt %>% 
                                             select(ends_with("dir")), id_noscale_trait)

d_selvec2_noscale_mc <- inner_join(id_noscale_mc, d_selvec, 
                           by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim_noscale_gc <- GetCosineSimilarity(h2_pd_noscale_mc, d_selvec2_noscale_mc %>% 
                                             select(ends_with("dir")), id_noscale_mc)

d_selvec2_m <- inner_join(id_m_valid, d_selvec, 
                          by = c("timePoint", "seed", "modelindex", "dataset", "model", "r"))

d_cossim_m <- GetCosineSimilarity(m_matrices_valid_pd, d_selvec2_m %>% select(ends_with("dir")), id_m_valid)

# Combine
d_cossim_noscale_gt$matrix <- TeX("$G_T$", output = "character")
d_cossim_noscale_gc$matrix <- TeX("$G_C$", output = "character")
d_cossim_m$matrix <- TeX("$M$", output = "character")

d_cossim_noscale <- rbind(d_cossim_noscale_gt, d_cossim_noscale_gc, d_cossim_m)

d_cossim_noscale <- AddCombosToDF(d_cossim_noscale)

d_cossim_noscale$model <- factor(d_cossim_noscale$model, levels = model_names, labels = model_names_noquote)

d_cossim_noscale_sum <- d_cossim_noscale %>%
  filter(isAdapted == T, log10(r) == -1) %>%
  group_by(model, matrix, dataset) %>%
  summarise(meanLogbTMb = mean(log10(bTMb)))



ggplot(d_cossim_noscale %>% 
         filter(isAdapted == T, log10(r) == -1),
       aes(x = model, y = log10(bTMb), colour = model)) +
  facet_nested("Matrix" + matrix ~ "Trait/selection alignment" + dataset, 
               labeller = labeller(matrix = label_parsed)) +
  geom_quasirandom(show.legend = F) +
  geom_point(data = d_cossim_noscale_sum, aes(y = meanLogbTMb), colour = "black", 
             fill = "white", stroke = 1,
             shape = 21, inherit.aes = T, show.legend = F) +
  scale_colour_manual(values = pal) +
  labs(x = "Model", y = TeX("$log_{10}(\\beta^T X \\beta)$")) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_btxb_noscale.png", device = png, bg = "white", width = 10, height = 8)

# PCASim
## Measure how similar the covariance within/between components/traits is between replicates
krz_in_gc_noscale <- id_noscale_mc %>%
  mutate(g = h2_pd_noscale_mc,
         group = interaction(model, dataset, log10(r)))


# Save krz_in: run this part on HPC
saveRDS(krz_in_gc_noscale, "/mnt/d/SLiMTests/tests/newMotifs/paper1/pca_in_gc_noscale.RDS")

# Bootstrap in ten parts for RAM reasons
# This is slow: uncomment to run, otherwise read in precalculated data
# Generate seeds
# newseed <- sample(1:.Machine$integer.max, 10)
# [1]  314145285 2009911717  267335506  231424073 1190700402 1189454198  395819651  848071181 1762114410
# [10] 1739509036
newseed <- c(314145285L, 2009911717L, 267335506L, 231424073L, 1190700402L, 
             1189454198L, 395819651L, 848071181L, 1762114410L, 1739509036L)
bootPCASim_noscale <- vector(mode = "list", length = length(newseed))

# Per model inputs
krz_in_gc_NAR_noscale <- krz_in_gc_noscale %>% filter(model == "NAR", log10(r) == -1)
krz_in_gc_PAR_noscale <- krz_in_gc_noscale %>% filter(model == "PAR", log10(r) == -1)
krz_in_gc_FFLC1_noscale <- krz_in_gc_noscale %>% filter(model == "FFLC1", log10(r) == -1)
krz_in_gc_FFLI1_noscale <- krz_in_gc_noscale %>% filter(model == "FFLI1", log10(r) == -1)
krz_in_gc_FFBH_noscale <- krz_in_gc_noscale %>% filter(model == "FFBH", log10(r) == -1)


for (i in seq_along(newseed)) {
  # Set seed
  set.seed(newseed[i])
  # Run replicate but only within models
  res_NAR <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_NAR_noscale, "group", T))
  res_PAR <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_PAR_noscale, "group", T))
  res_FFLC1 <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_FFLC1_noscale, "group", T))
  res_FFLI1 <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_FFLI1_noscale, "group", T))
  res_FFBH <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_gc_FFBH_noscale, "group", T))
  
  # To data.frame
  res_NAR <- unnest(as.data.frame(t(res_NAR)), cols = everything()) %>%
    mutate(model = "NAR")
  res_PAR <- unnest(as.data.frame(t(res_PAR)), cols = everything()) %>%
    mutate(model = "PAR")
  res_FFLC1 <- unnest(as.data.frame(t(res_FFLC1)), cols = everything()) %>%
    mutate(model = "FFLC1")
  res_FFLI1 <- unnest(as.data.frame(t(res_FFLI1)), cols = everything()) %>%
    mutate(model = "FFLI1")
  res_FFBH <- unnest(as.data.frame(t(res_FFBH)), cols = everything()) %>%
    mutate(model = "FFBH")
  
  # Combine to output
  bootPCASim_noscale[[i]] <- rbind(res_NAR, res_PAR, res_FFLC1, res_FFLI1, res_FFBH)
}

# Output list into combined df
bootPCASim2 <- bind_rows(bootPCASim_noscale)
bootPCASim_noscale <- bootPCASim2 %>%
  separate(group1, c("model1", "dataset1", "r1"), "\\.",
           extra = "merge") %>%
  separate(group2, c("model2", "dataset2", "r2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = as.numeric(r1),
         r2 = as.numeric(r2),
         dataset1 = factor(dataset1, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised")),
         dataset2 = factor(dataset2, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised"))) %>%
  rename(PCASim = krzCor)

saveRDS(bootPCASim_noscale, paste0("/mnt/d/SLiMTests/tests/newMotifs/paper1/d_bootPCASim_noscale.RDS"))
bootPCASim_noscale <- readRDS(paste0("/mnt/d/SLiMTests/tests/newMotifs/paper1/d_bootPCASim_noscale.RDS"))

# Plot

# Get model comparisons for labelling
bootPCASim_noscale <- bootPCASim_noscale %>%
  mutate(datasetCombo = GetModelComparison(dataset1, dataset2, c("Parallel",
                                                                 "Orthogonal",
                                                                 "Randomised")),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)),
         model = factor(model, levels = model_names_noquote),
         nMCs = unlist(lapply(molComp_names[as.character(model)], length)))

# recomb by modelCombo - we don't have all the recombination levels for the
# non-randomised datasets
bootPCASim_noscale_sum <- bootPCASim_noscale %>%
  group_by(model, dataset1, dataset2) %>%
  summarise(meanPCASim = mean(PCASim),
            ciPCASim = CI(PCASim),
            nMCs = length(molComp_names[as.character(model[1])][[1]]))

# Facet design
design <- c(
  "
  AABB
  CCDD
  #EE#
  "
)


ggplot(bootPCASim_noscale_sum, aes(
  x = (dataset1), y = (dataset2)
)) +
  #facet_nested("Model" + model ~ .) + 
  facet_manual(.~model, design = design, axes = T) + 
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim_noscale, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = seq(0.0, 1.0, by = 0.25),
                       labels = seq(0.0, 1.0, by = 0.25),
                       limits = c(0.0, 1)) +
  # scale_fill_gradientn(colours = paletteer_c("ggthemes::Blue-Green Sequential", n = 6),
  #                      values = c(0.0, 0.7, 0.75, 0.8, 0.9, 1.0)) +
  labs(x = "Trait/selection alignment (Matrix 1)", y = "Trait selection alignment (Matrix 2)", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 12), 
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "bottom",
        legend.key.width = unit(3.5, "lines"))
ggsave("PCASim_dataset_noscale.png", device = png, width = 8, height = 10)


