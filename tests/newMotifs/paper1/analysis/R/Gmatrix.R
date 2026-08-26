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
id$model <- factor(id$model, levels = model_names)

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
saveRDS(id, "id_trait")

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


DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar.csv"

d_m <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                         paste0("mean_", 1:4),
                                         paste0("var_", 1:4),
                                         paste0("cov_", c(12, 13, 14, 23, 24, 34))))

d_m <- d_m %>%
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex))




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

m_matrices_valid <- d_m_mean %>%
  rowwise() %>%
  group_map(~ row_to_m(.x))


id_m_valid <- d_m_mean %>%
  select(1:5)


# First convert to nearest positive definite matrix
m_matrices_valid_pd <- lapply(m_matrices_valid, function(x) {
  if (!matrixcalc::is.positive.definite(x)) {return (as.matrix(Matrix::nearPD(x)$mat))}
  return(x)
})


eig_mean_m_pd <- purrr::map(seq_along(m_matrices_valid_pd), function(i){
  eig <- eigen(m_matrices_valid_pd[[i]])
  
  circle <- cbind(cos(seq(0, 2*pi, length.out = 100)),
                  sin(seq(0, 2*pi, length.out = 100)))
  
  ellipse_pc <- circle %*% diag(sqrt(eig$values[1:2])) %*% t(eig$vectors[,1:2])
  
  # Project onto ellipse for the major/minor axes
  V12 <- eig$vectors[1:2,1:2]
  
  Sigma_proj <- V12 %*%
    diag(eig$values[1:2]) %*%
    t(V12)
  
  eig_proj <- eigen(Sigma_proj)
  
  result <- id_m_valid[rep(i, times = 100),]
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

# bTGb, bTCb, bTMb
## Measure how much variance is along the selection gradient

# PCASim
## Measure how similar the covariance within/between components/traits is between replicates

# Total variance explained
## Measure how similar the total 

