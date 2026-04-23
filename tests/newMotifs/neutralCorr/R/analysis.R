library(tidyverse)
library(ggh4x)
library(paletteer)
library(emmeans)
library(patchwork)
library(mvtnorm)
library(ggforce) # geom_ellipse
library(Rcpp)


# Helper functions
ModelFromIndex <- function(id) {
  motifs <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
  return(motifs[id])
}

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}


CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi / 180)}

cppFunction('IntegerVector double2int(NumericVector x) { 
              int n = x.size();
              IntegerVector result(n);
              
              for (int i = 0; i < n; ++i) {
                int64_t x64 = (int64_t)x[i];
                result[i] = static_cast<int32_t>(x64);
              }
              return result;
            }')


pivot_pairwise <- function(data, 
                           pivot_cols, 
                           other_cols = !pivot_cols,
                           names_to = "name",
                           values_to = "value",
                           pair_label = c("x", "y"),
                           pair_label_sep = "_",
                           row_id = "row_id") {
  
  # construct variable names
  x_value <- paste(pair_label[1], values_to, sep = pair_label_sep)
  y_value <- paste(pair_label[2], values_to, sep = pair_label_sep)
  x_name  <- paste(pair_label[1], names_to, sep = pair_label_sep)
  y_name  <- paste(pair_label[2], names_to, sep = pair_label_sep)
  
  # create an id column
  base <- data |> 
    dplyr::mutate({{row_id}} := dplyr::row_number())
  
  # variables to be retained but not pairwise-pivoted
  fixed_data <- base |> 
    dplyr::select(
      {{other_cols}}, 
      tidyselect::all_of(row_id)
    )
  
  # select pivoting columns, pivot to long, and relabel as x-var 
  long_x <- base |>
    dplyr::select(
      {{pivot_cols}},            
      tidyselect::all_of(row_id)
    ) |>
    tidyr::pivot_longer(
      cols = {{pivot_cols}},
      names_to = {{x_name}},
      values_to = {{x_value}}
    )
  
  # same data frame, but with new variable names for pivoted vars
  long_y <- long_x |> 
    dplyr::rename(
      {{y_name}} := {{x_name}}, 
      {{y_value}} := {{x_value}}
    )
  
  # full join with many-to-many gives all pairs; then restore other columns
  pairs <- dplyr::full_join(
    x = long_x,
    y = long_y,
    by = row_id,
    relationship = "many-to-many"
  ) |>
    dplyr::relocate({{y_name}}, .after = {{x_name}}) |>
    dplyr::left_join(fixed_data, by = row_id)
  
  return(pairs)
}

generate_vector_at_angle <- function(u, angle) {
  # https://stackoverflow.com/questions/72374034/generating-two-vectors-with-a-given-angle-between-them
  ab <- runif(2)
  a <- ab[1]
  b <- ab[2]
  h <- (b - a) - (b - a) * (u*u)
  
  cos(angle) * u + sin(angle) * h
}

test_u <- c(0.1, 0.2)
d_test_gen <- as.data.frame(t(replicate(100, generate_vector_at_angle(test_u, deg2rad(90)))))


ggplot(d_test_gen, 
       aes(x = V1, y = V2)) +
  geom_abline(slope = 0.2/0.1, colour = "red", linetype = "dashed") +
  geom_point() +
  theme_bw()

# Load in data
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

# Check if different from zero
lm.traitcor.nar <- lm(cor ~ traitCombo, d_qg_traitCor %>% filter(model == "NAR"))
lm.traitcor.par <- lm(cor ~ traitCombo, d_qg_traitCor %>% filter(model == "PAR"))
lm.traitcor.fflc1 <- lm(cor ~ traitCombo, d_qg_traitCor %>% filter(model == "FFLC1"))
lm.traitcor.ffli1 <- lm(cor ~ traitCombo, d_qg_traitCor %>% filter(model == "FFLI1"))
lm.traitcor.ffbh <- lm(cor ~ traitCombo, d_qg_traitCor %>% filter(model == "FFBH"))
summary(lm.traitcor.nar)
summary(lm.traitcor.par)
summary(lm.traitcor.fflc1)
summary(lm.traitcor.ffli1)
summary(lm.traitcor.ffbh)

em.nar <- emmeans(lm.traitcor.nar, ~ traitCombo)
pairs(em.nar, simple = "traitCombo")
plot(em.nar, comparisons = T)
joint_tests(em.nar)

em.par <- emmeans(lm.traitcor.par, ~ traitCombo)
pairs(em.par, simple = "traitCombo")
plot(em.par, comparisons = T)
joint_tests(em.par)

em.fflc1 <- emmeans(lm.traitcor.fflc1, ~ traitCombo)
pairs(em.fflc1, simple = "traitCombo")
plot(em.fflc1, comparisons = T)
joint_tests(em.fflc1)

em.ffli1 <- emmeans(lm.traitcor.ffli1, ~ traitCombo)
pairs(em.ffli1, simple = "traitCombo")
plot(em.ffli1, comparisons = T)
joint_tests(em.ffli1)

em.ffbh <- emmeans(lm.traitcor.ffbh, ~ traitCombo)
pairs(em.ffbh, simple = "traitCombo")
plot(em.ffbh, comparisons = T)
joint_tests(em.ffbh)

 
d_qg_traitCor_sum <- d_qg_traitCor %>%
group_by(model, traitCombo) %>%
  mutate(angle = rad2deg(acos(cor))) %>%
    summarise(mean = mean(cor),
              var = var(cor),
              angle = mean(angle)) %>%
    ungroup()


trait_comp_names_nar <- c(
  "Response time vs steady state" # 12
)

trait_comp_names_c1 <- c(
  "Response time vs response delay", # 12
  "Response time vs steady state", # 13
  "Response delay vs steady state" # 23
)

trait_comp_names_i1 <- c(
  "Time to half max expression vs max expression", # 12
  "Time to half max expression vs time above half max expression", # 13
  "Max expression vs time above half max expression" # 23
)

trait_comp_names_bh <- c(
  "Time to max expression vs max expression", # 12
  "Time to max expression vs time to final steady state", # 13
  "Time to max expression vs final steady state", # 14
  "Max expression vs time to final steady state", # 23
  "Max expression vs final steady state", # 24
  "Time to final steady state vs final steady state" # 34
)

trait_comp_names <- data.frame(model = c(rep("NAR", times = length(trait_comp_names_nar)),
                                         rep("PAR", times = length(trait_comp_names_nar)),
                                         rep("FFLC1", times = length(trait_comp_names_c1)),
                                         rep("FFLI1", times = length(trait_comp_names_i1)),
                                         rep("FFBH", times = length(trait_comp_names_bh))),
                               label = c(trait_comp_names_nar, trait_comp_names_nar,
                                         trait_comp_names_c1, trait_comp_names_i1,
                                         trait_comp_names_bh))

traitCorBoilerplate <- function(plt) {
  plt + 
    facet_nested(. ~ model) +
    geom_point(position = position_dodge(0.9), show.legend = F) +
    geom_errorbar(aes(ymin = mean - var, ymax = mean + var), width = 0.2,
                  position = position_dodge(0.9), show.legend = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    labs(x = "Trait combination", y = "Mean correlation", colour = "Model") +
    theme_bw() +
    theme(text = element_text(size = 12),
          legend.position = "bottom")
}


ggplot(d_qg_traitCor_sum %>%
         filter(model == "NAR") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
    aes(x = as.factor(traitCombo), y = mean, colour = model)) +
    scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "NAR",2],
                     guide = guide_axis(n.dodge = 2)) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                             5, direction = -1)[4]) -> plt_traitcor_nar
plt_traitcor_nar <- traitCorBoilerplate(plt_traitcor_nar)

ggplot(d_qg_traitCor_sum %>%
         filter(model == "PAR") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "PAR",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[5]) -> plt_traitcor_par
plt_traitcor_par <- traitCorBoilerplate(plt_traitcor_par)

ggplot(d_qg_traitCor_sum %>%
         filter(model == "FFLC1") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "FFLC1",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[2]) -> plt_traitcor_fflc1
plt_traitcor_fflc1 <- traitCorBoilerplate(plt_traitcor_fflc1)

ggplot(d_qg_traitCor_sum %>%
         filter(model == "FFLI1") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "FFLI1",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[3]) -> plt_traitcor_ffli1
plt_traitcor_ffli1 <- traitCorBoilerplate(plt_traitcor_ffli1)


ggplot(d_qg_traitCor_sum %>%
         filter(model == "FFBH") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "FFBH",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[1]) -> plt_traitcor_ffbh
plt_traitcor_ffbh <- traitCorBoilerplate(plt_traitcor_ffbh)


layout <-
"
AABB
CCCC
DDDD
EEEE
"

plt_traitcor <- plt_traitcor_nar +
plt_traitcor_par +
plt_traitcor_fflc1 +
plt_traitcor_ffli1 + 
plt_traitcor_ffbh

plt_traitcor <- plt_traitcor + plot_layout(design = layout)
plt_traitcor

ggsave("plt_trait_cor.png", device = png, bg = "white", width = 12, height = 8)

# Print correct correlations
d_qg_traitCor_sum %>%
  filter(!(model == "NAR" & grepl("[3-4]", traitCombo)),
         !(model == "PAR" & grepl("[3-4]", traitCombo)),
         !(model == "FFLC1" & grepl("[4]", traitCombo)),
         !(model == "FFLI1" & grepl("[4]", traitCombo))) %>%
  select(-var) -> d_traitCor_tab
d_traitCor_tab

# Get correlation matrices
C_NAR <- matrix(rep(1, times = 4), nrow = 2)
C_NAR[!diag(nrow = 2, ncol = 2)] <- d_traitCor_tab[d_traitCor_tab$model == "NAR",]$mean

C_PAR <- matrix(rep(1, times = 4), nrow = 2)
C_PAR[!diag(nrow = 2, ncol = 2)] <- d_traitCor_tab[d_traitCor_tab$model == "PAR",]$mean

C_FFLC1 <- matrix(rep(1, times = 9), nrow = 3)
C_FFLC1[upper.tri(C_FFLC1)] <- d_traitCor_tab[d_traitCor_tab$model == "FFLC1",]$mean
C_FFLC1[lower.tri(C_FFLC1)] <- d_traitCor_tab[d_traitCor_tab$model == "FFLC1",]$mean

C_FFLI1 <- matrix(rep(1, times = 9), nrow = 3)
C_FFLI1[upper.tri(C_FFLI1)] <- d_traitCor_tab[d_traitCor_tab$model == "FFLI1",]$mean
C_FFLI1[lower.tri(C_FFLI1)] <- d_traitCor_tab[d_traitCor_tab$model == "FFLI1",]$mean

C_FFBH <- matrix(rep(1, times = 16), nrow = 4)
C_FFBH[lower.tri(C_FFBH)] <- d_traitCor_tab[d_traitCor_tab$model == "FFBH",]$mean
C_FFBH[upper.tri(C_FFBH)] <- t(C_FFBH)[upper.tri(C_FFBH)]

C_NAR; C_PAR; C_FFLC1; C_FFLI1; C_FFBH

# Save direction vector to file
# SLiM can't do eigen(), will need to load them in
# scale eigenvectors to unit length so we can scale by selection strength in SLiM

# The parallel choice is the first eigenvector - need a value per model
cor_matrices <- list(C_NAR,
                     C_PAR,
                     C_FFLC1,
                     C_FFLI1,
                     C_FFBH)

d_parallel_motifs <- data.frame(t1 = double(length(cor_matrices)),
                                t2 = double(length(cor_matrices)),
                                t3 = double(length(cor_matrices)),
                                t4 = double(length(cor_matrices)),
                                l1 = double(length(cor_matrices)))

d_orth_motifs <- d_parallel_motifs

for (i in seq_along(cor_matrices)) {
  m <- cor_matrices[[i]]
  
  eig <- eigen(m)
  # Parallel direction
  v1 <- eig$vectors[,1]
  v1 <- v1 / sum(abs(v1)) # normalise
  
  # Minimal variance/most misaligned with correlations
  vlast <- eig$vectors[,ncol(eig$vectors)]
  vlast <- vlast / sum(abs(vlast)) # normalise
  
  if (length(v1) < 4) {
    v1 <- c(v1, double(4 - length(v1)))
    vlast <- c(vlast, double(4 - length(vlast)))
  }
  
  l1 <- eig$values[1]
  llast <- eig$values[length(eig$values)]
  
  d_parallel_motifs[i,] <- c(v1, l1)
  d_orth_motifs[i,] <- c(vlast, llast)
}

write_csv(d_parallel_motifs, "./parallel_traitdir.csv", col_names = F)
write_csv(d_orth_motifs, "./orth_traitdir.csv", col_names = F)

# Test direction is right
d_dir_nar <- data.frame(point = c("start", "end"),
                    dir = rep(c("parallel", "orthogonal"), each = 2),
                    x = c(0.052335, 0.0575675, 0.052335, 0.0471026),
                    y = c(2.58808, 2.59331, 2.58808, 2.59331))

ggplot(d_dir_nar,
       aes(x = x, y = y, colour = point)) +
  facet_grid(dir~.) +
  geom_point() +
  geom_line(aes(group = NA), 
            arrow = arrow(length = unit(0.03, "npc")), colour = "black") + 
  theme_bw()

d_dir_fflc1 <- data.frame(point = c("start", "end"),
                        dir = rep(c("parallel", "orthogonal"), each = 2),
                        x = c(0.15665, 0.17048, 0.15665, 0.167308),
                        y = c(0.2, 0.20938, 0.2, 0.185343),
                        z = c(7.41879, 7.4363, 7.41879, 7.41822))

ggplot(d_dir_fflc1,
       aes(x = y, y = z, colour = x, shape = point)) +
  facet_grid(dir~.) +
  geom_point() +
  geom_line(aes(group = NA), 
            arrow = arrow(length = unit(0.03, "npc"),
                          ends = "last"), colour = "black") + 
  theme_bw()


test_FFLC1 <- rmvnorm(1000, c(0, 0, 0), sigma = C_FFLC1)
plot(test_FFLC1[,1], test_FFLC1[,2])
plot(test_FFLC1[,1], test_FFLC1[,3])
plot(test_FFLC1[,2], test_FFLC1[,3])
plot(test_FFLC1 %*% (eig_fflc1$vectors))
plot(test_FFLC1 %*% (eig_fflc1$vectors / sum(abs(eig_fflc1$vectors))))


# trait correlations between ~0 and ~0.333
# Save in a table for the adjusted selection direction experiments
write_csv(d_traitCor_tab, "./traitCor.csv", col_names = F)

d_qg_means <- d_qg %>%
  mutate(model = ModelFromIndex(modelindex)) %>%
  group_by(gen, seed, model) %>%
  pivot_longer(cols = starts_with("phenomean"),
               names_to = "trait",
               values_to = "mean", names_prefix = "phenomean_") %>%
  select(gen, seed, model, trait, mean)

d_qg_vars <- d_qg %>%
  mutate(model = ModelFromIndex(modelindex)) %>%
  group_by(gen, seed, model) %>%
  pivot_longer(cols = starts_with("phenovar"),
               names_to = "trait",
               values_to = "var", names_prefix = "phenovar_") %>%
  select(gen, seed, model, trait, var)

d_qg_sum <- d_qg_means
d_qg_sum$var <- d_qg_vars$var

d_qg_sum <- d_qg_sum %>%
  group_by(model) %>%
  # Standardise
  mutate(mean = scale(mean),
         var = scale(var)) %>%
  group_by(gen, model, trait) %>%
  summarise(meanMean = mean(mean),
            meanVar = mean(var),
            CIMean = CI(mean),
            CIVar = CI(var))

ggplot(d_qg_sum, 
       aes(x = gen, y = meanMean, colour = as.factor(trait))) +
  facet_nested("Model" + model ~., scales = "free") +
  geom_line() +
  geom_ribbon(aes(ymin = meanMean - CIMean, ymax = meanMean + CIMean, fill = trait),
              colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generation", y = "Mean of means", colour = "Trait") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, direction = -1)) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, direction = -1)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

ggplot(d_qg_sum, 
       aes(x = gen, y = meanVar, colour = as.factor(trait))) +
  facet_nested("Model" + model ~., scales = "free") +
  geom_line() +
  geom_ribbon(aes(ymin = meanVar - CIVar, ymax = meanVar + CIVar, fill = trait),
              colour = NA, alpha = 0.2, show.legend = F) +
  labs(x = "Generation", y = "Mean of variances", colour = "Trait") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, direction = -1)) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, direction = -1)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# Plot mean over generations 
d_qg_sum <- d_qg_means
d_qg_sum$var <- d_qg_vars$var

d_qg_sum <- d_qg_sum %>%
  group_by(model) %>%
  # Standardise
  mutate(mean = scale(mean),
         var = scale(var)) %>%
  group_by(model, trait) %>%
  summarise(meanMean = mean(mean),
            meanVar = mean(var),
            CIMean = CI(mean),
            CIVar = CI(var))

ggplot(d_qg_sum %>%
         filter(!(model != "FFBH" & as.numeric(trait) > 3), # invalid traits
                !(nchar(model) == 3 & as.numeric(trait) > 2)), 
       aes(x = as.factor(trait), y = meanMean)) +
  facet_nested("Model" + model ~., scales = "free") +
  geom_point() +
  geom_errorbar(aes(ymin = meanMean - CIMean, ymax = meanMean + CIMean),
                width = 0.2) +
  labs(x = "Trait", y = "Mean of mean trait values") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
