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

# Calculate leading eigenvector for each
GetDirAtAngle <- function(c, angle, r, x0) {
  eig <- eigen(c)
  Q <- eig$vectors
  L <- diag(eig$values)
  
  # Whitened space
  W <- solve(sqrt(L)) %*% t(Q)
  Winv <- Q %*% sqrt(L)
  
  # Transform Q by whitened space and normalise
  vz <- W %*% Q[,1]
  vz <- vz / sqrt(sum(vz*vz))
  
  u <- rnorm(ncol(c))
  u <- u - sum(u * vz) * vz
  u <- u / sqrt(sum(u*u))
  
  d <- cos(angle) * vz + sin(angle) * u
  
  # Starting point in whitened space
  z0 <- W %*% x0
  
  # End point in whitened space
  z_new <- z0 + r * d
  
  # End point in original space
  x_new <- Winv %*% z_new
  return(x_new)
}

GetEllipse <- function(mat, trait_data) {
  ev <- eigen(mat, symmetric = T)
  
  # We need to find the angle to rotate the ellipse back onto trait space: 
  # this is using inverse tangent between the distance between lambda_gmax 
  # (ev$values[1]) and the first trait's variance (diag(g)[1]) and the covariance 
  # (g[upper.tri(g)][1]) between traits
  theta <- atan2((ev$values[1] - diag(mat)[1]), mat[upper.tri(mat)][1]) * (180/pi)
  
  # Get the length of the major and minor axes: these are the square root of eigenvalues
  # multiplied by the confidence interval for the ellipse (here it's a 95% CI)
  major_len <- qnorm(0.975) * sqrt(ev$values[1])
  minor_len <- qnorm(0.975) * sqrt(ev$values[2])
  
  # Center of ellipse is the sum of the traits' eigenvectors times their mean
  dplot_ellipse <- data.frame(vert_x = cos(theta * pi/180) * major_len,
                              vert_y = sin(theta * pi/180) * major_len,
                              covert_x = cos((theta - 90) * pi/180) * minor_len,
                              covert_y = sin((theta - 90) * pi/180) * minor_len,
                              mean_t1 = mean(trait_data[,1]),
                              mean_t2 = mean(trait_data[,2]),
                              theta = theta,
                              major_len = major_len,
                              minor_len = minor_len)
  return(dplot_ellipse)
}

Q_NAR <- GetDirAtAngle(C_NAR, 0, 1, rep(0, 2))
Q_PAR <- GetDirAtAngle(C_PAR, 0, 1, rep(0, 2))
Q_FFLC1 <- GetDirAtAngle(C_FFLC1, 0, 1, rep(0, 3))
Q_FFLI1 <- GetDirAtAngle(C_FFLI1, 0, 1, rep(0, 3))
Q_FFBH <- GetDirAtAngle(C_FFBH, 0, 1, rep(0, 4))

# Test
d_nar_test <- as.data.frame(rmvnorm(10000, sigma = C_NAR))
d_nar_ellipse <- GetEllipse(C_NAR, d_nar_test)

# Plot the ellipses
# geom_ellipse() expects angle in radians
ggplot(d_nar_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(t(Q_NAR)), aes(x = V1, y = V2),
             colour = "red") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 

d_par_test <- as.data.frame(rmvnorm(10000, sigma = C_PAR))
d_par_ellipse <- GetEllipse(C_PAR, d_par_test)

# Plot the ellipses
# geom_ellipse() expects angle in radians
ggplot(d_par_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(t(Q_PAR)), aes(x = V1, y = V2),
             colour = "red") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 


d_fflc1_test <- as.data.frame(rmvnorm(10000, sigma = C_FFLC1))
d_fflc1_ellipse <- GetEllipse(C_FFLC1, d_fflc1_test)

# Plot the ellipses
# geom_ellipse() expects angle in radians
ggplot(d_fflc1_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(t(Q_FFLC1)), aes(x = V1, y = V2),
             colour = "red") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 

fflc1_pca <- FactoMineR::PCA(d_fflc1_test, scale. = T)
fviz_pca_ind(fflc1_pca)

eigen(C_FFLC1)


d_ffli1_test <- as.data.frame(rmvnorm(10000, sigma = C_FFLI1))
d_ffli1_ellipse <- GetEllipse(C_FFLI1, d_ffli1_test)

# Plot the ellipses
# geom_ellipse() expects angle in radians
ggplot(d_ffli1_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(t(Q_FFLI1)), aes(x = V1, y = V2),
             colour = "red") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 


d_ffbh_test <- as.data.frame(rmvnorm(10000, sigma = C_FFBH))
d_ffbh_ellipse <- GetEllipse(C_FFBH, d_ffbh_test)

# Plot the ellipses
# geom_ellipse() expects angle in radians
ggplot(d_ffbh_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(t(Q_FFBH)), aes(x = V1, y = V2),
             colour = "red") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 

# Interpolate along C_max
eigen(C_FFLC1)$vectors[,1]



lin.interp <- function(p0, p1, s) {
  # Calculates a point s units away from the origin of a point across a direction vector 
  d <- p1 - p0
  norm_d <- sqrt(sum(d^2))
  d_unit <- d / norm_d
  
  outer(s, d_unit) + matrix(p0, nrow = length(s), ncol = length(p0), byrow = T)
}

interpol <- lin.interp(c(0,0,0),
                       eigen(C_FFLC1)$vectors[,1],
                       3)

known <- data.frame(
  x = c(0, eigen(C_FFLC1)$vectors[1,1], interpol[1,1]),
  y = c(0, eigen(C_FFLC1)$vectors[2,1], interpol[1,2]),
  z = c(0, eigen(C_FFLC1)$vectors[3,1], interpol[1,3])
)

GGally::ggpairs(known)

ggplot(d_fflc1_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(interpol), aes(x = V1, y = V2),
             colour = "red") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 

interpol <- lin.interp(c(0,0),
                       eigen(C_PAR)$vectors[,1],
                       3)


# Gets a random normalised direction orthogonal to another direction 
RandomDirectionOrthogonalTo <- function(dir) {
  u <- rnorm(length(dir))
  u <- u - sum(u * dir) * dir
  return (u / sqrt(sum(u^2)))
}

orth_point <- lin.interp(c(0,0),
                         RandomDirectionOrthogonalTo(eigen(C_PAR)$vectors[,1]),
                        3)

RandomValueOrthogonalTo(eigen(C_PAR)$vectors[,1])

known <- data.frame(
  x = c(0, eigen(C_PAR)$vectors[1,1], interpol[1,1], orth_point[1,1]),
  y = c(0, eigen(C_PAR)$vectors[2,1], interpol[1,2], orth_point[1,2])
)

GGally::ggpairs(known)

ggplot(d_par_ellipse, aes(x = mean_t1, y = mean_t2)) +
  #geom_point(data = as.data.frame(d_nar_test), aes(x = V1, y = V2)) +
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_point(data = as.data.frame(interpol), aes(x = V1, y = V2),
             colour = "red") +
  geom_point(data = as.data.frame(orth_point), aes(x = V1, y = V2),
             colour = "blue") +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 


interpol <- lin.interp(c(0,0,0),
                       eigen(C_FFLC1)$vectors[,1],
                       3)

orth_point <- lin.interp(c(0,0,0),
                         RandomDirectionOrthogonalTo(eigen(C_FFLC1)$vectors[,1]),
                         3)

known <- data.frame(
  x = c(0, eigen(C_FFLC1)$vectors[1,1], interpol[1,1], orth_point[1,1]),
  y = c(0, eigen(C_FFLC1)$vectors[2,1], interpol[1,2], orth_point[1,2]),
  z = c(0, eigen(C_FFLC1)$vectors[3,1], interpol[1,3], orth_point[1,3])
)

GGally::ggpairs(known)


# Save direction vector to file
# SLiM can't do eigen(), will need to load them in
# scale eigenvectors to unit length so we can scale by selection strength in SLiM

# The parallel choice is the first eigenvector - need a value per model

c1_NAR_u <- c(eigen(C_NAR)$vectors[,1], rep(0, 2))
c1_PAR_u <- c(eigen(C_PAR)$vectors[,1], rep(0, 2))
c1_FFLC1_u <- c(eigen(C_FFLC1)$vectors[,1], rep(0, 1))
c1_FFLI1_u <- c(eigen(C_FFLI1)$vectors[,1], rep(0, 1))
c1_FFBH_u <- eigen(C_FFBH)$vectors[,1]

d_parallel_motifs <- data.frame(t1 = double(length(cor_matrices)),
                                t2 = double(length(cor_matrices)),
                                t3 = double(length(cor_matrices)),
                                t4 = double(length(cor_matrices)),
                                l1 = double(length(cor_matrices)))

cor_matrices <- list(C_NAR,
                     C_PAR,
                     C_FFLC1,
                     C_FFLI1,
                     C_FFBH)

for (i in seq_along(cor_matrices)) {
  m <- cor_matrices[[i]]
  
  eig <- eigen(m)
  v <- eig$vectors[,1]
  v <- v / sum(abs(v)) # normalise
  
  if (length(v) < 4) {
    v <- c(v, double(4 - length(v)))
  }
  
  l <- eig$values[1]
  
  d_parallel_motifs[i,] <- c(v, l)
}


write_csv(d_parallel_motifs, "./parallel_traitdir.csv", col_names = F)

# Test direction is right

d_dir <- data.frame(point = c("start", "end"),
                    x = c(0.00473648, 0.00521012),
                    y = c(1.69566, 1.69614))


ggplot(d_dir,
       aes(x = x, y = y, colour = point)) +
  geom_point() +
  stat_ellipse(data = d_sim, aes(x = V1, y = V2), inherit.aes = F) +
 # geom_line(aes(group = NA), arrow = arrow(), colour = "black") + 
  theme_bw()


c2_NAR_u <- c(eigen(C_NAR)$vectors[,2], rep(0, 2))
c2_PAR_u <- c(eigen(C_PAR)$vectors[,2] , rep(0, 2))
c3_FFLC1_u <- c(eigen(C_FFLC1)$vectors[,3], rep(0, 1))
c3_FFLI1_u <- c(eigen(C_FFLI1)$vectors[,3], rep(0, 1))
c4_FFBH_u <- eigen(C_FFBH)$vectors[,4]

d_orth_motifs <- as.data.frame(matrix(c(c2_NAR_u, 
                                        c2_PAR_u, 
                                        c3_FFLC1_u, 
                                        c3_FFLI1_u, 
                                        c4_FFBH_u), nrow  = 5, byrow = T))

write_csv(d_orth_motifs, "./orth_traitdir.csv", col_names = F)


# The orthogonal is a randomly sampled other eigenvector - because randomly sampled,
# we need a value for each model and seed
# Load model seeds
seeds <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/R/newMotifs_seeds.csv", 
                  col_names = F)

# These are 64 bits but R only supports 32 so we will ignore the last 32 bits
seeds_t <- double2int(seeds$X1)

# Result matrix - number of seeds x max number of traits x number motifs 
result <- matrix(double(length(seeds_t) * ncol(C_FFBH) * 5), ncol = ncol(C_FFBH))

eig_motifs <- list(v_NAR_u, 
                v_PAR_u, 
                v_FFLC1_u, 
                v_FFLI1_u, 
                v_FFBH_u)

# Fill result matrix 
for (i in seq_along(seeds_t)) {
  seed <- seeds_t[i]
  for (j in seq_along(eig_motifs)) {
    motif <- eig_motifs[[j]]
    
    # Number of traits
    n <- length(motif)
    
    # 1D index
    index <- (i-1) * length(eig_motifs) + j
    
    # Set seed
    set.seed(seed)
    
    # Sample and normalise direction
    dir <- RandomDirectionOrthogonalTo(motif)
    #dir <- dir / (sum(abs(dir)))
    result[index, 1:n] <- dir
  }
}

# Test that they are orthogonal in some directions
known <- data.frame(
  t1 = c(0, eig_motifs[[3]][1], result[8,1]),
  t2 = c(0, eig_motifs[[3]][2], result[8,2]),
  t3 = c(0, eig_motifs[[3]][3], result[8,3])
)

sum(eig_fflc1$vectors[1,] * eig_fflc1$vectors[2,])

known <- known %>%
  pivot_pairwise(pivot_cols = starts_with("t"))

ggplot(known, 
       aes(x = x_value, y = y_value)) +
  ggh4x::facet_grid2(y_name ~ x_name, scales = "free",
                     render_empty = F) +
  geom_point() +
  geom_line() +
  labs(x = NULL, y = NULL) +
  theme_bw()

GGally::ggpairs(known)

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
