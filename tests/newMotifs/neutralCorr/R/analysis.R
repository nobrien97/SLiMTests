library(tidyverse)
library(ggh4x)
library(paletteer)
library(emmeans)
library(patchwork)
library(mvtnorm)
library(ggforce) # geom_ellipse


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

FactoMineR::PCA(d_fflc1_test, scale. = T)
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
