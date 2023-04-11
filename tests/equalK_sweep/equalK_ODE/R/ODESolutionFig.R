library(tidyverse)
library(gganimate)
library(gghighlight)
library(FactoMineR)
library(factoextra)
library(latex2exp)
library(lattice)
library(latticeExtra)
library(paletteer)
library(foreach)
library(doParallel)
library(future)

# plot ODE solution

d_ode <- read_csv("/mnt/d/SLiMTests/tests/equalK_sweep/equalK_ODE/data/ode_solutions.csv", 
                  col_names = F)

colnames(d_ode) <- c("gen", "seed", "modelindex", "nloci", "sigma", "phenomean", "time", "X", "Z",
                     "diff_bottom", "diff_top", "len_left", "len_right", 
                     "bottom_left_angle", "bottom_right_angle", "top_right_angle", "top_left_angle")


d_ode$nloci_cat <- cut(d_ode$nloci,
                       breaks=c(-Inf, 10, 50, 250, Inf))

d_ode$sigma_cat <- cut(d_ode$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))



d_ode <- d_ode %>% 
  mutate(id = paste(seed, modelindex, sep = "_"))

# measure lengths/angles of sides
seed <- sample(1:.Machine$integer.max, 1)
set.seed(seed)
# 18799215
set.seed(18799215)

sampled_seeds <- d_ode %>%
  group_by(nloci, sigma) %>% 
  select(nloci, sigma, seed, modelindex, id) %>%
  sample_n(1)

d_sample <- d_ode %>% filter(seed %in% sampled_seeds$seed[1], 
                                 modelindex %in% sampled_seeds$modelindex[1])

ggplot(d_sample %>% mutate(gen = gen - 50000), aes(x = time, y = Z)) +
  geom_line() +
  labs(x = "Developmental time", y = "Z expression") +
  theme_bw() +
  theme(text = element_text(size=20)) -> plt_ode

plt_ode <- plt_ode + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_ode, nframes = 42, duration = 10, width = 720, height = 720, 
        renderer = ffmpeg_renderer())
anim_save("ode_time.mp4", last_animation())

# plot phase diagram
Xstart <- 1
Xstop <- 6

plt_phase <- ggplot(d_sample %>% filter(gen == 50000) %>%
                      mutate(gen = gen - 50000, time = time/10), 
                    aes(x = X, y = Z)) + 
  geom_segment(aes(x = lag(X), y = lag(Z), xend = X, yend = Z), 
               arrow = arrow(length = unit(0.25, "cm"))) +
  xlab("X") + 
  ylab("Z") + 
  theme_bw() +
  theme(text = element_text(size=36), panel.spacing.x = unit(1, "lines"))

plt_phase <- plt_phase + transition_states(gen, wrap = F) +
  labs(title = "Generation: {closest_state}")

anim <- animate(plt_phase, nframes = 42, duration = 10, width = 1920, height = 1920,
                renderer = ffmpeg_renderer())
anim_save("ode_meanphase_test.mp4", anim)

# Analyse phase diagram differences as PCA

# Combine ODE measurements with phenotype combo data
d_com <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/slim_qg.csv")

d_com %>%
  group_by(seed, nloci, sigma) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  mutate(phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  ungroup() -> d_com

d_com %>% mutate(id = as_factor(paste(seed, modelindex, sep = "_"))) -> d_com

d_ode %>% mutate(modelindex = as_factor(modelindex),
                 seed = as_factor(seed)) -> d_ode

d_com_ode <- inner_join(d_com, d_ode %>% select(-c(phenomean, nloci, sigma)), 
                        by = c("gen", "modelindex", "seed"))

res.pca <- PCA(d_ode %>% select(10:17), scale.unit = T, graph = F)
fviz <- fviz_eig(res.pca, addlabels = TRUE)
ggsave("fviz_ode_shape.png", fviz, bg = "white")
var <- get_pca_var(res.pca)
head(var$contrib)

d_com_ode$pc1 <- res.pca$ind$coord[, 1]
d_com_ode$pc2 <- res.pca$ind$coord[, 2]

pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)

d_com_ode_short <- d_com_ode %>% filter(gen == 51950)

d_com_ode_short <- d_com_ode_short %>% 
  pivot_longer(cols = c(aZ, bZ, KZ, KXZ), names_to = "molTraitName", values_to = "molTraitValue")

saveRDS(d_com_ode, "d_com_ode.RDS")

ggplot(d_com_ode_short %>% filter(molTraitValue < 2.5) %>%
         mutate(gen = gen - 50000), aes(x = pc1, y = pc2, colour = molTraitValue)) +
  facet_grid(molTraitName~.) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(shape = 1, size = 2) +
  scale_colour_viridis_c() +
  labs(x = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), y = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), colour = "Population adapted") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_molTrait_pca

ggsave("moltrait_ode_pca.png", plt_molTrait_pca, width = 10, height = 10)
ggsave("moltrait_pca_adapted_facet.png", plt_isAdapted_pca + 
         facet_grid(nloci_cat~sigma_cat) +
         scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                breaks = NULL, labels = NULL)) +
         scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                breaks = NULL, labels = NULL)), 
       width = 10, height = 10)

plt_isAdapted_pca <- plt_isAdapted_pca + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_isAdapted_pca, nframes = 41, duration = 10, width = 1920, height = 1920, renderer = ffmpeg_renderer())
anim_save("moltrait_pca.mp4", last_animation())



# sample 10 from each nloci/sigma group and plot
sampled_seeds <- d_ode %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(nloci, sigma, seed, modelindex, id) %>%
  sample_n(10)



d_sample <- d_ode %>% filter(seed %in% sampled_seeds$seed, 
                   modelindex %in% sampled_seeds$modelindex)

highlight_ids <- sampled_seeds %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(id) %>%
  sample_n(1)

cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")


plt_ode_total <- ggplot(d_sample %>%
                          mutate(gen = gen - 50000, time = time/10),
                        aes(x = time, y = Z, colour = id, group = id)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line(linewidth = 2) +
  scale_colour_manual(values = rep(cc_ibm[3], 11), guide = "none") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1), 
                     sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  gghighlight(id %in% highlight_ids$id,
              calculate_per_facet = T, use_direct_label = F) +
  labs(x = "Developmental time", y = "Mean Z expression") +
  theme_bw() +
  theme(text = element_text(size=36), panel.spacing.x = unit(1, "lines"))

plt_ode_total <- plt_ode_total + transition_states(gen, wrap = F) +
  labs(title = "Generation: {closest_state}")

animate(plt_ode_total, nframes = 42, duration = 10, width = 1920, height = 1920, 
        renderer = ffmpeg_renderer())
anim_save("ode_time_total.mp4", last_animation())


# Make table of phase diagram with molecular trait and phenotype values
# Note: this is for population mean molecular trait/phenomean values
d_ode_phasemeasures <- d_ode %>% select(gen, seed, modelindex, nloci, sigma, phenomean,
                                        len_left, len_right, bottom_left_angle, bottom_right_angle,
                                        top_right_angle, top_left_angle) %>% 
  distinct(gen, seed, modelindex, nloci, sigma, .keep_all = T) %>%
  mutate(seed = as_factor(seed), modelindex = as_factor(modelindex))

d_qg <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/checkpoint/d_qg.RDS")

d_qg <- d_qg %>%
  select(gen, seed, modelindex, nloci, sigma, aZ, bZ, KZ, KXZ)

d_ode_phasemeasures <- inner_join(d_ode_phasemeasures, d_qg, 
                                  by = c("gen", "seed", "modelindex", "nloci", "sigma"))

d_ode_phasemeasures <- d_ode_phasemeasures %>%
  select(gen, seed, modelindex, nloci, sigma, aZ, bZ, KZ, KXZ, len_left, 
         len_right, bottom_left_angle, bottom_right_angle,
         top_right_angle, top_left_angle, phenomean) %>%
  rename(L1 = len_left, L2 = len_right, 
         thetaBL = bottom_left_angle, thetaBR = bottom_right_angle,
         thetaTR = top_right_angle, thetaTL = top_left_angle)

saveRDS(d_ode_phasemeasures, "d_phaseDiagramMeasures.RDS")
write_csv(d_ode_phasemeasures, "d_phaseDiagramMeasures.csv")
d_ode_phasemeasures <- readRDS("d_phaseDiagramMeasures.RDS")

# PCA of all the measures
# First filter out massive molecular trait values
d_ode_phasemeasures <- d_ode_phasemeasures %>%
  mutate(aZ99 = quantile(aZ, 0.99),
         bZ99 = quantile(bZ, 0.99),
         KZ99 = quantile(KZ, 0.99),
         KXZ99 = quantile(KXZ, 0.99)) %>%
  filter(aZ <= aZ99 & bZ <= bZ99, KZ <= KZ99, KXZ <= KXZ99)

res.pca <- PCA(d_ode_phasemeasures %>% select(10:15), scale.unit = T, graph = F)
fviz <- fviz_eig(res.pca, addlabels = TRUE)
fviz
ggsave("fviz_ode_shape.png", fviz, bg = "white")
var <- get_pca_var(res.pca)
var$contrib

# make dataframe to plot pca
d_pca <- d_ode_phasemeasures %>% select(-c(10:15, 17:20))

d_pca$pc1 <- res.pca$ind$coord[, 1]
d_pca$pc2 <- res.pca$ind$coord[, 2]

pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars

d_pca$nloci_cat <- cut(d_pca$nloci,
                             breaks=c(-Inf, 10, 50, 250, Inf))

d_pca$sigma_cat <- cut(d_pca$sigma,
                             breaks=c(-Inf, 0.05, 0.5, 1, Inf))

# Find range of data for axis limits and filtering massive molecular component values
pc1_lim <- quantile(d_pca$pc1, c(0.1, 0.99))
pc2_lim <- quantile(d_pca$pc2, c(0.1, 0.99))



# plot each timepoint for each molecular trait
dir.create("./anim")

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
foreach (i = seq_along(unique(d_pca$gen))) %dopar% {
  library(lattice); library(latticeExtra); library(dplyr); library(latex2exp)
  library(paletteer)
  this_gen <- unique(d_pca$gen)[i]
  plt_level_aZ <- levelplot(aZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
                              select(pc1, pc2, aZ), 
                            col.regions = paletteer_c("viridis::viridis", 100, -1),
                            panel = panel.levelplot.points, 
                            xlab = list(label = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), cex = 1.2), 
                            ylab = list(label = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), cex = 1.2), 
                            colorkey = list(space = "bottom",
                                            title = TeX("$\\alpha_Z$"),
                                            labels = list(cex = 1.2)),
                            main = paste("Generations post-optimum shift: ", this_gen),
                            par.settings = list(layout.heights = list(xlab.key.padding = 4),
                                                par.main.text = list(just = "left",
                                                                     x = grid::unit(5, "mm"))),
                            scales = list(tck = c(1, 0),
                                          cex = 1.2),
                            xlim = seq(from = pc1_lim[1], to = pc1_lim[2], length.out = 6),
                            ylim = seq(from = pc2_lim[1], to = pc2_lim[2], length.out = 5),
                            pretty = T) + 
    layer_(panel.2dsmoother(..., n = 200))
  
  plt_contour_aZ <- 
    contourplot(
      aZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
        select(pc1, pc2, aZ), 
      panel=panel.2dsmoother, col = "white",
      labels = list(col = "white",
                    cex = 1.2))
  
  plt_aZ <- plt_level_aZ + plt_contour_aZ
  # save the plot
  png(sprintf("./anim/pca_phase_aZ_gen%04d.png", i),
      type = "cairo",
      units = "in", 
      pointsize = 12,
      res = 300, 
      bg = "white",
      width = 7, 
      height = 7)
  print(plt_aZ)
  dev.off()
  
  rm(plt_aZ, plt_contour_aZ, plt_level_aZ)
  gc()

  # beta
  plt_level_bZ <- levelplot(bZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
                              select(pc1, pc2, bZ), 
                            col.regions = paletteer_c("viridis::viridis", 100, -1),
                            panel = panel.levelplot.points, 
                            xlab = list(label = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), cex = 1.2), 
                            ylab = list(label = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), cex = 1.2), 
                            colorkey = list(space = "bottom",
                                            title = TeX("$\\beta_Z$"),
                                            labels = list(cex = 1.2)),
                            main = paste("Generations post-optimum shift: ", this_gen),
                            par.settings = list(layout.heights = list(xlab.key.padding = 4),
                                                par.main.text = list(just = "left",
                                                                     x = grid::unit(5, "mm"))),
                            scales = list(tck = c(1, 0),
                                          cex = 1.2),
                            xlim = seq(from = pc1_lim[1], to = pc1_lim[2], length.out = 6),
                            ylim = seq(from = pc2_lim[1], to = pc2_lim[2], length.out = 5),
                            pretty = T) + 
    layer_(panel.2dsmoother(..., n = 200))
  
  plt_contour_bZ <- 
    contourplot(
      bZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
        select(pc1, pc2, bZ), 
      panel=panel.2dsmoother, col = "white",
      labels = list(col = "white",
                    cex = 1.2))
  
  plt_bZ <- plt_level_bZ + plt_contour_bZ
  # save the plot
  png(sprintf("./anim/pca_phase_bZ_gen%04d.png", i),
      type = "cairo",
      units = "in", 
      pointsize = 12,
      res = 300, 
      bg = "white",
      width = 7, 
      height = 7)
  print(plt_bZ)
  dev.off()
  
  rm(plt_bZ, plt_contour_bZ, plt_level_bZ)
  gc()

  # KZ
  plt_level_KZ <- levelplot(KZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
                              select(pc1, pc2, KZ), 
                            col.regions = paletteer_c("viridis::viridis", 100, -1),
                            panel = panel.levelplot.points, 
                            xlab = list(label = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), cex = 1.2), 
                            ylab = list(label = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), cex = 1.2), 
                            colorkey = list(space = "bottom",
                                            title = TeX("$K_Z$"),
                                            labels = list(cex = 1.2)),
                            main = paste("Generations post-optimum shift: ", this_gen),
                            par.settings = list(layout.heights = list(xlab.key.padding = 4),
                                                par.main.text = list(just = "left",
                                                                     x = grid::unit(5, "mm"))),
                            scales = list(tck = c(1, 0),
                                          cex = 1.2),
                            xlim = seq(from = pc1_lim[1], to = pc1_lim[2], length.out = 6),
                            ylim = seq(from = pc2_lim[1], to = pc2_lim[2], length.out = 5),
                            pretty = T) + 
    layer_(panel.2dsmoother(..., n = 200))
  
  plt_contour_KZ <- 
    contourplot(
      KZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
        select(pc1, pc2, KZ), 
      panel=panel.2dsmoother, col = "white",
      labels = list(col = "white",
                    cex = 1.2))
  
  plt_KZ <- plt_level_KZ + plt_contour_KZ
  # save the plot
  png(sprintf("./anim/pca_phase_KZ_gen%04d.png", i),
      type = "cairo",
      units = "in", 
      pointsize = 12,
      res = 300, 
      bg = "white",
      width = 7, 
      height = 7)
  print(plt_KZ)
  dev.off()
  
  rm(plt_KZ, plt_contour_KZ, plt_level_KZ)
  gc()
  
  # beta
  plt_level_KXZ <- levelplot(KXZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
                              select(pc1, pc2, KXZ), 
                            col.regions = paletteer_c("viridis::viridis", 100, -1),
                            panel = panel.levelplot.points, 
                            xlab = list(label = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), cex = 1.2), 
                            ylab = list(label = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), cex = 1.2), 
                            colorkey = list(space = "bottom",
                                            title = TeX("$K_{XZ}$"),
                                            labels = list(cex = 1.2)),
                            main = paste("Generations post-optimum shift: ", this_gen),
                            par.settings = list(layout.heights = list(xlab.key.padding = 4),
                                                par.main.text = list(just = "left",
                                                                     x = grid::unit(5, "mm"))),
                            scales = list(tck = c(1, 0),
                                          cex = 1.2),
                            xlim = seq(from = pc1_lim[1], to = pc1_lim[2], length.out = 6),
                            ylim = seq(from = pc2_lim[1], to = pc2_lim[2], length.out = 5),
                            pretty = T) + 
    layer_(panel.2dsmoother(..., n = 200))
  
  plt_contour_KXZ <- 
    contourplot(
      KXZ ~ pc1*pc2, d_pca %>% filter(gen == this_gen) %>%
        select(pc1, pc2, KXZ), 
      panel=panel.2dsmoother, col = "white",
      labels = list(col = "white",
                    cex = 1.2))
  
  plt_KXZ <- plt_level_KXZ + plt_contour_KXZ
  # save the plot
  png(sprintf("./anim/pca_phase_KXZ_gen%04d.png", i),
      type = "cairo",
      units = "in", 
      pointsize = 12,
      res = 300, 
      bg = "white",
      width = 7, 
      height = 7)
  print(plt_KXZ)
  dev.off()
  
  rm(plt_KXZ, plt_contour_KXZ, plt_level_KXZ)
  gc()
  
}
parallel::stopCluster(cl)

