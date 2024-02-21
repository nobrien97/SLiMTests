library(tidyverse)
library(latex2exp)

d_epistasis_150 <- read_csv("../../calcMutationStats/R/d_epistasis_150.csv", 
                            col_names = c("gen", "seed", "modelindex", "mutType_ab",
                                          "meanEP", "meanEW", "sdEP", "sdEW"))

d_epistasis_freq_150 <- read_csv("../../calcMutationStats/R/d_epistasis_freqweight_150.csv", 
                                 col_names = c("gen", "seed", "modelindex", "mutType_ab",
                                               "meanEP", "meanEW", "sdEP", "sdEW"))


ggplot(d_epistasis_150, aes(x = meanEP, fill = mutType_ab)) +
  geom_histogram(bins = 50) +
  coord_cartesian(xlim = c(-0.75, 1.5)) +
  ggtitle("Uniformly sampled epistasis") +
  labs(x = TeX("Mean trait epistasis $(\\bar{\\epsilon}_Z)$"))


ggplot(d_epistasis_freq_150, aes(x = meanEP)) +
  geom_histogram(bins = 50) +
  coord_cartesian(xlim = c(-0.75, 1.5)) +
  ggtitle("Frequency-weighted epistasis") +
  labs(x = TeX("Mean trait epistasis $(\\bar{\\epsilon}_Z)$"))


ggplot(d_epistasis_150, aes(x = meanEW)) +
  geom_histogram(bins = 20) +
  coord_cartesian(xlim = c(0.03, 0.07)) +
  ggtitle("Uniformly sampled epistasis") +
  labs(x = TeX("Mean fitness epistasis $(\\bar{\\epsilon}_w)$"))


ggplot(d_epistasis_freq_150, aes(x = meanEW)) +
  geom_histogram(bins = 20) +
  coord_cartesian(xlim = c(0.03, 0.07)) +
  ggtitle("Frequency-weighted epistasis") +
  labs(x = TeX("Mean fitness epistasis $(\\bar{\\epsilon}_w)$"))



ggplot(d_epistasis_150, aes(x = sdEP^2)) +
  geom_histogram(bins = 50) +
  coord_cartesian(xlim = c(0, 0.2)) +
  ggtitle("Uniformly sampled epistasis") +
  labs(x = TeX("Trait epistasis variance $(Var{\\epsilon}_Z)$"))


ggplot(d_epistasis_freq_150, aes(x = sdEP^2)) +
  geom_histogram(bins = 50) +
  coord_cartesian(xlim = c(0, 0.2)) +
  ggtitle("Frequency-weighted epistasis") +
  labs(x = TeX("Trait epistasis variance $(Var{\\epsilon}_Z)$"))


ggplot(d_epistasis_150, aes(x = sdEW^2)) +
  geom_histogram(bins = 20) +
  coord_cartesian(xlim = c(0, 7e-04)) +
  ggtitle("Uniformly sampled epistasis") +
  labs(x = TeX("Fitness epistasis variance $(Var{\\epsilon}_w)$"))


ggplot(d_epistasis_freq_150, aes(x = sdEW^2)) +
  geom_histogram(bins = 20) +
  coord_cartesian(xlim = c(0, 7e-04)) +
  ggtitle("Frequency-weighted epistasis") +
  labs(x = TeX("Fitness epistasis variance $(Var{\\epsilon}_w)$"))

# Density plots for epistasis
d_epistasis_density <- read_csv("/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/d_epi_density194.csv", 
                                col_names = c("optPerc", "modelindex", "mutType_ab",
                                              "wa_x","wa_y","wb_x","wb_y",
                                              "wab_x", "wab_y", "Pa_x", "Pa_y",
                                              "Pb_x", "Pb_y", "Pab_x", "Pab_y",
                                              "ew_x", "ew_y", "ep_x", "ep_y"))

# mock data
# d_epistasis <- data.frame(
#   optPerc = rep(dpdt$V1, each = 4800/5),
#   wa = runif(4800),
#   wb = runif(4800),
#   wab = runif(4800),
#   Pwt = rnorm(4800),
#   Pa = rnorm(4800),
#   Pb = rnorm(4800),
#   Pab = rnorm(4800)
# )
# d_epistasis$ew <- log(d_epistasis$wab) - (log(d_epistasis$wa) + log(d_epistasis$wb))
# d_epistasis$ep <- (d_epistasis$Pab - d_epistasis$Pwt) - 
#   ((d_epistasis$Pa - d_epistasis$Pwt) + (d_epistasis$Pb - d_epistasis$Pwt))

ggplot(d_epistasis_density, 
       aes(x = ep_x, y = optPerc, height = ep_y)) +
  facet_grid(.~mutType_ab) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# weighted by frequency
d_epistasis_density <- read_csv("/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/d_epi_freqweight_density194.csv", 
                                col_names = c("optPerc", "modelindex", "mutType_ab",
                                              "wa_x","wa_y","wb_x","wb_y",
                                              "wab_x", "wab_y", "Pa_x", "Pa_y",
                                              "Pb_x", "Pb_y", "Pab_x", "Pab_y",
                                              "ew_x", "ew_y", "ep_x", "ep_y"))

ggplot(d_epistasis_density, 
       aes(x = Pab_x, y = optPerc, height = Pab_y)) +
  facet_grid(.~mutType_ab) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)
