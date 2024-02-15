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
