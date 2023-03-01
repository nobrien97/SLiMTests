# Find example populations that hae and haven't adapted

library(tidyverse)
library(grid)
library(gridExtra)
library(latex2exp)
library(cowplot)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

# Filepath
path_add <- "/mnt/d/SLiMTests/tests/newTestCross/additive/getH2_newTestCross/data/"
path_net <- "/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/"

# Colour palette
cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")


# Data for groups that have/haven't adapted
path_both <- "/mnt/d/SLiMTests/tests/newTestCross/moreReps2/getH2_newTestCross/data/"

d_com_adapted_eg <- readRDS(paste0(path_both, "d_com_prefiltered_adapted_eg.RDS"))
d_com_maladapted_eg <- readRDS(paste0(path_both, "d_com_prefiltered_maladapted_eg.RDS"))
d_com_wasadapted_eg <- readRDS(paste0(path_both, "d_com_prefiltered_wasadapted_eg.RDS"))

d_com_adapted_eg %>%
  filter(gen >= 49500) %>%
  group_by(model, nloci, sigma) %>%
  mutate(seed = paste(model, nloci, sigma, seed, sep = "_")) %>%
  ggplot(aes(x = gen, y = phenomean, color = nloci, group = seed)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 1.9, ymax = 2.1, alpha = 0.2) +
  facet_grid(model~sigma) +
  #coord_cartesian(ylim = c(0, 2.5)) +
  geom_line() +
  scale_color_manual(values = cc_ibm) +
  geom_point(size = 0.25) +
  labs(x = "Generations after optimum shift", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

d_com_maladapted_eg %>%
  filter(gen >= 49500) %>%
  mutate(seed = paste(model, nloci, sigma, seed, sep = "_")) %>%
  ggplot(aes(x = gen, y = phenomean, color = nloci, group = seed)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 1.9, ymax = 2.1, alpha = 0.2) +
  facet_grid(model~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm) +
  #coord_cartesian(ylim = c(0, 2.5)) +
  geom_point(size = 0.25) +
  labs(x = "Generations after optimum shift", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

d_com_wasadapted_eg %>%
  filter(gen >= 49500) %>%
  group_by(model, nloci, sigma) %>%
  mutate(seed = paste(model, nloci, sigma, seed, sep = "_")) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = phenomean, color = nloci, group = seed)) +
  facet_grid(model~sigma) +
  geom_line() +
  geom_point(size = 0.25) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 1.9, ymax = 2.1, alpha = 0.2) +
  scale_color_manual(values = cc_ibm) +
  labs(x = "Generations after optimum shift", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

View(d_com_adapted_eg %>% filter (gen >= 49500, model == "NAR") %>% 
       distinct(gen, seed, model, nloci, sigma, .keep_all = T))

d_com_adapted <- readRDS(paste0(path_both, "d_com_prefiltered_adapted.RDS"))
View(d_com_adapted %>% filter(model == "NAR"))

d_com <- readRDS(paste0(path_both, "d_com_prefiltered.RDS"))

# Proportions of populations that adapted vs those that didn't
View(d_com %>%
       filter(any(gen > 51800)) %>%
       group_by(seed, model, nloci, sigma) %>%
       filter(row_number() == n()) %>%
       mutate(isAdapted = between(phenomean, 1.9, 2.1)) %>%
       ungroup(seed) %>%
       summarise(isAdapted = mean(isAdapted)))

# Get seeds for each group
path_ind <- "/mnt/c/GitHub/SLiMTests/tests/indTrack/R/"

d_com_adapted_eg %>%
  distinct(seed, model, nloci, sigma, .keep_all = T) %>%
  group_by(model, nloci, sigma) %>%
  select(seed) %>%
  arrange(model, nloci, sigma) -> adapted_seeds

write_csv(adapted_seeds, paste0(path_ind, "adapted_seeds.csv"))

d_com_maladapted_eg %>%
  distinct(seed, model, nloci, sigma, .keep_all = T) %>%
  group_by(model, nloci, sigma) %>%
  select(seed) %>%
  arrange(model, nloci, sigma) -> maladapted_seeds

write_csv(maladapted_seeds, paste0(path_ind, "maladapted_seeds.csv"))

d_com_wasadapted_eg %>%
  distinct(seed, model, nloci, sigma, .keep_all = T) %>%
  group_by(model, nloci, sigma) %>%
  select(seed) %>%
  arrange(model, nloci, sigma) -> wasadapted_seeds

write_csv(wasadapted_seeds, paste0(path_ind, "wasadapted_seeds.csv"))
