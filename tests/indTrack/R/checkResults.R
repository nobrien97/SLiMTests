# Compare the QG results with those from newTestCross to make sure they are the same
library(tidyverse)
library(gghighlight)
library(latex2exp)

cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")

cc_10cols <- c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#012749")

d_match <- readRDS("/mnt/d/SLiMTests/tests/indTrack/d_com_match.RDS")

ggplot(d_match, aes(x = gen, y = phenomean, colour = model, group = seed)) +
  facet_grid(nloci~sigma) +
  geom_line()

d_new <- read_csv("/mnt/d/SLiMTests/tests/indTrack/slim_qg.csv", col_names = F)

names(d_new) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                  "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")

d_com_eg <- read_table("./combos.csv", col_names = F)
names(d_com_eg) <- c("model", "nloci", "sigma", "seed")
d_com_eg %>% 
  mutate(model = as_factor(model), 
         model = recode_factor(model, "\"Additive\"" = "Additive", "\"NAR\"" = "NAR")) -> d_com_eg

d_com_eg %>%
  ungroup() %>%
  mutate(seed = as_factor(seed),
         nloci = as_factor(nloci),
         sigma = as_factor(sigma),
         id = fct_cross(seed, model, nloci, sigma)) -> d_com_eg

d_new %>% mutate(model = d_com_eg$model[.$modelindex],
                 nloci = d_com_eg$nloci[.$modelindex],
                 sigma = d_com_eg$sigma[.$modelindex]) -> d_new

d_new %>%
  group_by(seed, model, nloci, sigma) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  mutate(phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  ungroup() -> d_new_adapted


d_new %>%
  group_by(seed, model, nloci, sigma) %>%
  filter(all(phenomean < 1.9 | phenomean > 2.1)) %>%
  ungroup() -> d_new_maladapted

sampled_seeds <- d_new_adapted %>% filter(gen > 49500, phenomean < 5, model == "NAR") %>%
  group_by(nloci, sigma) %>% 
  select(seed, modelindex) %>%
  sample_n(1)

d_new_adapted %>% filter(gen > 49500, phenomean < 5, model == "NAR") %>%
  mutate(gen = gen - 50000) %>% 
ggplot(aes(x = gen, y = phenomean, group = seed, colour = seed)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  gghighlight(modelindex %in% sampled_seeds$modelindex, calculate_per_facet = T, use_direct_label = F) +
  scale_colour_manual(values = rep(cc_ibm[3], 4), guide = "none") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adapted
plt_adapted

ggsave("pheno_adapted.png", plt_adapted, width = 8, height = 8)


View(d_new_adapted %>% filter(sigma == 0.1, model == "NAR", nloci == 100))

ggplot(d_new_maladapted %>% filter(gen > 49500, phenomean < 5), aes(x = gen, y = phenomean, colour = model, group = seed)) +
  facet_grid(nloci~sigma) +
  geom_line()


d_indPheno <- read_csv("/mnt/d/SLiMTests/tests/indTrack/slim_indPheno.csv", col_names = F)
names(d_indPheno) <- c("gen", "seed", "modelindex", "ind", "pheno", "aZ", "bZ", "KZ", "KXZ")

d_indPheno %>% mutate(model = d_com_eg$model[.$modelindex],
                 nloci = d_com_eg$nloci[.$modelindex],
                 sigma = d_com_eg$sigma[.$modelindex]) -> d_indPheno

d_indPheno %>%
  group_by(gen, seed, modelindex) %>%
  mutate(ind = 1:n()) %>%
  ungroup() %>%
  mutate(ind = as_factor(ind)) -> d_indPheno

ggplot(d_indPheno %>% filter(gen > 49500, modelindex == 69) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = pheno, colour = ind, group = ind)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_colour_manual(values = cc_10cols) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Phenotypic value", colour = "Individual") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_indpheno
plt_indpheno

ggsave("pheno_inds_adapted.png", plt_indpheno, width = 8, height = 8)

ggplot(d_indPheno %>% filter(gen > 49500, modelindex == 69) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = aZ, colour = ind, group = ind)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_colour_manual(values = cc_10cols) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = TeX("$\\alpha_Z$"), colour = "Individual") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_indaZ
plt_indaZ

ggsave("aZ_inds_adapted.png", plt_indaZ, width = 8, height = 8)

# Mutation data
d_indMuts <- read_csv("/mnt/d/SLiMTests/tests/indTrack/slim_indMut.csv", col_names = F)

names(d_indMuts) <- c("gen", "seed", "modelindex", "ind", "mutType", 
                      "mutid", "pos", "constraint", "originGen", "value", "chi", "fixGen")

d_indMuts %>% mutate(model = d_com_eg$model[.$modelindex],
                      nloci = d_com_eg$nloci[.$modelindex],
                      sigma = d_com_eg$sigma[.$modelindex]) -> d_indMuts

View(d_indMuts %>% filter(modelindex == 69, mutType > 1))

d_indMuts %>% 
  group_by(gen, seed, modelindex) %>%
  mutate(ind_id = match(ind, unique(ind))) %>%
  ungroup() %>%
  mutate(ind = as_factor(ind_id)) %>% select(!ind_id) -> d_indMuts

d_ind_com <- inner_join(d_indMuts, d_indPheno, by = c("gen", "seed", "modelindex", "model", "nloci", "sigma", "ind"))

saveRDS(d_ind_com, "/mnt/d/SLiMTests/tests/indTrack/d_ind_com.RDS")

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$"),
  TeX("$K_Z$"),
  TeX("$K_{XZ}$")
)


ggplot(d_ind_com %>% filter(gen >= 49500, modelindex == 69, mutType > 1, ind == 8) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = log(value), colour = as.factor(mutType))) +
  facet_grid(nloci~sigma) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  geom_point() +
  #geom_boxplot(notch = FALSE, alpha = 0.2, outlier.size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Molecular effect size", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts
plt_indmuts

ggsave("pheno_indMuts_adapted.png", plt_indmuts, width = 8, height = 8)

# look at the other highlighted examples
ggplot(d_indPheno %>% filter(gen >= 49500, modelindex %in% sampled_seeds$modelindex) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = pheno, colour = ind, group = ind)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_colour_manual(values = cc_10cols) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Phenotypic value", colour = "Individual") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_indpheno
plt_indpheno

ggsave("pheno_inds_adapted.png", plt_indpheno, width = 8, height = 8)

# Heritability
ggplot(d_new_adapted %>% filter(gen > 49500, modelindex %in% sampled_seeds$modelindex) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = meanH)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Phenotypic value", colour = "Individual") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_meanH
plt_meanH

ggplot(d_indPheno %>% filter(gen > 49500, modelindex %in% sampled_seeds$modelindex) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = aZ, colour = ind, group = ind)) +
  facet_grid(nloci~sigma, scales = "free") +
  geom_line() +
  scale_colour_manual(values = cc_10cols) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = TeX("$\\alpha_Z$"), colour = "Individual") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_indaZ
plt_indaZ

ggsave("aZ_inds_adapted.png", plt_indaZ, width = 8, height = 8)

ggplot(d_ind_com %>% filter(gen >= 49500, modelindex == 28, mutType > 1, ind == 5) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = log(value), colour = as.factor(mutType))) +
  facet_grid(nloci~sigma) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  geom_point() +
  #geom_boxplot(notch = FALSE, alpha = 0.2, outlier.size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Molecular effect size", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts
plt_indmuts

ggsave("pheno_indMuts_adapted.png", plt_indmuts, width = 8, height = 8)
