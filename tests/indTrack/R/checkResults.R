# Compare the QG results with those from newTestCross to make sure they are the same
library(tidyverse)
library(gghighlight)
library(latex2exp)
library(ggh4x)
library(cowplot)
library(grid)
library(gridExtra)

cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")

cc_10cols <- c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#012749")

d_match <- readRDS("/mnt/d/SLiMTests/tests/indTrack/d_com_match.RDS")

ggplot(d_match, aes(x = gen, y = phenomean, colour = model, group = seed)) +
  facet_grid(nloci~sigma) +
  geom_line()

d_new <- read_csv("/mnt/d/SLiMTests/tests/indTrack/slim_qg.csv", col_names = F)
d_new <- read_csv("/mnt/c/GitHub/SLiMTests/tests/indTrack/R/slim_qg.csv", col_names = F)

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

seed <- sample(1:.Machine$integer.max, 1)
set.seed(seed)
# 1263390450
#set.seed(1263390450)

sampled_seeds <- d_new_adapted %>% filter(gen > 49500, phenomean < 5, model == "NAR") %>%
  group_by(nloci, sigma) %>% 
  select(nloci, sigma, seed, modelindex) %>%
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
d_indPheno <- read_csv("/mnt/c/GitHub/SLiMTests/tests/indTrack/R/slim_indPheno.csv", col_names = F)

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
d_indMuts <- read_csv("/mnt/c/GitHub/SLiMTests/tests/indTrack/R/slim_indMut.csv", col_names = F)



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
saveRDS(d_ind_com, "/mnt/c/GitHub/SLiMTests/tests/indTrack/R/d_ind_com.RDS")

d_ind_com <- readRDS("/mnt/d/SLiMTests/tests/indTrack/d_ind_com.RDS")

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


d_ind_com %>% mutate(ind_mod = paste(modelindex, ind, sep = "_")) -> d_ind_com
d_indPheno %>% mutate(ind_mod = paste(modelindex, ind, sep = "_")) -> d_indPheno

chosen_inds <- c("66_6", "69_8", "71_9", "73_7")

# look at the other highlighted examples
ggplot(d_indPheno %>% filter(gen >= 49500, modelindex %in% sampled_seeds$modelindex) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = pheno, group = as.factor(gen))) +
  facet_grid(nloci~sigma) +
  geom_boxplot() +
  #gghighlight(ind_mod %in% chosen_inds, calculate_per_facet = T) +
  scale_colour_manual(values = cc_10cols) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Phenotypic value", colour = "Individual") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_indpheno
plt_indpheno

ggsave("pheno_inds_adapted_total.png", plt_indpheno, width = 8, height = 8)

# Heterozygosity
ggplot(d_new_adapted %>% filter(gen > 49500, modelindex %in% sampled_seeds$modelindex) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = meanH)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean population heterozygosity") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_meanH
plt_meanH

ggsave("meanH_adapted.png", plt_meanH, width = 8, height = 8)


ggplot(d_indPheno %>% filter(gen >= 49500, modelindex %in% sampled_seeds$modelindex) %>%
         pivot_longer(cols = c(aZ, bZ, KZ, KXZ), names_to = "molTrait", values_to = "molTraitVal") %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = molTraitVal, colour = molTrait, group = interaction(as.factor(gen), molTrait))) +
  facet_grid(nloci~sigma, scales = "free") +
  geom_boxplot(outlier.size = 0.5) +
  #coord_cartesian(ylim = c(0, 11)) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Molecular trait value", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indMol
plt_indMol

ggsave("molTrait_inds_adapted.png", plt_indMol, width = 8, height = 8)
ggsave("molTrait_inds_adapted_zoom.png", plt_indMol + coord_cartesian(ylim = c(0, 11)), width = 8, height = 8)


# Mutations for all individuals in the sampled models (red highlight in pheno figure)
ggplot(d_ind_com %>% filter(gen >= 49500, mutType != 1) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = log(value), colour = as.factor(mutType), 
           group = interaction(as.factor(gen), as.factor(mutType)))) +
  facet_grid(nloci~sigma) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  #geom_point() +
  geom_boxplot(outlier.size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Molecular effect size", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts
plt_indmuts

ggsave("pheno_indMuts_adapted_allInds.png", plt_indmuts, width = 8, height = 8)

# Mutations in one individual only - the maximum outlier
d_ind_com %>% 
  filter(gen >= 49500, mutType != 1) %>%
  group_by(gen, modelindex, mutType) %>%
  mutate(valMed = median(log(value)), 
         devVal = (log(value) - valMed)) %>%
  filter(devVal == max(devVal)) -> d_ind_com_filter


ggplot(d_ind_com %>% filter(gen >= 49500, mutType != 1, ind_mod %in% chosen_inds, sigma == 0.1) %>% 
         mutate(gen = gen - 50000, sigma = fct_drop(sigma)), 
       aes(x = gen, y = value, colour = as.factor(mutType))) +
  facet_grid(nloci~sigma, drop = F) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  geom_point(size = 0.5) +
  #geom_boxplot(outlier.size = 0.1) +
  #scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
  #                                       breaks = NULL, labels = NULL)) +
  #scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
  #                                       breaks = NULL, labels = NULL)) +
  labs(x = NULL, y = "Molecular effect size", colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts_s01
plt_indmuts_s01

ggplot(d_ind_com %>% filter(gen >= 49500, mutType != 1, ind_mod %in% chosen_inds, sigma == 1) %>% 
         mutate(gen = gen - 50000, sigma = fct_drop(sigma)), 
       aes(x = gen, y = value, colour = as.factor(mutType))) +
  facet_grid(nloci~sigma, drop = F) +
  scale_colour_manual(values = cc_ibm, labels = mutType_names) +
  geom_point(size = 0.5) +
  #geom_boxplot(outlier.size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  #scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
  #                                       breaks = NULL, labels = NULL)) +
  labs(x = NULL, y = NULL, colour = "Molecular trait") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_indmuts_s1
plt_indmuts_s1

indMuts_leg <- get_legend(plt_indmuts_s01)
y.grob <- textGrob("Number of QTLs", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Mutational effect size variance", 
                     gp = gpar(fontsize = 16))
x_lab <- textGrob("Generations post-optimum shift", 
                  gp = gpar(fontsize = 16))

plt_indmuts <- plot_grid(plt_indmuts_s01 + theme(strip.text.y = element_blank(), legend.position = "none"), 
                         plt_indmuts_s1 + theme(axis.title.y.left = element_blank(), legend.position = "none"), 
                         ncol = 2)
plt_indmuts <- plot_grid(top.grob, plt_indmuts, x_lab, indMuts_leg, ncol = 1, rel_heights = c(0.03, 1, 0.03, 0.1))
plt_indmuts

ggsave("pheno_indMuts_adapted_selectedInds.png", plt_indmuts, width = 8, height = 8, bg = "white")


# Number adapted
d_new %>%
  filter(model != "Additive") %>%
  group_by(seed, nloci, sigma) %>%
  mutate(isAdapted = any(gen >= 51800 & between(phenomean, 1.9, 2.1)),
         phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  distinct(.keep_all = T) %>%
  select(seed, nloci, sigma, isAdapted) %>%
  ungroup(seed) %>%
  summarise(pAdapted = mean(isAdapted)) %>%
  ungroup() -> d_isAdapted

ggplot(d_isAdapted, aes(x = nloci, y = pAdapted, fill = sigma)) +
  geom_col(position = position_dodge(0.9)) +
  labs(x = "Number of loci", y = "Probability of adaptation", fill = "Mutational effect size variance") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_pAdapted

ggsave("probAdapted.png", plt_pAdapted, width = 5, height = 5)

# Probability of adaptation by molecular trait

d_new %>%
  filter(model != "Additive", gen >= 49500) %>%
  group_by(seed, nloci, sigma) %>%
  mutate(isAdapted = any(gen >= 51800 & between(phenomean, 1.9, 2.1)),
         seed = as_factor(seed)) %>%
  ungroup() -> d_isAdapted


library(factoextra)
library(FactoMineR)

res.pca <- PCA(d_isAdapted %>% select(aZ, bZ, KZ, KXZ), scale.unit = T, graph = F)
fviz_eig(res.pca, addlabels = TRUE) -> fig
ggsave("fviz.png", fig, bg = "white")
var <- get_pca_var(res.pca)
head(var$contrib)


# https://tem11010.github.io/Plotting-PCAs/
d_isAdapted$pc1 <- res.pca$ind$coord[, 1]
d_isAdapted$pc2 <- res.pca$ind$coord[, 2]
pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars 

ggplot(d_isAdapted %>%
         mutate(gen = gen - 50000), aes(x = pc1, y = pc2, colour = isAdapted)) +
  #facet_grid2(nloci~sigma, scales = "free", independent = T) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(shape = 1) +
  labs(x = "PC 1 (36.1%)", y = "PC 2 (25.1%)", colour = "Population adapted") +
  theme_bw() +
  theme(text = element_text(size = 10)) -> plt_isAdapted_pca

ggsave("moltrait_pca_adapted.png", plt_isAdapted_pca, width = 10, height = 10)
