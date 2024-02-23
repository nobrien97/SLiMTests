# Have a look at the epistasis data to see which models differ the most

library(tidyverse)
library(latex2exp)
library(paletteer)

# functions
source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

# combos
d_combos <- read.table("../../R/combos.csv", header = F,
                            col.names = c("nloci", "tau", "r", "model"))

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/epistasisDensity/"
# data
d_epi_means <- read.table(paste0(DATA_PATH, "d_epi_mean.csv"), header = F, sep = ",",
                       col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                     "meanEW", "sdEW", "count"))


# Get pairwise differences between in mean epistasis between models 
d_epi_means_sbst <- d_epi_means %>% filter(optPerc == "[0.9, Inf)") %>% distinct()
d_epi_diff <- expand.grid(modelindex1 = d_epi_means_sbst$modelindex,
                          modelindex2 = d_epi_means_sbst$modelindex)

d_epi_diff <- d_epi_diff %>%
  group_by(modelindex1, modelindex2) %>%
  mutate(meanEP_1 = d_epi_means_sbst$meanEP[d_epi_means_sbst$modelindex == modelindex1],
         meanEP_2 = d_epi_means_sbst$meanEP[d_epi_means_sbst$modelindex == modelindex2],
         meanEW_1 = d_epi_means_sbst$meanEW[d_epi_means_sbst$modelindex == modelindex1],
         meanEW_2 = d_epi_means_sbst$meanEW[d_epi_means_sbst$modelindex == modelindex2],
         diffEP = meanEP_2 - meanEP_1,
         diffEW = meanEW_2 - meanEW_1)

# Outliers: model 396 generated some huge phenotypes...
boxplot(d_epi_diff$diffEP)

d_epi_diff %>% 
  filter(modelindex1 != 396 & modelindex2 != 396) -> d_epi_diff

boxplot(d_epi_diff$diffEP)
hist(d_epi_diff$diffEP, breaks = 100)

# Attach combos
# Adds the parameter combination to a dataframe
AddCombosToDiffDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model_1 = d_combos$model[modelindex1],
           nloci_1 = d_combos$nloci[modelindex1],
           tau_1 = d_combos$tau[modelindex1],
           r_1 = d_combos$r[modelindex1],
           model_2 = d_combos$model[modelindex2],
           nloci_2 = d_combos$nloci[modelindex2],
           tau_2 = d_combos$tau[modelindex2],
           r_2 = d_combos$r[modelindex2],)
}

d_epi_diff <- AddCombosToDiffDF(d_epi_diff)

# Plot distributions of differences between models
# heatmap: model1 vs model2, EP?

ggplot(d_epi_diff,
       aes(x = modelindex1, y = modelindex2, fill = diffEP)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Model 1", y = "Model 2", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Difference in variables
d_epi_diff <- d_epi_diff %>%
  mutate(r_diff = r_2 - r_1,
         tau_diff = tau_2 - tau_1,
         nloci_diff = nloci_2 - nloci_1,
         model_diff = paste(model_1, model_2, sep = "_"))

d_epi_diff_plt <- d_epi_diff %>%
  filter(r_diff >= 0, tau_diff >= 0, nloci_diff >= 0,
         model_diff != "K_ODE", 
         model_diff != "Add_ODE", 
         model_diff != "Add_K") %>%
  mutate(model_diff = ordered(model_diff, levels = c("ODE_ODE", "K_K", "Add_Add", 
                                                    "ODE_K", "ODE_Add", "K_Add")))

ggplot(d_epi_diff_plt  %>%
         group_by(r_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = r_diff, y = diffEP, colour = model_diff)) +
  geom_jitter(size = 1) +
  scale_colour_paletteer_d("awtools::a_palette") +
  labs(x = "Difference in recombination rate", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# The above plot shows that the difference in recombination rate doesn't
# particularly matter for the change in trait epistasis between models:
# small changes in rate can still lead to big differences in trait epistasis

ggplot(d_epi_diff_plt %>% 
         group_by(tau_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = tau_diff, y = diffEP, colour = model_diff)) +
  geom_jitter(size = 2) +
  scale_colour_paletteer_d("awtools::a_palette") +
  labs(x = "Difference in mutational effect variance", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# Mutational effect size variance has a couple of groups: the difference matters
# more than recombination

ggplot(d_epi_diff_plt %>% 
         group_by(nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = nloci_diff, y = diffEP, colour = model_diff)) +
  geom_jitter(size = 2) +
  scale_colour_paletteer_d("awtools::a_palette") +
  labs(x = "Difference in number of loci", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# ODE comparisons tend to cluster together, K has little difference


ggplot(d_epi_diff_plt %>% 
         group_by(r_1, r_2, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)) %>%
         filter(r_2 > 1e-03),
       aes(x = nloci_diff, y = diffEP, colour = model_diff)) +
  facet_grid(r_1~r_2) +
  geom_jitter(size = 2, width = 100, height = 0.15) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate 2", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate 1", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("awtools::a_palette") +
  labs(x = "Difference in number of loci", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# There is elevated average trait epistasis comparing models with high recombination
# 


ggplot(d_epi_diff %>% filter(r_diff >= 0, tau_diff >= 0, nloci_diff >= 0) %>%
         group_by(r_diff, tau_diff, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = r_diff, y = diffEP, colour = tau_diff)) +
  facet_grid(model_diff~nloci_diff) +
  geom_line() +
  labs(x = "Difference in recombination rate", y = "Difference in trait epistasis",
       colour = "Difference in mutational\neffect size variance") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))
