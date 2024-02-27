# Have a look at the epistasis data to see which models differ the most

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)

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

d_epi_diff <- crossing(modelindex1 = d_epi_means_sbst$modelindex, 
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

# Identify combos which are the same but inverted and remove the inverted one
d_epi_diff %>%
  filter(modelindex1 != modelindex2) %>%
  mutate(absDiffEP = abs(diffEP)) %>%
  arrange(absDiffEP) %>%
  ungroup() %>%
  mutate(absGroup = row_number() - (row_number() %% 2 != 1)) %>%
  group_by(absGroup) %>%
  filter(diffEP == absDiffEP) -> d_epi_diff


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
  filter(r_diff >= 0, tau_diff >= 0, nloci_diff >= 0) %>%
  mutate(model_diff = ordered(model_diff, levels = c("ODE_ODE", "K_K", "Add_Add", 
                                                    "ODE_K", "ODE_Add", "K_Add",
                                                    "K_ODE", "Add_ODE", "Add_K")))

ggplot(d_epi_diff_plt  %>%
         group_by(r_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = r_diff, y = diffEP, colour = model_diff)) +
  geom_point() +
  #geom_jitter(size = 2, width = 0.005) +
  scale_colour_paletteer_d("tvthemes::simpsons") +
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
  scale_colour_paletteer_d("tvthemes::simpsons") +
  labs(x = "Difference in mutational effect variance", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# Mutational effect size variance is important: big differences between models
# in mutational variance leads to big differences in epistasis

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

# Differences between models in number of loci doesn't have any effect on epistasis
# probably because of the standardised mutation rate


ggplot(d_epi_diff_plt %>% 
         group_by(r_1, r_2, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
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
# with any other level of recombination: recombination has an effect, difference
# doesn't matter, but the magnitude. Pattern flips with 0.01 and 0.1 according to
# number of loci
# With a big difference in nloci (one model has high, the other has low), 
# 0.1 recombination rate produces the highest trait epistasis
# with 0.01 recombination rate, smaller differences in nloci produce the greatest
# difference in epistasis

ggplot(d_epi_diff_plt %>%
         group_by(r_1, r_2, nloci_1, nloci_2) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = as.factor(log10(r_1)), y = as.factor(log10(r_2)), fill = diffEP)) +
  facet_grid(nloci_1~nloci_2) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Differences in trait epistasis between models correlate with very 
# high recombination rates but only in cases where at least one model has
# many loci

ggplot(d_epi_diff_plt %>%
         group_by(r_1, r_2, tau_1, tau_2) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = as.factor(log10(r_1)), y = as.factor(log10(r_2)), fill = diffEP)) +
  facet_grid(tau_1~tau_2) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16))

# So the trait epistasis difference requires one of the models to have high
# mutational variance and high recombination

d_epi_means_plt <- AddCombosToDF(d_epi_means_sbst %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier
d_epi_means_plt <- d_epi_means_plt %>%
  filter(modelindex != 396)

ggplot(d_epi_means_plt %>% filter(meanEP < 3),
       aes(x = as.factor(log10(r)), y = as.factor(nloci), fill = meanEP)) +
  facet_grid(as.factor(tau)~as.factor(model)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate (log10)", y = "Number of loci", 
       fill = "Trait epistasis") +
  guides(fill = guide_colourbar(barwidth = 10)) +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Two models with high trait epistasis: ODE, many loci, high mutational variance

ggplot(d_epi_means_plt %>% filter(sdEP < 15),
       aes(x = as.factor(log10(r)), y = as.factor(nloci), fill = sdEP^2)) +
  facet_grid(as.factor(tau)~as.factor(model)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate (log10)", y = "Number of loci", 
       fill = "Fitness epistasis variance") +
  guides(fill = guide_colourbar(barwidth = 10)) +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Find the distributions of epistasis components for these outliers,
# compare across models

# Find the models of interest
# tau = 1.25
# r = 0.01, 0.1, 1e-10
# nloci = 1024


# Load in density data
d_epi_dens <- fread(paste0(DATA_PATH, "d_epi_density.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "mutType_ab", 
                                        "wa_x", "wa_y", "wb_x", "wb_y", "wab_x",
                                        "wab_y", "Pa_x", "Pa_y", "Pb_x", "Pb_y",
                                        "Pab_x", "Pab_y",
                                        "ew_x", "ew_y", "ep_x", "ep_y" 
                                        ))

# Remove duplicates
d_epi_dens <- d_epi_dens %>%
  distinct() 

d_epi_dens <- AddCombosToDF(d_epi_dens %>% 
                                   mutate(modelindex = as.factor(modelindex)))

d_epi_dens_sbst <- d_epi_dens %>%
  filter(optPerc == "[0.9, Inf)", 
         tau == 1.25, nloci == 1024)

ggplot(d_epi_dens_sbst %>% filter(model == "ODE"), 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(.~mutType_ab) +
  coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

ggplot(d_epi_dens_sbst %>% filter(model == "K"), 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(.~mutType_ab) +
  coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# frequency weighted sampling