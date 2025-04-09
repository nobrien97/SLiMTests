# Have a look at the epistasis data to see which models differ the most

library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)

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
d_epi_means_sbst <- d_epi_means %>% 
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)") %>% 
  distinct()

r_subsample <- c(1e-10, 1e-5, 1e-1)

d_epi_diff <- crossing(optPerc = d_epi_means_sbst$optPerc,
                       modelindex1 = d_epi_means_sbst$modelindex, 
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
  facet_grid(.~optPerc) + 
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
  mutate(model_diff = ordered(model_diff, levels = c("K-_K-", "K+_K+", "Add_Add", 
                                                    "K-_K+", "K-_Add", "K+_Add",
                                                    "K+_K-", "Add_K-", "Add_K+")))

ggplot(d_epi_diff_plt  %>%
         group_by(optPerc, r_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = r_diff, y = diffEP, colour = model_diff)) +
  facet_grid(.~optPerc) +
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
         group_by(optPerc, tau_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = tau_diff, y = diffEP, colour = model_diff)) +
  facet_grid(.~optPerc) +
  geom_point(size = 2) +
  scale_colour_paletteer_d("tvthemes::simpsons") +
  labs(x = "Difference in mutational effect variance", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# Mutational effect size variance is important: big differences between models
# in mutational variance leads to big differences in epistasis

ggplot(d_epi_diff_plt %>% 
         group_by(optPerc, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = nloci_diff, y = diffEP, colour = model_diff)) +
  facet_grid(.~optPerc) +
  geom_point(size = 2) +
  scale_colour_paletteer_d("awtools::a_palette") +
  labs(x = "Difference in number of loci", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# Differences between models in number of loci doesn't have any effect on epistasis
# probably because of the standardised mutation rate

ggplot(d_epi_diff_plt %>% 
         group_by(optPerc, r_1, r_2, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = nloci_diff, y = diffEP, colour = model_diff)) +
  facet_nested(log10(r_1)~optPerc+log10(r_2)) +
  geom_point(shape = 5) +
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
         group_by(optPerc, r_1, r_2, nloci_1, nloci_2) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = as.factor(log10(r_1)), y = as.factor(log10(r_2)), fill = diffEP)) +
  facet_nested(nloci_1~optPerc+nloci_2) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 20))

# Differences in trait epistasis between models correlate with very 
# high recombination rates but only in cases where at least one model has
# many loci

ggplot(d_epi_diff_plt %>%
         group_by(optPerc, r_1, r_2, tau_1, tau_2) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = as.factor(log10(r_1)), y = as.factor(log10(r_2)), fill = diffEP)) +
  facet_nested(tau_1~optPerc+tau_2) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 20))


# So the trait epistasis difference requires one of the models to have high
# mutational variance and high recombination

d_epi_means_plt <- AddCombosToDF(d_epi_means %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier
d_epi_means_plt <- d_epi_means_plt %>%
  filter(modelindex != 396)


d_epi_means_plt <- d_epi_means_plt %>%
  filter(modelindex != 396, tau == 0.0125, r %in% r_subsample)

d_epi_means_plt_sum <- d_epi_means_plt %>%
  group_by(model, r) %>%
  summarise(meanEWBar = mean(meanEW),
            CIEWBar = CI(meanEW),
            varEWBar = var(meanEW),
            n = n())

ggplot(d_epi_means_plt %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = model, y = meanEW, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_epi_means_plt_sum %>% 
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanEWBar, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", y = "Average fitness epistasis", colour = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave("plt_ew_sml.png", width = 4, height = 7, device = png)

# For preso
d_epi_means_plt <- d_epi_means_plt %>%
  filter(modelindex != 396, tau == 0.0125, r %in% c(1e-10, 0.1))

d_epi_means_plt_sum <- d_epi_means_plt %>%
  group_by(model, r) %>%
  summarise(meanEWBar = mean(meanEW),
            CIEWBar = CI(meanEW),
            n = n())

ggplot(d_epi_means_plt %>% 
         mutate(r_title = "Recombination rate",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"), 
       aes(x = model, y = meanEW, colour = model)) +
  facet_nested(r_title + r ~ .,
               labeller = labeller(r = as_labeller(c(`1e-10` = "Low",
                                                     `0.1` = "High")))) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_epi_means_plt_sum %>% 
               mutate(r_title = "Recombination rate",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanEWBar, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      guide = "none") +
  labs(x = "Model", y = "Average fitness epistasis", colour = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave("plt_ew_talk.png", width = 4, height = 4, device = png)


# Box plots
ggplot(d_epi_means_plt, aes(x = model, y = meanEP)) +
  facet_grid(.~optPerc) +
  geom_boxplot() +
  labs(x = "Model", y = "Average trait epistasis") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))

# Zoom in, ignore big outlier
ggplot(d_epi_means_plt %>% filter(meanEP < 3), aes(x = model, y = meanEP)) +
  facet_grid(.~optPerc) +
  geom_boxplot() +
  labs(x = "Model", y = "Average trait epistasis") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))


ggplot(d_epi_means_plt, aes(x = model, y = meanEW)) +
  facet_grid(.~optPerc) +
  geom_boxplot() +
  labs(x = "Model", y = "Average fitness epistasis") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))

d_epi_means_plt$model <- as.factor(d_epi_means_plt$model)
d_epi_means_plt$tau <- factor(d_epi_means_plt$tau, 
                              levels = c("0.0125", "0.125", "1.25"))
d_epi_means_plt$nloci_title <- "Number of loci"
d_epi_means_plt$tau_title <- "Mutational effect variance"
d_epi_means_plt$r_title <- "Recombination rate (log10)"

ggplot(d_epi_means_plt %>% filter(meanEP < 3), 
       aes(x = model, 
           y = meanEP, colour = tau)) +
  facet_nested(r_title + log10(r)~optPerc + nloci_title + nloci) +
  geom_point() +
  geom_errorbar(aes(ymin = meanEP - sdEP/sqrt(count), 
                    ymax = meanEP + sdEP/sqrt(count))) +
  labs(x = "Model", y = "Average trait epistasis",
       colour = "Mutational effect variance") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_epi_means_plt, 
       aes(x = model, 
           y = meanEW, colour = tau)) +
  facet_nested(r_title + log10(r)~optPerc + nloci_title + nloci) +
  geom_point() +
  geom_errorbar(aes(ymin = meanEW - sdEW/sqrt(count), 
                    ymax = meanEW + sdEW/sqrt(count))) +
  labs(x = "Model", y = "Average fitness epistasis", 
       colour = "Mutational effect variance") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


ggplot(d_epi_means_plt %>% filter(meanEP < 3),
       aes(x = as.factor(log10(r)), y = as.factor(nloci), fill = meanEP)) +
  facet_nested(as.factor(tau)~optPerc+as.factor(model)) +
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
  facet_nested(as.factor(tau)~optPerc+as.factor(model)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate (log10)", y = "Number of loci", 
       fill = "Fitness epistasis variance") +
  guides(fill = guide_colourbar(barwidth = 10)) +
  theme(legend.position = "bottom", text = element_text(size = 16))


# compare VA to VI
d_epi_means_plt <- d_epi_means_plt %>% 
  mutate(VI = sdEP^2,
         optPerc = as.factor(optPerc))

source("h2Figures.R")

d_VA_VI <- inner_join(d_h2 %>% 
                        group_by(modelindex, model, 
                                 nloci, tau, r, width) %>%
                        summarise(VA = mean(VA_Z)), 
                      d_epi_means_plt %>%
                        group_by(modelindex, model, 
                                 nloci, tau, r, width) %>%
                        summarise(VI = mean(VI)),  # mean across time points
                      by = c("modelindex", 
                             "model", "nloci", "tau", "r", "width"))

# Proportions of VA vs VI
ggplot(d_VA_VI %>% 
         rowwise() %>%
         mutate(propVA = VA / ( VI),
                propVI = VI / ( VA)) %>%
         mutate(r_title = "Recombination rate",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance"),
       aes(x = as.factor(r), y = propVI, colour = model)) +
  facet_nested(tau_title + tau ~ width_title + width, scales = "free") +
  geom_point(position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  #coord_cartesian(ylim = c(0, 0.0001)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom")



# Now the frequency-weighted sampling data:
# Not much difference to the non-weighted results
{
  
d_epi_freq_means <- read.table(paste0(DATA_PATH, "d_epi_freqweight_mean.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "count"))

# Get pairwise differences between in mean epistasis between models 
d_epi_freq_means_sbst <- d_epi_freq_means %>% 
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)") %>% distinct()

d_epi_freq_diff <- crossing(optPerc = d_epi_freq_means_sbst$optPerc,
                            modelindex1 = d_epi_freq_means_sbst$modelindex, 
                       modelindex2 = d_epi_freq_means_sbst$modelindex)

d_epi_freq_diff <- d_epi_freq_diff %>%
  group_by(modelindex1, modelindex2) %>%
  mutate(meanEP_1 = d_epi_means_sbst$meanEP[d_epi_means_sbst$modelindex == modelindex1],
         meanEP_2 = d_epi_means_sbst$meanEP[d_epi_means_sbst$modelindex == modelindex2],
         meanEW_1 = d_epi_means_sbst$meanEW[d_epi_means_sbst$modelindex == modelindex1],
         meanEW_2 = d_epi_means_sbst$meanEW[d_epi_means_sbst$modelindex == modelindex2],
         diffEP = meanEP_2 - meanEP_1,
         diffEW = meanEW_2 - meanEW_1)

# Outliers: model 396 generated some huge phenotypes...
boxplot(d_epi_freq_diff$diffEP)

d_epi_freq_diff %>% 
  filter(modelindex1 != 396 & modelindex2 != 396) -> d_epi_freq_diff

boxplot(d_epi_freq_diff$diffEP)
hist(d_epi_freq_diff$diffEP, breaks = 100)

# Identify combos which are the same but inverted and remove the inverted one
d_epi_freq_diff %>%
  filter(modelindex1 != modelindex2) %>%
  mutate(absDiffEP = abs(diffEP)) %>%
  arrange(absDiffEP) %>%
  ungroup() %>%
  mutate(absGroup = row_number() - (row_number() %% 2 != 1)) %>%
  group_by(absGroup) %>%
  filter(diffEP == absDiffEP) -> d_epi_freq_diff


# Attach combos
d_epi_freq_diff <- AddCombosToDiffDF(d_epi_freq_diff)


# Plot distributions of differences between models
# heatmap: model1 vs model2, EP?

ggplot(d_epi_freq_diff,
       aes(x = modelindex1, y = modelindex2, fill = diffEP)) +
  facet_grid(.~optPerc) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Model 1", y = "Model 2", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Difference in variables
d_epi_freq_diff <- d_epi_freq_diff %>%
  mutate(r_diff = r_2 - r_1,
         tau_diff = tau_2 - tau_1,
         nloci_diff = nloci_2 - nloci_1,
         model_diff = paste(model_1, model_2, sep = "_"))

d_epi_freq_diff_plt <- d_epi_freq_diff %>%
  filter(r_diff >= 0, tau_diff >= 0, nloci_diff >= 0) %>%
  mutate(model_diff = ordered(model_diff, levels = c("ODE_ODE", "K_K", "Add_Add", 
                                                     "ODE_K", "ODE_Add", "K_Add",
                                                     "K_ODE", "Add_ODE", "Add_K")))

ggplot(d_epi_freq_diff_plt  %>%
         group_by(optPerc, r_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = r_diff, y = diffEP, colour = model_diff)) +
  facet_grid(.~optPerc) +
  geom_point() +
  scale_colour_paletteer_d("tvthemes::simpsons") +
  labs(x = "Difference in recombination rate", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

  # Same as before, difference in recombination rate between models doesn't
  # relate to difference in trait epistasis between model pairs

ggplot(d_epi_freq_diff_plt %>% 
         group_by(optPerc, tau_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = tau_diff, y = diffEP, colour = model_diff)) +
  facet_grid(.~optPerc) +
  geom_point(size = 2) +
  scale_colour_paletteer_d("tvthemes::simpsons") +
  labs(x = "Difference in mutational effect variance", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# And same as before, when model pairs differ in mutational effect variance, there
# is a big difference in trait epistasis for ODE model comparisons

ggplot(d_epi_freq_diff_plt %>% 
         group_by(optPerc, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = nloci_diff, y = diffEP, colour = model_diff)) +
  facet_grid(.~optPerc) +
  geom_point(size = 2) +
  scale_colour_paletteer_d("awtools::a_palette") +
  labs(x = "Difference in number of loci", y = "Difference in trait epistasis",
       colour = "Model difference") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, byrow = T))

# Same as before, nloci doesn't particularly matter for trait epistasis

ggplot(d_epi_freq_diff_plt %>% 
         group_by(optPerc, r_1, r_2, nloci_diff, model_diff) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = nloci_diff, y = diffEP, colour = model_diff)) +
  facet_nested(r_1~optPerc + r_2) +
  geom_point(size = 2, shape = 5) +
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

# Same result as before: high recombination in combination with differences in 
# number of loci produces those models with higher trait epistasis

ggplot(d_epi_freq_diff_plt %>%
         group_by(optPerc, r_1, r_2, nloci_1, nloci_2) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = as.factor(log10(r_1)), y = as.factor(log10(r_2)), fill = diffEP)) +
  facet_nested(nloci_1~optPerc + nloci_2) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Same as before: when one model has 1024 loci, and the other has 0.1 recombination,
# there is a big difference between the models in trait epistasis, regardless
# of the level of recombination in the first model and the number of loci of the
# second model
# when one model has 256 loci, and the other has 0.01 recombination, there is a
# similar difference between models as above.

ggplot(d_epi_freq_diff_plt %>%
         group_by(optPerc, r_1, r_2, tau_1, tau_2) %>%
         summarise(diffEP = mean(diffEP)),
       aes(x = as.factor(log10(r_1)), y = as.factor(log10(r_2)), fill = diffEP)) +
  facet_nested(tau_1~optPerc + tau_2) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "Difference in trait epistasis\nbetween models") +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Same as before: When one model has high mutational variance, and recombination
# is high in the other model, there is a big difference in trait epistasis regardless
# of the mutational variance of the second model and recombination rate of the 
# first model

d_epi_freq_means_plt <- AddCombosToDF(d_epi_freq_means_sbst %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Remove outlier
d_epi_freq_means_plt <- d_epi_freq_means_plt %>%
  filter(modelindex != 396)

d_epi_freq_means_plt$model <- as.factor(d_epi_freq_means_plt$model)
d_epi_freq_means_plt$tau <- factor(d_epi_freq_means_plt$tau, 
                              levels = c("0.0125", "0.125", "1.25"))
d_epi_freq_means_plt$nloci_title <- "Number of loci"
d_epi_freq_means_plt$tau_title <- "Mutational effect variance"
d_epi_freq_means_plt$r_title <- "Recombination rate (log10)"

ggplot(d_epi_freq_means_plt %>% filter(meanEP < 3), 
       aes(x = model, 
           y = meanEP, colour = tau)) +
  facet_nested(r_title + log10(r)~optPerc+nloci_title + nloci) +
  geom_point() +
  geom_errorbar(aes(ymin = meanEP - sdEP/sqrt(count), 
                    ymax = meanEP + sdEP/sqrt(count))) +
  labs(x = "Model", y = "Average trait epistasis",
       colour = "Mutational effect variance") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggplot(d_epi_freq_means_plt, 
       aes(x = model, 
           y = meanEW, colour = tau)) +
  facet_nested(r_title + log10(r)~optPerc + nloci_title + nloci) +
  geom_point() +
  geom_errorbar(aes(ymin = meanEW - sdEW/sqrt(count), 
                    ymax = meanEW + sdEW/sqrt(count))) +
  labs(x = "Model", y = "Average fitness epistasis", 
       colour = "Mutational effect variance") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_colour_paletteer_d("nationalparkcolors::Badlands") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")



ggplot(d_epi_freq_means_plt %>% filter(meanEP < 3),
       aes(x = as.factor(log10(r)), y = as.factor(nloci), fill = meanEP)) +
  facet_nested(as.factor(tau)~optPerc + as.factor(model)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate (log10)", y = "Number of loci", 
       fill = "Trait epistasis") +
  guides(fill = guide_colourbar(barwidth = 20)) +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Quite different to the un-frequency-weighted model: much smaller mean epistasis
# still variance in the same ODE model, but lower magnitude.


ggplot(d_epi_freq_means_plt %>% filter(sdEP < 15),
       aes(x = as.factor(log10(r)), y = as.factor(nloci), fill = sdEP^2)) +
  facet_nested(as.factor(tau)~optPerc+as.factor(model)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Recombination rate (log10)", y = "Number of loci", 
       fill = "Fitness epistasis variance") +
  guides(fill = guide_colourbar(barwidth = 10)) +
  theme(legend.position = "bottom", text = element_text(size = 16))

# Again, lower variance overall, but variance distributed in the same two models
# as non-freq-weighted

  }

# Find the distributions of epistasis components for outliers in both freq and
# non-freq weighted data,
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
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 1.25, nloci == 1024)

# Recombination 1e-10 - 1e-2 not really any effect, shifts e_p positive
# more synergistic interactions segregating in population with very high
# recombination: with this level of recombination, might be more realistic
# to separate these overshooting combinations so they aren't as deleterious
# as in a more linked background
ggplot(d_epi_dens_sbst %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# Not the case in the K environment, no difference at all with recombination:
# maybe with more moving parts, it is more difficult to navigate the interactions
# between molecular components, even with high recombination?
# need a more complex example to see if this holds up
ggplot(d_epi_dens_sbst %>% filter(model == "K"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  coord_cartesian(xlim = c(-5, 10)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)



# frequency weighted sampling: is it any different?
d_epi_freq_dens <- fread(paste0(DATA_PATH, "d_epi_freqweight_density.csv"), header = F, sep = ",",
                    col.names = c("optPerc", "modelindex", "mutType_ab", 
                                  "wa_x", "wa_y", "wb_x", "wb_y", "wab_x",
                                  "wab_y", "Pa_x", "Pa_y", "Pb_x", "Pb_y",
                                  "Pab_x", "Pab_y",
                                  "ew_x", "ew_y", "ep_x", "ep_y" 
                    ))

# Remove duplicates
d_epi_freq_dens <- d_epi_freq_dens %>%
  distinct() 

d_epi_freq_dens <- AddCombosToDF(d_epi_freq_dens %>% 
                              mutate(modelindex = as.factor(modelindex)))

d_epi_freq_dens_sbst <- d_epi_freq_dens %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 1.25, nloci == 1024) %>%
  mutate(ep_y = as.numeric(ep_y))

# Same result for recombination with the frequency adjusted model, but
# the distributions are all narrower: rare alleles contribute to large
# trait epistasis
ggplot(d_epi_freq_dens_sbst %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# Not the case in the K environment, no difference at all with recombination:
# maybe with more moving parts, it is more difficult to navigate the interactions
# between molecular components, even with high recombination?
# need a more complex example to see if this holds up
ggplot(d_epi_freq_dens_sbst %>% filter(model == "K"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  coord_cartesian(xlim = c(-5, 10)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)


# load data averaged over mutType
d_epi_dens_nomut <- fread(paste0(DATA_PATH, "d_epi_density_nomut.csv"), header = F, sep = ",",
                         col.names = c("optPerc", "modelindex", 
                                       "wa_x", "wa_y", "wb_x", "wb_y", "wab_x",
                                       "wab_y", "Pa_x", "Pa_y", "Pb_x", "Pb_y",
                                       "Pab_x", "Pab_y",
                                       "ew_x", "ew_y", "ep_x", "ep_y" 
                         ))

d_epi_dens_nomut <- d_epi_dens_nomut %>%
  distinct() 

d_epi_dens_nomut <- AddCombosToDF(d_epi_dens_nomut %>% 
                                   mutate(modelindex = as.factor(modelindex)))

d_epi_dens_nomut_sbst <- d_epi_dens_nomut %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 1.25, nloci == 1024) %>%
  mutate(ep_y = as.numeric(ep_y))


ggplot(d_epi_dens_nomut_sbst %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(.~optPerc) +
  #coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# load data averaged over mutType
d_epi_freqweight_dens_nomut <- fread(paste0(DATA_PATH, "d_epi_freqweight_density_nomut.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", 
                                        "wa_x", "wa_y", "wb_x", "wb_y", "wab_x",
                                        "wab_y", "Pa_x", "Pa_y", "Pb_x", "Pb_y",
                                        "Pab_x", "Pab_y",
                                        "ew_x", "ew_y", "ep_x", "ep_y" 
                          ))

d_epi_freqweight_dens_nomut <- d_epi_freqweight_dens_nomut %>%
  distinct() 

d_epi_freqweight_dens_nomut <- AddCombosToDF(d_epi_freqweight_dens_nomut %>% 
                                    mutate(modelindex = as.factor(modelindex)))

d_epi_freqweight_dens_nomut_sbst <- d_epi_freqweight_dens_nomut %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 1.25, nloci == 1024) %>%
  mutate(ep_y = as.numeric(ep_y))

d_epi_dens_nomut_sbst <- d_epi_dens_nomut_sbst %>%
  mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE"))

ggplot(d_epi_dens_nomut_sbst, 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-5, 25)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 1024") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Trait epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

d_epi_freqweight_dens_nomut_sbst <- d_epi_freqweight_dens_nomut_sbst %>%
  mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE"))

ggplot(d_epi_freqweight_dens_nomut_sbst, 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~model) +
  coord_cartesian(xlim = c(-5, 25)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 1024") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Trait epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))


############

ggplot(d_epi_dens_nomut_sbst, 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-0.05, 0.15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 1024") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Fitness epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_freqweight_dens_nomut_sbst, 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(optPerc~model) +
  coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 1024") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Fitness epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

#####

ggplot(d_epi_dens_nomut_sbst, 
       aes(x = Pab_x, y = as.factor(log10(r)), height = Pab_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 1024") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Phenotypic effect of double mutant (Pab)", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_freqweight_dens_nomut_sbst, 
       aes(x = Pab_x, y = as.factor(log10(r)), height = Pab_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 1024") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Phenotypic effect of double mutant (Pab)", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))


################
# A different case: 4 loci, 0.0125 mut var
d_epi_freq_dens_sbst2 <- d_epi_freq_dens %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 0.0125, nloci == 4) %>%
  mutate(ep_y = as.numeric(ep_y))

# Not much difference between mutation type combinations in K- models
ggplot(d_epi_freq_dens_sbst2 %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# Most trait epistasis in K+ models comes from comparisons between two KXZ mutations
ggplot(d_epi_freq_dens_sbst2 %>% filter(model == "K"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

{
# frequency weighted sampling: is it any different?
d_epi_freq_dens_sbst2 <- d_epi_freq_dens %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 0.0125, nloci == 4) %>%
  mutate(ep_y = as.numeric(ep_y))


ggplot(d_epi_freq_dens_sbst2 %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  #coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# More or less the same as non-frequency weighted
ggplot(d_epi_freq_dens_sbst2 %>% filter(model == "K"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  #coord_cartesian(xlim = c(-5, 10)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)
}

# Average over mutType
d_epi_dens_nomut_sbst2 <- d_epi_dens_nomut %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 0.0125, nloci == 4) %>%
  mutate(ep_y = as.numeric(ep_y))


ggplot(d_epi_dens_nomut_sbst2 %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(.~optPerc) +
  #coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

d_epi_freqweight_dens_nomut_sbst2 <- d_epi_freqweight_dens_nomut %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 0.0125, nloci == 4) %>%
  mutate(ep_y = as.numeric(ep_y))

d_epi_dens_nomut_sbst2 <- d_epi_dens_nomut_sbst2 %>%
  mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE"))

ggplot(d_epi_dens_nomut_sbst2, 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~model) +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  ggtitle("Mutational variance = 0.0125; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Trait epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

d_epi_freqweight_dens_nomut_sbst2 <- d_epi_freqweight_dens_nomut_sbst2 %>%
  mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE"))

ggplot(d_epi_freqweight_dens_nomut_sbst2, 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~model) +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  ggtitle("Mutational variance = 0.0125; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Trait epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_dens_nomut_sbst2, 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-0.05, 0.15)) +
  ggtitle("Mutational variance = 0.0125; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Fitness epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_freqweight_dens_nomut_sbst2, 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 0.0125; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Fitness epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

#####

ggplot(d_epi_dens_nomut_sbst2, 
       aes(x = Pab_x, y = as.factor(log10(r)), height = Pab_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 0.0125; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Phenotypic effect of double mutant (Pab)", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_freqweight_dens_nomut_sbst2, 
       aes(x = Pab_x, y = as.factor(log10(r)), height = Pab_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 0.0125; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Phenotypic effect of double mutant (Pab)", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

{
# Do we get different distributions for a larger mutational variance even with
# a similar mean?

d_epi_dens_nomut_sbst3 <- d_epi_dens_nomut %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 1.25, nloci == 4) %>%
  mutate(ep_y = as.numeric(ep_y))


ggplot(d_epi_dens_nomut_sbst3 %>% filter(model == "ODE"), 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(.~optPerc) +
  #coord_cartesian(xlim = c(-5, 25)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

d_epi_freqweight_dens_nomut_sbst3 <- d_epi_freqweight_dens_nomut %>%
  filter(optPerc == "[0.9, Inf)" | optPerc == "[-Inf,0.25)", 
         tau == 1.25, nloci == 4) %>%
  mutate(ep_y = as.numeric(ep_y))

d_epi_dens_nomut_sbst3 <- d_epi_dens_nomut_sbst3 %>%
  mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE"))

ggplot(d_epi_dens_nomut_sbst3, 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-0.1, 0.1)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Trait epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

d_epi_freqweight_dens_nomut_sbst3 <- d_epi_freqweight_dens_nomut_sbst3 %>%
  mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE"))

ggplot(d_epi_freqweight_dens_nomut_sbst3, 
       aes(x = ep_x, y = as.factor(log10(r)), height = ep_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-0.1, 0.1)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Trait epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_dens_nomut_sbst3, 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-0.05, 0.15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Fitness epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_freqweight_dens_nomut_sbst3, 
       aes(x = ew_x, y = as.factor(log10(r)), height = ew_y)) +
  facet_grid(optPerc~model) +
  #coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Fitness epistasis", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

#####

ggplot(d_epi_dens_nomut_sbst3, 
       aes(x = Pab_x, y = as.factor(log10(r)), height = Pab_y)) +
  facet_grid(optPerc~model) +
  coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Phenotypic effect of double mutant (Pab)", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_epi_freqweight_dens_nomut_sbst3, 
       aes(x = Pab_x, y = as.factor(log10(r)), height = Pab_y)) +
  facet_grid(optPerc~model) +
  coord_cartesian(xlim = c(-15, 15)) +
  ggtitle("Mutational variance = 1.25; Number of loci = 4") +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4) +
  labs(x = "Phenotypic effect of double mutant (Pab)", y = "Recombination rate (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14))
  }

# Combine these figures to compare across model combinations
plot_grid()
