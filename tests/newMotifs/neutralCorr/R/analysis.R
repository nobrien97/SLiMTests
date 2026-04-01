library(tidyverse)
library(ggh4x)
library(paletteer)
library(emmeans)

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
    summarise(mean = mean(cor),
              var = var(cor)) %>%
    ungroup()

ggplot(d_qg_traitCor_sum %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
    aes(x = as.factor(traitCombo), y = mean, colour = model)) +
    geom_point(position = position_dodge(0.9)) +
    geom_errorbar(aes(ymin = mean - var, ymax = mean + var), width = 0.2,
                  position = position_dodge(0.9)) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1)) +
    labs(x = "Trait combination", y = "Mean correlation", colour = "Model") +
    theme_bw() +
    theme(text = element_text(size = 12),
          legend.position = "bottom")
ggsave("plt_trait_cor.png", device = png, bg = "white", width = 7, height = 8)

# trait correlations between ~0 and ~0.333

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
  labs(x = "Trait", y = "Mean of means") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
