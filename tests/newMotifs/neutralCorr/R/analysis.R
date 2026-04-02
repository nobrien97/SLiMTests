library(tidyverse)
library(ggh4x)
library(paletteer)
library(emmeans)
library(patchwork)

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


trait_comp_names_nar <- c(
  "Response time vs steady state" # 12
)

trait_comp_names_c1 <- c(
  "Response time vs response delay", # 12
  "Response time vs steady state", # 13
  "Response delay vs steady state" # 23
)

trait_comp_names_i1 <- c(
  "Time to half max expression vs max expression", # 12
  "Time to half max expression vs time above half max expression", # 13
  "Max expression vs time above half max expression" # 23
)

trait_comp_names_bh <- c(
  "Time to max expression vs max expression", # 12
  "Time to max expression vs time to final steady state", # 13
  "Time to max expression vs final steady state", # 14
  "Max expression vs time to final steady state", # 23
  "Max expression vs final steady state", # 24
  "Time to final steady state vs final steady state" # 34
)

trait_comp_names <- data.frame(model = c(rep("NAR", times = length(trait_comp_names_nar)),
                                         rep("PAR", times = length(trait_comp_names_nar)),
                                         rep("FFLC1", times = length(trait_comp_names_c1)),
                                         rep("FFLI1", times = length(trait_comp_names_i1)),
                                         rep("FFBH", times = length(trait_comp_names_bh))),
                               label = c(trait_comp_names_nar, trait_comp_names_nar,
                                         trait_comp_names_c1, trait_comp_names_i1,
                                         trait_comp_names_bh))

traitCorBoilerplate <- function(plt) {
  plt + 
    facet_nested(. ~ model) +
    geom_point(position = position_dodge(0.9), show.legend = F) +
    geom_errorbar(aes(ymin = mean - var, ymax = mean + var), width = 0.2,
                  position = position_dodge(0.9), show.legend = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    labs(x = "Trait combination", y = "Mean correlation", colour = "Model") +
    theme_bw() +
    theme(text = element_text(size = 12),
          legend.position = "bottom")
}


ggplot(d_qg_traitCor_sum %>%
         filter(model == "NAR") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
    aes(x = as.factor(traitCombo), y = mean, colour = model)) +
    scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "NAR",2],
                     guide = guide_axis(n.dodge = 2)) +
    scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                             5, direction = -1)[4]) -> plt_traitcor_nar
plt_traitcor_nar <- traitCorBoilerplate(plt_traitcor_nar)

ggplot(d_qg_traitCor_sum %>%
         filter(model == "PAR") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "PAR",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[5]) -> plt_traitcor_par
plt_traitcor_par <- traitCorBoilerplate(plt_traitcor_par)

ggplot(d_qg_traitCor_sum %>%
         filter(model == "FFLC1") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "FFLC1",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[2]) -> plt_traitcor_fflc1
plt_traitcor_fflc1 <- traitCorBoilerplate(plt_traitcor_fflc1)

ggplot(d_qg_traitCor_sum %>%
         filter(model == "FFLI1") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "FFLI1",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[3]) -> plt_traitcor_ffli1
plt_traitcor_ffli1 <- traitCorBoilerplate(plt_traitcor_ffli1)


ggplot(d_qg_traitCor_sum %>%
         filter(model == "FFBH") %>%
         filter(!(model != "FFBH" & as.numeric(traitCombo) %% 10 > 3), # illegal comparisons
                !(nchar(model) == 3 & as.numeric(traitCombo) %% 10 > 2)),
       aes(x = as.factor(traitCombo), y = mean, colour = model)) +
  scale_x_discrete(labels = trait_comp_names[trait_comp_names$model == "FFBH",2],
                   guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = -1)[1]) -> plt_traitcor_ffbh
plt_traitcor_ffbh <- traitCorBoilerplate(plt_traitcor_ffbh)


layout <-
"
AABB
CCCC
DDDD
EEEE
"

plt_traitcor <- plt_traitcor_nar +
plt_traitcor_par +
plt_traitcor_fflc1 +
plt_traitcor_ffli1 + 
plt_traitcor_ffbh

plt_traitcor <- plt_traitcor + plot_layout(design = layout)
plt_traitcor

ggsave("plt_trait_cor.png", device = png, bg = "white", width = 12, height = 8)

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
  labs(x = "Trait", y = "Mean of mean trait values") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
