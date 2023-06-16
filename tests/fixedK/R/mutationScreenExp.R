# From each adaptive step, generate some mutations and get the distribution of 
# fitness effects
library(tidyverse)
library(future)
library(doParallel)
library(foreach)
source("wrangle_data.R")

.createEffectDataframeAdd <- function(dat, sampled_effects) {
  dat <- dat %>% filter(modelindex == 1)
  
  dat$rowID <- as.integer(rownames(dat))
  
  dat_out <- dat[rep(seq_len(nrow(dat)), each = length(sampled_effects)),]
  rownames(dat_out) <- NULL # reset row names
  
  
  Aa <- calcAddFitness(dat_out$fixEffectSum + sampled_effects, 2, 0.05)
  AA <- calcAddFitness(dat_out$fixEffectSum + sampled_effects * 2, 2, 0.05)
  aa <- calcAddFitness(dat_out$fixEffectSum, 2, 0.05)
  
  dat <- dat_out
  dat$avFit <- Aa - aa
  dat$avFit_AA <- AA - aa
  dat$value_AA <- dat$value * 2
  dat$wAA <- AA
  dat$wAa <- Aa
  dat$waa <- aa
  return(dat)
}

.createEffectDataframe <- function(dat, sampled_effects) {
  dat <- dat %>% filter(modelindex == 2)
  
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = d_popfx$fixEffectSum + 2 * value
  # Aa = d_popfx$fixEffectSum + value
  # aa = d_popfx$fixEffectSum
  
  # Calculate base phenotype/fitness with only fixations (reference)
  dat$rowID <- as.integer(rownames(dat))
  
  write.table(dat %>% ungroup() %>% mutate(KZ = 1, KXZ = 1) %>%
                dplyr::select(rowID, fixEffectSum_aZ, fixEffectSum_bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_base <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  
  # Duplicate rows so that each has the sampled additional effect
  dat_out <- dat[rep(seq_len(nrow(dat)), each = length(sampled_effects)),]
  rownames(dat_out) <- NULL # reset row names
  
  d_popfx <- d_base[rep(seq_len(nrow(d_base)), each = length(sampled_effects)*2),]
  rownames(d_popfx) <- NULL # reset row names
  
  # We need to add to aZ/bZ separately
  d_dat_withFX_aZ <- dat_out
  d_dat_withFX_bZ <- dat_out

  # Now add a random effect: calculate seperately for alpha/beta
  d_dat_withFX_aZ$aZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + sampled_effects)
  d_dat_withFX_aZ$bZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + sampled_effects)

  # Homozygous estimation - for dominance calculation
  d_dat_withFX_aZ$aZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + 2 * sampled_effects)
  d_dat_withFX_aZ$bZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + 2 * sampled_effects)
  
  d_dat_withFX <- rbind(d_dat_withFX_aZ, d_dat_withFX_bZ)
  
  # Add KZ/KXZ values
  d_dat_withFX$KZ <- 1
  d_dat_withFX$KXZ <- 1
  
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ_AA, bZ_AA, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  AA <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Get the effect size by taking away the phenotype missing that fixation
  # Ensure that the tables are aligned by id before we join them
  dat <- d_dat_withFX %>% arrange(rowID)
  Aa <- Aa %>% arrange(id)
  AA <- AA %>% arrange(id)
  
  dat$avFX <- Aa$pheno - d_popfx$pheno
  dat$avFit <- Aa$fitness - d_popfx$fitness
  dat$avFX_AA <- AA$pheno - d_popfx$pheno
  dat$avFit_AA <- AA$fitness - d_popfx$fitness
  dat$wAA <- AA$fitness
  dat$wAa <- Aa$fitness
  dat$waa <- d_popfx$fitness
  dat <- dat %>% mutate(value = ifelse((aZ - fixEffectSum_aZ) > 0, aZ - fixEffectSum_aZ, bZ - fixEffectSum_bZ))
  dat <- dat %>% select(gen, rank, seed, modelindex, value, phenomean, 
                        w, aZ, bZ, fixEffectSum_aZ, fixEffectSum_bZ,
                        avFX, avFX_AA, avFit, avFit_AA, wAA, wAa, waa)
  return(dat)
}

MutationScreenExp <- function(fixed, n, isAdditive = F) {
  # at each step in the walk, sample n mutations for each molecular component
  # and add them to the ancestor (the previous step) - then add that to a dataframe
  
  # sample effects
  fx <- rnorm(n)
  
  # Run mutation sweep
  if (isAdditive) {
    return(.createEffectDataframeAdd(fixed, fx))
  }
  
  .createEffectDataframe(fixed, fx)
}

test <- MutationScreenExp(d_fix_ranked, 100)
test <- CalcDominance(test)

test_add <- MutationScreenExp(d_fix_ranked_add, 100, T)
test_add <- CalcDominance(test_add)


# plot 
ggplot(test, aes(x = as.factor(rank), y = s)) +
  geom_boxplot() +
  labs(x = "Adaptive step", y = "Fitness effect (s)") +
  theme_bw()

ggplot(test_add %>% filter(rank > 0), aes(x = as.factor(rank), y = s)) +
  geom_boxplot() +
  labs(x = "Adaptive step", y = "Fitness effect (s)") +
  theme_bw()


test %>% 
  group_by(rank) %>%
  summarise(percBeneficial = mean(as.integer(s > 0)),
            CIperc = CI(as.integer(s > 0)),
            meanEffect = mean(s),
            CIEffect = CI(s)) -> test_sum

test %>%
  group_by(rank, seed) %>%
  summarise(percBeneficial = mean(as.integer(s > 0))) -> test_percBeneficial

test_add %>% 
  group_by(rank) %>%
  summarise(percBeneficial = mean(as.integer(s > 0)),
            CIperc = CI(as.integer(s > 0)),
            meanEffect = mean(s),
            CIEffect = CI(s)) -> test_add_sum

test_add %>%
  group_by(rank, seed) %>%
  summarise(percBeneficial = mean(as.integer(s > 0))) -> test_add_percBeneficial


test_add_sum$model <- "Additive"
test_sum$model <- "NAR"

test_add_percBeneficial$model <- "Additive"
test_percBeneficial$model <- "NAR"


test_sum_combined <- rbind(test_add_sum, test_sum)
test_perc_combined <- rbind(test_add_percBeneficial, test_percBeneficial)

ggplot(test_sum_combined, aes(x = as.factor(rank), y = percBeneficial)) +
  facet_grid(.~model) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = percBeneficial - CIperc, 
                              ymax = percBeneficial + CIperc),
                width = 0.2) +
  geom_line(group = 1) +
  labs(x = "Adaptive step", y = "Proportion of mutations where s > 0") +
  theme_bw() +
  theme(text = element_text(size = 16))
ggsave("propbeneficial.png", device = png, bg = "white")

ggplot(test_perc_combined, aes(x = as.factor(rank), y = percBeneficial)) +
  facet_grid(.~model) +
  geom_boxplot(alpha = 0.1) +
  geom_point(data = test_sum_combined, aes(x = as.factor(rank), y = percBeneficial),
            size = 2) +
  geom_line(data = test_sum_combined, group = 1,
           linetype = "dashed", linewidth = 1) +
  geom_errorbar(data = test_sum_combined, mapping = aes(ymin = percBeneficial - CIperc, 
                              ymax = percBeneficial + CIperc),
                width = 0.2) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Adaptive step", y = "Beneficial mutation\nprobability (s > 0)") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_percBeneficial
plt_percBeneficial

ggplot(test_sum_combined, aes(x = as.factor(rank), y = meanEffect, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = meanEffect - CIEffect, 
                              ymax = meanEffect + CIEffect),
                width = 0.2) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  geom_line(aes(group=model)) +
  labs(x = "Adaptive step", y = "Mean effect", colour = "Model") +
  theme_bw() + 
  theme(text = element_text(size = 12))
ggsave("meaneffectstep.png", device = png, bg = "white")
test_combined <- rbind(test, test_add)

test_combined$model <- ifelse(test_combined$modelindex == 1, "Additive", "NAR")

ggplot(test_combined, aes(x = as.factor(rank), y = s)) +
  facet_grid(.~model) +
  geom_boxplot() +
  geom_point(data = test_sum_combined, mapping = aes(y = meanEffect)) +
  geom_errorbar(data = test_sum_combined, 
                mapping = aes(y = meanEffect, ymin = meanEffect - CIEffect, 
                              ymax = meanEffect + CIEffect),
                width = 0.2) +
  geom_line(data = test_sum_combined, 
            mapping = aes(y = meanEffect), group = 1) +
  labs(x = "Adaptive step", y = "Fitness effect (s)") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_effectsizerandom

plot_grid(plt_percBeneficial, plt_effectsizerandom,
          nrow = 2, labels = "AUTO")
ggsave("effectstepdist.png", width = 8, height = 6, device = png, bg = "white")
