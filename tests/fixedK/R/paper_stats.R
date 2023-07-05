library(tidyverse)
library(fitdistrplus)
library(multimode)
library(boot)

setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("wrangle_data.R")
source("mutationScreenExp.R")

# Fisher test - number of pops adapted
fisher.test(table((d_qg %>% distinct(seed, modelindex, .keep_all = T))$modelindex, 
                  (d_qg %>% distinct(seed, modelindex, .keep_all = T))$isAdapted))
# Fisher's Exact Test for Count Data
# 
# data:  table((d_qg %>% distinct(seed, modelindex, .keep_all = T))$modelindex, (d_qg %>% distinct(seed, modelindex, .keep_all = T))$isAdapted)
# p-value = 4.987e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.6125830 0.7565588
# sample estimates:
# odds ratio 
#  0.6808441 

# Fisher test - number of steps taken in walk
d_fix_ranked_combined %>% 
  group_by(model, rank) %>%
  filter(rank != 0) %>%
  summarise(n = n())

fisher.test(table((d_fix_ranked_combined %>% filter(rank > 0))$model, 
      (d_fix_ranked_combined %>% filter(rank > 0))$rank))
# Fisher's Exact Test for Count Data
# 
# data:  table((d_fix_ranked_combined %>% filter(rank > 0))$model, (d_fix_ranked_combined %>% filter(rank > 0))$rank)
# p-value = 0.5898
# alternative hypothesis: two.sided

# mean number of steps
d_fix_ranked_combined %>% 
  group_by(model, seed) %>%
  summarise(nSteps = max(rank)) %>%
  ungroup() %>%
  summarise(meanSteps = mean(nSteps),
            CISteps = CI(nSteps))

d_fix_ranked_combined %>% 
  group_by(model, seed) %>%
  summarise(nSteps = max(rank)) %>%
  group_by(nSteps) %>%
  summarise(n_nSteps = n(),
            nRows = nrow(.),
            perc = n_nSteps/nRows)
  
# timing of steps among populations that haven't made it yet
d_fix_ranked_combined %>% filter(rank > 0, phenomean < 1.9) %>%
  mutate(gen = gen - 50000) %>%
  group_by(model, rank) %>%
  summarise(meanGen = mean(gen), CIGen = CI(gen)) -> d_meanGenTiming

# bootstrap difference between fixation gen sampled from models
fixGenDiff <- function(data, n) {
  samples_add <- sample((data %>% filter(modelindex == 1))$gen, n, replace = T)
  samples_net <- sample((data %>% filter(modelindex == 2))$gen, n, replace = T)
  
  return(samples_net - samples_add)
}

# Difference with non adapted populations
d_nonadapted <- fixGenDiff(d_fix_ranked_combined %>% 
                             filter(rank > 0, phenomean < 1.9), 100000)

ggplot(data.frame(diff = d_nonadapted), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_nonadapted) %>%
  filter(abs(diff) < 1000) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))

# Difference with all populations
d_all <- fixGenDiff(d_fix_ranked_combined %>% filter(rank > 0), 10000)

ggplot(data.frame(diff = d_all) %>% filter(abs(diff) < 1000), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_all) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))


# Difference with adapted populations
d_adapteddiff <- fixGenDiff(d_fix_ranked_combined %>% filter(rank > 0, phenomean > 1.9, 
                                            phenomean < 2.1), 10000)

ggplot(data.frame(diff = d_adapteddiff) %>% filter(abs(diff) < 1000), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_adapteddiff) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))



mdl <- lm(gen ~ model + as.factor(rank), data = d_fix_ranked_combined %>% filter(rank > 0, phenomean < 1.9) %>%
      mutate(gen = gen - 50000))
plot(mdl)
summary(mdl)

# fit gamma distribution to fixed effects
fit_nar <- fitdist((d_fix_ranked %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")
fit_add <- fitdist((d_fix_ranked_add %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")
summary(fit_nar)
summary(fit_add)
plot(fit_nar)
plot(fit_add)

# Look at distribution of fitness effects of fixed deleterious alleles in the context
# of the burn-in environment


ggplot(d_fix_del %>% filter(rank > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") 

d_fix_del %>%
  group_by(model) %>%
  mutate(isBen = (s > 0)) %>%
  summarise(n = n(),
            nBen = sum(isBen),
            propBen = nBen/n)

# Fit multimodal distribution to beneficial effects

ggplot(mutExp, aes(s)) + 
  geom_density() + 
  xlim((min(mutExp$s)-1),(max(mutExp$s)+1) )


# Compare means - proportion of beneficial muts, waiting time to beneficial mut
aov(s ~ model * rank, mutExp_combined)


# alpha/beta optimum ratio
(plt_aZbZratio$data$pheno)




# is there a difference between the number of populations adapting with only fixations
# in additive vs network models?
d_fix_ranked_combined %>% filter(rank > 0) %>%
  group_by(model, seed) %>% filter(rank == max(rank)) %>%
  ungroup() %>%
  mutate(FixedOnly = (AA_pheno/phenomean > 0.999 & 
                        AA_pheno/phenomean < 1.001)) -> d_fixedOnlyTable

d_fixedOnlyTable %>%
  group_by(model) %>%
  summarise(nFixedOnly = sum(FixedOnly),
            percFixedOnly = nFixedOnly/n())

fisher.test(table(d_fixedOnlyTable$model, 
                   d_fixedOnlyTable$FixedOnly))
