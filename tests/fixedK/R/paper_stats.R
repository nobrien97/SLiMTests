library(tidyverse)
library(fitdistrplus)
library(multimode)
library(boot)
library(moments)


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

# bootstrap difference in time to fixation
fixTimeDiff <- function(data, n) {
  data_add <- data %>% filter(modelindex == 1)
  data_net <- data %>% filter(modelindex == 2)
  samples_add <- data_add[sample(nrow(data_add), n, replace = T),]
  samples_net <- data_net[sample(nrow(data_net), n, replace = T),]

  ttf_add <- samples_add$gen - samples_add$originGen
  ttf_net <- samples_net$gen - samples_net$originGen
  
  return(ttf_net - ttf_add)
}

# Difference add vs NAR with populations not yet at the optimum
d_nonadapted <- fixGenDiff(d_fix_ranked_combined %>% 
                             filter(rank > 0, phenomean < 1.9 | phenomean > 2.1), 
                           100000)

ggplot(data.frame(diff = d_nonadapted), 
       aes(x = diff)) +
  geom_histogram() +
  theme_bw()

ggplot()

data.frame(diff = d_nonadapted) %>%
  filter(abs(diff) < 1000) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))

# Difference add vs NAR with all populations (adapted or not)
d_all <- fixGenDiff(d_fix_ranked_combined %>% filter(rank > 0), 10000)

ggplot(data.frame(diff = d_all) %>% filter(abs(diff) < 1000), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_all) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))


# Difference add vs NAR with only populations that have already adapted
d_adapteddiff <- fixGenDiff(d_fix_ranked_combined %>% 
                              filter(rank > 0, 
                                     between(phenomean, 1.9, 2.1)), 10000)

ggplot(data.frame(diff = d_adapteddiff) %>% filter(abs(diff) < 1000), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_adapteddiff) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))


# again but with time to fixation
# Difference add vs NAR with populations not yet at the optimum
d_nonadapted <- fixTimeDiff(d_fix_ranked_combined %>% 
                             filter(rank > 0, phenomean < 1.9 | phenomean > 2.1), 
                           100000)

ggplot(data.frame(diff = d_nonadapted), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_nonadapted) %>%
  filter(abs(diff) < 1000) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))

# Difference add vs NAR with all populations (adapted or not)
d_all <- fixTimeDiff(d_fix_ranked_combined %>% filter(rank > 0), 10000)

ggplot(data.frame(diff = d_all) %>% filter(abs(diff) < 1000), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_all) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))


# Difference add vs NAR with only populations that have already adapted
d_adapteddiff <- fixTimeDiff(d_fix_ranked_combined %>% 
                              filter(rank > 0, 
                                     between(phenomean, 1.9, 2.1)), 10000)

ggplot(data.frame(diff = d_adapteddiff) %>% filter(abs(diff) < 1000), 
       aes(x = diff)) +
  geom_density() +
  theme_bw()

data.frame(diff = d_adapteddiff) %>%
  summarise(meanDiff = mean(diff),
            CIDiff = CI(diff))


# linear model instead of bootstrap - not a **super** good fit, but 
# sample size is large enough it shouldn't be a problem 
mdl <- lm(gen ~ model, data = d_fix_ranked_combined %>% 
            filter(rank > 1, phenomean < 1.9 | phenomean > 2.1) %>%
      mutate(gen = gen - 50000))
plot(mdl)
summary(mdl)

# 95% CI
sqrt(diag(vcov(mdl))) * qnorm(0.975)

plot(fitted(mdl), mdl$residuals)
abline(0, 0)

step_labs <- paste0("$", levels(d_fix_ranked_combined$rankFactor), "$")

ggplot(d_fix_ranked_combined %>% 
         filter(rank > 0, phenomean < 1.9 | phenomean > 2.1) %>%
         mutate(gen = gen - 50000),
       aes(y = rankFactor, x = gen, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_y_discrete(labels = parse(text=TeX(step_labs[2:4]))) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(x = "Fixation generation (post-optimum shift)", y = "Adaptive step", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

ggsave("s_fig_fixgendist.png", device = png)

# fit gamma distribution to fixed effects
fit_nar <- fitdist((d_fix_ranked %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")
fit_add <- fitdist((d_fix_ranked_add %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")

# 95% CI
fit_nar$sd * qnorm(0.975)
fit_add$sd * qnorm(0.975)
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

# waiting time to beneficial mut, compare means
mutExp_combined %>% 
  group_by(seed, model, rankFactor) %>%
  summarise(percBeneficial = mean(as.integer(s > 0)),
            waitingTime = 1/(10000 * (9.1528*10^-6) * percBeneficial)) -> mutExp_wt

percBen_aov <- lm(percBeneficial ~ model * rankFactor, 
                   mutExp_wt %>% filter(is.finite(percBeneficial)))
summary(percBen_aov)
sqrt(diag(vcov(percBen_aov))) * qnorm(0.975)

# Can't do lm for waitingTime, residuals are too wack (strongly right-heavy tailed)
# Instead let's bootstrap 
waitingTimeDiff <- function(data, n) {
  data <- data %>% filter(is.finite(waitingTime))
  samples_add <- sample((data %>% filter(model == "Additive"))$waitingTime, n, replace = T)
  samples_net <- sample((data %>% filter(model == "NAR"))$waitingTime, n, replace = T)
  
  result <- data.frame(
    diff = samples_net - samples_add,
    sample_net = samples_net,
    sample_add = samples_add
  )
  return(result)
}

d_waitingTime <- waitingTimeDiff(mutExp_wt, 100000)

ggplot(d_waitingTime, aes(x = diff)) +
  geom_histogram(bins = 100)

d_waitingTime %>%
  summarise(meanNetWT = mean(sample_net),
            meanAddWT = mean(sample_add),
            CINetWT = CI(sample_net),
            CIAddWT = CI(sample_add),
            meanDiff = mean(diff),
            CIDiff = CI(diff)) -> res

ggplot(mutExp_wt %>% filter(is.finite(waitingTime)), 
       aes(y = rankFactor, x = percBeneficial, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_y_discrete(labels = parse(text=TeX(step_labs))) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(x = "Waiting time to beneficial mutation", y = "Adaptive step", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")


mdl_wt <- anova(aov(percBeneficial ~ model * rankFactor, 
                    mutExp_wt %>% filter(is.finite(waitingTime))))
summary(mdl_wt)
plot(mdl_wt)


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

# usage of molecular components
mean(d_molCompDiff$molCompDiff)
CI(d_molCompDiff$molCompDiff)

# percentage of models with more than 1 step that used only alpha, only beta, or both
d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(evoBybZ = all(value_aZ == 0, na.rm = T),
         evoByaZ = all(value_bZ == 0, na.rm = T)) %>%
  ungroup() %>%
  distinct(seed, .keep_all = T) %>% 
  summarise(percEvoByaZ = mean(evoByaZ),
            percEvoBybZ = mean(evoBybZ),
            percEvoByBoth = 1 - (percEvoByaZ + percEvoBybZ),
            countEvoByaZ = sum(evoByaZ),
            countEvoBybZ = sum(evoBybZ),
            countEvoByBoth = n() - (countEvoByaZ + countEvoBybZ))

# ratio of seg/fixed effects
View(d_fix_ranked_combined %>% filter(rank > 0) %>%
  mutate(rat = AA_pheno/phenomean) %>%
  group_by(model) %>%
  summarise(meanRat = mean(rat),
            CIRat = CI(rat)))
