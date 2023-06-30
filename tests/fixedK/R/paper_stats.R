library(tidyverse)
library(fitdistrplus)

setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("wrangle_data.R")

# Fisher test - number of pops adapted
fisher.test(table((d_qg %>% distinct(seed, modelindex, .keep_all = T))$modelindex, 
                  (d_qg %>% distinct(seed, modelindex, .keep_all = T))$isAdapted))
# data:  table((d_qg %>% distinct(seed, modelindex, .keep_all = T))$modelindex, (d_qg %>% distinct(seed, modelindex, .keep_all = T))$isAdapted)
# p-value = 6.279e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5904316 0.7972828
# sample estimates:
#   odds ratio 
# 0.6862194 

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
# p-value = 0.3891
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
  
# fit gamma distribution to fixed effects
fit_nar <- fitdist((d_fix_ranked %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")
fit_add <- fitdist((d_fix_ranked_add %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")
summary(fit_nar)
summary(fit_add)
plot(fit_nar)
plot(fit_add)

# Compare means - proportion of beneficial muts, waiting time to beneficial mut
aov(s ~ model * rank, mutExp_combined)


# alpha/beta optimum ratio
(plt_aZbZratio$data$pheno)