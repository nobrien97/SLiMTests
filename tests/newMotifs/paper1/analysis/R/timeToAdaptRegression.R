library(tidyverse)
library(paletteer)
library(Rcpp)
library(ggh4x)
library(latex2exp)
library(randomForest)
library(pROC)
library(patchwork)
library(ggalt)
library(cowplot)
library(ggbeeswarm)

source("helperFn.R")

d_btgb_Malign_rf <- readRDS("d_btgb_Malign_rf.RDS")


# Filter out r < 0.1, makes analysis simpler
d_btgb_Malign_rf_nor <- d_btgb_Malign_rf %>% filter(r == -1, isAdapted == "Adapted") %>%
  select(-r, -isAdapted)

# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215
set.seed(seed)
# Sample per group to avoid unbalanced groups
adapted_counts <- table(d_btgb_Malign_rf_nor$isAdapted)
total_counts <- sum(adapted_counts)
num_responses <- length(adapted_counts)
adapted_weights <- total_counts / (num_responses * adapted_counts)
names(adapted_weights) <- levels(d_btgb_Malign_rf_nor$isAdapted)


idx <- sample(2, nrow(d_btgb_Malign_rf_nor), replace = T, prob = c(0.7, 0.3))
train_mbeta_adapted <- d_btgb_Malign_rf_nor[idx == 1,]
test_mbeta_adapted <- d_btgb_Malign_rf_nor[idx == 2,]


d_timeToAdapt <- d_btgb_Malign_rf_nor %>%
  group_split(model)

rf_timeToAdapt <- lapply(d_timeToAdapt, function(x) 
  {
  set.seed(seed)
  idx <- sample(2, nrow(x), replace = T, prob = c(0.7, 0.3))
  train_mbeta_adapted <- x[idx == 1,]
  test_mbeta_adapted <- x[idx == 2,]
  
  rf_model <- randomForest(formula = timeToAdapt ~ .,
               data = train_mbeta_adapted,
               strata = train_mbeta_adapted$isAdapted,
               classwt = adapted_weights,
               ntree = 500,
               proximity = T,
               importance = T,
               type = "regression")
  
  
})

seed <- 18799215

rf_result_tta <- RunRandomForestPerMotifTimeToAdapt(d_btgb_Malign_rf_nor, seed)
saveRDS(rf_result_tta, "rf_result_tta.RDS")
rf_result_tta <- readRDS("rf_result_tta.RDS")

rf_result_tta[["NAR"]]

