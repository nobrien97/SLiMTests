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
library(iml)
library(emmeans)

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


# Try a linear model
lm_tta <- lm(log(timeToAdapt) ~ model * dataset + ., data = d_btgb_Malign_rf_nor[-1707,])
summary(lm_tta)
plot(lm_tta)

em_tta <- emmeans(lm_tta, ~ model * dataset, type = "response")
test(em_tta)
pwpp(em_tta, by = "model")

prop_adapted <- d_btgb_Malign_rf %>%
  group_by(model, dataset) %>%
  summarise(sum(isAdapted == "Adapted") / n())

emmip(em_tta, model ~ dataset, CIs = T) +
  theme_bw() +
  labs(x = "Model", y= TeX("Predicted adaptation time (generations)"),
       colour = "Model") +
  scale_colour_manual(values = pal) +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
ggsave("plt_pred_tta.png", device = png, bg = "white",
       width = 8, height = 6)

