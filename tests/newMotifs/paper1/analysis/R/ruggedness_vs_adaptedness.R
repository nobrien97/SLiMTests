library(tidyverse)
library(paletteer)
library(ggbeeswarm)
library(xtable)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

model_levels <- c("NAR", "PAR", "FFLC1", 
                  "FFLI1", "FFBH")
model_labels <- c("NAR", "PAR", "FFL-C1", "FFL-I1", "FFBH")

DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/ruggedness/R/"

RUG_DATA_PATH <- "/mnt/i/SLiMTests/tests/newMotifs/ruggedness/"


# Load data
RugRes <- data.table::fread(paste0(RUG_DATA_PATH, "d_ruggedness_permolcomp.csv"), header = F,
                            col.names = c("model", "startW", "endW", "netChangeW",
                                          "sumChangeW", "numFitnessHoles", "molComp",
                                          "bkg"))
RugRes$ruggedness <- RugRes$netChangeW - RugRes$sumChangeW



# Average over genetic backgrounds
d_ruggedness <- RugRes %>%
  mutate(model = factor(model, levels = model_levels)) %>%
  group_by(model, molComp, bkg) %>%
  mutate(replicate = row_number()) %>% # create replicate row
  ungroup() %>%
  # Filter out irrelevant rows
  filter( ( model == "NAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) |
            ( model == "PAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) | 
            ( model == "FFLC1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                                "KXZ", "KY", "Hilln", "base", "XMult") ) |
            ( model == "FFLI1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                                "KXZ", "KY", "Hilln", "base", "XMult") ) |
            model == "FFBH" )

numSteps <- 10
d_ruggedness <- d_ruggedness %>%
  group_by(model, bkg, replicate) %>%
  mutate(nHolesAcrossMolComps = sum(numFitnessHoles)) %>%
  ungroup() %>%
  mutate(propHoles = numFitnessHoles / nHolesAcrossMolComps)


# Per model (average over molcomps and over genetic backgrounds)
d_ruggedness_sum <- d_ruggedness %>%
  group_by(model) %>% # replicate: walk direction is the same, averaged across backgrounds
  summarise(meanRuggedness = mean(ruggedness),
            CIRuggedness = CI(ruggedness),
            meanFitnessHoles = mean(numFitnessHoles),
            CIFitnessHoles = CI(numFitnessHoles))

# Load in QG data
QG_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/analysis/R/"
d_btgb_Malign_rf <- readRDS(paste0(QG_PATH, "d_btgb_Malign_rf.RDS"))


# Fit RF model - use median ruggedness/fitness values
d_btgb_Malign_rf <- inner_join(d_btgb_Malign_rf, d_ruggedness_sum %>% 
                                 select(-CIRuggedness, -CIFitnessHoles), 
                               by = c("model"))


seed <- 18799215
set.seed(seed)
# Sample per group to avoid unbalanced groups
adapted_counts <- table(d_btgb_Malign_rf$isAdapted)
total_counts <- sum(adapted_counts)
num_responses <- length(adapted_counts)
adapted_weights <- total_counts / (num_responses * adapted_counts)
names(adapted_weights) <- levels(d_btgb_Malign_rf$isAdapted)


idx <- sample(2, nrow(d_btgb_Malign_rf), replace = T, prob = c(0.7, 0.3))
train_mbeta_adapted <- d_btgb_Malign_rf[idx == 1,]
test_mbeta_adapted <- d_btgb_Malign_rf[idx == 2,]

# adjust for class imbalance
rf_mbeta_adapted_bal <- randomForest(formula = isAdapted ~ model * dataset * absCS_Gb * absCS_Mb * 
                                       bTGb * bTMb * vrel_g * vrel_m * meanRuggedness * meanFitnessHoles,
                                     data = train_mbeta_adapted,
                                     strata = train_mbeta_adapted$isAdapted,
                                     classwt = adapted_weights,
                                     ntree = 500,
                                     proximity = T,
                                     importance = T,
                                     type = "classification")

p_train_mbeta_adapted_bal <- predict(rf_mbeta_adapted_bal, train_mbeta_adapted)
caret::confusionMatrix(p_train_mbeta_adapted_bal, train_mbeta_adapted$isAdapted)

p_test_mbeta_adapted_bal <- predict(rf_mbeta_adapted_bal, test_mbeta_adapted)
p_test_mbeta_adapted_bal_probs <- predict(rf_mbeta_adapted_bal, test_mbeta_adapted,
                                          type = "prob")[,1]
caret::confusionMatrix(p_test_mbeta_adapted_bal, test_mbeta_adapted$isAdapted)

d_roc_bal <- roc(response = test_mbeta_adapted$isAdapted,
                 predictor = p_test_mbeta_adapted_bal_probs)

roc_aucs <- pROC::auc(d_roc_bal)
d_rocs <- data.frame(model = c(rep("Weighted", times = length(d_roc_bal$sensitivities))),
                     sens = c(d_roc_bal$sensitivities),
                     spec = c(d_roc_bal$specificities))


ggplot(d_rocs,
       aes(x = 1 - spec, y = sens, colour = model)) +
  geom_line() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  annotate("text", x = c(0.75, 0.75), y = c(0.375, 0.25), 
           label = paste("AUC:", round(roc_aucs, digits = 3)),
           colour = pal[1:2]) +
  theme_bw() +
  scale_colour_manual(values = pal) +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = "RF Model") +
  theme(legend.position = "bottom",
        text = element_text(size = 12))

plot(rf_mbeta_adapted_bal)

# Boruta importance of ruggedness and holeyness
d_btgb_Malign_rf_model <- d_btgb_Malign_rf %>%
  group_split(model)

bor_rug <- lapply(d_btgb_Malign_rf_model, function(x) {
  Boruta::Boruta(isAdapted ~ ., data = d_btgb_Malign_rf %>% select(-model))
})

pltlst_bor_rug <- lapply(bor_rug, function(x) {
  plot(x)
})

plot_grid(plotlist = pltlst_bor_rug)


summary(lm(meanRuggedness ~ model, d_btgb_Malign_rf))

cor(d_btgb_Malign_rf$meanRuggedness, d_btgb_Malign_rf$model)


seed <- 538108254

rf_mc_rug <- RunRandomForestPerMotif(d_btgb_Malign_rf, seed)
saveRDS(rf_mc, "rf_mc_rug.RDS")
rf_result <- readRDS("rf_result.RDS")
