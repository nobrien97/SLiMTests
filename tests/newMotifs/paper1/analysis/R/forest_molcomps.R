# Fits a random forest for adapted outcome vs covariance between molecular components and 
# mutational variance in each trait
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

source("helperFn.R")

# Load in data from the selection experiments

# Combos
COMBO_PATH <- '/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
COMBO_PATH <- '/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/combos.csv'
d_combos <- read_delim(COMBO_PATH, 
                       delim = " ", col_names = F)
names(d_combos) <- c("model", "r")

# Load in total QG data from other forest analysis
PATH_QG <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_qg.csv"
PATH_QG <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_qg.csv"
PATH_QG <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_qg.csv"

d_qg <- data.table::fread(PATH_QG, header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 29)), 
                          col.names = c("gen", "seed", "modelindex", "meanH",
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var", 
                                        "trait3_var", "trait4_var", "dist", 
                                        "dist1", "dist2", "dist3", "dist4", "mean_w",
                                        "var_w", "deltaPheno", "deltaW", 
                                        "meanMC1", "meanMC2", "meanMC3", "meanMC4", 
                                        "meanMC5", "meanMC6", "meanMC7", "meanMC8", 
                                        "meanMC9", "meanMC10", "meanMC11"), 
                          fill = T)

# Summarise adapted/maladapted at end of sim
d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98)) %>%
  mutate(model = factor(model, levels = model_names)) %>%
  ungroup() -> d_qg




# Load in data
PATH_QG_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"
PATH_QG_ORTH <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/orthSel/slim_qg.csv"

d_qg_orth <- data.table::fread(PATH_QG_ORTH, header = F, 
                               sep = ",", colClasses = c("integer", "factor", "factor", 
                                                         rep("numeric", times = 29)), 
                               col.names = c("gen", "seed", "modelindex", "meanH",
                                             "trait1_mean", "trait2_mean", "trait3_mean",
                                             "trait4_mean", "trait1_var", "trait2_var", 
                                             "trait3_var", "trait4_var", "dist", 
                                             "dist1", "dist2", "dist3", "dist4", "mean_w",
                                             "var_w", "deltaPheno", "deltaW", 
                                             "meanMC1", "meanMC2", "meanMC3", "meanMC4", 
                                             "meanMC5", "meanMC6", "meanMC7", "meanMC8", 
                                             "meanMC9", "meanMC10", "meanMC11"), 
                               fill = T)

# Summarise adapted/maladapted at end of sim
d_qg_orth <- AddCombosToDF(d_qg_orth) 

d_qg_orth %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98)) %>%
  mutate(model = factor(model, levels = model_names),
         modelindex = factor(modelindex, levels = 1:15)) %>%
  ungroup() -> d_qg_orth


PATH_QG_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_qg.csv"
PATH_QG_PAR <- "/mnt/j/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_qg.csv"

d_qg_par <- data.table::fread(PATH_QG_PAR, header = F, 
                              sep = ",", colClasses = c("integer", "factor", "factor", 
                                                        rep("numeric", times = 29)), 
                              col.names = c("gen", "seed", "modelindex", "meanH",
                                            "trait1_mean", "trait2_mean", "trait3_mean",
                                            "trait4_mean", "trait1_var", "trait2_var", 
                                            "trait3_var", "trait4_var", "dist", 
                                            "dist1", "dist2", "dist3", "dist4", "mean_w",
                                            "var_w", "deltaPheno", "deltaW", 
                                            "meanMC1", "meanMC2", "meanMC3", "meanMC4", 
                                            "meanMC5", "meanMC6", "meanMC7", "meanMC8", 
                                            "meanMC9", "meanMC10", "meanMC11"), 
                              fill = T)

# Summarise adapted/maladapted at end of sim
d_qg_par <- AddCombosToDF(d_qg_par) 

d_qg_par %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.98)) %>%
  mutate(model = factor(model, levels = model_names),
         modelindex = factor(modelindex, levels = 1:15)) %>%
  ungroup() -> d_qg_par

# Join qg together
d_qg$dataset <- "Randomised"
d_qg_orth$dataset <- "Orthogonal"
d_qg_par$dataset <- "Parallel"

d_qg_tot <- rbind(d_qg, d_qg_orth, d_qg_par)

# Load in mutational covariance with components
# Does the correlation structure among traits reflect correlations among molecular 
# components?
## Can look at covariance of trait M with mol comps?
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar_percomp.csv"
DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/R/slim_mutvar_percomp.csv"
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/slim_mutvar_percomp.csv"

DATA_PATH_ORTH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/orthSel/R/slim_mutvar_percomp.csv"
DATA_PATH_ORTH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/orthSel/R/slim_mutvar_percomp.csv"
DATA_PATH_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/slim_mutvar_percomp.csv"

DATA_PATH_PAR <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/parallelSel/R/slim_mutvar_percomp.csv"
DATA_PATH_PAR <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/parallelSel/R/slim_mutvar_percomp.csv"
DATA_PATH_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/slim_mutvar_percomp.csv"


t_mc_combos <- expand.grid(1:11, 1:4)

d_m_molcomp <- read_csv(DATA_PATH, col_names = c("gen", "seed", "modelindex",
                                                 paste0("cov_", t_mc_combos$Var2, "_", t_mc_combos$Var1)))

d_m_molcomp <- d_m_molcomp %>%
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex),
         dataset = "Randomised"
  )

# Split cov_ wider
d_m_molcomp_long <- d_m_molcomp %>% filter(gen == 60000) %>%
  pivot_longer(cols = starts_with("cov"),
               names_to = c("misc", "trait", "component"),
               names_sep = "_",
               values_to = "cov") %>%
  select(-misc) %>%
  mutate(dataset = "Randomised")

# orthogonal set
d_m_molcomp_orth <- data.table::fread(DATA_PATH_ORTH, header = F, sep = ",",
                               col.names = c("gen", "seed", "modelindex",
                                                 paste0("cov_", t_mc_combos$Var2, "_", t_mc_combos$Var1)),
                             colClasses = c("integer", "character", "integer",
                                            rep("numeric", times = 44)), fill = 47)

data.table::setnafill(d_m_molcomp_orth, type = "const", fill = 0, cols = 4:47)

d_m_molcomp_orth <- d_m_molcomp_orth %>%
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex),
         dataset = "Orthogonal"
  ) 

# Split cov_ wider
d_m_molcomp_orth_long <- d_m_molcomp_orth %>% filter(gen == 60000) %>%
  pivot_longer(cols = starts_with("cov"),
               names_to = c("misc", "trait", "component"),
               names_sep = "_",
               values_to = "cov") %>%
  select(-misc) %>%
  mutate(dataset = "Orthogonal")

# Parallel set
d_m_molcomp_par <- data.table::fread(DATA_PATH_PAR,  header = F, sep = ",",
                                     fill = 47,
                                     col.names = c("gen", "seed", "modelindex",
                                                 paste0("cov_", t_mc_combos$Var2, "_", t_mc_combos$Var1)),
                                     colClasses = c("integer", "character", "integer",
                                                    rep("numeric", times = 44)))

data.table::setnafill(d_m_molcomp_par, type = "const", fill = 0, cols = 4:47)

d_m_molcomp_par <- d_m_molcomp_par %>%
  mutate(model = ModelFromIndexWithR(modelindex),
         r = RFromIndex(modelindex),
         dataset = "Parallel"
)

# Split cov_ wider
d_m_molcomp_par_long <- d_m_molcomp_par %>% filter(gen == 60000) %>%
  pivot_longer(cols = starts_with("cov"),
               names_to = c("misc", "trait", "component"),
               names_sep = "_",
               values_to = "cov") %>%
  select(-misc) %>%
  mutate(dataset = "Parallel")

# bind
d_m_molcomp_long_tot <- rbind(d_m_molcomp_long, d_m_molcomp_orth_long, d_m_molcomp_par_long)
d_m_molcomp_tot <- rbind(d_m_molcomp, d_m_molcomp_orth, d_m_molcomp_par)


d_m_molcomp_long_tot <- left_join(d_m_molcomp_long_tot %>% mutate(seed = factor(seed),
                                                modelindex = factor(modelindex),
                                                r = RFromIndex(modelindex)),
                         d_qg_tot %>% select(gen, seed, modelindex, model, dataset, r, isAdapted) %>%
                           mutate(model = str_remove_all(model, "'")),
                         by = c("gen", "seed", "modelindex", "dataset", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote))


d_m_molcomp_tot <- left_join(d_m_molcomp_tot %>% mutate(seed = factor(seed),
                                                                  modelindex = factor(modelindex),
                                                                  r = RFromIndex(modelindex)),
                                  d_qg_tot %>% select(gen, seed, modelindex, model, dataset, r, isAdapted) %>%
                                    mutate(model = str_remove_all(model, "'")),
                                  by = c("gen", "seed", "modelindex", "dataset", "model", "r")) %>%
  mutate(model = factor(model, levels = model_names_noquote))



d_m_molcomp_long_sum <- d_m_molcomp_long_tot %>%
  filter(abs(cov) < 1e2) %>%
  group_by(model, isAdapted, trait, component) %>%
  summarise(meanCov = mean(abs(cov)),
            SECov = se(abs(cov)))

# Plot covariances as heatmap
ggplot(d_m_molcomp_long_sum, 
       aes(x = trait, y = component, fill = meanCov)) +
  facet_nested("Model" + model ~ "Population adapted" + isAdapted) +
  geom_tile() +
  theme_bw() +
  labs(x = "Trait", y = "Component", fill = "Covariance")


# Random forest for each motif
d_m_molcomp_rf <- d_m_molcomp_tot %>%
  filter(log10(r) == -1) %>%
  select(-r, -seed, -gen, -modelindex) %>%
  mutate(isAdapted = factor(isAdapted, levels = c("TRUE", "FALSE"), 
                            labels = c("Adapted", "Maladapted")),
         dataset = factor(dataset, levels = c("Parallel", "Orthogonal", "Randomised")))

# Model
mc_permodel <- list("NAR" = t_mc_combos[(t_mc_combos$Var1 < 8) & !(t_mc_combos$Var2 %in% c(3,4)),],
                    "PAR" = t_mc_combos[(t_mc_combos$Var1 < 8) & !(t_mc_combos$Var2 %in% c(3,4)),],
                    "FFLC1" = t_mc_combos[(t_mc_combos$Var1 < 10) & !(t_mc_combos$Var2 %in% c(4)),],
                    "FFLI1" = t_mc_combos[(t_mc_combos$Var1 < 10) & !(t_mc_combos$Var2 %in% c(4)),],
                    "FFBH" = t_mc_combos)


d_m_molcomp_sbst <- d_m_molcomp_rf %>% filter(model == "FFLC1") %>% 
  mutate(isAdapted = as.integer(isAdapted) - 1)

x <- model.matrix(as.formula(paste0("isAdapted ~ dataset +", 
                                    paste("cov_", mc_permodel[["FFLC1"]]$Var2, "_", mc_permodel[["FFLC1"]]$Var1, 
                                          sep="", collapse = "+"))),
                  d_m_molcomp_sbst
                  )
y <- d_m_molcomp_sbst$isAdapted


# Plot ROC for the different GLM models 
lm.molcomp <- cv.glmnet(x, y,
                  family = "binomial", type.measure = "auc",
                  keep = T)
lm.molcomp
plot(lm.molcomp)
lm.molcomp$lambda.min
coef(lm.molcomp, s = "lambda.min")
best <- lm.molcomp$index["min",]
rocs <- roc.glmnet(lm.molcomp$fit.preval, newy = y)
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col = "grey"))
lines(rocs[[best]], lwd = 2, col = "red")

# Confusion matrix
cnf <- confusion.glmnet(lm.molcomp, newx = x, newy = y)
print(cnf)

seed <- 123
dataset <- d_m_molcomp_rf
train.test = c(0.7, 0.3)
motifs <- levels(dataset$model)
motif <- motifs[1]
RunRandomForestMolCompPerMotif <- function(dataset, seed = NULL, train.test = c(0.7, 0.3)) {
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }
  
  motifs <- levels(dataset$model)
  result <- vector(mode = "list")
  
  for (motif in motifs) {
    set.seed(seed)
    d_rf <- dataset %>%
      filter(model == motif) %>%
      select(-model)
    
    adapted_counts <- table(d_rf$isAdapted)
    total_counts <- sum(adapted_counts)
    num_responses <- length(adapted_counts)
    adapted_weights <- total_counts / (num_responses * adapted_counts)
    names(adapted_weights) <- levels(d_rf$isAdapted)
    
    
    idx <- sample(2, nrow(d_rf), replace = T, prob = train.test)
    train <- d_rf[idx == 1,]
    test <- d_rf[idx == 2,]
    
    
    # no balancing
    rf_nobal <- randomForest(formula = isAdapted ~ .,
                                 data = train,
                                 ntree = 500,
                                 proximity = T,
                                 importance = T,
                                 type = "classification")
    
    print(rf_nobal)
    
    # With balancing (class weights)
    rf_bal <- randomForest(formula = isAdapted ~ .,
                               data = train,
                               strata = train$isAdapted,
                               classwt = adapted_weights,
                               ntree = 500,
                               proximity = T,
                               importance = T,
                               type = "classification")
    
    print(rf_bal)
    
    # Training data
    p_train_bal <- predict(rf_bal, train)
    caret::confusionMatrix(p_train_bal, train$isAdapted)
    
    p_train_nobal <- predict(rf_nobal, train)
    caret::confusionMatrix(p_train_nobal, train$isAdapted)
    
    
    # Test data
    p_test_bal <- predict(rf_bal, test)
    p_test_bal_probs <- predict(rf_bal, test,
                                            type = "prob")[,1]
    
    
    
    result[[motif]][["cMat"]] <- caret::confusionMatrix(p_test_bal, test$isAdapted)
    
    p_test_nobal <- predict(rf_nobal, test)
    p_test_nobal_probs <- predict(rf_nobal, test,
                                              type = "prob")[,1]
    
    caret::confusionMatrix(p_test_nobal, test$isAdapted)
    
    
    # roc
    d_roc_bal <- roc(response = test$isAdapted,
                         predictor = p_test_bal_probs)
    
    d_roc_nobal <- roc(response = test$isAdapted,
                           predictor = p_test_nobal_probs,
                           levels = rev(levels(test$isAdapted)))
    
    d_rocs <- data.frame(model = c(rep("Weighted", times = length(d_roc_bal$sensitivities)),
                                       rep("Unbalanced", times = length(d_roc_nobal$sensitivities))),
                             sens = c(d_roc_bal$sensitivities, d_roc_nobal$sensitivities),
                             spec = c(d_roc_bal$specificities, d_roc_nobal$specificities))
    
    
    roc_aucs <- c(pROC::auc(d_roc_nobal), pROC::auc(d_roc_bal))
    
    
    ggplot(d_rocs,
           aes(x = 1 - spec, y = sens, colour = model)) +
      geom_line() + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
      annotate("text", x = c(0.75, 0.75), y = c(0.375, 0.25), 
               label = paste("AUC:", round(roc_aucs, digits = 3)),
               colour = pal[1:2]) +
      theme_bw() +
      scale_colour_manual(values = pal) +
      ggtitle(motif) +
      labs(x = "1 - Specificity", y = "Sensitivity", colour = "RF Model") +
      theme(legend.position = "bottom",
            text = element_text(size = 12))
    ggsave(paste0("plt_RF_ROC_", motif, ".png"), 
           device = png, width = 4, height = 4, bg = "white")
    
    
    # Importance measures
    ## Boruta, permutation importance, sobol MDA
    bor <- Boruta::Boruta(isAdapted ~ ., data = d_rf)
    bor
    plot(bor)
    
    d_bor <- process_the_Boruta_data(bor)

    result[[motif]][["bor"]] <- d_bor
    
    # Permutation
    predictor <- iml::Predictor$new(rf_bal, 
                                        data = test[, 1:4], 
                                        y = test$isAdapted,
                                        type = "prob")
    
    # Need to set the option future globals maxsize
    options(future.globals.maxSize = 3221225472)
    imp <- iml::FeatureImp$new(predictor,
                                   loss = "ce",
                                   n.repetitions = 100)
    
    result[[motif]][["FeatImp"]] <- imp
    
      ggplot(imp$results,
           aes(x = feature, y = importance)) +
      geom_point() +
      geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                    width = 0.2) +
      scale_x_discrete(
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Permutation Importance") +
      theme_bw() +
      theme(text = element_text(size = 12)) -> plt_perm_imp
    plt_perm_imp
    ggsave(paste0("plt_perm_feat_imp_align_", motif,".png"), 
           device = png, width = 9, height = 5, bg = "white")
    
    # Sobol MDA
    rf_sob <- sobolMDA::ranger(isAdapted ~ .,
                                             data = train, num.trees = 500, 
                                             importance = "sobolMDA")
    sob <- rf_sob$variable.importance
    d_sob <- data.frame(feature = names(sob),
                                      sobelMDA = sob)
    
    d_sob$feature <- factor(d_sob$feature)
    
    ggplot(d_sob,
           aes(x = feature, y = sobelMDA)) +
      geom_point() +
      geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
                   linewidth = 0.5) +
      theme_bw() +
      scale_x_discrete(
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Sobel MDA") +
      theme(text = element_text(size = 12)) -> plt_sob
    plt_sob
    
    
    layout <- "
    AAAA
    AAAA
    AAAA
    BBCC
    BBCC
    "
    # plt_featimp <- plt_boruta_imp +
    #   plt_perm_imp +
    #   plt_sob +
    #   plot_layout(design = layout) +
    #   plot_annotation(tag_levels = 'A') &
    #   theme(plot.tag = element_text(face = "bold"))
    
    result[[motif]][["sob"]] <- d_sob

    # ggsave(paste0("plt_feat_imp_align", motif, ".png"), 
    #        device = png, width = 12, height = 10, bg = "white")
    
    
    # Accumulated local effects
    ale <- FeatureEffects$new(predictor, grid.size = 10)
    ale$plot()
    
    result[[motif]][["ale"]] <- ale
  }
  
  return(result)
}


seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 538108254
seed <- 538108254

rf_mc <- RunRandomForestMolCompPerMotif(d_m_molcomp_tot, seed)


d_m_molcomp_NAR_rf <- d_m_molcomp_tot %>%
  filter(model == "NAR", log10(r) == -1) %>%
  select(-r)

seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215
set.seed(seed)
# Sample per group to avoid unbalanced groups
adapted_counts <- table(d_m_molcomp_NAR_rf$isAdapted)
total_counts <- sum(adapted_counts)
num_responses <- length(adapted_counts)
adapted_weights <- total_counts / (num_responses * adapted_counts)
names(adapted_weights) <- levels(d_m_molcomp_NAR_rf$isAdapted)


idx <- sample(2, nrow(d_m_molcomp_NAR_rf), replace = T, prob = c(0.7, 0.3))
train_nar <- d_m_molcomp_NAR_rf[idx == 1,]
test_nar <- d_m_molcomp_NAR_rf[idx == 2,]

# no balancing
rf_nar_nobal <- randomForest(formula = isAdapted ~ .,
                                       data = train_nar,
                                       ntree = 500,
                                       proximity = T,
                                       importance = T,
                                       type = "classification")

print(rf_nar_nobal)

# With balancing (class weights)
rf_nar_bal <- randomForest(formula = isAdapted ~ .,
                                     data = train_nar,
                                     strata = train_nar$isAdapted,
                                     classwt = adapted_weights,
                                     ntree = 500,
                                     proximity = T,
                                     importance = T,
                                     type = "classification")

print(rf_nar_bal)

# Training data
p_train_adapted_nar_bal <- predict(rf_nar_bal, train_nar)
caret::confusionMatrix(p_train_adapted_nar_bal, train_nar$isAdapted)

p_train_adapted_nar_nobal <- predict(rf_nar_nobal, train_nar)
caret::confusionMatrix(p_train_adapted_nar_nobal, train_nar$isAdapted)


# Test data
p_test_adapted_nar_bal <- predict(rf_nar_bal, test_nar)
p_test_adapted_nar_bal_probs <- predict(rf_nar_bal, test_nar,
                                          type = "prob")[,1]

caret::confusionMatrix(p_test_adapted_nar_bal, test_nar$isAdapted)

p_test_adapted_nar_nobal <- predict(rf_nar_nobal, test_nar)
p_test_adapted_nar_nobal_probs <- predict(rf_nar_nobal, test_nar,
                                            type = "prob")[,1]

caret::confusionMatrix(p_test_adapted_nar_nobal, test_nar$isAdapted)


# roc
d_roc_nar_bal <- roc(response = test_nar$isAdapted,
                 predictor = p_test_adapted_nar_bal_probs)

d_roc_nar_nobal <- roc(response = test_nar$isAdapted,
                   predictor = p_test_adapted_nar_nobal_probs,
                   levels = rev(levels(test_nar$isAdapted)))

d_rocs_nar <- data.frame(model = c(rep("Weighted", times = length(d_roc_nar_bal$sensitivities)),
                               rep("Unbalanced", times = length(d_roc_nar_nobal$sensitivities))),
                     sens = c(d_roc_nar_bal$sensitivities, d_roc_nar_nobal$sensitivities),
                     spec = c(d_roc_nar_bal$specificities, d_roc_nar_nobal$specificities))


roc_aucs_nar <- c(pROC::auc(d_roc_nar_nobal), pROC::auc(d_roc_nar_bal))


ggplot(d_rocs_nar,
       aes(x = 1 - spec, y = sens, colour = model)) +
  geom_line() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  annotate("text", x = c(0.75, 0.75), y = c(0.375, 0.25), 
           label = paste("AUC:", round(roc_aucs, digits = 3)),
           colour = pal[1:2]) +
  theme_bw() +
  scale_colour_manual(values = pal) +
  ggtitle("NAR") +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = "RF Model") +
  theme(legend.position = "bottom",
        text = element_text(size = 12))
ggsave("plt_RF_ROC_NAR.png", device = png, width = 4, height = 4, bg = "white")


# Importance measures
## Boruta, permutation importance, sobol MDA
bor_nar <- Boruta::Boruta(isAdapted ~ ., data = d_m_molcomp_NAR_rf)
bor_nar
plot(bor_nar)

d_bor_nar <- process_the_Boruta_data(bor_nar)

feature_names <- c(TeX("Shadow Variable (min)", output = "character"),
                   TeX("Shadow Variable (mean)", output = "character"),
                   TeX("Shadow Variable (max)", output = "character"),
                   TeX("$G/\\beta$ alignment",
                       output = "character"),
                   TeX("G Evolvability", output = "character"),
                   TeX("$V_{rel}$ (G)", output = "character"),
                   TeX("$R_{max} / \\beta$ alignment", output = "character"),
                   TeX("$M/\\beta$ alignment",
                       output = "character"),
                   TeX("$V_{rel}$ (M)", output = "character"),
                   TeX("M Evolvability", output = "character"),
                   TeX("Motif", output = "character"))


pal_boruta <- c(rep("#00C0EA", times = 3),
                rep("#00A000", times = 8))

ggplot(d_bor_nar %>% pivot_longer(everything()) %>%
         mutate(x = fct_reorder(name, value, median)),
       aes(x = x, y = value, fill = x)) +
  geom_boxplot(show.legend = F, linewidth = 0.25) +
  theme_bw() +
  scale_fill_manual(values = pal_boruta) +
  scale_x_discrete(labels = parse(text = feature_names),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Boruta Importance") +
  theme(text = element_text(size = 12)) -> plt_boruta_imp_nar
plt_boruta_imp_nar
ggsave("plt_boruta_import_nar.png", device = png, bg = "white",
       width = 12, height = 8)


# Permutation
predictor_nar <- iml::Predictor$new(rf_nar_bal, 
                                data = test_nar[, 2:9], 
                                y = test_nar$isAdapted,
                                type = "prob")

# Need to set the option future globals maxsize
options(future.globals.maxSize = 3221225472)
imp_nar <- iml::FeatureImp$new(predictor_nar,
                           loss = "ce",
                           n.repetitions = 100)

ggplot(imp_nar$results %>% 
         mutate(feature = factor(feature, levels = c("absCS_Gb", "bTGb", "vrel_g", "dataset",  
                                                     "absCS_Mb", "vrel_m", "bTMb", "model"))),
       aes(x = feature, y = importance)) +
  geom_point() +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                width = 0.2) +
  scale_x_discrete(labels = parse(text = feature_names[4:11]),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Permutation Importance") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_perm_imp_nar
plt_perm_imp_nar
ggsave("plt_perm_feat_imp_align.png", device = png, width = 9, height = 5, bg = "white")

# Sobol MDA
rf_sob_mbeta_adapted <- sobolMDA::ranger(isAdapted ~ .,
                                         data = train_mbeta_adapted, num.trees = 500, 
                                         importance = "sobolMDA")
sob_mbeta_adapted <- rf_sob_mbeta_adapted$variable.importance
d_sob_mbeta_adapted <- data.frame(feature = names(sob_mbeta_adapted),
                                  sobelMDA = sob_mbeta_adapted)

d_sob_mbeta_adapted$feature <- factor(d_sob_mbeta_adapted$feature,
                                      levels = c("absCS_Gb",
                                                 "bTGb",
                                                 "vrel_g",
                                                 "dataset",
                                                 "absCS_Mb",
                                                 "vrel_m",
                                                 "bTMb",
                                                 "model"))

ggplot(d_sob_mbeta_adapted,
       aes(x = feature, y = sobelMDA)) +
  geom_point() +
  geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
               linewidth = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = parse(text = feature_names[4:11]),
                   guide = guide_axis(n.dodge = 2)) +
  labs(x = "Feature", y = "Sobel MDA") +
  theme(text = element_text(size = 12)) -> plt_sob_mbeta_adapted
plt_sob_mbeta_adapted


layout <- "
AAAA
AAAA
AAAA
BBCC
BBCC
"
plt_boruta_imp +
  plt_perm_imp +
  plt_sob_mbeta_adapted +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

ggsave("plt_feat_imp_align.png", device = png, width = 12, height = 10, bg = "white")


# Accumulated local effects
ale <- FeatureEffects$new(predictor, grid.size = 10)
ale$plot()

ale_plots <- vector(mode = "list", length = 8)

ale_labels <- feature_names[c(11, 7, 8, 4, 5, 10, 6, 9)]
names(ale_labels) <- names(ale$results)

for (i in seq_along(ale_plots)) {
  d_ale <- ale$results[[i]]
  d_ale <- d_ale %>% filter(.class == "Adapted")
  x_label <- ale_labels[d_ale$.feature[1]] 
  
  if (d_ale$.feature[1] == "model") {
    d_ale$.borders <- factor(d_ale$.borders, 
                             levels = model_names_noquote,
                             labels = model_names_noquote)
  }
  
  if (!is.numeric(d_ale$.borders[1])) {
    geom_fn <- geom_lollipop
    scale_fn <- scale_x_discrete
  } else {
    geom_fn <- geom_line
    scale_fn <- scale_x_continuous
    # Remove outliers
    d_ale <- d_ale %>% filter(.borders < 1000)
  }
  ale_plots[[i]] <- ggplot(d_ale,
                           aes(x = .borders, y = .value)) +
    geom_fn() +
    scale_fn() +
    theme_bw() +
    labs(x = parse(text = x_label), y = "ALE of adaptation probability")
}

plot_grid(plotlist = ale_plots,
          labels= "AUTO")
ggsave("plt_ale_align.png", device = png, bg = "white",
       width = 12, height = 9)
