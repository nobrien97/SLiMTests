# Calculates H max for a fitness landscape, a measure of ruggedness
# Refs: Munoz et al. 2015, Zhu et al. 2022, Vassilev et al. 2000

library(tidyverse)
library(DoE.wrapper)
library(deSolve)
library(mvtnorm)
library(broom)
library(latex2exp)
library(paletteer)
library(ggbeeswarm)

library(future)
library(doParallel)
library(foreach)


source("./fitnesslandscapefunctions.R")


# Nosil et al. 2020 method: based on Poursoltan and Neumann 2015
CalculateRuggedness <- function(g, model, optima, sigma, n = 10, 
                                width = 0.004, 
                                seed = sample(1:.Machine$integer.max, n)) {
  # g = genotypes (molecular components). Replicate starting points for the walk
  # w = fitnesses of the starting points
  # n = number of steps in the walk
  # seed = replicate seed for the run
  result <- data.frame(netChangeW = numeric(nrow(g)),
                       sumChangeW = numeric(nrow(g)),
                       startW = numeric(nrow(g)),
                       endW = numeric(nrow(g)),
                       numFitnessHoles = numeric(nrow(g)))
  nComps <- ncol(g)
  rollingGenotypes <- g[1:(n+1), ]
  rollingFitnesses <- numeric(n+1)
  
  for (row_index in seq_len(nrow(g))) {
    # Set the seed for each walk
    set.seed(seed)
    # Sample n steps per genotype per a normal distribution with a given width
    # Assume width is split evenly across the components
    mutations <- rmvnorm(n, sigma = diag(nComps) * ( width / nComps ))
    mutations <- rbind(rep(0.0, nComps), mutations)
    
    # cumulative sum each column to add it to rollingGenotypes
    mutations <- apply(mutations, 2, cumsum)
    rollingGenotypes <- g[rep(row_index, times = n+1),] + mutations
    for (j in seq_len(n+1)) {
      rollingFitnesses[j] <- CalcTraitAndFitness(rollingGenotypes[j,], 
                                                 model,
                                                 optima, 
                                                 sigma)
    }
    # Calculate results - add in original fitness
    # remove invalid fitnesses from bad solutions
    changeFitnesses <- rollingFitnesses[rollingFitnesses >= 0.0]
    result[row_index,]$startW <- rollingFitnesses[1]
    result[row_index,]$endW <- rollingFitnesses[n+1]
    
    result[row_index,]$netChangeW <- changeFitnesses[length(changeFitnesses)] - changeFitnesses[1]
    result[row_index,]$sumChangeW <- sum(abs(diff(changeFitnesses)))
    # Number of fitness holes: where there was no solution at a point in the walk
    result[row_index,]$numFitnessHoles <- sum(rollingFitnesses <= 0.0)
  }
  return(result)
}

CalculateRuggednessParallel <- function(g, model, optima, sigma, n = 10, nCores = availableCores(),
                                        width = 0.004,
                                        seed = sample(1:.Machine$integer.max, nrow(g))) {
  # g = genotypes (molecular components). Replicate starting points for the walk
  # w = fitnesses of the starting points
  # n = number of steps in the walk
  # seed = replicate seed for the run

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  
  df_result <- foreach (row_index = seq_len(nrow(g)), .combine = rbind) %dopar% {
    require(tidyverse)
    require(deSolve)
    require(mvtnorm)
    
    setwd("/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/R/")
    source("./fitnesslandscapefunctions.R")
    
    nComps <- ncol(g)
    rollingGenotypes <- g[1:(n+1), ]
    rollingFitnesses <- numeric(n+1)

    # Set the seed for each walk
    set.seed(seed[row_index])
    # Sample n steps per genotype per a normal distribution with a given width
    # Assume width is split evenly across the components
    mutations <- rmvnorm(n, sigma = diag(nComps) * ( width / nComps ))
    mutations <- rbind(rep(0.0, nComps), mutations)
    
    # cumulative sum each column to add it to rollingGenotypes
    mutations <- apply(mutations, 2, cumsum)
    rollingGenotypes <- exp(log(g[rep(row_index, times = n+1),]) + mutations)
    for (j in seq_len(n+1)) {
      rollingFitnesses[j] <- CalcTraitAndFitness(rollingGenotypes[j,], 
                                                 model,
                                                 optima, 
                                                 sigma)
    }
    # Calculate results - add in original fitness
    # remove invalid fitnesses from bad solutions
    changeFitnesses <- rollingFitnesses[rollingFitnesses >= 0.0]
    
    result <- data.frame(model = model,
                         startW = rollingFitnesses[1],
                         endW = rollingFitnesses[n+1],
                         netChangeW = changeFitnesses[length(changeFitnesses)] - changeFitnesses[1],
                         sumChangeW = sum(abs(diff(changeFitnesses))),
                         numFitnessHoles = sum(rollingFitnesses <= 0.0))
    return(result)
  }
  
  stopCluster(cl)
  return(df_result)
}
################################################################################
# Script
################################################################################

models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")
NUM_RUNS <- 1000
MAX_COMP_SIZE <- 3

# Generate hypercube of parameters
pars <- lhs.design(NUM_RUNS, 12, type = "random", factor.names = comps)
pars <- pars * MAX_COMP_SIZE

GGally::ggpairs(pars)
saveRDS(pars, "pars_lhc.RDS")

pars <- readRDS("pars_lhc.RDS")

# Test
model <- models[1]
# randomly sample an optimum
parsMasked <- ParsMask(pars, model)
startSolution <- SolveModel(exp(parsMasked[1,]), model)
startTraits <- GetTraitValues(startSolution, model, exp(parsMasked[1,]))
sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
opt <- CalcOptima(startTraits, sigma, 0.9)

S <- numeric(nrow(pars))
for (i in seq_len(nrow(pars))) {
  S[i] <- CalcTraitAndFitness(cbind(parsMasked[i,]), model, opt, sigma)
}

df_test <- data.frame(epsilon = seq(0, 1, length.out = 1000),
                      H = numeric(1000),
                      M = numeric(1000))

for (i in 1:nrow(df_test)) {
  psi <- GeneratePsi(parsMasked, S, df_test$epsilon[i])
  #psi <- GeneratePsiVassilev(parsMasked, S, df_test$epsilon[i])
  df_test$H[i] <- CalculateInformationContent(psi)
  df_test$M[i] <- CalculatePartialInformationContent(psi)
}

df_test <- df_test %>% 
  mutate(HM = H / M) %>%
  pivot_longer(cols = c(H, M, HM), names_to = "stat", values_to = "value")

ggplot(df_test, aes(x = log10(epsilon), y = value, colour = stat)) +
  geom_line() +
  theme_bw() +
  scale_colour_paletteer_d("nationalparkcolors::Badlands") +
  labs(x = TeX("$log_{10}(\\epsilon)$"), y = "Information curve") +
  theme(text = element_text(size = 14))
    
# Nosil method
# Generate Latin hypercube starting points
NUM_RUNS <- 10000
MAX_COMP_SIZE <- log(3)

#seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 18799215
seed <- 18799215

pars <- lhs.design(NUM_RUNS, 12, type = "random", factor.names = comps)
pars <- pars * MAX_COMP_SIZE

# # Test
# model <- models[1]
# # randomly sample an optimum
# parsMasked <- ParsMask(pars, model)
# startSolution <- SolveModel(exp(parsMasked[1,]), model)
# startTraits <- GetTraitValues(startSolution, model, exp(parsMasked[1,]))
# sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
# opt <- CalcOptima(startTraits, sigma, 0.9)
# fitnesses <- numeric(nrow(parsMasked)) 
# 
# for (i in seq_len(nrow(parsMasked))) {
#   fitnesses[i] <- CalcTraitAndFitness(cbind(parsMasked[i,]), model, opt, sigma)
# }
# 
# RugRes <- CalculateRuggedness(parsMasked, fitnesses, model, opt, sigma, 
#                               width = 0.01)
# 
# RugRes$ruggedness <- RugRes$netChangeW - RugRes$sumChangeW

# saveRDS(RugRes, "rugres.RDS")

# Plot ruggedness ratio
ggplot(RugRes %>% filter(netChangeW != 0), 
       aes(x = netChangeW, y = sumChangeW, colour = ruggedness)) +
  geom_point() +
  theme_bw() +
  scale_colour_paletteer_c("grDevices::Viridis") +
  labs(x = "Net change in fitness", y = "Total absolute change in fitness",
       colour = "Landscape Ruggedness") +
  theme(text = element_text(size = 14), legend.position = "bottom")

mean((RugRes %>% filter(netChangeW != 0))$ruggedness)

# Run across all models
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

d_ruggedness <- foreach (model_index = seq_along(models), .combine = rbind) %dopar% {
  require(tidyverse)
  require(deSolve)
  require(mvtnorm)
  
  model <- models[model_index]
  # randomly sample an optimum
  parsMasked <- ParsMask(pars, model)
  optMolComps <- as.data.frame(t(runif(ncol(parsMasked), 0, MAX_COMP_SIZE)))
  colnames(optMolComps) <- colnames(parsMasked)
  startSolution <- SolveModel(exp(optMolComps), model)
  startTraits <- GetTraitValues(startSolution, model, exp(optMolComps))
  sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
  opt <- CalcOptima(startTraits, sigma, 0.9)

  RugRes <- CalculateRuggedness(parsMasked, model, opt, sigma, 
                                n = 20, width = 0.01)
  
  RugRes$ruggedness <- RugRes$netChangeW - RugRes$sumChangeW
  RugRes$model <- model
  
  # Fill result data frame
  #output_index <- match(model, models)
  return(RugRes)
  #d_ruggedness[(NUM_RUNS * (output_index - 1) + 1):(NUM_RUNS * output_index),] <- RugRes
}

stopCluster(cl)

# Run in parallel
#seeds <- sample(1:.Machine$integer.max, NUM_RUNS)
#saveRDS(seeds, "seeds.RDS")
seeds <- readRDS("seeds.RDS")

d_ruggedness <- data.frame(
  model = character(NUM_RUNS * length(models)),
  startW = numeric(NUM_RUNS * length(models)),
  endW = numeric(NUM_RUNS * length(models)),
  netChangeW = numeric(NUM_RUNS * length(models)),
  sumChangeW = numeric(NUM_RUNS * length(models)),
  numFitnessHoles = integer(NUM_RUNS * length(models))
)

for (model in models) {
  # randomly sample an optimum
  parsMasked <- ParsMask(pars, model)
  optMolComps <- as.data.frame(t(runif(ncol(parsMasked), 0, MAX_COMP_SIZE)))
  colnames(optMolComps) <- colnames(parsMasked)
  startSolution <- SolveModel(exp(optMolComps), model)
  startTraits <- GetTraitValues(startSolution, model, exp(optMolComps))
  sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
  opt <- CalcOptima(startTraits, sigma, 0.9)
  
  RugRes <- CalculateRuggednessParallel(parsMasked, model, opt, sigma, seed = seeds)
  
  output_index <- match(model, models)
  d_ruggedness[(NUM_RUNS * (output_index - 1) + 1):(NUM_RUNS * output_index),] <- RugRes
  
}

d_ruggedness$ruggedness <- d_ruggedness$netChangeW - d_ruggedness$sumChangeW
saveRDS(d_ruggedness, "nosil_ruggedness2.RDS")

# Plot ruggedness ratio
ggplot(d_ruggedness %>% filter(netChangeW != 0 & !is.na(netChangeW)) %>%
         mutate(model = as_factor(model),
                ruggedness = abs(ruggedness)), 
       aes(x = netChangeW, y = sumChangeW, colour = ruggedness)) +
  facet_grid(.~model) +
  geom_point() +
  theme_bw() +
  scale_colour_paletteer_c("grDevices::Viridis") +
  labs(x = "Net change in fitness", y = "Total absolute change in fitness",
       colour = "Landscape Ruggedness") +
  theme(text = element_text(size = 14), legend.position = "bottom",
        panel.spacing.x = unit(1.5, "lines"))
ggsave("plt_landscaperuggedness.png")

# Plot number of fitness holes in the walk
ggplot(d_ruggedness %>% filter(numFitnessHoles > 0) %>%
         mutate(model = as_factor(model)), 
       aes(x = numFitnessHoles)) +
  facet_grid(.~model) +
  geom_histogram(bins = 11) +
  theme_bw() +
  scale_colour_paletteer_c("grDevices::Viridis") +
  labs(x = "Number of fitness holes during the random walk", y = "Count",
       colour = "Landscape Ruggedness") +
  theme(text = element_text(size = 14), legend.position = "bottom",
        panel.spacing.x = unit(1.5, "lines"))
ggsave("plt_landscapeholeyness.png")

d_ruggedness %>%
  group_by(model) %>%
  summarise(MeanFitnessHoles = mean(numFitnessHoles),
            MeanRuggedness = mean(ruggedness, na.rm = T))
