# Calculates H max for a fitness landscape, a measure of ruggedness
# Refs: Munoz et al. 2015, Zhu et al. 2022, Vassilev et al. 2000

library(tidyverse)
library(DoE.wrapper)
library(deSolve)
library(mvtnorm)
library(broom)
library(latex2exp)
library(paletteer)



################################################################################
# Function definitions
################################################################################
# clamp function
clamp <- function(x, lower=-Inf, upper = Inf) {
  pmax(pmin(x, upper), lower)
}

# ODE system for feedback autoregulation:
ODEs_NAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dZ <- base + bZ * (t > Xstart && t <= Xstop) * ((X^Hilln)/(KXZ^Hilln + X^Hilln)) * ((KZ^Hilln)/(KZ^Hilln + Z^Hilln)) - aZ*Z
    return(list(c(dZ)))
  })
}

# ODE system for positive feedback autoregulation:
# When Z = 0, this is a stable point
# Will need to give a basal expression level of Z to activate the circuit
ODEs_PAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    # Add small amount of Z for basal expression
    # addBase <- ((t > Xstart) && (t < Xstart + 0.01))
    
    dZ <- base + bZ * (X^Hilln/(KXZ^Hilln + X^Hilln)) * ((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- base + X * bZ *((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- bZ * (Z^Hilln)/((KZ^Hilln)+(Z^Hilln)) - aZ * Z
    return(list(c(dZ)))
  })
}

# ODE system for C1 FFL:
ODEs_C1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dY <- bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base + bZ * ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# I1 FFL
ODEs_I1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    
    dY <- bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base + bZ * ((X * KY)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# ODE system for Feed forward/back hybrid:
ODEs_FFBH <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    
    # X is step function augmented by Z product
    # baseline X given by environmental cue in Xstart -> Xstop
    # change in X at Xstart is 1, change in X at Xstop is -1
    # Z adds a bit more signal to that
    
    # Manually set X
    X <- XMult * (t >= Xstart && t <= Xstop)
    X <- X + XH
    
    # Hill function component of X, XH
    dXH <- ( Z^Hilln / (KZX^Hilln + Z^Hilln) ) - aX*XH
    
    # Update X
    X <- X + dXH
    
    dY <- bY * X^Hilln/( KY^Hilln + X^Hilln ) - aY*Y
    dZ <- base + bZ *  ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dXH, dY, dZ)))
  })
}

# Solution functions
solve_NAR <- function(Xstart = 1,
                      Xstop = 6,
                      tmax = 10,
                      dt = 0.1,
                      pars) {
  
  #pars <- c(base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  iniState <- c(Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_NAR, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Z)
  return(solution)
}

solve_PAR <- function(Xstart = 1,
                      Xstop = 6,
                      tmax = 10,
                      dt = 0.1,
                      pars) {
  
  #pars <- c(base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  iniState <- c(Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_PAR, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Z)
  return(solution)
}

solve_FFLC1 <- function(Xstart = 1,
                        Xstop = 6,
                        tmax = 10,
                        dt = 0.1,
                        pars) {
  
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  #pars <- c(aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
  iniState <- c(Y=0, Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_C1_FFL, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Y, Z)
  return(solution)
}

solve_FFLI1 <- function(Xstart = 1,
                        Xstop = 6,
                        tmax = 10,
                        dt = 0.1,
                        pars) {
  
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  #pars <- c(aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
  iniState <- c(Y=0, Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_I1_FFL, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Y, Z)
  return(solution)
}

solve_FFBH <- function(Xstart = 1,
                       Xstop = 6,
                       tmax = 10,
                       dt = 0.1,
                       pars) {
  
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  #pars <- c(aX = 1, KZX = 1, aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
  iniState <- c(XH=0, Y=0, Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_FFBH, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Y, Z)
  return(solution)
}


Interpolate <- function(x1, y1, x2, y2, y_target) {
  return(x1 + (y_target - y1) * (x2 - x1) / (y2 - y1));
}

# Returns response time [1], steady state [2], time to steady state [3] 
SteadyState <- function(df, startTime, stopTime, solutionIndex) {
  epsilon = 0.001;
  # Output 
  result <- c(0, 0, 0)
  # Count of entries with stable phenotype
  steadyCount <- 0
  
  # Start index, multiplied by 10 because 10 entries per unit time
  start <- startTime * 10 + 1
  
  stop <- stopTime * 10
  
  if (start > nrow(df)) {
    start <- nrow(df)
  }
  
  if (stop > nrow(df)) {
    stop <- nrow(df)
  }
  
  # Iterate over dataframe, find concentrations and change over time
  for (i in start:stop) {
    c1 = df[i-1, solutionIndex]
    c2 = df[i, solutionIndex]
    
    # If we're stable, increment the count and store the result
    if (abs(c2 - c1) < epsilon) {
      steadyCount = steadyCount + 1
      
      # We're done finding steady state if we've been there a while
      if (steadyCount >= 3) {
        result[2] <- c2
        result[3] <- df[i, 1]
        break
      }
      next
    }
    steadyCount = 0
  }
  
  # Response time is time to halfway the steady state
  half = result[2] * 0.5;
  
  # Figure out where the halfway point is
  for (i in start:nrow(df)) {
    c1 = df[i-1, solutionIndex]
    c2 = df[i, solutionIndex]
    t1 = df[i-1, 1]
    t2 = df[i, 1]
    
    if ((c1 < half & c2 >= half) | (c1 > half & c2 <= half)) {
      result[1] = Interpolate(t1, c1, t2, c2, half) - startTime
      break
    }
  }
  return(result)
}

SecondSteadyState <- function(df, prevSteadyState, prevSteadyStateTime, solutionIndex) {
  result <- c(0.0, 0.0, 0.0)
  half = 0.0;
  epsilon = 0.001;
  steadyCount = 0;
  maxSteadyCount = 5;
  
  # Find start index in solution history
  startIndex = prevSteadyStateTime * 10 + 1
  
  if (startIndex > nrow(df)) {
    startIndex = nrow(df);
  }
  
  # Iterate over dataframe, find concentrations and change over time
  for (i in startIndex:(nrow(df))) {
    c1 = df[i-1, solutionIndex]
    c2 = df[i, solutionIndex]
    
    # If we're stable, increment the count and store the result
    if (abs(c2 - c1) < epsilon) {
      steadyCount = steadyCount + 1
      
      # We're done finding steady state if we've been there a while
      if (steadyCount >= maxSteadyCount) {
        result[2] <- c2
        result[3] <- df[i, 1]
        break
      }
      next
    }
    steadyCount = 0
  }
  
  #Find the response time: taken from difference between old and new steady state
  half = abs(prevSteadyState - result[2]) * 0.5;
  
  #Figure out where the halfway point is
  for (i in startIndex:nrow(df))
  {
    t1 = df[i-1,1];
    t2 = df[i,1];
    c1 = df[i-1, solutionIndex];
    c2 = df[i, solutionIndex];
    if ((c1 < half & c2 >= half) | (c1 > half & c2 <= half)) {
      result[1] = Interpolate(t1, c1, t2, c2, half) - prevSteadyStateTime
      break;
    }
  }
  return(result)
}

# Find max expression [1] and time to half max expression [2]
MaxExpression <- function(df, startTime, solutionIndex) {
  result <- c(0.0, 0.0)
  curMax = 0.0
  curTime = 0.0
  
  startIndex = startTime * 10 + 1
  
  if (startIndex > nrow(df)) {
    startIndex = nrow(df)
  }
  
  for (i in startIndex:nrow(df))
  {
    curVal = df[i, solutionIndex]
    if (curVal > (curMax + 0.0001)) {
      curMax = curVal;
    }
  }
  
  # Response time is time to halfway the steady state
  half = curMax * 0.5;
  
  # Figure out where the halfway point is
  for (i in 2:nrow(df)) {
    c1 = df[i-1, solutionIndex]
    c2 = df[i, solutionIndex]
    t1 = df[i-1, 1]
    t2 = df[i, 1]
    
    if ((c1 < half & c2 >= half) | (c1 > half & c2 <= half)) {
      curTime = Interpolate(t1, c1, t2, c2, half)
      break
    }
  }
  
  
  result[1] = curMax;
  result[2] = curTime - startTime;
  return(result)
}

# Find time above a threshold expression level
TimeAboveThreshold <- function(df, threshold, solutionIndex) {
  DELTA = 0.1; # delta between two measurements
  timeAboveThreshold = 0.0;
  
  for (i in seq_len(nrow(df)))
  {
    curVal = df[i,solutionIndex];
    if (curVal >= threshold) {
      timeAboveThreshold = timeAboveThreshold + DELTA;
    }
  }
  
  return(timeAboveThreshold);
  
}

# Sign sensitive delay: rate of change shift > threshold
# measure where inflection point is: where there is a large shift in the roc
DelayTime <- function(df, startTime, stopTime, solutionIndex, baseline, aZ) {
  result = 0.0
  
  # Threshold measured relative to the rate of change for the simple reg case before XStart, base - aZ * Z
  threshold = max(diff(SimpleRegulation(seq(0, 1, by = 0.1), aZ, baseline)))
  
  if (threshold < 0.001) {
    threshold = 0.001
  }
  
  startIndex = 2 #startTime * 10 + 1
  stopIndex = stopTime * 10
  
  if (startIndex >= nrow(df)) {
    startIndex = nrow(df);
  }
  
  if (stopIndex > nrow(df)) {
    stopIndex = nrow(df)
  }
  
  for (i in startIndex:stopIndex)
  {
    t = df[i,1]
    c1 = df[i-1, solutionIndex]
    c2 = df[i, solutionIndex]
    
    newDiff = abs(c2 - c1)
    
    
    # Check if we are increasing faster than under simple regulation
    if (newDiff > threshold + 0.0001) {
      result = t - startTime
      break;
    }
    
    # if (newDiff < prevGreatestSlope) {
    #   result = t - startTime;
    #   break;
    # }
  }
  
  return(result)
}

# Simple regulation function: steady state is bZ/aZ
SimpleRegulation <- function(t, aZ, baseline) {
  return((baseline / aZ) * (1 - exp(-aZ*t)))
}

# Get width of fitness function
CalcSelectionSigmas <- function(optTraits, fitnessCost, fitnessCostDistance, minimumSigma = 0.1) {
  # Calculate the strength of selection per trait, assuming equal contributions per trait
  
  # Get correct row in optimum matrix for the model
  n = length(optTraits)
  
  # If any of the traits are zero, we need to set a minimum width for that trait
  optTraits[optTraits == 0] = minimumSigma
  
  # change in fitnessCostDistance% of trait optimum over the distance for fitnessCost% change in fitness is sigma	
  cost = sqrt(-(1 / (2 * log(1 - fitnessCost))))
  result = (fitnessCostDistance * optTraits * cost) # Note this is a standard deviation NOT a variance
  
  return(result)
  
}

# Fitness function
FitnessFunction <- function(traits, optima, sigma) {
  # Calculate fitness from a vector of inds' traits
  nTraits = length(traits)
  
  # Calculate fitness width based on the given distance
  sigma = diag(sigma * sigma);
  
  # Normalise fitness
  fitnessNorm = dmvnorm(c(optima), c(optima), sigma)
  
  traits = matrix(traits, ncol = nTraits, byrow = T)
  fitnesses = dmvnorm(traits, c(optima), sigma) / fitnessNorm
  
  return(fitnesses)
}

# Get trait values from expression curve
GetTraitValues <- function(solution, model, p) {
  solution <- as.data.frame(solution)
  
  # Steady state
  if (model == "NAR") {
    data <- (SteadyState(solution, 1.0, 6.0, 3))[1:2]
    return(data)
  }
  
  if (model == "PAR") {
    data <- double(3)
    data[1:2] <- (SteadyState(solution, 1.0, 6.0, 3)[1:2])
    data[3] <- DelayTime(solution, 1.0, 6.0, 3, p["base"], p["aZ"])
    return(data)
  }
  
  if (model == "FFLC1") {
    data <- double(3)
    data[3] <- DelayTime(solution, 1.0, 6.0, 4, p["base"], p["aZ"])
    data[1:2] <- (SteadyState(solution, 1 + data[3], 6.0, 4)[1:2]) # Start at the delay time to avoid identifying that as steady state
    
    return(data)
  }
  
  if (model == "FFLI1") {
    data <- double(3)
    
    data[1:2] <- MaxExpression(solution, 1.0, 4)
    data[3] <- TimeAboveThreshold(solution, data[1] * 0.5, 4)
    return(data)
  }
  
  if (model == "FFBH") {
    data <- double(4)
    data[1:2] <- MaxExpression(solution, 1.0, 4)
    data[3:4] <- SecondSteadyState(solution, data[1], 6, 4)[1:2]
    return(data)
  }
}

# Solve an ODE
SolveModel <- function(p, model) {
  switch (model,
          "NAR"   = { solution <-   solve_NAR(pars = p) },
          "PAR"   = { solution <-   solve_PAR(pars = p) },
          "FFLC1" = { solution <- solve_FFLC1(pars = p) },
          "FFLI1" = { solution <- solve_FFLI1(pars = p) },
          "FFBH"  = { solution <-  solve_FFBH(pars = p) }
  )
  return(solution)
}

# Calculate trait values and fitness
CalcTraitAndFitness <- function(p, model, optima, sigma) {
  # Solve model, catch any warnings from lsoda
  solution <- tryCatch(SolveModel(p, model),
                       warning = function(w) { 
                         # If we have a lsoda warning, return NA
                         if (as.character(w$call[[1]]) == "lsoda")
                           return(NA) 
                       }
  )
  
  if (all(is.na(solution))) {
    return(-1)
  }
  
  traits <- GetTraitValues(solution, model, p)
  
  return(FitnessFunction(traits, optima, sigma))
}

# Calculate optima
CalcOptima <- function(oldOptimum, sigma, desiredFitnessAtShift) {
  n = length(oldOptimum);
  optNew = sqrt(-2 * sigma * sigma * log(desiredFitnessAtShift)/n) + oldOptimum;
  return(optNew);
}

# Only feed in appropriate parameters to the right model
ParsMask <- function(pars, model) {
  switch (model,
          "NAR" = { return(pars[c("aZ", "bZ", "KZ", "KXZ", "Hilln", "XMult", "base")]) },
          "PAR" = { return(pars[c("aZ", "bZ", "KZ", "KXZ", "Hilln", "XMult", "base")]) },
          "FFLC1" = { return(pars[c("aY", "bY", "KY", "KXZ", "aZ", "bZ", "Hilln", 
                                    "XMult", "base")]) },
          "FFLI1" = { return(pars[c("aY", "bY", "KY", "KXZ", "aZ", "bZ", "Hilln", 
                                    "XMult", "base")]) },
          "FFBH" = { return(pars) }
  )
}

# Maximum entropy functions

# Generate the symbol sequence Psi from fitness values y, adjusted for distance
# between component hypercube samples x for a given epsilon value

#e.g. data
x <- matrix(c(1, 20, 10, 2, 50, 20), nrow = 3)
y <- c(0.1, 0.6, 0.3)

GeneratePsi <- function(x, y, epsilon) {
  result <- integer(length(y))
  
  # Sort X by nearest neighbour
  matX <- as.matrix(x)
  n <- nrow(matX)
  orderX <- integer(n)
  orderX[1] <- 1 # start with the first row
  deltaX <- integer(n - 1)
  
  # ignore already sampled rows (i.e. 1)
  available <- setdiff(seq_len(n), 1)
  
  for (i in 2:n) {
    last_idx <- orderX[i - 1]
    
    # get distances between this row and all others
    rowDiffs <- (matX[available,] - matrix(matX[last_idx, ], 
                                           nrow = length(available), 
                                           ncol = ncol(matX), byrow = TRUE))^2

    # Only two rows -> keep matrix form
    if (is.null(dim(rowDiffs))) {
      rowDiffs <- matrix(rowDiffs, nrow = 1, byrow = T)
    }
    
    dX <- sqrt(rowSums(rowDiffs))
    
    # Minimum is the deltaX we want
    deltaX[i-1] <- min(dX)
    
    # Now get index of dX in original matrix (minus first entry)
    closest.idx <- which.min(dX)
    orderX[i] <- available[closest.idx]
    
    # Remove this index from available
    available <- available[-closest.idx]
  }
  
  # reorder y to calculate dX: indices are offset by 1 in the y vector
  y <- y[orderX]
  deltaY <- diff(y)
  
  # Adjust dY by dX
  diffXY <- (deltaY/deltaX)
  
  # If we try to divide by 0, fix
  diffXY[is.nan(diffXY)] <- 0
  
  result[diffXY < -epsilon] <- -1
  result[abs(diffXY) <= epsilon] <- 0
  result[diffXY > epsilon] <- 1
  
  return(result)
}

# Calculate H(epsilon) for a given symbol sequence psi
CalculateInformationContent <- function(x) {
  # Ensure x has at least two elements
  if (length(x) < 2) {
    stop("Input vector must have at least two elements.")
  }
  
  # Generate all adjacent pairs
  pairs <- paste(head(x, -1), tail(x, -1), sep = ",")
  
  # Compute frequencies of each pair
  pair_counts <- table(pairs)
  pair_probs <- pair_counts / sum(pair_counts)
  
  # Compute entropy H using base-6 logarithm
  H <- -sum(pair_probs * log(pair_probs, base = 6))
  
  return(H)
}

# Calculate M(e) for a given symbol sequence psi
CalculatePartialInformationContent <- function(x) {
  # Remove all repeated symbols
  marked <- logical(length(x))
  
  # mark first value, then subsequent ones
  marked[1] <- x[1] == 0
  marked <- diff(x) == 0 | (x[-1] == 0)
  
  # for (i in 2:length(x)) {
  #   if (x[i] == x[i - 1]) {
  #     marked[i] <- T
  #   }
  # }
  
  mu <- sum(!marked)
  return(mu / length(x))
  
}

################################################################################
# Script
################################################################################

models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")
NUM_RUNS <- 10000
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

df_test <- data.frame(epsilon = seq(0, 10, length.out = 1000),
                      H = numeric(1000),
                      M = numeric(1000))

for (i in 1:nrow(df_test)) {
  psi <- GeneratePsi(parsMasked, S, df_test$epsilon[i])
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
    

