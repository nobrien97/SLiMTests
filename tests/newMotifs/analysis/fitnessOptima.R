# Calculates ruggedness of fitness landscape by number of unique fitness peaks
library(deSolve)
library(tidyverse)
library(ggraph)
library(igraph)
library(mvtnorm)
library(raster)
library(foreach)
library(doParallel)
library(future)
library(testit)
library(broom)

select <- dplyr::select

################################################################################
# Function definitions
################################################################################
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

################################################################################
# Main Script
################################################################################

# Optimise fitness function from different starting points for each motif, 
# moving in random directions

#pars <- c(base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
#pars <- c(base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
#pars <- c(aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
#pars <- c(aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
pars <- c(aX = 0, KZX = 0, aY = 0, bY = 0, KY = 0, KZ = 0, KXZ = 0, 
          aZ = 0, bZ = 0, Hilln = 0, XMult = 0, base = -Inf)
#pars <- c(base = 0, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
model <- "NAR"

# Get starting trait values
parsMasked <- ParsMask(pars, model)
startSolution <- SolveModel(exp(parsMasked), model)
startTraits <- GetTraitValues(startSolution, model, exp(parsMasked))
sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)

# Get optimum
opt <- CalcOptima(startTraits, sigma, 0.9)

# Should = 0.9 if all is working
assert("Fitness incorrect - something wrong in fitness calculation", 
       CalcTraitAndFitness(exp(parsMasked), model, opt, sigma) == 0.9)

# Optimise using differential evolution (in case of discontinuities)
#optimSolution <- JDEoptim()

# Set up parallel processing
cl <- parallel::makeCluster(future::availableCores())
doParallel::registerDoParallel(cl)
N_SAMPLES <- 100

# Output list
result <- vector(mode = "list", length = 5L)
names(result) <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")

# General pars for all models - masked by modeltype in the loop
pars <- c(aX = 0, KZX = 0, aY = 0, bY = 0, KY = 0, KZ = 0, KXZ = 0, 
          aZ = 0, bZ = 0, Hilln = 0, XMult = 0, base = -Inf)

# Repeat for each model
for (model in "FFBH") {#names(result)) {
  # Setup trait optima etc.
  parsMasked <- ParsMask(pars, model)
  startSolution <- SolveModel(exp(parsMasked), model)
  startTraits <- GetTraitValues(startSolution, model, exp(parsMasked))
  sigma <- CalcSelectionSigmas(startTraits, 0.1, 0.1, 0.1)
  opt <- CalcOptima(startTraits, sigma, 0.9)
  
  # Iterate with random samples
  optima <- foreach (i = seq_len(N_SAMPLES), .combine = rbind) %dopar% {
    require(raster)
    require(tidyverse)
    require(deSolve)
    require(mvtnorm)
    # Randomly sample a new starting point and try to find a new optimum
    newPars <- clamp(exp(parsMasked + rnorm(length(parsMasked), 0, 0.01)), 0.0, 3)
    
    # Solve
    optimSolution <- optim(newPars, CalcTraitAndFitness, method = "L-BFGS-B", 
                           control = list(fnscale = -1), lower = 0.001, upper = 3,
                           model = model, optima = opt, sigma = sigma
                           )
    
    # If we have converged, fill in the table
    if (optimSolution$convergence == 0) {
      return(c(optimSolution$value, optimSolution$par))
    }
    
    # If we haven't solved, return an NA row
    return(rep(NA, length(parsMasked) + 1))
  }
  
  # Save results
  result[[model]] <- optima 
}

parallel::stopCluster(cl)

test_res <- as.data.frame(t(rep(NA, length(parsMasked) + 1)))
names(test_res) <- c("w", names(parsMasked))
test_res[N_SAMPLES,] <- NA

# single thread debug
for (i in seq_len(N_SAMPLES)) {
  # Randomly sample a new starting point and try to find a new optimum
  newPars <- clamp(exp(parsMasked + runif(length(parsMasked), -3, 3)), 0.0, 3)
  
  optimSolution <- tryCatch(optim(newPars, CalcTraitAndFitness, method = "L-BFGS-B", 
                                  control = list(fnscale = -1), lower = 0.001, upper = 3,
                                  model = model, optima = opt, sigma = sigma
                    ),
                    error = function(e) {
                      return(rep(NA, length(parsMasked) + 1))
                    },
                    warning = function(w) {
                      return(rep(NA, length(parsMasked) + 1))
                    }
              )
  
  # First check if NA
  if (sum(is.na(optimSolution)) == length(parsMasked) + 1) {
    return(optimSolution)
  }
  
  # If we have converged, fill in the table
  if (optimSolution$convergence == 0) {
    test_res[i,] <- (c(optimSolution$value, optimSolution$par))
    next
  }
  
  # If we haven't solved, return an NA row
  test_res[i,] <- (rep(NA, length(parsMasked) + 1))
}

# Test for bad solutions 
newPars <- exp(parsMasked + rnorm(length(parsMasked), 0, 10))
test_sol <- SolveModel(newPars, model)
test_sol <- CalcTraitAndFitness(newPars, model, opt, sigma)

# read in fitness optima
fitnessOptima <- readRDS("./fitnessOptima.RDS")

# Remove invalid solutions
for (i in seq_along(fitnessOptima)) {
  fitnessOptima[[i]] <- as.data.frame(fitnessOptima[[i]]) %>% 
    drop_na() %>%
    rename(w = V1) %>%
    filter(w >= 0) 
}

# Look for unique optima: those which are within a certain distance of each other
norm_vec <- function(x) sqrt(crossprod(x))
dists <- vector(mode = "list", length = length(fitnessOptima))
for (i in seq_along(fitnessOptima)) {
  curData <- fitnessOptima[[i]]
  dists[[i]] <- dist(as.matrix(curData))
}

# Plot distances
dist_NAR <- dists[[1]]

df_nar <- tidy(dist_NAR)

ggplot(df_nar, aes(x = distance)) +
  geom_density() +
  theme_bw() +
  labs(x = "Distance", y = "Density") +
  theme(text = element_text(size = 12))

df_par <- tidy(dists[[2]])
ggplot(df_par, aes(x = distance)) +
  geom_density() +
  theme_bw() +
  labs(x = "Distance", y = "Density") +
  theme(text = element_text(size = 12))

df_fflc1 <- tidy(dists[[3]])
ggplot(df_fflc1, aes(x = distance)) +
  geom_density() +
  theme_bw() +
  labs(x = "Distance", y = "Density") +
  theme(text = element_text(size = 12))

df_ffli1 <- tidy(dists[[4]])
ggplot(df_ffli1, aes(x = distance)) +
  geom_density() +
  theme_bw() +
  labs(x = "Distance", y = "Density") +
  theme(text = element_text(size = 12))

df_ffbh <- tidy(dists[[5]])
ggplot(df_ffbh, aes(x = distance)) +
  geom_density() +
  theme_bw() +
  labs(x = "Distance", y = "Density") +
  theme(text = element_text(size = 12))

# Mark similar optima
uniqueOptima <- vector(mode = "list", length(dists))
for (i in seq_along(dists)) {
 uniqueOptima[[i]] <- tidy(dists[[i]]) %>% 
   mutate(largeDist = distance > 0.5) %>%
   group_by(item1) %>%
   summarise(isUnique = all(largeDist))
}

# Number of optima
df_opt <- data.frame(model = "",
                     totalOptima = 0,
                     countUniqueOptima = 0,
                     propUniqueOptima = 0.0)
df_opt[length(uniqueOptima),] <- NA
models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
for (i in seq_along(uniqueOptima)) {
  df_opt[i,] <- c(models[i],
                  nrow(uniqueOptima[[i]]),
                  sum(uniqueOptima[[i]]$isUnique),
                  sum(uniqueOptima[[i]]$isUnique) / nrow(uniqueOptima[[i]])
                  )
}


