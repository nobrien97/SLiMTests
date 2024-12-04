library(tidyverse)
library(deSolve)
library(DescTools)

Hilln <- 1 # constant for Hill function

# ODE system for feedback autoregulation:
ODEs_FBA <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dZ <- base + bZ * (t > Xstart && t <= Xstop) * ((X^Hilln)/(KXZ^Hilln + X^Hilln)) * ((KZ^Hilln)/(KZ^Hilln + Z^Hilln)) - aZ*Z
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, base = 0, aZ = 2, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
iniState <- c(Z=0)
times <- seq(0,10,by=0.1)
solution <- ode(iniState, times, ODEs_FBA, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Z)

plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotZ

# AUC
d_auc <- data.frame(
  time = numeric(nrow(solution)),
  Z = numeric(nrow(solution)))

d_auc$time <- solution$time
d_auc$Z <- solution$Z

d_auc %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))

# Steady state + response time
Interpolate <- function(x1, y1, x2, y2, y_target) {
  return(x1 + (y_target - y1) * (x2 - x1) / (y2 - y1));
}

# Returns response time [1], steady state [2], time to steady state [3] 
SteadyState <- function(df, startTime, solutionIndex) {
  epsilon = 0.001;
  # Output 
  result <- c(0, 0, 0)
  # Count of entries with stable phenotype
  steadyCount <- 0
  
  # Start index, multiplied by 10 because 10 entries per unit time
  start <- startTime * 10 + 1
  
  if (start > nrow(df)) {
    start <- nrow(df)
  }
    # Iterate over dataframe, find concentrations and change over time
    for (i in start:(nrow(df))) {
      c1 = df[i-1, solutionIndex]
      c2 = df[i, solutionIndex]
      
      # If we're stable, increment the count and store the result
      if (abs(c2 - c1) < epsilon) {
        steadyCount = steadyCount + 1
        
        # We're done finding steady state if we've been there a while
        if (steadyCount >= 5) {
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
    for (i in 2:nrow(df)) {
      c1 = df[i-1, solutionIndex]
      c2 = df[i, solutionIndex]
      t1 = df[i-1, 1]
      t2 = df[i, 1]
      
      if ((c1 < half & c2 >= half) | (c1 > half & c2 <= half)) {
        result[1] = Interpolate(t1, c1, t2, c2, half)
        break
      }
    }
    return(result)
}

# Returns response time to steady state [1], second steady state [2], time to steady state [3]
SecondSteadyState <- function(df, prevSteadyState, prevSteadyStateTime, solutionIndex) {
  result <- c(0.0, 0.0, 0.0)
  half = 0.0;
  epsilon = 0.001;
  steadyCount = 0;
  maxSteadyCount = 4;
  
  # Find start index in solution history
  startIndex = prevSteadyStateTime * 10
  
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
    print(prevSteadyState)
    print(result[2])
    print(half)
    #Figure out where the halfway point is
    for (i in 2:nrow(df))
    {
        t1 = df[i-1,1];
        t2 = df[i,1];
        c1 = df[i-1, solutionIndex];
        c2 = df[i, solutionIndex];
        if ((c1 < half & c2 >= half) | (c1 > half & c2 <= half)) {
            result[1] = Interpolate(t1, c1, t2, c2, half)
            break;
        }
    }
    return(result)
}

# Find max expression and time to max expression
MaxExpression <- function(df, solutionIndex) {
  result <- c(0.0, 0.0)
  curMax = 0.0
  curTime = 0.0
  for (i in 1:nrow(df))
  {
    curVal = df[i, solutionIndex]
    if (curVal > curMax) {
      curMax = curVal;
      curTime = df[i, 1];
    }
  }
  
  result[1] = curMax;
  result[2] = curTime;
  return(result)
}

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

# Sign sensitive delay
DelayTime <- function(df, startTime, solutionIndex) {
  result = 0.0;
  epsilon = 0.001;
  
  startIndex = startTime * 10 + 1
  
  if (startIndex >= nrow(df)) {
    startIndex = nrow(df);
  }
  
  for (i in startIndex:nrow(df))
  {
    t = df[i,1]
    c1 = df[i-1, solutionIndex]
    c2 = df[i, solutionIndex]
    
    if (abs(c2 - c1) > epsilon) {
      result = t - startTime;
      break;
    }
  }
  
  return(result)
}


NAR_steadyState <- SteadyState(d_auc, 1.0, 2)
DelayTime(d_auc, 1.0, 2)
MaxExpression(d_auc, 2)

SecondSteadyState(d_auc, NAR_steadyState[2], 6, 2)

# Plot the steady state and halfway point
ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  geom_point(x = result[2], y = half, colour = "cornflowerblue") +
  geom_point(x = result[2] * 3, y = result[1], colour = "forestgreen") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")



# ODE system for positive feedback autoregulation:
# When Z = 0, this is a stable point
# Will need to give a basal expression level of Z to activate the circuit
ODEs_FBA_pos <- function(t, state, parameters) {
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

params <- c(Xstart = 1, Xstop = 6, base = 1.2214 - 0.99, aZ = 1, bZ = 1, 
            KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)

iniState <- c(Z=0)
times <- seq(0,10,by=0.1)
solution <- ode(iniState, times, ODEs_FBA_pos, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Z)

# plot
plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotZ

# AUC
solution %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))

d_auc <- data.frame(
  time = numeric(nrow(solution)),
  Z = numeric(nrow(solution)))

d_auc$time <- solution$time
d_auc$Z <- solution$Z



PAR_steadyState <- SteadyState(d_auc, 1.0, 2)
DelayTime(d_auc, 1.0, 2)


# ODE system for C1 FFL:
ODEs_C1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dY <- base + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base + bZ * ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}
log(2)
params <- c(Xstart = 1, Xstop = 6, aY = 1, bY = 1, KY = 0.5, KXZ = 2, 
            aZ = 1, bZ = 1, Hilln = 2, XMult = 1, base = 0)

iniState <- c(Y=0, Z=0)
times <- seq(0,10,by=0.1)
solution <- ode(iniState, times, ODEs_C1_FFL, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Y, Z)

# plot
plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotZ

# AUC
solution %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))


d_auc <- data.frame(
  time = numeric(nrow(solution)),
  Z = numeric(nrow(solution)))

d_auc$time <- solution$time
d_auc$Z <- solution$Z

C1FFL_steadyState <- SteadyState(d_auc, 1.0, 2)
C1FFLDelay <- DelayTime(d_auc, 1.0, 2)

c(C1FFL_steadyState[1], C1FFLDelay, C1FFL_steadyState[2])



# I1 FFL
ODEs_I1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    
    dY <- base + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base + bZ *  ((X * KY)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}


params <- c(Xstart = 1, Xstop = 6, 
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 0)

iniState <- c(Y=0, Z=0)
times <- seq(0,10,by=0.1)
solution <- ode(iniState, times, ODEs_I1_FFL, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Y, Z)

# plot
plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotZ

# AUC
solution %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))

d_auc <- data.frame(
  time = numeric(nrow(solution)),
  Z = numeric(nrow(solution)))

d_auc$time <- solution$time
d_auc$Z <- solution$Z

I1FFL_maxExpression <- MaxExpression(d_auc, 2)

I1FFL_timeAboveThreshold <- TimeAboveThreshold(d_auc, I1FFL_maxExpression[1]/2, 2)

c(I1FFL_maxExpression[2], I1FFL_maxExpression[1], I1FFL_timeAboveThreshold)


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

params <- c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 0)

iniState <- c(XH=0,Y=0, Z=0)
times <- seq(0,10,by=0.1)
solution <- ode(iniState, times, ODEs_FFBH, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, XH, Y, Z)

# plot
plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotZ

solution %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))

d_auc <- data.frame(
  time = numeric(nrow(solution)),
  Z = numeric(nrow(solution)))

d_auc$time <- solution$time
d_auc$Z <- solution$Z

FFBH_maxExpression <- MaxExpression(d_auc, 2)
FFBH_secondState <- SecondSteadyState(d_auc, FFBH_maxExpression[1], 6, 2)

c(FFBH_maxExpression[2], FFBH_maxExpression[1], FFBH_secondState[1], FFBH_secondState[2])