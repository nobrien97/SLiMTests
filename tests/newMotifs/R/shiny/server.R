library(deSolve)
library(tidyverse)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)
library(shiny)
library(RColorBrewer)
library(shinyjs)

# Colour scheme:
colX   <- brewer.pal(6, "Paired")[2]
colXbg <- brewer.pal(6, "Paired")[1]
colY   <- brewer.pal(6, "Paired")[6]
colYbg <- brewer.pal(6, "Paired")[5]
colZ   <- brewer.pal(6, "Paired")[4]
colZbg <- brewer.pal(6, "Paired")[3]

colX   <- "#44546A"
colXOutline <- "#222A35"


# plot graphs:

genes <- tibble(name = c(expression("X"), expression("Z")))
activationNAR <- tibble(from = c(1, 2),
                        to =   c(2, 2))
gNAR <- tbl_graph(edges = activationNAR, directed = TRUE, nodes = genes)
plt_NAR <- ggraph(gNAR, layout = "manual", x = c(0,0), y=c(1,0)) +
  geom_node_point(size = 16, color = c(colXOutline)) +
  geom_node_point(size = 15, color = c(colX)) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = 0.0) +
  geom_edge_loop(aes(start_cap = circle(8, unit = "mm"),
                     end_cap = circle(8, unit = "mm"),
                     span = 90, direction = 0, strength = 0.7),
                 arrow = arrow(type = "open", angle = 90, length = unit(3, 'mm'))) +
  scale_x_continuous(limits = c(-0.5,0.5)) +
  scale_y_continuous(limits = c(-0.2,1.1)) +
  coord_fixed() +
  theme_graph()

# PAR
plt_PAR <- ggraph(gNAR, layout = "manual", x = c(0,0), y=c(1,0)) +
  geom_node_point(size = 16, color = c(colXOutline)) +
  geom_node_point(size = 15, color = c(colX)) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = 0.0) +
  geom_edge_loop(aes(start_cap = circle(8, unit = "mm"),
                     end_cap = circle(8, unit = "mm"),
                     span = 90, direction = 0, strength = 0.7),
                 arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm'))) +
  scale_x_continuous(limits = c(-0.5,0.5)) +
  scale_y_continuous(limits = c(-0.2,1.1)) +
  coord_fixed() +
  theme_graph()

genes <- tibble(name = c(expression("X"), expression("Y"), expression("Z")))
activationFFL <- tibble(from = c(1, 1, 2),
                        to =   c(2, 3, 3))
gFFL <- tbl_graph(edges = activationFFL, directed = TRUE, nodes = genes)
plt_FFLC1 <- ggraph(gFFL, layout = "manual", x = c(0,-0.7,0.7), y=c(1,0,0)) +
  geom_node_point(size = 16, color = colXOutline) +
  geom_node_point(size = 15, color = colX) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = c(-0.2, 0.2, -0.2)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
  coord_fixed() +
  theme_graph()

# I1 FFL
plt_FFLI1 <- ggraph(gFFL, layout = "manual", x = c(0,-0.7,0.7), y=c(1,0,0)) +
  geom_node_point(size = 16, color = colXOutline) +
  geom_node_point(size = 15, color = colX) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = c(30, 30, 90), length = unit(3, 'mm')),
                strength = c(-0.2, 0.2, -0.2)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
  coord_fixed() +
  theme_graph()

# Arabidopsis FFL with feedback
genes <- tibble(name = c(expression("X"), expression("Y"), expression("Z")))
activationFFLA <- tibble(from = c(1, 1, 2, 3),
                         to =   c(2, 3, 3, 1))
gFFLA <- tbl_graph(edges = activationFFLA, directed = TRUE, nodes = genes)
plt_FFBH <- ggraph(gFFLA, layout = "manual", x = c(0,-0.7,0.7), y=c(1,0,0)) +
  geom_node_point(size = 16, color = colXOutline) +
  geom_node_point(size = 15, color = colX) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = c(-0.2, 0.2, -0.2, 0.2)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
  coord_fixed() +
  theme_graph()


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

plotDynamics_NAR <- function(solution, Xstart, Xstop) {
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_NAR, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
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

plotDynamics_PAR <- function(solution, Xstart = 1, Xstop = 6) {
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    scale_x_continuous(breaks = seq(0, 10, by = 0.5)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_PAR, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
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

plotDynamics_FFLC1 <- function(solution, Xstart = 1, Xstop = 6) {
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    scale_x_continuous(breaks = seq(0, 10, by = 0.5)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_FFLC1, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
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

plotDynamics_FFLI1 <- function(solution, Xstart = 1, Xstop = 6) {
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_FFLI1, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
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

plotDynamics_FFBH <- function(solution, Xstart = 1, Xstop = 6) {
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_FFBH, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
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
  
  if (startTime > nrow(df)) {
    startTime = nrow(df)
  }
  
  for (i in startTime:nrow(df))
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

getStats <- function(solution, model, pars = NULL) {
  # Steady state
  if (model == "NAR") {
    data <- (SteadyState(solution, 1.0, 6.0, 3))[1:2]
    d <- as.data.frame(as.character(round(data, 3)))
    row.names(d) <- c("Response time", "Steady state concentration")
    return(d)
  }
  
  if (model == "PAR") {
    data <- double(3)
    # if (pars["base"] == 0.1) {
    #   browser()
    # }
    data[1:2] <- (SteadyState(solution, 1.0, 6.0, 3)[1:2])
    data[3] <- DelayTime(solution, 1.0, 6.0, 3, pars["base"], pars["aZ"])
    d <- as.data.frame(as.character(round(data, 3)))
    row.names(d) <- c("Response time", "Steady state concentration", "Sign-Sensitive Delay Time")
    return(d)
  }
  
  if (model == "FFL-C1") {
    data <- double(3)
    prevDT <- solution[11, "Z"] - solution[10, "Z"]
    data[3] <- DelayTime(solution, 1.0, 6.0, 4, pars["base"], pars["aZ"])
    data[1:2] <- (SteadyState(solution, 1 + data[3], 6.0, 4)[1:2]) # Start at the delay time to avoid identifying that as steady state
    
    if (pars["Hilln"] == 4.1) {
      browser()
    }
    d <- as.data.frame(as.character(round(data, 3)))
    row.names(d) <- c("Response time", "Steady state concentration", "Sign-Sensitive Delay Time")
    return(d)
  }
  
  if (model == "FFL-I1") {
    data <- double(3)
    
    data[1:2] <- MaxExpression(solution, 1.0, 4)
    data[3] <- TimeAboveThreshold(solution, data[1] * 0.5, 4)
    d <- as.data.frame(as.character(round(data, 3)))
    row.names(d) <- c("Maximum expression", "Time to half maximum expression", "Time above half maximum expression")
    return(d)
  }
  
  if (model == "FFBH") {
    data <- double(4)
    data[1:2] <- MaxExpression(solution, 1.0, 4)
    data[3:4] <- SecondSteadyState(solution, data[1], 6, 4)[1:2]
    d <- as.data.frame(as.character(round(data, 3)))
    row.names(d) <- c("Maximum expression", "Time to half maximum expression", "Response time to final steady state", 
                      "Final steady state concentration")
    return(d)
  }
}



motif_strings <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
server<-function(input, output) {
  
  pars_NAR <- reactive({
    c(base = input$base_NAR, aZ = input$aZ_NAR, bZ = input$bZ_NAR, KXZ = input$KXZ_NAR,
      KZ = input$KZ_NAR, Hilln = input$Hilln_NAR, XMult = input$XMult_NAR)
  })
  
  pars_PAR <- reactive({
    c(base = input$base_PAR, aZ = input$aZ_PAR, bZ = input$bZ_PAR, KXZ = input$KXZ_PAR,
      KZ = input$KZ_PAR, Hilln = input$Hilln_PAR, XMult = input$XMult_PAR)
  })
  
  pars_FFLC1 <- reactive({
    c(base = input$base_FFLC1, aY = input$aY_FFLC1, bY = input$bY_FFLC1, KY = input$KY_FFLC1,
      aZ = input$aZ_FFLC1, bZ = input$bZ_FFLC1,  KXZ = input$KXZ_FFLC1, Hilln = input$Hilln_FFLC1, 
      XMult = input$XMult_FFLC1)
  })
  
  pars_FFLI1 <- reactive({
    c(base = input$base_FFLI1, aY = input$aY_FFLI1, bY = input$bY_FFLI1, KY = input$KY_FFLI1,
      aZ = input$aZ_FFLI1, bZ = input$bZ_FFLI1, KXZ = input$KXZ_FFLI1, Hilln = input$Hilln_FFLI1, 
      XMult = input$XMult_FFLI1)
  })
  
  pars_FFBH <- reactive({
    c(base = input$base_FFBH, aX = input$aX_FFBH, KZX = input$KZX_FFBH, aY = input$aY_FFBH, 
      bY = input$bY_FFBH, KY = input$KY_FFBH, aZ = input$aZ_FFBH, bZ = input$bZ_FFBH, 
      KXZ = input$KXZ_FFBH, Hilln = input$Hilln_FFBH, XMult = input$XMult_FFBH)
  })
  
  solution_NAR <- reactive({
    solve_NAR(Xstart = input$Xstart_NAR, Xstop = input$Xstop_NAR,
              tmax = input$tmax_NAR, pars = pars_NAR())
  })
  
  solution_PAR <- reactive({
    solve_PAR(Xstart = input$Xstart_PAR, Xstop = input$Xstop_PAR,
              tmax = input$tmax_PAR, pars = pars_PAR())
  })
  
  solution_FFLC1 <- reactive({
    solve_FFLC1(Xstart = input$Xstart_FFLC1, Xstop = input$Xstop_FFLC1,
                tmax = input$tmax_FFLC1, pars = pars_FFLC1())
  })
  
  solution_FFLI1 <- reactive({
    solve_FFLI1(Xstart = input$Xstart_FFLI1, Xstop = input$Xstop_FFLI1,
                tmax = input$tmax_FFLI1, pars = pars_FFLI1())
  })
  
  solution_FFBH <- reactive({
    solve_FFBH(Xstart = input$Xstart_FFBH, Xstop = input$Xstop_FFBH, 
               tmax = input$tmax_FFBH, pars = pars_FFBH())
  })
  
  output$main_plot_NAR <- renderPlot({
    if (input$motif != 1)
      return()
      plotDynamics_NAR(solution_NAR(), Xstart = input$Xstart_NAR, Xstop = input$Xstop_NAR)
    
  })
  
  output$main_plot_PAR <- renderPlot({
    if (input$motif != 2)
      return()
    
      plotDynamics_PAR(solution_PAR(), Xstart = input$Xstart_PAR, Xstop = input$Xstop_PAR)

  })
  
  output$main_plot_FFLC1 <- renderPlot({
    if (input$motif != 3)
      return()
    
      plotDynamics_FFLC1(solution_FFLC1(), Xstart = input$Xstart_FFLC1, Xstop = input$Xstop_FFLC1)

  })
  
  output$main_plot_FFLI1 <- renderPlot({
    if (input$motif != 4)
      return()
    
      plotDynamics_FFLI1(solution_FFLI1(), Xstart = input$Xstart_FFLI1, Xstop = input$Xstop_FFLI1)

  })
  
  output$main_plot_FFBH <- renderPlot({
    if (input$motif != 5)
      return()
    
      plotDynamics_FFBH(solution_FFBH(), Xstart = input$Xstart_FFBH, Xstop = input$Xstop_FFBH)

  })
  

  output$NARTabOutput <- renderTable({
    if (input$motif != 1)
      return()
    
    getStats(as.data.frame(solution_NAR()), "NAR")
  }, colnames = F, rownames = T, width = "100%", align = "lr")
  
  output$PARTabOutput <- renderTable({
    if (input$motif != 2)
      return()
    
    getStats(as.data.frame(solution_PAR()), "PAR", pars_PAR())
  }, colnames = F, rownames = T, width = "100%", align = "lr")
  
  output$FFLC1TabOutput <- renderTable({
    if (input$motif != 3)
      return()
    
    getStats(as.data.frame(solution_FFLC1()), "FFL-C1", pars_FFLC1())
  }, colnames = F, rownames = T, width = "100%", align = "lr")
  
  output$FFLI1TabOutput <- renderTable({
    if (input$motif != 4)
      return()
    
    getStats(as.data.frame(solution_FFLI1()), "FFL-I1")
  }, colnames = F, rownames = T, width = "100%", align = "lr")
  
  output$FFBHTabOutput <- renderTable({
    if (input$motif != 5)
      return()
    
    getStats(as.data.frame(solution_FFBH()), "FFBH")
  }, colnames = F, rownames = T, width = "100%", align = "lr")
  
  
}