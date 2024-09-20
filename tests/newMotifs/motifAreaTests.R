library(tidyverse)
library(deSolve)
library(DescTools)

Hilln <- 8 # constant for Hill function

# ODE system for feedback autoregulation:
ODEs_FBA <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * (t > Xstart && t <= Xstop) * (X^Hilln)/(KXZ^Hilln + X^Hilln) * (KZ^Hilln)/(KZ^Hilln + Z^Hilln) - aZ*Z
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 8)
iniState <- c(Z=0)
times <- seq(0,10,by=0.01)
solution <- ode(iniState, times, ODEs_FBA, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Z)

# AUC
d_auc <- data.frame(
  time = numeric(nrow(solution)),
  Z = numeric(nrow(solution)))

d_auc$time <- solution$time
d_auc$Z <- solution$Z

d_auc %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))


# ODE system for positive feedback autoregulation:
ODEs_FBA_pos <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- X * 0.01 + bZ * X^Hilln/(KXZ^Hilln + X^Hilln) * (Z^Hilln)/((KZ^Hilln)+(Z^Hilln)) - aZ*Z
    #dZ <- bZ * (t > Xstart && t <= Xstop) * (X^Hilln)/(KXZ^Hilln + X^Hilln) * (KZ^Hilln)/(KZ^Hilln + Z^Hilln) - aZ*Z
    
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 8)

iniState <- c(Z=0)
times <- seq(0,10,by=0.01)
solution <- ode(iniState, times, ODEs_FBA_pos, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Z)


# AUC
d_auc$time <- solution$time
d_auc$Z <- solution$Z

d_auc %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))
