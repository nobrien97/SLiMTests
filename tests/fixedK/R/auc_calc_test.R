library(tidyverse)
library(deSolve)
library(DescTools)


Hilln <- 8 # constant for Hill function

ODEs_FBA <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * (t > Xstart && t <= Xstop) * (X^Hilln)/(KYZ^Hilln + X^Hilln) * (KZ^Hilln)/(KZ^Hilln + Z^Hilln) - aZ*Z
    dZnoFB <- aZ * (t > Xstart && t <= Xstop) - aZ*ZnoFB
    return(list(c(dZ, dZnoFB)))
  })
}


plotDynamics_FBA <- function(Xstart = 1,
                             Xstop = 6,
                             tmax = 10,
                             dt = 0.1,
                             pars = list(aZ = 20,
                                         bZ1 = 29,
                                         KYZ = 1,
                                         KZ = 1)) {
  params <- c(Xstart = Xstart, Xstop = Xstop, aZ = pars$aZ, bZ = pars$bZ1, KYZ = pars$KYZ, KZ = pars$KZ)
  iniState <- c(Z=0, ZnoFB = 0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_FBA, params) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0),
           aZ = as.factor(pars$aZ), bZ = as.factor(pars$bZ1), 
           KXZ = as.factor(pars$KYZ), KZ = as.factor(pars$KZ)) %>%
    select(time, aZ, bZ, KXZ, KZ, X, Z, ZnoFB)
}

dat <- plotDynamics_FBA(pars = list(aZ = 20, bZ1 = 29, KYZ = 1, KZ = 1))
AUC(dat$time, dat$Z, absolutearea = T)

dat <- plotDynamics_FBA(pars = list(aZ = 11, bZ1 = 16, KYZ = 1, KZ = 1))
AUC(dat$time, dat$Z, absolutearea = T)
