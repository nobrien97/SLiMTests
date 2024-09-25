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
    dZ <- bZ * (t > Xstart && t <= Xstop) * ((X^Hilln)/(KXZ^Hilln + X^Hilln)) * ((KZ^Hilln)/(KZ^Hilln + Z^Hilln)) - aZ*Z
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 25, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 8)
iniState <- c(Z=0)
times <- seq(0,30,by=0.01)
solution <- ode(iniState, times, ODEs_FBA, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Z)

plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 25, ymin = 0, ymax = 1.05,
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


# ODE system for positive feedback autoregulation:
# When Z = 0, this is a stable point
# Will need to give a basal expression level of Z to activate the circuit
ODEs_FBA_pos <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    # Add small amount of Z for basal expression
    addBase <- ((t > Xstart) && (t < Xstart + 0.01))
    dZ <- addBase * base + X * bZ * (X^Hilln/(KXZ^Hilln + X^Hilln)) * ((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- bZ * (Z^Hilln)/((KZ^Hilln)+(Z^Hilln)) - aZ * Z
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 25, base = 0.01, aZ = 1, bZ = 2, 
            KXZ = 0.1, KZ = 1, Hilln = 1)

iniState <- c(Z=0)
times <- seq(0,30,by=0.01)
solution <- ode(iniState, times, ODEs_FBA_pos, params) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
  select(time, X, Z)

# plot
plotZ <- ggplot(solution, aes(x = time, y = Z)) +
  annotate("rect", xmin = 1, xmax = 25, ymin = 0, ymax = 1.05,
           alpha = .2, fill = "#11241E") +
  geom_line(linewidth = 1.5, colour = "#AA2610") +
  #scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotZ



# AUC
d_auc$time <- solution$time
d_auc$Z <- solution$Z

d_auc %>%
  summarise(auc = AUC(time, Z, absolutearea = T, method = "trapezoid"))
