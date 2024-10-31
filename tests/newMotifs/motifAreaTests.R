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
    dZ <- base * X + bZ * (t > Xstart && t <= Xstop) * ((X^Hilln)/(KXZ^Hilln + X^Hilln)) * ((KZ^Hilln)/(KZ^Hilln + Z^Hilln)) - aZ*Z
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1)
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
    
    dZ <- base * X + bZ * (X^Hilln/(KXZ^Hilln + X^Hilln)) * ((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- base + X * bZ *((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- bZ * (Z^Hilln)/((KZ^Hilln)+(Z^Hilln)) - aZ * Z
    return(list(c(dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, base = 1, aZ = 1, bZ = 1, 
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

# ODE system for C1 FFL:
ODEs_C1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dY <- base * X + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base * X + bZ * ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, aY = 1, bY = 1, KY = 1, KXZ = 1, 
            aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)

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

# I1 FFL
ODEs_I1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    
    dY <- base * X + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base * X + bZ *  ((X * KY)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}


params <- c(Xstart = 1, Xstop = 6, 
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)

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
    
    dY <- base * X + bY * X^Hilln/( KY^Hilln + X^Hilln ) - aY*Y
    dZ <- base * X + bZ *  ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dXH, dY, dZ)))
  })
}

params <- c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)

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
