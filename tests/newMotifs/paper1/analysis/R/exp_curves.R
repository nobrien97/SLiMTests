library(tidyverse)
library(tidygraph)
library(ggraph)
library(deSolve)
library(cowplot)
library(deSolve)
library(cowplot)
library(latex2exp)
library(paletteer)
library(patchwork)

source("helperFn.R")


# Colour palettes
colX   <- "#44546A"
colXOutline <- "#222A35"


pal_change <- c("#222", "#B00", "#0BB")

# ODE system for feedback autoregulation:
ODEs_FBA <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    dZ <- bZ * (t > Xstart && t <= Xstop) * 1/(1 + Z^h) - aZ*Z
    return(list(c(dZ)))
  })
}

# ODE system for positive feedback autoregulation:
ODEs_FBA_pos <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- X * 0.01 + bZ * X^h/(KXZ^h + X^h) * (Z^h)/((KZ^h)+(Z^h)) - aZ*Z
    #dZ <- bZ * (t > Xstart && t <= Xstop) * (X^h)/(KXZ^h + X^h) * (KZ^h)/(KZ^h + Z^h) - aZ*Z
    
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
    dY <- base * X + bY * X^h/(KY^h + X^h) - aY*Y
    dZ <- base * X + bZ * ((X * Y)^h)/((KXZ^h + X^h) * (KY^h + Y^h)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# ODE system for I1 FFL:
ODEs_I1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    
    dY <- base * X + bY * X^h/(KY^h + X^h) - aY*Y
    dZ <- base * X + bZ *  ((X * KY)^h)/((KXZ^h + X^h) * (KY^h + Y^h)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

ODEs_FFBH <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    
    # X is step function augmented by Z product
    # baseline X given by environmental cue in Xstart -> Xstop
    # change in X at Xstart is 1, change in X at Xstop is -1
    # Z adds a bit more signal to that
    
    # Manually set X
    X <- XMult * (t > Xstart && t <= Xstop)
    X <- X + XH
    
    # Hill function component of X, XH
    dXH <- ( Z^h / (KZX^h + Z^h) ) - aX*XH
    
    # Update X
    X <- X + dXH
    
    dY <- base * X + bY * X^h/( KY^h + X^h ) - aY*Y
    dZ <- base * X + bZ *  ((X * Y)^h)/((KXZ^h + X^h) * (KY^h + Y^h)) - aZ*Z
    
    return(list(c(dXH, dY, dZ)))
  })
}


math_labels <- c("aZ" = "$\\alpha_Z$",
                 "bZ" = "$\\beta_Z",
                 "h" = "$h$")

solveODE <- function(ode_fn, params, iniState, times = NULL, param_treatment = NULL) {
  solution <- NULL
  
  if (is.null(times)) {
    times <- seq(0,10,by=0.01)
  }
  
  if (is.null(param_treatment)) {
    param_treatment <- 1:(length(params[[1]]))
  }
  
  
  for (i in seq_along(params)) {
    solutionCur <- ode(iniState, times, ode_fn, params[[i]]) %>%
      as_tibble() %>%
      mutate(X = ifelse(time > params[[i]]["Xstart"] & time <= params[[i]]["Xstop"], 1, 0)) %>%
      mutate(treatment = TeX(paste(math_labels[names(params[[i]][param_treatment])], "=", 
                               params[[i]][param_treatment], collapse = ","), output = "character")) #%>%
      #select(time, treatment, X, Z)
    solution <- rbind(solution, solutionCur)
  }
  
  param_values <- sort(unlist(lapply(params, function(x){
    x[[param_treatment]]
  })))
  
  solution <- solution %>%
    mutate(treatment = factor(treatment,
                              levels = TeX(paste(math_labels[param_treatment], "=", param_values), 
                                           output = "character")))
  return(solution)
}

params_h <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 30, XMult = 1, base = 0.0)
)

solution_nar <- solveODE(ODEs_FBA, params_h, c(Z = 0), param_treatment = "h")
plotNAR_h <- ggplot(solution_nar, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_nar$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotNAR_h


params_h_par <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 30, XMult = 1, base = 0.0)
)

solution_par <- solveODE(ODEs_FBA_pos, params_h_par, c(Z = 0), param_treatment = "h")
plotPAR_h <- ggplot(solution_par, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_par$Z),
           alpha = .2, fill = colX) +
  geom_line(aes(linetype = treatment), linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  guides(linetype = "none") +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotPAR_h


params_h_fflc1 <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 2, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 2, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 2, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 30, XMult = 1, base = 0.0)
)

solution_fflc1 <- solveODE(ODEs_C1_FFL, params_h_fflc1, c(Y = 0, Z = 0), param_treatment = "h")
plotFFLC1_h <- ggplot(solution_fflc1, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_fflc1$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFLC1_h


params_h_ffli1 <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 30, XMult = 1, base = 0.0)
)

solution_ffli1 <- solveODE(ODEs_I1_FFL, params_h_ffli1, c(Y = 0, Z = 0), param_treatment = "h")
plotFFLI1_h <- ggplot(solution_ffli1, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_ffli1$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFLI1_h


params_h_ffbh <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 2, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 2, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 2, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 2, h = 30, XMult = 1, base = 0.0)
)

solution_ffbh <- solveODE(ODEs_FFBH, params_h_ffbh, c(XH = 0, Y = 0, Z = 0), param_treatment = "h")
plotFFBH_h <- ggplot(solution_ffbh, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_ffbh$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFBH_h

leg <- get_legend(plotFFBH_h)

layout <- "
AABBCC
#DDEE#
"

free(plotNAR_h + ggtitle("NAR")) +
  plotPAR_h + ggtitle("PAR") +
  plotFFLC1_h + ggtitle("FFLC1") +
  plotFFLI1_h + ggtitle("FFLI1") +
  plotFFBH_h + ggtitle("FFBH") +
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "bottom")

ggsave("plt_h_eg.png", device = png, bg = "white", width = 12, height = 7)



params_b <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 2, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 5, h = 1, XMult = 1, base = 0.0)
)

solution_nar <- solveODE(ODEs_FBA, params_b, c(Z = 0), param_treatment = "bZ")
plotNAR_b <- ggplot(solution_nar, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_nar$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotNAR_b



solution_par <- solveODE(ODEs_FBA_pos, params_b, c(Z = 0), param_treatment = "bZ")
plotPAR_b <- ggplot(solution_par, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_par$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotPAR_b

solution_fflc1 <- solveODE(ODEs_C1_FFL, params_b, c(Y = 0, Z = 0), param_treatment = "bZ")
plotFFLC1_b <- ggplot(solution_fflc1, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_fflc1$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFLC1_b

params_b_ffli1 <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 2, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 5, h = 2, XMult = 1, base = 0.0)
)

solution_ffli1 <- solveODE(ODEs_I1_FFL, params_b_ffli1, c(Y = 0, Z = 0), param_treatment = "bZ")
plotFFLI1_b <- ggplot(solution_ffli1, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_ffli1$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFLI1_b


params_b_ffbh <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 2, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 2, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 5, h = 1, XMult = 1, base = 0.0)
)


solution_ffbh <- solveODE(ODEs_FFBH, params_b_ffbh, c(XH = 0, Y = 0, Z = 0), param_treatment = "bZ")
plotFFBH_b <- ggplot(solution_ffbh, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_ffbh$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFBH_b

free(plotNAR_b + ggtitle("NAR")) +
  plotPAR_b + ggtitle("PAR") +
  plotFFLC1_b + ggtitle("FFLC1") +
  plotFFLI1_b + ggtitle("FFLI1") +
  plotFFBH_b + ggtitle("FFBH") +
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "bottom")

ggsave("plt_b_eg.png", device = png, bg = "white", width = 12, height = 7)




params_a <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 4, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 2, bZ = 1, h = 4, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 1, KXZ = 1,
    aZ = 5, bZ = 1, h = 4, XMult = 1, base = 0.0)
)

solution_nar <- solveODE(ODEs_FBA, params_a, c(Z = 0), param_treatment = "aZ")
plotNAR_a <- ggplot(solution_nar, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_nar$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotNAR_a

solution_par <- solveODE(ODEs_FBA_pos, params_a, c(Z = 0), param_treatment = "aZ")
plotPAR_a <- ggplot(solution_par, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_par$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotPAR_a

solution_fflc1 <- solveODE(ODEs_C1_FFL, params_a, c(Y = 0, Z = 0), param_treatment = "aZ")
plotFFLC1_a <- ggplot(solution_fflc1, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_fflc1$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFLC1_a


params_a_ffli1 <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 1, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 2, bZ = 1, h = 2, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 4, KY = 1, KXZ = 1,
    aZ = 5, bZ = 1, h = 2, XMult = 1, base = 0.0)
)

solution_ffli1 <- solveODE(ODEs_I1_FFL, params_a_ffli1, c(Y = 0, Z = 0), param_treatment = "aZ")
plotFFLI1_a <- ggplot(solution_ffli1, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_ffli1$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFLI1_a


params_a_ffbh <- list(
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 1, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 2, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 2, bZ = 1, h = 1, XMult = 1, base = 0.0),
  c(Xstart = 1, Xstop = 6, 
    aX = 1, KZX = 1, KZ = 1,
    aY = 1, bY = 1, KY = 0.1, KXZ = 0.5,
    aZ = 5, bZ = 1, h = 1, XMult = 1, base = 0.0)
)

solution_ffbh <- solveODE(ODEs_FFBH, params_a_ffbh, c(XH = 0, Y = 0, Z = 0), param_treatment = "aZ")
plotFFBH_a <- ggplot(solution_ffbh, aes(x = time, y = Z, colour = treatment)) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = max(solution_ffbh$Z),
           alpha = .2, fill = colX) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = pal_change, labels = function(l) parse(text = l)) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z expression", colour = "Treatment") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
plotFFBH_a

free(plotNAR_a + ggtitle("NAR")) +
  plotPAR_a + ggtitle("PAR") +
  plotFFLC1_a + ggtitle("FFLC1") +
  plotFFLI1_a + ggtitle("FFLI1") +
  plotFFBH_a + ggtitle("FFBH") +
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "bottom")

ggsave("plt_a_eg.png", device = png, bg = "white", width = 12, height = 7)
