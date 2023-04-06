##############################################################
# Author: Daniel Ortiz-Barrientos
# Goal: To understand the NAR model
##############################################################
# NAR
# Load libraries
library(ggplot2) 
library(deSolve) 
# Function
nar_ode <- function(t, state, parameters) { 
  X <- state[1] 
  Z <- state[2] 
  
  with(as.list(parameters), { 
    if (t > t_stop_production) { 
      r_X <- 0 
    } 
    
    dX_dt <- r_X - alpha_X * X 
    dZ_dt <- beta_Z * (X ^ n_H) / (K_XZ ^ n_H + X ^ n_H) * (K_Z ^ n_H) / (K_Z ^ n_H + Z ^ n_H) - alpha_Z * Z 
    return(list(c(dX_dt, dZ_dt))) 
  }) 
} 

# Set the time at which you want to stop the production of X:
t_stop_production <- 50

# Create the initial genotype:
genotype <- c(r_X = 1.0, 
              alpha_X = 1, 
              beta_Z = 0.2, 
              alpha_Z = 0.5) 
n_H <- 10
K_XZ <- 1
K_Z <- 1
initial_state <- c(X = 0, Z = 0) 
time_span <- seq(0, 120, by = 1) 

parameter_names <- 
  c("r_X", "alpha_X", "beta_Z", "alpha_Z", "n_H", "K_XZ", "K_Z") 
named_genotype <- 
  setNames(c(genotype, n_H, K_XZ, K_Z), parameter_names) 

nar_ode(0, c(X = 1, Z = 0), named_genotype) 

ode_result <- 
  ode( 
    y = initial_state, 
    times = time_span, 
    func = nar_ode, 
    parms = named_genotype 
  ) 
# Find the maximum values of X and Z:
max_X <- max(ode_result[, "X"]) 
max_Z <- max(ode_result[, "Z"]) 

# Determine the upper limit for the Y-axis by taking the maximum value of max_X and max_Z:
y_upper_limit <- max(max_X, max_Z) 

# Plot the results with the adaptive Y-axis:
plot( 
  ode_result[, "time"], 
  ode_result[, "X"], 
  type = "l", 
  xlab = "Time", 
  ylab = "Concentration", 
  main = "ODE Results", 
  ylim = c(0, y_upper_limit) 
) 
lines(ode_result[, "time"], ode_result[, "Z"], col = "blue") 
legend( 
  "topright", 
  legend = c("X", "Z"), 
  col = c("black", "blue"), 
  lty = 1
) 

# Function to calculate the derivatives of X and Z
nar_ode <- function(t, state, parameters) { 
  X <- state[1] 
  Z <- state[2] 
  
  with(as.list(parameters), { 
    if (t > t_stop_production) { 
      r_X <- 0 
    } 
    dX_dt <- r_X - alpha_X * X 
    dZ_dt <- beta_Z * (X ^ n_H) / (K_XZ ^ n_H + X ^ n_H) * (K_Z ^ n_H) / (K_Z ^ n_H + Z ^ n_H) - alpha_Z * Z 
    return(list(c(dX_dt, dZ_dt))) 
  }) 
} 
# Set the parameter values
genotype <- c(r_X = 1.0, 
              alpha_X = 1, 
              beta_Z = 0.2, 
              alpha_Z = 0.5) 
n_H <- 10
K_XZ <- 1
K_Z <- 1

# Set the time at which you want to stop the production of X
t_stop_production <- 70

# Set the initial values of X and Z
initial_state <- c(X = 0, Z = 0) 

# Set the time span for the simulation
time_span <- seq(0, 120, by = 1) 

# Set up the parameters for the ode solver
parameter_names <- c("r_X", "alpha_X", "beta_Z", "alpha_Z", "n_H", "K_XZ", "K_Z") 
named_genotype <- setNames(c(genotype, n_H, K_XZ, K_Z), parameter_names) 

# Solve the ODE system
ode_result <- ode(y = initial_state, times = time_span, func = nar_ode, parms = named_genotype) 
ode_df <- data.frame(time = ode_result[, "time"], X = ode_result[, "X"], Z = ode_result[, "Z"]) 

# Convert ode_result into a data frame
ode_df <- as.data.frame(ode_result) 
# Plot the phase diagram
library(dplyr) 
ggplot(data = ode_df, aes(x = X, y = Z)) + 
  geom_segment(aes(x = lag(X), y = lag(Z), xend = X, yend = Z), 
               arrow = arrow(length = unit(0.25, "cm")), color = "black") + 
  xlab("X") + 
  ylab("Z") + 
  ggtitle("NAR Phase Diagram") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        legend.position = "none")

#################################################################
# With step function for X
# Function to calculate the derivative of Z
nar_ode <- function(t, state, parameters) { 
  Z <- state[1]
  
  with(as.list(parameters), { 
    X <- (t > Xstart && t <= Xstop)
    dZ_dt <- beta_Z * (X ^ n_H) / (K_XZ ^ n_H + X ^ n_H) * (K_Z ^ n_H) / (K_Z ^ n_H + Z ^ n_H) - alpha_Z * Z 
    return(list(dZ_dt)) 
  }) 
} 
# Set the parameter values
genotype <- c( 
              beta_Z = 0.2, 
              alpha_Z = 0.5
              ) 
n_H <- 10
K_XZ <- 1
K_Z <- 1

# Set the time at which you want to stop the production of X
Xstart <- 10
Xstop <- 60

# Set the initial values of Z
initial_state <- c(Z = 0) 

# Set the time span for the simulation
time_span <- seq(0, 120, by = 1) 

# Set up the parameters for the ode solver
parameter_names <- c("beta_Z", "alpha_Z", "n_H", "K_XZ", "K_Z") 
named_genotype <- setNames(c(genotype, n_H, K_XZ, K_Z), parameter_names) 

# Solve the ODE system
ode_result <- ode(y = initial_state, times = time_span, func = nar_ode, parms = named_genotype) 
ode_df <- data.frame(time = ode_result[, "time"], 
                     X = as.integer((ode_result[, "time"] > Xstart & ode_result[, "time"] <= Xstop)), 
                     Z = ode_result[, "Z"]) 

# Plot the phase diagram
library(dplyr) 
ggplot(data = ode_df, aes(x = X, y = Z)) + 
  geom_segment(aes(x = lag(X), y = lag(Z), xend = X, yend = Z), 
               arrow = arrow(length = unit(0.25, "cm")), color = "black") + 
  xlab("X") + 
  ylab("Z") + 
  ggtitle("NAR Phase Diagram") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        legend.position = "none")


########################################################################
# Multiple Phase Diagrams
########################################################################
library(ggplot2)

parameters <- expand.grid(beta_Z = c(0.2, 0.5), alpha_Z = c(0.5, 1),
                          n_H = c(10), K_XZ = c(1), K_Z = c(1))
parameter_names <- c("beta_Z", "alpha_Z", "n_H", "K_XZ", "K_Z") 

# Set the time span for the simulation
time_span <- seq(0, 120, by = 1)

# Function to calculate the derivatives of X and Z
nar_ode <- function(t, state, parameters) {
  Z <- state[1]
  
  with(as.list(parameters), {
    X <- (t > Xstart && t <= Xstop)
    dZ_dt <-
      beta_Z * (X ^ n_H) / (K_XZ ^ n_H + X ^ n_H) * (K_Z ^ n_H) / (K_Z ^ n_H + Z ^ n_H) - alpha_Z * Z
    return(list(dZ_dt))
  })
}

# Set the time at which you want to stop the production of X
Xstart <- 10
Xstop <- 60

# Set the initial values of X and Z
initial_state <- c(Z = 0)

# Set the time span for the simulation
time_span <- seq(0, 120, by = 1)

# Create an empty list to store the results
ode_results <- vector("list", nrow(parameters))

# Run the simulation for each set of parameters
for (i in 1:nrow(parameters)) {
  named_genotype <- setNames(as.list(parameters[i, ]), parameter_names)
  ode_result <- ode(y = initial_state, times = time_span, func = nar_ode, parms = named_genotype)
  ode_results[[i]] <- data.frame(X = as.integer((ode_result[, "time"] > Xstart & ode_result[, "time"] <= Xstop)), 
                                 Z = ode_result[, "Z"])
}

# Combine the results into a single data frame
ode_df <- data.frame()

for (i in seq_along(ode_results)) {
  tmp_df <- data.frame(X = ode_results[[i]]$X,
                       Z = ode_results[[i]]$Z,
                       ParameterCombination = rep(i, times = length(ode_results[[i]]$X)))
  ode_df <- rbind(ode_df, tmp_df)
  
}

ode_df$ParameterCombination <- factor(ode_df$ParameterCombination)

# Phase diagram plot
library(viridis)

# Create a color vector using the viridis palette
color_vector <- viridis(length(ode_results), alpha = 0.1)

ode_df$ParameterCombination <- rep(paste("beta_Z =", parameters$beta_Z,
                                         "alpha_Z =", parameters$alpha_Z),
                                   each = length(time_span))

ggplot(data = ode_df, aes(x = X, y = Z)) +
  facet_wrap(.~ParameterCombination) +
  geom_segment(aes(x = lag(X), y = lag(Z), xend = X, yend = Z), 
               arrow = arrow(length = unit(0.25, "cm"))) + 
  xlab("X") +
  ylab("Z") +
  scale_color_viridis(name = "Parameter Combination", discrete = TRUE, begin = 0, end = 1, direction = 1) +
  ggtitle("Phase Diagram for Multiple Models") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        legend.key = element_rect(colour = "black", fill = "white"))
