library(deSolve)
library(phaseR)
library(tidyverse)

# Function to calculate the derivatives of X and Z
nar_ode <- function(t, state, parameters) {
  X <- state[1]
  Z <- state[2]
  
  with(as.list(parameters), {
    if (t > t_stop_production) {
      r_X <- 0
    }
    dX_dt <- r_X - alpha_X * X
    dZ_dt <-
      beta_Z * (X ^ n_H) / (K_XZ ^ n_H + X ^ n_H) * (K_Z ^ n_H) / (K_Z ^ n_H + Z ^ n_H) - alpha_Z * Z
    return(list(c(dX_dt, dZ_dt)))
  })
}

# Set the parameter values
genotype <- c(r_X = 0.2347692,
              alpha_X = 0.16215139,
              beta_Z = 0.55608044,
              alpha_Z = 0.82151065)
n_H <- 0.4
K_XZ <- 0.5
K_Z <- 0.7

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
ggplot(data = ode_df, aes(x = X, y = Z)) +
  geom_path() +
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
