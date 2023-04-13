# KZ seems to drift to really large values, and even with small locisigma,
# it consistently larger than other molecular trait values
# There could be diminishing returns with modifying it - at a certain point, 
# the removal becomes so efficient that it doesn't matter what the value is

# To test this, I plot KZ values vs Z for fixed levels of the other molecular
# traits


library(tidyverse)
library(deSolve)
library(DescTools)
NARODE <-function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * X^nXZ/(KXZ^nXZ + X^nXZ) * (KZ^nZ/(KZ^nZ+Z^nZ)) - aZ*Z
    list(c(dZ))
  })
}

state <- c(Z = 0)
times <- seq(0, 10, by = 0.1)

KZ_vals = 10^c(seq(from = -5, to = 1, by = 1))
out <- data.frame(
  KZ = numeric(length(KZ_vals)),
  Z = numeric(length(KZ_vals)))

for (i in seq_along(KZ_vals)) {
  params <- c(Xstart = 1, Xstop = 6, 			
              aZ = 3.62633, bZ = 0.706658,
              nXZ = 8, nZ = 8, KXZ = 0.930639,
              KZ = KZ_vals[i])
  solution <- ode(state, times, Freya, params, "rk4") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
    select(time, X, Z)
  
  out[i,] <- cbind(KZ_vals[i],
               AUC(solution$time, solution$Z, absolutearea = T, method = "trapezoid"))
}

plot(out$KZ, out$Z, type = "b")
