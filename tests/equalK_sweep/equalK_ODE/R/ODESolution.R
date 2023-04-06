library(dplyr)
library(deSolve)
library(DescTools)

args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])

d_new <- readRDS("/g/data/ht96/nb9894/equalK_sweep/d_new_adapted.RDS")

# Get appropriate subset - 42 timepoints
d_new <- d_new[((run-1)*42+1):(run*42),]

# Setup NAR
nar <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * X^nXZ/(KXZ^nXZ + X^nXZ) * (KZ^nZ/(KZ^nZ+Z^nZ)) - aZ*Z
    list(c(dZ))
  })
}

state <- c(Z = 0)
times <- seq(0, 10, by = 0.1)
len_row <- nrow(d_new) * length(times)

solution <- data.frame(gen = integer(len_row),
                       seed = integer(len_row),
                       modelindex = integer(len_row),
                       nloci = integer(len_row), 
                       sigma = integer(len_row), 
                       phenomean = integer(len_row), 
                       time = integer(len_row),
                       X = integer(len_row), 
                       Z = integer(len_row))

for(i in 1:nrow(d_new)) {
  params <- c(Xstart = 1, Xstop = 6, 			
              aZ = d_new$aZ[i], bZ = d_new$bZ[i],
              nXZ = 8, nZ = 8, KXZ = d_new$KXZ[i],
              KZ = d_new$KZ[i])
  res <- ode(state, times, nar, params, "rk4") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(gen = d_new$gen[i], seed = d_new$seed[i], modelindex = d_new$modelindex[i],
           nloci = d_new$nloci[i], sigma = d_new$sigma[i],
           phenomean = d_new$phenomean[i],
           X = as.integer((time > Xstart & time <= Xstop))) %>%
    select(gen, seed, modelindex, nloci, sigma, phenomean, time, X, Z)
  
  solution[((i-1)*101+1):(i*101),] <- res
}

solution <- solution %>%
  group_by(gen, seed, modelindex) %>%
  mutate(diff_bottom = min(.$Z[.$X == 1]) - min(.$Z[.$X == 0]),
         diff_top = max(.$Z[.$X == 1]) - max(.$Z[.$X == 0]),
         len_left = max(.$Z[.$X == 0]) - min(.$Z[.$X == 0]),
         len_right = max(.$Z[.$X == 1]) - min(.$Z[.$X == 1]),
         bottom_left_angle = 90 - atan2(diff_bottom, 1) * (180 / pi),
         top_right_angle = atan2(1, diff_top) * (180 / pi),
         bottom_right_angle = 90 + atan2(diff_bottom, 1) * (180 / pi),
         top_left_angle = 90 + (90 - top_right_angle)
         )

write.table(solution, paste0("ode_", run, ".csv"), sep = ",", col.names = F, row.names = F)
