library(tidyverse)
library(gganimate)
library(deSolve)
library(DescTools)
library(future)
library(doParallel)
library(foreach)



d_new <- readRDS(paste0(local_path, "checkpoint/d_qg.RDS"))

d_new %>% mutate(id = as_factor(paste(seed, modelindex, sep = "_"))) -> d_new

d_new$nloci_cat <- cut(d_new$nloci,
                       breaks=c(-Inf, 10, 50, 250, Inf))

d_new$sigma_cat <- cut(d_new$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))


d_new %>%
  group_by(seed, nloci, sigma) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  mutate(phenoCI = qnorm(0.975) * phenovar/sqrt(5000),
         seed = as_factor(seed)) %>%
  ungroup() -> d_new_adapted

seed <- sample(1:.Machine$integer.max, 1)
set.seed(seed)
# 563655642
set.seed(563655642)

sampled_seeds <- d_new_adapted %>% filter(gen > 49500, phenomean < 5) %>%
  group_by(nloci, sigma) %>% 
  select(nloci, sigma, seed, modelindex, id) %>%
  sample_n(1)

sample_num <- 1
d_sample <- d_new %>% filter(gen > 49000, seed %in% sampled_seeds$seed[sample_num],
                             modelindex %in% sampled_seeds$modelindex[sample_num])
                            

nar <-function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * X^nXZ/(KXZ^nXZ + X^nXZ) * (KZ^nZ/(KZ^nZ+Z^nZ)) - aZ*Z
    list(c(dZ))
  })
}

#Values for the ODE
state <- c(Z = 0)
times <- seq(0, 10, by = 0.1)
len_row <- nrow(d_new_adapted)

solution <- data.frame(gen = integer(len_row), 
                       nloci = integer(len_row), 
                       sigma = integer(len_row), 
                       nloci_cat = integer(len_row), 
                       sigma_cat = integer(len_row), 
                       id = integer(len_row), 
                       time = integer(len_row), 
                       Z = integer(len_row))

cl <- parallel::makeForkCluster(future::availableCores())
doParallel::registerDoParallel(cl)

foreach(i = 1:nrow(d_new_adapted)) %dopar% {
  params <- c(Xstart = 1, Xstop = 6, 			
              aZ = d_new_adapted$aZ[i], bZ = d_new_adapted$bZ[i],
              nXZ = 8, nZ = 8, KXZ = d_new_adapted$KXZ[i],
              KZ = d_new_adapted$KZ[i])
  res <- ode(state, times, nar, params, "rk4") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0),
           gen = d_new_adapted$gen[i],
           nloci = d_new_adapted$nloci[i], sigma = d_new_adapted$sigma[i],
           nloci_cat = d_new_adapted$nloci_cat[i], sigma_cat = d_new_adapted$sigma_cat[i],
           id = d_new_adapted$id[i]) %>%
    select(gen, nloci, sigma, nloci_cat, sigma_cat, id, time, Z)
  
  solution[((i-1)*140+1):(i*140),] <- res
}

parallel::stopCluster(cl)

saveRDS(solution, "ode_solutions.RDS")

ggplot(solution %>% mutate(gen = gen - 50000), aes(x = time, y = Z)) +
  geom_line() +
  labs(x = "Developmental time", y = "Z expression") +
  theme_bw() +
  theme(text = element_text(size=20)) -> plt_ode

plt_ode <- plt_ode + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_ode, nframes = 42, duration = 10, width = 720, height = 720, 
        renderer = ffmpeg_renderer())
anim_save("ode_time.mp4", last_animation())
