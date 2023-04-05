library(tidyverse)
library(gganimate)
library(gghighlight)

# plot ODE solution

d_ode <- read_csv("/mnt/d/SLiMTests/tests/equalK_sweep/equalK_ODE/data/ode_solutions.csv", 
                  col_names = F)

colnames(d_ode) <- c("gen", "seed", "modelindex", "nloci", "sigma", "phenomean", "time", "Z")


d_ode$nloci_cat <- cut(d_ode$nloci,
                       breaks=c(-Inf, 10, 50, 250, Inf))

d_ode$sigma_cat <- cut(d_ode$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))

d_ode <- d_ode %>% mutate(id = paste(seed, modelindex, sep = "_"))

seed <- sample(1:.Machine$integer.max, 1)
set.seed(seed)
# 18799215
set.seed(18799215)

sampled_seeds <- d_ode %>%
  group_by(nloci, sigma) %>% 
  select(nloci, sigma, seed, modelindex, id) %>%
  sample_n(1)

d_sample <- d_ode %>% filter(seed %in% sampled_seeds$seed[1], 
                                 modelindex %in% sampled_seeds$modelindex[1])


ggplot(d_sample %>% mutate(gen = gen - 50000), aes(x = time, y = Z)) +
  geom_line() +
  labs(x = "Developmental time", y = "Z expression") +
  theme_bw() +
  theme(text = element_text(size=20)) -> plt_ode

plt_ode <- plt_ode + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_ode, nframes = 42, duration = 10, width = 720, height = 720, 
        renderer = ffmpeg_renderer())
anim_save("ode_time.mp4", last_animation())

# sample 10 from each nloci/sigma group and plot

sampled_seeds <- d_ode %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(nloci, sigma, seed, modelindex, id) %>%
  sample_n(10)



d_sample <- d_ode %>% filter(seed %in% sampled_seeds$seed, 
                   modelindex %in% sampled_seeds$modelindex)

highlight_ids <- sampled_seeds %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(id) %>%
  sample_n(1)

cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")


plt_ode_total <- ggplot(d_sample %>%
                          mutate(gen = gen - 50000, time = time/10),
                        aes(x = time, y = Z, colour = id, group = id)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line(linewidth = 2) +
  scale_colour_manual(values = rep(cc_ibm[3], 11), guide = "none") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1), 
                     sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  gghighlight(id %in% highlight_ids$id,
              calculate_per_facet = T, use_direct_label = F) +
  labs(x = "Developmental time", y = "Mean Z expression") +
  theme_bw() +
  theme(text = element_text(size=36), panel.spacing.x = unit(1, "lines"))

plt_ode_total <- plt_ode_total + transition_states(gen, wrap = F) +
  labs(title = "Generation: {closest_state}")

animate(plt_ode_total, nframes = 42, duration = 10, width = 1920, height = 1920, 
        renderer = ffmpeg_renderer())
anim_save("ode_time_total.mp4", last_animation())
