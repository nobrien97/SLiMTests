library(tidyverse)
library(gganimate)
library(gghighlight)

setwd("/g/data/ht96/nb9894/equalK_sweep/equalK_ODE/indPheno/")

# plot ODE solution
d_ode <- read_csv("/g/data/ht96/nb9894/equalK_sweep/equalK_ODE/indPheno/ode_solutions_inds.csv", 
                  col_names = F)

colnames(d_ode) <- c("gen", "seed", "modelindex", "nloci", "sigma", "ind", "phenotype", "time", "Z")


d_ode$nloci_cat <- cut(d_ode$nloci,
                       breaks=c(-Inf, 10, 50, 250, Inf))

d_ode$sigma_cat <- cut(d_ode$sigma,
                       breaks=c(-Inf, 0.05, 0.5, 1, Inf))

d_ode <- d_ode %>% mutate(id = paste(seed, modelindex, sep = "_"))

seed <- sample(1:.Machine$integer.max, 1)
set.seed(seed)
# 1158538046
set.seed(1121489353)

# sample 1 individual and seed from each group
sampled_seeds <- d_ode %>%
  group_by(nloci_cat, sigma_cat) %>% 
  select(nloci, sigma, seed, modelindex, id, ind) %>%
  sample_n(1)

d_sample <- d_ode %>% filter(seed %in% sampled_seeds$seed[1], 
                             modelindex %in% sampled_seeds$modelindex[1]) %>%
                        mutate(idx = interaction(ind, id))

highlight_ids <- sampled_seeds %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(ind, id) %>%
  sample_n(1) %>% mutate(idx = interaction(ind, id))

cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")


plt_ode <- ggplot(d_sample %>%
                    mutate(gen = gen - 50000, time = time/10),
                  aes(x = time, y = Z, colour = idx, group = idx)) +
  geom_line() +
  scale_colour_manual(values = rep(cc_ibm[3], 11), guide = "none") +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  gghighlight(idx %in% highlight_ids$idx,
              calculate_per_facet = T, use_direct_label = F) +
  labs(x = "Developmental time", y = "Mean Z expression") +
  theme_bw() +
  theme(text = element_text(size=36), panel.spacing.x = unit(1, "lines"))

saveRDS(plt_ode, "testpltode.RDS")
plt_ode <- readRDS("testpltode.RDS")

plt_ode <- plt_ode + transition_states(gen, wrap = F) +
  labs(title = "Generation: {closest_state}")

anim <- animate(plt_ode, nframes = 42, duration = 10, width = 1920, height = 1920,
        renderer = ffmpeg_renderer())
anim_save("ode_time_test.mp4", anim)

# sample 10 from each nloci/sigma group and plot
sampled_seeds <- d_ode %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(nloci, sigma, seed, modelindex, id, ind) %>%
  sample_n(2)

d_sample <- d_ode %>% filter(seed %in% sampled_seeds$seed, 
                             modelindex %in% sampled_seeds$modelindex) %>%
  mutate(idx = interaction(ind, id))

highlight_ids <- sampled_seeds %>%
  group_by(nloci_cat, sigma_cat) %>%
  select(ind, id) %>%
  sample_n(1) %>% mutate(idx = interaction(ind, id))

cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")

# Since we have huge values, we're going to cap the y-axis so we ignore those
max_y <- quantile(d_sample$Z, probs = 0.99, na.rm = T) + 0.1


plt_ode_total <- ggplot(d_sample %>%
                          mutate(gen = gen - 50000, time = time/10),
                        aes(x = time, y = Z, colour = idx, group = idx)) +
  facet_grid(nloci_cat~sigma_cat) +
  geom_line() +
  scale_colour_manual(values = rep(cc_ibm[3], 11), guide = "none") +
  coord_cartesian(ylim = c(0, max_y)) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  gghighlight(idx %in% highlight_ids$idx,
              calculate_per_facet = T, use_direct_label = F) +
  labs(x = "Developmental time", y = "Z expression") +
  theme_bw() +
  theme(text = element_text(size=36), panel.spacing.x = unit(1, "lines"))

plt_ode_total <- plt_ode_total + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

anim <- animate(plt_ode_total, nframes = 42, duration = 10, width = 1920, height = 1920,
        renderer = ffmpeg_renderer())
anim_save("ode_time_total.mp4", anim)
