# Plot proportion of beneficial/neutral/deleterious mutations
library(tidyverse)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)

source("helperFn.R")


# Load in RDS files
PATH_FX_ORTH <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/orthSel/calcMutationStats/"
PATH_FX_PAR <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/parallelSel/calcMutationStats/"
PATH_FX_RAND <- "/mnt/d/SLiMTests/tests/newMotifs/paper1/randomisedStartsM/calcMutationStats/"

d_fx_orth <- readRDS(paste0(PATH_FX_ORTH, "d_fx_sum.RDS"))
d_fx_par <- readRDS(paste0(PATH_FX_PAR, "d_fx_sum.RDS"))
d_fx_rand <- readRDS(paste0(PATH_FX_RAND, "d_fx_sum.RDS"))

d_fx_orth$dataset <- "Orthogonal"
d_fx_par$dataset <- "Parallel"
d_fx_rand$dataset <- "Randomised"

d_fx <- rbind(d_fx_rand, d_fx_orth, d_fx_par)

d_fx$dataset <- factor(d_fx$dataset, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised"))

# Create secondary axis
d_fx <- d_fx %>%
  mutate(gen_k = (gen - 50000) / 1000)  # Convert to "k" units

# Assign numeric x-axis positions
gen_levels <- sort(unique(d_fx$gen_k))
gap <- 1  # space between model blocks

d_fx <- d_fx %>%
  group_by(model, dataset) %>%
  mutate(gen_index = as.integer(factor(gen_k, levels = gen_levels))) %>%
  ungroup() %>%
  mutate(model_index = as.integer(model),
         dataset_index = as.integer(dataset),
         x_numeric = (model_index - 1) * (length(gen_levels) + gap) + gen_index)

gen_labels <- d_fx %>% filter(gen_k != 0) %>%
  group_by(model, dataset, gen_k) %>%
  summarise(x = mean(x_numeric), .groups = "drop")

fx_model_labels <- d_fx %>% filter(gen_k != 0) %>%
  group_by(model, dataset, isAdapted) %>%
  summarise(x = mean(x_numeric), .groups = "drop")

ggplot(d_fx %>% filter(gen_k > 0, isAdapted == T), 
       aes(x = x_numeric, y = meanProp, fill = mutClass)) +
  facet_nested("Trait/selection alignment" + dataset ~ .) +
  geom_bar(stat = "identity", width = 1) +
  coord_cartesian(ylim = c(0, 1), expand = F, clip = "off") +
  scale_fill_manual(values = paletteer_d("ggprism::viridis", 
                                         5, direction = 1)[c(5, 1, 3)],
                    labels = c("Beneficial", "Deleterious", "Neutral")) +
  labs(x = TeX("Generations post-optimum shift ($x10^3$) / Model"), 
       y = "Mean proportion of mutations",
       fill = "Mutation class") +
  scale_x_continuous(breaks = gen_labels$x, labels = gen_labels$gen_k) +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  geom_text(data = fx_model_labels %>% filter(isAdapted == T), aes(x = x, y = -0.15 * 2.5,
                                                             label = model), inherit.aes = F) +
  geom_segment(data = fx_model_labels %>% filter(isAdapted == T),
               aes(x = x - 5,
                   xend = x + 5,
                   y = -0.085 * 2.5, yend = -0.085 * 2.5),
               inherit.aes = FALSE,
               color = "black", size = 0.4) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines"),
        legend.position = "bottom",
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = 28))) -> plt_prop
plt_prop
ggsave("plt_propFX.png", plt_prop, device = png, bg = "white",
       width = 8, height = 5)


# Distributions
d_fx_dens_orth <- readRDS(paste0(PATH_FX_ORTH, "d_fx_density.RDS"))
d_fx_dens_par <- readRDS(paste0(PATH_FX_PAR, "d_fx_density.RDS"))
d_fx_dens_rand <- readRDS(paste0(PATH_FX_RAND, "d_fx_density.RDS"))

d_fx_dens_orth$dataset <- "Orthogonal"
d_fx_dens_par$dataset <- "Parallel"
d_fx_dens_rand$dataset <- "Randomised"

d_fx_dens <- rbind(d_fx_dens_rand, d_fx_dens_orth, d_fx_dens_par)

d_fx_dens$dataset <- factor(d_fx_dens$dataset, levels = c("Parallel",
                                                "Orthogonal",
                                                "Randomised"))

# Plot
ggplot(d_fx_dens %>% filter(isAdapted == T, gen == 50000 | gen == 55000 | gen == 59000) %>%
         group_by(model, dataset) %>%
         mutate(dens = dens / sum(dens)),
       aes(x = s, y = gen - 50000, height = dens, 
           group = gen - 50000, colour = model, fill = model)) +
  facet_nested("Model" + model ~ "Trait/selection alignment" + dataset) +
  geom_ridgeline(alpha = 0.2, scale = 10000, show.legend = F) +
  coord_cartesian(xlim = c(-1, 0.5)) +
  scale_y_continuous(breaks = seq(0, 10000, by = 2000), labels = scales::comma) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "Selection coefficient", y = "Generations post-optimum shift") +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plt_fx_dist.png", device = png, width = 6, height = 9.5)

# Peak of each density curve
print(xtable::xtable(d_fx_dens %>% filter(isAdapted == T, gen == 59000) %>%
  group_by(model, dataset) %>%
  mutate(dens = dens / sum(dens)) %>%
  slice_max(dens, n = 1) %>%
    select(model, dataset, s, dens), digits = 3), include.rownames = F)


