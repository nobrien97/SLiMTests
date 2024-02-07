

# Pairwise epistasis
breaks <- seq(min(d_add_e$gen - 50000), max(d_add_e$gen - 50000), by = 1000)
d_add_e$gen_group <- breaks[findInterval(d_add_e$gen - 50000, breaks, rightmost.closed = TRUE)]

ggplot(d_add_e, aes(x = meanEP, y = as.factor(gen_group))) +
  geom_density_ridges() +
  labs(y = "Generations post-optimum shift", x = "Mean pairwise trait epistasis") +
  theme_bw()

breaks <- seq(min(d_net_e$gen - 50000), max(d_net_e$gen - 50000), by = 1000)
d_net_e$gen_group <- breaks[findInterval(d_net_e$gen - 50000, breaks, rightmost.closed = TRUE)]

ggplot(d_net_e, aes(x = meanEW, y = as.factor(gen_group), fill = mutType_ab)) +
  geom_density_ridges(alpha = 0.4) +
  labs(y = "Generations post-optimum shift", x = "Mean pairwise trait epistasis") +
  theme_bw()

breaks <- seq(min(d_epistasis$gen - 50000), max(d_epistasis$gen - 50000), by = 1000)
d_epistasis$gen_group <- breaks[findInterval(d_epistasis$gen - 50000, breaks, rightmost.closed = TRUE)]

ggplot(d_epistasis, aes(x = meanEW, y = as.factor(gen_group), fill = mutType_ab)) +
  geom_density_ridges(alpha = 0.4) +
  labs(y = "Generations post-optimum shift", x = "Mean pairwise trait epistasis") +
  theme_bw()
