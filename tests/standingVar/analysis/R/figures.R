

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

# dpdt
ggplot(d_dpdt,
       aes(x = dPdT, y = optPerc, fill = model)) +
  facet_wrap(vars(nloci)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Change in phenotype per generation ($\\frac{\\delta P}{\\delta t}$)"), 
       y = "Closeness to optimum (%)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

# SFS
ggplot(d_SFS,
       aes(x = freq, y = optPerc, height = after_stat(density), fill = model)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Allele frequency"), 
       y = "Closeness to optimum (%)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

