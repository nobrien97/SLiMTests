# plot molecular components vs phenotype w/ samples overlaid 
# as points connected by arrows (for time)
library(tidyverse)
library(gganimate)
library(latex2exp)
library(lattice)
library(latticeExtra)
library(paletteer)
library(reshape2)
library(scales)

# Load data
d_ode_phasemeasures <- readRDS("d_phaseDiagramMeasures.RDS")

d_ode_phasemeasures$id <- as.factor(paste(d_ode_phasemeasures$seed, d_ode_phasemeasures$modelindex, sep = "_"))

# Get range of each mol trait
molTraitRange <- list(aZ = quantile(d_ode_phasemeasures$aZ, c(0.01, 0.99)),
                      bZ = quantile(d_ode_phasemeasures$bZ, c(0.01, 0.99)),
                      KZ = quantile(d_ode_phasemeasures$KZ, c(0.01, 0.99)),
                      KXZ = quantile(d_ode_phasemeasures$KXZ, c(0.01, 0.99)))

# Filter data
d_ode_phasemeasures <- d_ode_phasemeasures %>% 
  filter(aZ >= molTraitRange$aZ[1] & aZ <= molTraitRange$aZ[2],
         bZ >= molTraitRange$bZ[1] & bZ <= molTraitRange$bZ[2],
         KZ >= molTraitRange$KZ[1] & KZ <= molTraitRange$KZ[2],
         KXZ >= molTraitRange$KXZ[1] & KXZ <= molTraitRange$KXZ[2])

# Get sequences between that range
molTraitSeq <- lapply(molTraitRange, function(x) {
  seq(from = x[1], to = x[2], length.out = 30)
})

# Generate data frame of mol trait combos
d_molTrait <- expand.grid(molTraitSeq)
write.table(d_molTrait, "./samples.csv", sep = ",", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# Use ODE Landscaper
width <- 0.05
optimum <- 2
threads <- 12
system(sprintf("ODELandscaper -i ./samples.csv -o ./landscape.csv -w %f -p %f -t %i",
               width, optimum, threads))

d_landscape <- read_csv("./landscape.csv", col_names = F, col_types = "d")
names(d_landscape) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")

d_landscape$pheno_rescaled <- ifelse(d_landscape$pheno > 3, 4, d_landscape$pheno)
d_ode_phasemeasures$pheno_rescaled <- ifelse(d_ode_phasemeasures$phenomean > 3, 4, 
                                             d_ode_phasemeasures$phenomean)

# Plot the landscape
ggplot(d_landscape, aes(x = aZbZ, y = KZKXZ, fill = pheno_rescaled)) +
  geom_tile() +
  geom_segment(data = d_ode_phasemeasures %>% filter(id %in% sample(.$id, 1)), 
               mapping = aes(x = lag(aZbZ), y = lag(KZKXZ), xend = aZbZ, yend = KZKXZ, 
                                        group = id), 
               arrow = arrow(length = unit(0.25, "cm")), color = "white") + 
  scale_fill_gradientn(colors = c("#4B0055", "#FDE333", "#009F94"), values = rescale(c(0, 2, 4)),
                         limits = c(0, 4),
                       labels = c(seq(0, 3, 1), TeX("$\\geq 3$")),
                       breaks = seq(0, 4, 1)
                       ) +
  labs(x = TeX("$\\alpha_Z + \\beta_Z$"), y = TeX("$K_Z + K_{XZ}$"), fill = "Phenotype (Z)") +
  theme_bw()

library(GGally)
library(ggarrow)

# Line width for time - scale gens
d_ode_phasemeasures$gen_width <- rescale(d_ode_phasemeasures$gen, to = c(0.001, 1))

sampled_id <- sample(d_ode_phasemeasures$id, 1)

cc <- paletteer_c("ggthemes::Orange-Blue Diverging", 3)
cc <- c(cc[1], cc[3], cc[2])

d_ode_phasemeasures2 <- d_ode_phasemeasures
d_ode_phasemeasures <- d_ode_phasemeasures %>% filter(id %in% sampled_id)

aZbZScatter <- ggplot(d_ode_phasemeasures, aes(x = aZ, y = bZ, colour = pheno_rescaled)) +
  geom_point() +
  geom_arrow_segment(data = d_ode_phasemeasures %>% filter(id %in% sampled_id), 
               mapping = aes(x = lag(aZ), y = lag(bZ), xend = aZ, yend = bZ, 
                             group = id, linewidth_head = gen_width, 
                             linewidth_fins = gen_width * 0.8), color = "black", 
               arrow_head = arrow_head_line()) +
  scale_linewidth(guide = "none") +
  scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                       limits = c(0, 4),
                       labels = c(seq(0, 2, 1), "3"),
                       breaks = seq(0, 3, 1)
  ) +
  labs(x = TeX("$\\alpha_Z$"), y = TeX("$\\beta_Z$"), colour = "Phenotype (Z)") +
  theme_bw() + 
  theme(legend.position = "bottom") +
  guides(colour=guide_colourbar(barwidth=20))

aZKZScatter <- ggplot(d_ode_phasemeasures, aes(x = aZ, y = KZ, colour = pheno_rescaled)) +
  geom_point() +
  geom_arrow_segment(data = d_ode_phasemeasures %>% filter(id %in% sampled_id), 
               mapping = aes(x = lag(aZ), y = lag(KZ), xend = aZ, yend = KZ, 
                             group = id, linewidth_head = gen_width, 
                             linewidth_fins = gen_width * 0.8), color = "black", 
               arrow_head = arrow_head_line()) +
  scale_linewidth(guide = "none") +
  scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                         limits = c(0, 4),
                         labels = c(seq(0, 2, 1), "3"),
                         breaks = seq(0, 3, 1)
                         ) +
  labs(x = TeX("$\\alpha_Z$"), y = TeX("$K_Z$"), colour = "Phenotype (Z)") +
  theme_bw()

aZKXZScatter <- ggplot(d_ode_phasemeasures, aes(x = aZ, y = KXZ, colour = pheno_rescaled)) +
  geom_point() +
  geom_arrow_segment(data = d_ode_phasemeasures %>% filter(id %in% sampled_id), 
               mapping = aes(x = lag(aZ), y = lag(KXZ), xend = aZ, yend = KXZ, 
                             group = id, linewidth_head = gen_width, 
                             linewidth_fins = gen_width * 0.8), color = "black", 
               arrow_head = arrow_head_line()) +
  scale_linewidth(guide = "none") +
  scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                         limits = c(0, 4),
                         labels = c(seq(0, 2, 1), "3"),
                         breaks = seq(0, 3, 1)
  ) +
  labs(x = TeX("$\\alpha_Z$"), y = TeX("$K_{XZ}$"), colour = "Phenotype (Z)") +
  theme_bw()

bZKZScatter <- ggplot(d_ode_phasemeasures, aes(x = bZ, y = KZ, colour = pheno_rescaled)) +
  geom_point() +
  geom_arrow_segment(data = d_ode_phasemeasures %>% filter(id %in% sampled_id), 
               mapping = aes(x = lag(bZ), y = lag(KZ), xend = bZ, yend = KZ, 
                             group = id, linewidth_head = gen_width, 
                             linewidth_fins = gen_width * 0.8), color = "black", 
               arrow_head = arrow_head_line()) + 
  scale_linewidth(guide = "none") +
  scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                         limits = c(0, 4),
                         labels = c(seq(0, 2, 1), "3"),
                         breaks = seq(0, 3, 1)
  ) +
  labs(x = TeX("$\\beta_Z$"), y = TeX("$K_Z$"), colour = "Phenotype (Z)") +
  theme_bw()

bZKXZScatter <- ggplot(d_ode_phasemeasures, aes(x = bZ, y = KXZ, colour = pheno_rescaled)) +
  geom_point() +
  geom_arrow_segment(data = d_ode_phasemeasures %>% filter(id %in% sampled_id), 
               mapping = aes(x = lag(bZ), y = lag(KXZ), xend = bZ, yend = KXZ, 
                             group = id, linewidth_head = gen_width, 
                             linewidth_fins = gen_width * 0.8), color = "black", 
               arrow_head = arrow_head_line()) + 
  scale_linewidth(guide = "none") +
  scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                         limits = c(0, 4),
                         labels = c(seq(0, 2, 1), "3"),
                         breaks = seq(0, 3, 1)
  ) +
  labs(x = TeX("$\\beta_Z$"), y = TeX("$K_{XZ}$"), colour = "Phenotype (Z)") +
  theme_bw()

KZKXZScatter <- ggplot(d_ode_phasemeasures, aes(x = KZ, y = KXZ, colour = pheno_rescaled)) +
  geom_point() +
  geom_arrow_segment(data = d_ode_phasemeasures %>% filter(id %in% sampled_id), 
               mapping = aes(x = lag(KZ), y = lag(KXZ), xend = KZ, yend = KXZ, 
                             group = id, linewidth_head = gen_width, 
                             linewidth_fins = gen_width * 0.8), color = "black", 
               arrow_head = arrow_head_line()) +
  scale_linewidth(guide = "none") +
  scale_colour_gradientn(colors = cc, values = rescale(c(0, 2, 4)),
                         limits = c(0, 4),
                         labels = c(seq(0, 2, 1), "3"),
                         breaks = seq(0, 3, 1)
  ) +
  labs(x = TeX("$K_Z$"), y = TeX("$K_{XZ}$"), colour = "Phenotype (Z)") +
  theme_bw()

plot_list <- list(ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = aZ)), 
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = aZ, y = bZ)),
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = aZ, y = KZ)),
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = aZ, y = KXZ)),
                  
                  aZbZScatter, 
                  ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = bZ)), 
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = bZ, y = KZ)),
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = bZ, y = KXZ)),
                  
                  aZKZScatter,
                  bZKZScatter,
                  ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = KZ)), 
                  ggally_cor(d_ode_phasemeasures, mapping = aes(x = KZ, y = KXZ)),
                  
                  aZKXZScatter,
                  bZKXZScatter,
                  KZKXZScatter,
                  ggally_densityDiag(d_ode_phasemeasures, mapping = aes(x = KXZ)))

xlabs <- c(TeX("$\\alpha_Z$", output = "character"), TeX("$\\beta_Z$", output = "character"), 
           TeX("$K_Z$", output = "character"), TeX("$K_{XZ}$", output = "character"))

ggmatrix(plot_list, nrow = 4, ncol = 4, xAxisLabels = xlabs, yAxisLabels = xlabs,
         progress = T, byrow = T, labeller = "label_parsed", 
         legend = grab_legend(plot_list[[5]])) + theme_bw() + 
  theme(legend.position = "bottom") -> pair_mat

ggsave("molTrait_landscape_singleWalk.png", pair_mat, width = 11, height = 8)

ggpairs(d_ode_phasemeasures, columns = c(6:9))


