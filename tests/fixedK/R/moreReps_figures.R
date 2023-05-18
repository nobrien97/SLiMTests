library(tidyverse)
library(paletteer)
library(latex2exp)
library(GGally)
library(cowplot)
library(ggnewscale)
library(ggalt)
library(ggarrow)

# Functions
se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}


data_path <- "/mnt/d/SLiMTests/tests/fixedK/moreReps/"

# load data
d_qg <- read.table(paste0(data_path, "slim_qg.csv"), header = F, 
                 sep = ",", colClasses = c("integer", "factor", "factor", 
                                           rep("numeric", times = 12)), 
                 col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                               "phenomean", "phenovar", "dist", "w", "deltaPheno",
                               "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                 fill = T)

d_qg %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg

d_qg %>% group_by(modelindex) %>%
  summarise(pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted))

d_adapted <- d_qg %>% filter(isAdapted)

d_adapted %>% 
  group_by(gen, modelindex) %>%
  summarise(meanPheno = mean(phenomean),
            CIPheno = CI(phenomean)) -> d_adapted_sum

d_muts <- read.table(paste0(data_path, "slim_muts.csv"), header = F, 
                     sep = ",", colClasses = c("integer", "factor", "factor", 
                                               "factor", rep("integer", times = 4),
                                               rep("numeric", times = 3),
                                               rep("integer", times = 2)), 
                     col.names = c("gen", "seed", "modelindex", "mutType", "mutID",
                                   "pos", "constraint", "originGen", "value", "chi",
                                   "Freq", "Count", "fixGen"), 
                     fill = T)

d_fix <- d_muts %>%
  filter(Freq == 1) %>%
  group_by(seed, modelindex, mutType) %>%
  distinct(mutID, .keep_all = T) 

d_fix_adapted <- d_fix %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_adapted$seed, d_adapted$modelindex))

# Fig 2 - phenotype mean
ggplot(d_adapted_sum %>% filter(gen > 49000) %>% mutate(gen = gen - 50000), 
       aes(x = gen, y = meanPheno, colour = modelindex)) +
  geom_line(size = 0.7) +
  geom_ribbon(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno, 
                  fill = modelindex), alpha = 0.2, colour = NA
              ) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  scale_fill_paletteer_d("ggsci::nrc_npg", guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean population phenotype",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenomean
plt_phenomean
ggsave("phenomean.png", plt_phenomean, png)
# Fig 3 - effect sizes
d_fix_adapted$fixTime <- d_fix_adapted$gen - d_fix_adapted$originGen

## Additive
ggplot(d_fix_adapted %>% filter(modelindex == 1), 
       aes(x = value)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Phenotypic effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_phenotype_additive
dfe_phenotype_additive

ggsave("dfe_phenotype_additive.png", dfe_phenotype_additive, png)

# Get fitness effect by subtracting fitness
d_fix_add <- d_fix_adapted %>% filter(modelindex == 1, gen >= 50000)

d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_add$gen, d_fix_add$seed, d_fix_add$modelindex)) %>%
  select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_add <- inner_join(d_fix_add, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

calcAddFitness <- function(phenotypes, optimum, width) {
  dists <- (phenotypes - optimum)^2
  return(exp(-(dists * width)))
}

d_absenceFitness <- d_fix_add %>% mutate(phenomean = phenomean - value)
d_absenceFitness$absenceW <- calcAddFitness(d_absenceFitness$phenomean, 2, 0.05)

d_fix_add$avFit <- d_fix_add$w - d_absenceFitness$absenceW

ggplot(d_fix_add %>% filter(modelindex == 1), 
       aes(x = avFit)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Fitness effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_fitness_additive
dfe_fitness_additive

ggsave("dfe_fitness_additive.png", dfe_fitness_additive, png)


runLandscaper <- function(df_path, output, width, optimum, threads) {
  system(sprintf("ODELandscaper -i %s -o ./%s -w %f -p %f -t %i",
                 df_path, output, width, optimum, threads))
  result <- read_csv(paste0("./", output), col_names = F, col_types = "d")
  names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  return(result)
}

# On average fitness effect is negative - so mutations require specific
# genetic backgrounds to be positive
# Need to measure fitness effect relative to actual background then rather than 
# average across the entire range experienced by every population
# so need to get the moltrait values at a fixed effect's given timepoint, and 
# subtract the fixed effect from it to measure the effect in context of the 
# population's background - can use mean pop phenotype and fitness to measure the
# background
d_fix_nar <- d_fix_adapted %>% filter(modelindex == 2, gen >= 50000)

# First get matched mol trait data for generations where we have fixations
d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_nar$gen, d_fix_nar$seed, d_fix_nar$modelindex)) %>%
  select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_nar2 <- inner_join(d_fix_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

d_fix_aZ <- d_fix_nar2 %>% filter(mutType == 3)
d_fix_bZ <- d_fix_nar2 %>% filter(mutType == 4)

# Calculate the mean phenotypes with the sampled mean aZ/bZ values
write.table(d_fix_aZ %>% ungroup() %>% select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_aZ <- runLandscaper("d_grid_aZ.csv", "data_popfx_aZ.csv", 0.05, 2, 8)

write.table(d_fix_bZ %>% ungroup() %>% select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_bZ <- runLandscaper("d_grid_bZ.csv", "data_popfx_bZ.csv", 0.05, 2, 8)


# Calculate the phenotypes when we take away the fixed effect in question from aZ/bZ
d_fix_aZ_diff <- d_fix_aZ
d_fix_aZ_diff$aZ <- exp(log(d_fix_aZ$aZ) - d_fix_aZ$value)

d_fix_bZ_diff <- d_fix_bZ
d_fix_bZ_diff$bZ <- exp(log(d_fix_bZ$bZ) - d_fix_bZ$value)

write.table(d_fix_aZ_diff %>% ungroup() %>% select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_aZ_diff <- runLandscaper("d_grid_aZ_diff.csv", "data_popfx_aZ_diff.csv", 0.05, 2, 8)

write.table(d_fix_bZ_diff %>% ungroup() %>% select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_bZ_diff <- runLandscaper("d_grid_bZ_diff.csv", "data_popfx_bZ_diff.csv", 0.05, 2, 8)

# Get the effect size by taking away the phenotype missing that fixation
d_fix_aZ$avFX <- d_fix_aZ$phenomean - d_popfx_aZ_diff$pheno
d_fix_aZ$avFit <- d_fix_aZ$w - d_popfx_aZ_diff$fitness

d_fix_bZ$avFX <- d_fix_bZ$phenomean - d_popfx_bZ_diff$pheno
d_fix_bZ$avFit <- d_fix_bZ$w - d_popfx_bZ_diff$fitness

d_fix_nar <- rbind(d_fix_aZ, d_fix_bZ)

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$")
)


ggplot(d_fix_nar, 
       aes(x = avFX, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Phenotypic effect", y = "Density",
       colour = "Molecular component") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5), 
        legend.position = "bottom") -> dfe_phenotype_nar
dfe_phenotype_nar
ggsave("dfe_phenotype_nar.png", dfe_phenotype_nar, png)

ggplot(d_fix_nar, 
       aes(x = avFit, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Fitness effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5),
        legend.position = "none") -> dfe_fitness_nar
dfe_fitness_nar
ggsave("dfe_fitness_nar.png", dfe_fitness_nar, png)

# Combine the above figures
leg <- get_legend(dfe_phenotype_nar)

fig3 <- plot_grid(dfe_phenotype_additive,
          dfe_phenotype_nar + theme(legend.position = "none"),
          dfe_fitness_additive,
          dfe_fitness_nar + theme(legend.position = "none"),
          ncol = 2, nrow = 2, labels = "AUTO")

fig3 <- plot_grid(fig3, get_legend(dfe_phenotype_nar), ncol = 1, rel_heights = c(1, 0.1))
fig3
ggsave("dfe_combined.png", fig3, png, bg = "white")

# Fig 4: Plot adaptive walks
cc <- paletteer_c("grDevices::Burg", 3)
cc2 <- paletteer_c("grDevices::Blues", 50, -1)

plotPairwiseScatter <- function(dat, x, y, labels) {
  GRID_RES <- 200
  d_grid <- expand.grid(seq(from = min(dat[[x]]), 
                                   to = max(dat[[x]]), length.out = GRID_RES), 
                        seq(from = min(dat[[y]]), 
                            to = max(dat[[y]]), length.out = GRID_RES), 1, 1)
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  mp <- sum(range(d_landscape$fitness))/2
  
  ggplot(dat %>% filter(modelindex == 2), aes(x = .data[[x]], y = .data[[y]])) +
    geom_raster(data = d_landscape, mapping = 
                     aes(x = .data[[x]], y = .data[[y]], fill = fitness)) +
    geom_point() +
    geom_encircle(colour = "black", mapping = aes(group = seed)) +
    scale_fill_gradient2(low = cc[1], mid = cc[2], high = cc[3], midpoint = mp) +
    labs(fill = "Relative fitness (w)") +
    
    
    new_scale_colour() +
    geom_arrow_segment(data = dat %>% filter(seed %in% sampled_seed),
                       mapping = aes(x = lag(.data[[x]]), y = lag(.data[[y]]),
                                     xend = .data[[x]], yend = .data[[y]],
                                     group = seed, linewidth_head = gen_width,
                                     linewidth_fins = gen_width * 0.8,
                                     colour = gen_width),
                       arrow_head = arrow_head_line()) +
    scale_colour_gradientn(colors = cc2, labels = scales::comma(c(0, 0.25*10000, 0.5*10000,
                                                    0.75*10000, 10000))) +
    scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2], colour = "Generation") +
    theme_bw() + 
    theme(legend.position = "bottom", text = element_text(size = 14)) +
    guides(colour=guide_colourbar(barwidth=10),
           fill = guide_colourbar(barwidth = 10))
}

d_qg_adapting <- d_adapted %>% filter(gen >= 50000)
d_qg_adapting$gen_width <- scales::rescale(d_qg_adapting$gen, to = c(0.1, 1))


sampled_seed <- sample(d_qg_adapting[d_qg_adapting$modelindex == 2,]$seed, 3)
walk <- plotPairwiseScatter(d_qg_adapting %>% filter(modelindex == 2, seed %in% sampled_seed), 
                    "aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
walk
ggsave("example_adaptiveWalk.png", walk, png)
