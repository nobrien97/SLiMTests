library(tidyverse)
## Network: need to get the average phenotype effect over a range of genetic backgrounds
## to do this, for a given alpha effect that was sampled, 
## we look at its effect on phenotype and fitness over a range of beta values
## each effect is a deviation from the standard rate (aZ = bZ = 1)
## need to measure the average effect on phenotype compared to this standard rate: just
## the difference between the phenotype values

# first need to generate the standard effect for alpha and beta 
# from the range of mutational effects - remember it needs to be exponentiated



runLandscaper <- function(df_path, output, width, optimum, threads) {
  system(sprintf("ODELandscaper -i %s -o ./%s -w %f -p %f -t %i",
                 df_path, output, width, optimum, threads))
  result <- read_csv(paste0("./", output), col_names = F, col_types = "d")
  names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  return(result)
}

d_fix_nar <- d_fix_adapted %>% filter(modelindex == 2, gen >= 50000)

mutRange <- d_fix_nar %>% 
  group_by(mutType) %>% 
  reframe(mutRange = range(value))

# Get the default average phenotype over the range of values for the other mol comp
GRID_RES <- 1000

d_aZ_defaultgrid <- expand.grid(0, seq(from = mutRange$mutRange[3], 
                                       to = mutRange$mutRange[4], length.out = GRID_RES), 0, 0)

d_aZ_defaultgrid <- d_aZ_defaultgrid %>% mutate(across(everything(), exp))

d_bZ_defaultgrid <- expand.grid(seq(from = mutRange$mutRange[1], 
                                    to = mutRange$mutRange[2], length.out = GRID_RES),
                                0, 0, 0)

d_bZ_defaultgrid <- d_bZ_defaultgrid %>% mutate(across(everything(), exp))

write.table(d_aZ_defaultgrid, "d_aZ_defaultgrid.csv", sep = ",", col.names = F, row.names = F)
write.table(d_bZ_defaultgrid, "d_bZ_defaultgrid.csv", sep = ",", col.names = F, row.names = F)


aZ_out <- runLandscaper("d_aZ_defaultgrid.csv", "default_aZ.csv", 0.05, 2, 8)
aZ_default_avg <- aZ_out %>% summarise(meanPheno = mean(pheno),
                                       meanFitness = mean(fitness))

bZ_out <- runLandscaper("d_bZ_defaultgrid.csv", "default_bZ.csv", 0.05, 2, 8)
bZ_default_avg <- bZ_out %>% summarise(meanPheno = mean(pheno),
                                       meanFitness = mean(fitness))

d_fix_aZ <- d_fix_nar[d_fix_nar$mutType == 3,]
d_fix_bZ <- d_fix_nar[d_fix_nar$mutType == 4,]

# For each alpha value, we need to calculate the average value by running the landscaper -
# we'll generate a list of all the effect size values to only run the landscaper once
d_aZ_datagrid <- expand.grid(seq(from = mutRange$mutRange[3], 
                                 to = mutRange$mutRange[4], length.out = GRID_RES), d_fix_aZ$value, 0, 0)
d_aZ_datagrid <- d_aZ_datagrid %>% mutate(across(everything(), exp)) %>% 
  select(c(2, 1, 3, 4))

d_bZ_datagrid <- expand.grid(seq(from = mutRange$mutRange[1], 
                                 to = mutRange$mutRange[2], length.out = GRID_RES), d_fix_bZ$value, 0, 0)
d_bZ_datagrid <- d_bZ_datagrid %>% mutate(across(everything(), exp))

write.table(d_aZ_datagrid, "d_aZ_grid.csv", sep = ",", col.names = F, row.names = F)
aZ_out <- runLandscaper("d_aZ_grid.csv", "data_aZ.csv", 0.05, 2, 8)

aZ_out <- aZ_out %>% mutate(aZ_value_id = rep(seq_along(d_fix_aZ$value), each = GRID_RES))

write.table(d_bZ_datagrid, "d_bZ_grid.csv", sep = ",", col.names = F, row.names = F)
bZ_out <- runLandscaper("d_bZ_grid.csv", "data_bZ.csv", 0.05, 2, 8)

bZ_out <- bZ_out %>% mutate(bZ_value_id = rep(seq_along(d_fix_bZ$value), each = GRID_RES))


d_mean_aZ_vals <- aZ_out %>%
  group_by(aZ_value_id) %>%
  summarise(meanPheno = mean(pheno),
            meanFitness = mean(fitness))

d_mean_bZ_vals <- bZ_out %>%
  group_by(bZ_value_id) %>%
  summarise(meanPheno = mean(pheno),
            meanFitness = mean(fitness))

d_fix_aZ$avFX <- d_mean_aZ_vals$meanPheno - aZ_default_avg$meanPheno
d_fix_aZ$avFit <- d_mean_aZ_vals$meanFitness - aZ_default_avg$meanFitness

d_fix_bZ$avFX <- d_mean_bZ_vals$meanPheno - bZ_default_avg$meanPheno
d_fix_bZ$avFit <- d_mean_bZ_vals$meanFitness - bZ_default_avg$meanFitness

d_fix_nar <- rbind(d_fix_aZ, d_fix_bZ)

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$")
)


ggplot(d_fix_nar, 
       aes(x = avFX, colour = mutType)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Average fixed phenotypic effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_fix_nar, 
       aes(x = avFit, colour = mutType)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Average fixed fitness effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_fix_nar, 
       aes(x = abs(value), colour = mutType)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Effect size on molecular component", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))
