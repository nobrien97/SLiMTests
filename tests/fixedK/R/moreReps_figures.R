library(tidyverse)
library(paletteer)


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
  group_by(gen, modelindex) %>%
  summarise(meanPheno = mean(phenomean),
            CIPheno = CI(phenomean)) -> d_qg_sum

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


# Fig 2 - phenotype mean
ggplot(d_qg_sum %>% filter(gen > 49000) %>% mutate(gen = gen - 50000), 
       aes(x = gen, y = meanPheno, colour = modelindex)) +
  geom_line() +
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
  theme(text = element_text(size = 16))

# Fig 3 - effect sizes
d_fix$fixTime <- d_fix$gen - d_fix$originGen

## Additive
ggplot(d_fix %>% filter(modelindex == 1), 
       aes(x = abs(value), colour = modelindex)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(x = "Fixation effect size", y = "Density",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16))

## Network: need to get the average phenotype effect over a range of genetic backgrounds
## to do this, for a given alpha effect that was sampled, 
## we look at its effect on phenotype and fitness over a range of beta values
## each effect is a deviation from the standard rate (aZ = bZ = 1)
## need to measure the average effect on phenotype compared to this standard rate: just
## the difference between the phenotype values

# first need to generate the standard effect for alpha and beta 
# from the range of mutational effects - remember it needs to be exponentiated

d_fix_nar <- d_fix %>% filter(modelindex == 2, gen >= 50000)

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


runLandscaper <- function(df_path, output, width, optimum, threads) {
  system(sprintf("ODELandscaper -i %s -o ./%s -w %f -p %f -t %i",
                 df_path, output, width, optimum, threads))
  result <- read_csv(paste0("./", output), col_names = F, col_types = "d")
  names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  return(result)
}

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

ggplot(d_fix_nar, 
       aes(x = avFX, colour = mutType)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("aZ", "bZ")) +
  labs(x = "Average fixed phenotypic effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_fix_nar, 
       aes(x = avFit, colour = mutType)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("aZ", "bZ")) +
  labs(x = "Average fixed fitness effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_fix_nar, 
       aes(x = abs(value), colour = mutType)) +
  geom_density() +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("aZ", "bZ")) +
  labs(x = "Effect size on molecular component", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))

# On average fitness effect is negative - so mutations require specific
# genetic backgrounds to be positive
# Need to measure fitness effect relative to background then rather than average
# so need to get the moltrait values at a fixed effect's given timepoint, and 
# subtract the fixed effect from it to measure the effect in context of the 
# population's background

# First get matched mol trait data for generations where we have fixations
d_qg_matched_fix <- d_qg %>% 
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

d_popfx_aZ <- d_popfx_aZ %>% 
  mutate(pheno_diff = d_popfx_aZ_diff$pheno - d_fix_aZ$phenomean,
         fitness_diff = d_popfx_aZ_diff$fitness - d_fix_aZ$w) 

mean(d_popfx_aZ$pheno_diff)
mean(d_popfx_aZ$fitness_diff)

d_popfx_bZ <- d_popfx_bZ %>% 
  mutate(pheno_diff = d_popfx_bZ_diff$pheno - d_fix_bZ$phenomean,
         fitness_diff = d_popfx_bZ_diff$fitness - d_fix_bZ$w) 

mean(d_popfx_bZ$pheno_diff)
mean(d_popfx_bZ$fitness_diff)
