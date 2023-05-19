library(tidyverse)

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

d_fix_adapted$fixTime <- d_fix_adapted$gen - d_fix_adapted$originGen
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
