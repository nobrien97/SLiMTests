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
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg

d_qg %>% group_by(modelindex) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
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

d_muts_adapted <- d_muts %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_adapted$seed, d_adapted$modelindex))

d_com_adapted <- inner_join(d_adapted, d_muts_adapted, by = c("gen", "seed", "modelindex"))

d_fix <- d_muts %>%
  filter(Freq == 1) %>%
  group_by(seed, modelindex, mutType) %>%
  distinct(mutID, .keep_all = T) 

d_fix_adapted <- d_fix %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_adapted$seed, d_adapted$modelindex))

d_fix_adapted$fixTime <- d_fix_adapted$gen - d_fix_adapted$originGen

d_indPheno <- read.table(paste0(data_path, "slim_indPheno.csv"), header = F, 
                     sep = ",", colClasses = c("integer", "factor", "factor",
                                               "integer", rep("numeric", times = 5)), 
                     col.names = c("gen", "seed", "modelindex", "index", "pheno",
                                   "aZ", "bZ", "KZ", "KXZ"), 
                     fill = T)

d_indPheno_adapted <- d_indPheno %>% filter(interaction(gen, seed, modelindex) %in% 
                                    interaction(d_fix_adapted$gen, d_fix_adapted$seed, d_fix_adapted$modelindex))

rm(d_indPheno)

# Get fitness effect by subtracting fitness
d_fix_add <- d_fix_adapted %>% filter(modelindex == 1, gen >= 50000)

d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_add$gen, d_fix_add$seed, d_fix_add$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_add <- inner_join(d_fix_add, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

calcAddFitness <- function(phenotypes, optimum, width) {
  dists <- (phenotypes - optimum)^2
  return(exp(-(dists * width)))
}

BASELINE_PHENO_ADD <- 0
BASELINE_PHENO <- 2.53531
BASELINE_FITNESS_ADD <- calcAddFitness(BASELINE_PHENO_ADD, 2, 0.05)
BASELINE_FITNESS <- calcAddFitness(BASELINE_PHENO, 2, 0.05)

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
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_nar2 <- inner_join(d_fix_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

d_fix_aZ <- d_fix_nar2 %>% filter(mutType == 3)
d_fix_bZ <- d_fix_nar2 %>% filter(mutType == 4)

# Calculate the mean phenotypes with the sampled mean aZ/bZ values
write.table(d_fix_aZ %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_aZ <- runLandscaper("d_grid_aZ.csv", "data_popfx_aZ.csv", 0.05, 2, 8)

write.table(d_fix_bZ %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_bZ <- runLandscaper("d_grid_bZ.csv", "data_popfx_bZ.csv", 0.05, 2, 8)

# measure difference relative to baseline - 1, 1, 1, 1 case
# so (baseline + mol comp value) - baseline 
# maybe also look over many backgrounds?

# Calculate the phenotypes when we take away the fixed effect in question from aZ/bZ
d_fix_aZ_diff <- d_fix_aZ
d_fix_aZ_diff$aZ <- exp(log(d_fix_aZ$aZ) - d_fix_aZ$value)

d_fix_bZ_diff <- d_fix_bZ
d_fix_bZ_diff$bZ <- exp(log(d_fix_bZ$bZ) - d_fix_bZ$value)

write.table(d_fix_aZ_diff %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_aZ_diff <- runLandscaper("d_grid_aZ_diff.csv", "data_popfx_aZ_diff.csv", 0.05, 2, 8)

write.table(d_fix_bZ_diff %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_bZ_diff <- runLandscaper("d_grid_bZ_diff.csv", "data_popfx_bZ_diff.csv", 0.05, 2, 8)

# Get the effect size by taking away the phenotype missing that fixation
d_fix_aZ$avFX <- d_fix_aZ$phenomean - d_popfx_aZ_diff$pheno
d_fix_aZ$avFit <- d_fix_aZ$w - d_popfx_aZ_diff$fitness

d_fix_bZ$avFX <- d_fix_bZ$phenomean - d_popfx_bZ_diff$pheno
d_fix_bZ$avFit <- d_fix_bZ$w - d_popfx_bZ_diff$fitness

d_fix_nar <- rbind(d_fix_aZ, d_fix_bZ)

# ranked mutations and average effects
d_fix_ranked <- d_fix_nar %>%
  group_by(seed, modelindex) %>%
  arrange(gen, .by_group = T) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(c(rank, seed, modelindex, mutType, value, aZ, bZ, phenomean, w, avFX, avFit))

step0_pheno <- d_adapted %>% 
  filter(modelindex == 2, gen == 49500, interaction(seed, modelindex) %in%
           interaction(d_fix_ranked$seed, d_fix_ranked$modelindex))

step0_pheno$rank <- 0
step0_pheno$value <- NA
step0_pheno$avFit <- NA

d_fix_ranked <- rbind(d_fix_ranked, step0_pheno %>% 
                            dplyr::select(rank, seed, modelindex, value, phenomean, w, avFit))


# additive: attach step 0 (phenomean from before the first step in the walk)
d_fix_ranked_add <- d_fix_add %>%
  group_by(seed, modelindex) %>%
  arrange(gen, .by_group = T) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(c(rank, seed, modelindex, value, phenomean, w, avFit)) %>%
  ungroup()

step0_pheno <- d_adapted %>% 
  filter(modelindex == 1, gen == 49500, interaction(seed, modelindex) %in%
           interaction(d_fix_ranked_add$seed, d_fix_ranked_add$modelindex))

step0_pheno$rank <- 0
step0_pheno$value <- NA
step0_pheno$avFit <- NA

d_fix_ranked_add <- rbind(d_fix_ranked_add, step0_pheno %>% 
             dplyr::select(rank, seed, modelindex, value, phenomean, w, avFit))

d_fix_ranked_add %>% filter(rank != 0) %>%
  group_by(rank) %>% 
  summarise(n())

d_fix_ranked %>% filter(rank != 0) %>%
  group_by(rank) %>% 
  summarise(n())

rbind(d_fix_ranked %>% mutate(model = "NAR"), 
      d_fix_ranked_add %>% mutate(model = "Additive")) %>% 
  group_by(model, rank) %>%
  filter(rank != 0) %>%
  summarise(n = n())


# Adapted and Maladapted
d_maladapted <- d_qg %>% filter(!isAdapted)

d_maladapted %>% 
  group_by(gen, modelindex) %>%
  summarise(meanPheno = mean(phenomean),
            CIPheno = CI(phenomean)) -> d_maladapted_sum

d_fix_maladapted <- d_fix %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_maladapted$seed, d_maladapted$modelindex))

# If we want to combine adapted and maladapted
d_fix_maladapted <- d_fix 


d_fix_maladapted$fixTime <- d_fix_maladapted$gen - d_fix_maladapted$originGen

d_fix_mal_add <- d_fix_maladapted %>% filter(modelindex == 1, gen >= 50000)

d_qg_matched_fix <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_mal_add$gen, d_fix_mal_add$seed, d_fix_mal_add$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w, isAdapted) %>% distinct()

d_fix_mal_add <- inner_join(d_fix_mal_add, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

d_absenceFitness <- d_fix_mal_add %>% mutate(phenomean = phenomean - value)
d_absenceFitness$absenceW <- calcAddFitness(d_absenceFitness$phenomean, 2, 0.05)
d_fix_mal_add$avFit <- d_fix_mal_add$w - d_absenceFitness$absenceW

d_fix_mal_nar <- d_fix_maladapted %>% filter(modelindex == 2, gen >= 50000)

# First get matched mol trait data for generations where we have fixations
d_qg_matched_fix <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_mal_nar$gen, d_fix_mal_nar$seed, d_fix_mal_nar$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w, isAdapted) %>% distinct()

d_fix_mal_nar <- inner_join(d_fix_mal_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

d_fix_mal_aZ <- d_fix_mal_nar %>% filter(mutType == 3)
d_fix_mal_bZ <- d_fix_mal_nar %>% filter(mutType == 4)

# Calculate the mean phenotypes with the sampled mean aZ/bZ values
write.table(d_fix_mal_aZ %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_mal_aZ <- runLandscaper("d_grid_aZ.csv", "data_popfx_aZ.csv", 0.05, 2, 8)

write.table(d_fix_mal_bZ %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_mal_bZ <- runLandscaper("d_grid_bZ.csv", "data_popfx_bZ.csv", 0.05, 2, 8)

# Calculate the phenotypes when we take away the fixed effect in question from aZ/bZ
d_fix_mal_aZ_diff <- d_fix_mal_aZ
d_fix_mal_aZ_diff$aZ <- exp(log(d_fix_mal_aZ$aZ) - d_fix_mal_aZ$value)

d_fix_mal_bZ_diff <- d_fix_mal_bZ
d_fix_mal_bZ_diff$bZ <- exp(log(d_fix_mal_bZ$bZ) - d_fix_mal_bZ$value)

write.table(d_fix_mal_aZ_diff %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_mal_aZ_diff <- runLandscaper("d_grid_aZ_diff.csv", "data_popfx_aZ_diff.csv", 0.05, 2, 8)

write.table(d_fix_mal_bZ_diff %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_mal_bZ_diff <- runLandscaper("d_grid_bZ_diff.csv", "data_popfx_bZ_diff.csv", 0.05, 2, 8)

# Get the effect size by taking away the phenotype missing that fixation
d_fix_mal_aZ$avFX <- d_fix_mal_aZ$phenomean - d_popfx_mal_aZ_diff$pheno
d_fix_mal_aZ$avFit <- d_fix_mal_aZ$w - d_popfx_mal_aZ_diff$fitness

d_fix_mal_bZ$avFX <- d_fix_mal_bZ$phenomean - d_popfx_mal_bZ_diff$pheno
d_fix_mal_bZ$avFit <- d_fix_mal_bZ$w - d_popfx_mal_bZ_diff$fitness

d_fix_mal_nar <- rbind(d_fix_mal_aZ, d_fix_mal_bZ)

# ranked mutations and average effects
d_fix_mal_ranked <- d_fix_mal_nar %>%
  group_by(seed, modelindex) %>%
  arrange(gen, .by_group = T) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(c(rank, seed, modelindex, mutType, value, aZ, bZ, 
           phenomean, w, avFX, avFit, isAdapted))

step0_pheno <- d_qg %>% 
  filter(modelindex == 2, gen == 49500, interaction(seed, modelindex) %in%
           interaction(d_fix_mal_ranked$seed, d_fix_mal_ranked$modelindex))

step0_pheno$rank <- 0
step0_pheno$value <- NA
step0_pheno$avFit <- NA

d_fix_mal_ranked <- rbind(d_fix_mal_ranked, step0_pheno %>% 
                        dplyr::select(rank, seed, modelindex, value, phenomean, w, avFit, isAdapted))


# additive: attach step 0 (phenomean from before the first step in the walk)
d_fix_ranked_add_mal <- d_fix_mal_add %>%
  group_by(seed, modelindex) %>%
  arrange(gen, .by_group = T) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(c(rank, seed, modelindex, value, phenomean, w, avFit, isAdapted)) %>%
  ungroup()

step0_pheno <- d_qg %>% 
  filter(modelindex == 1, gen == 49500, interaction(seed, modelindex) %in%
           interaction(d_fix_ranked_add_mal$seed, d_fix_ranked_add_mal$modelindex))

step0_pheno$rank <- 0
step0_pheno$value <- NA
step0_pheno$avFit <- NA

d_fix_ranked_add_mal <- rbind(d_fix_ranked_add_mal, step0_pheno %>% 
                            dplyr::select(rank, seed, modelindex, value, phenomean, w, 
                                   avFit, isAdapted))

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$")
)

seed <- sample(0:.Machine$integer.max, 1)
set.seed(seed)
#set.seed(1875704954)
sampled_seed <- sample(d_com_adapted[d_com_adapted$modelindex == 2,]$seed, 3)
d_com_adapted %>% 
  filter(modelindex == 2, seed %in% sampled_seed) %>% distinct() -> d_com_nar_sample

sampled_seed_add <- sample(d_com_adapted[d_com_adapted$modelindex == 1,]$seed, 3)
d_com_adapted %>% 
  filter(modelindex == 1, seed %in% sampled_seed_add) %>% distinct() -> d_com_add_sample

# Calculate fitness effects
## Additive
d_absenceFitness <- d_com_add_sample %>% mutate(phenomean = phenomean - value)
d_absenceFitness$absenceW <- calcAddFitness(d_absenceFitness$phenomean, 2, 0.05)
d_com_add_sample$avFit <- d_com_add_sample$w - d_absenceFitness$absenceW


## NAR
d_fix_aZ <- d_com_nar_sample %>% filter(mutType == 3)
d_fix_bZ <- d_com_nar_sample %>% filter(mutType == 4)

# Calculate the mean phenotypes with the sampled mean aZ/bZ values
write.table(d_fix_aZ %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_aZ <- runLandscaper("d_grid_aZ.csv", "data_popfx_aZ.csv", 0.05, 2, 8)

write.table(d_fix_bZ %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ.csv", sep = ",", col.names = F, row.names = F)
d_popfx_bZ <- runLandscaper("d_grid_bZ.csv", "data_popfx_bZ.csv", 0.05, 2, 8)

# Calculate the phenotypes when we take away the fixed effect in question from aZ/bZ
d_fix_aZ_diff <- d_fix_aZ
d_fix_aZ_diff$aZ <- exp(log(d_fix_aZ$aZ) - d_fix_aZ$value)

d_fix_bZ_diff <- d_fix_bZ
d_fix_bZ_diff$bZ <- exp(log(d_fix_bZ$bZ) - d_fix_bZ$value)

write.table(d_fix_aZ_diff %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_aZ_diff <- runLandscaper("d_grid_aZ_diff.csv", "data_popfx_aZ_diff.csv", 0.05, 2, 8)

write.table(d_fix_bZ_diff %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_bZ_diff.csv", sep = ",", col.names = F, row.names = F)
d_popfx_bZ_diff <- runLandscaper("d_grid_bZ_diff.csv", "data_popfx_bZ_diff.csv", 0.05, 2, 8)

d_fix_aZ$avFX <- d_fix_aZ$phenomean - d_popfx_aZ_diff$pheno
d_fix_aZ$avFit <- d_fix_aZ$w - d_popfx_aZ_diff$fitness

d_fix_bZ$avFX <- d_fix_bZ$phenomean - d_popfx_bZ_diff$pheno
d_fix_bZ$avFit <- d_fix_bZ$w - d_popfx_bZ_diff$fitness

d_com_nar_sample <- rbind(d_fix_aZ, d_fix_bZ)

# save data frames
d_fix_combined <- rbind(d_fix_add, d_fix_nar)
d_fix_combined <- d_fix_combined %>% dplyr::select(-c(fixGen, constraint, chi, Count))
write_csv(rbind(d_fix_add, d_fix_nar), "d_fix_combined.csv")

d_fix_mal_combined <- rbind(d_fix_mal_add, d_fix_mal_nar)
write_csv(d_fix_mal_combined, "d_fix_combined_all.csv")

d_ranked_combined <- rbind(d_fix_mal_ranked, d_fix_ranked_add_mal)
write_csv(d_ranked_combined, "d_ranked_combined_all.csv")

d_rank_combined_tbl <- d_ranked_combined %>% filter(rank != 0, modelindex == 2) %>%
  group_by(rank, modelindex, isAdapted, mutType) %>% 
  summarise(n = n()) %>% ungroup() %>% 
  complete(rank, isAdapted, mutType)
d_rank_combined_tbl <- d_rank_combined_tbl[-c(15, 16),]
d_rank_combined_tbl

d_ranked_combined %>% filter(modelindex == 2) %>%
  group_by(seed, isAdapted) %>%
  mutate(aZbZ = aZ/bZ) -> d_ratio
