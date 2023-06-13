library(tidyverse)
library(latex2exp)

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
d_fix_add <- d_fix_adapted %>% filter(modelindex == 1)

d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_add$gen, d_fix_add$seed, d_fix_add$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_add <- inner_join(d_fix_add, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

calcAddFitness <- function(phenotypes, optimum, width) {
  dists <- (phenotypes - optimum)^2
  return(exp(-(dists * width)))
}

CalcAddEffects <- function(dat, isFixed = T, dat_fixed = dat) {
  # If we are calculating fitness for segregating sites, need to evaluate fitness
  # vs the fixed effect background at the timepoint (i.e. all fixations at timepoint gen)
  # Fixation effect is multiplied by 2 because diploid
  dat_fixed <- dat_fixed %>% filter(modelindex == 1)
  dat <- dat %>% 
    group_by(gen, seed) %>%
    mutate(fixEffectSum = 2 * sum(dat_fixed[dat_fixed$gen <= cur_group()$gen &
                                              dat_fixed$seed == cur_group()$seed,]$value))
  
  if (isFixed) {
    # For fixed comparisons:
    # AA = 1+s; Aa = 1+hs; aa = 1
    # AA = fixEffectSum
    # Aa = fixEffectSum - value
    # aa = fixEffectSum - 2 * value
    Aa <- calcAddFitness(dat$fixEffectSum - dat$value, 2, 0.05)
    AA <- calcAddFitness(dat$fixEffectSum, 2, 0.05)
    aa <- calcAddFitness(dat$fixEffectSum - 2 * dat$value, 2, 0.05)
    dat$avFit <- Aa - aa
    dat$avFit_AA <- AA - aa
    dat$value_AA <- dat$value * 2
    dat$wAA <- AA
    dat$wAa <- Aa
    dat$waa <- aa
    return(dat)
  }
  
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = fixEffectSum + 2 * value
  # Aa = fixEffectSum + value
  # aa = fixEffectSum
  # Get effect
  Aa <- calcAddFitness(dat$fixEffectSum + dat$value, 2, 0.05)
  AA <- calcAddFitness(dat$fixEffectSum + dat$value * 2, 2, 0.05)
  aa <- calcAddFitness(dat$fixEffectSum, 2, 0.05)
  dat$avFit <- Aa - aa
  dat$avFit_AA <- AA - aa
  dat$value_AA <- dat$value * 2
  dat$wAA <- AA
  dat$wAa <- Aa
  dat$waa <- aa
  
  return(dat)
}

# Calculate dominance coefficient
CalcDominance <- function(dat) {
  dat <- dat %>%
    mutate(s = (wAA - waa),
           h = (waa - wAa) / s)
  return(dat)
}


d_fix_add <- CalcAddEffects(d_fix_add) %>% filter(gen >= 50000)
d_fix_add <- CalcDominance(d_fix_add)

runLandscaper <- function(df_path, output, width, optimum, threads, useID = FALSE) {
  command <- "ODELandscaper -i %s -o ./%s -w %f -p %f -t %i"
  if (useID) {
    command <- paste(command, "-I")
  }
  system(sprintf(command,
                 df_path, output, width, optimum, threads))
  result <- read_csv(paste0("./", output), col_names = F, col_types = "d")
  if (useID) {
    names(result) <- c("id", "fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  } else {
    names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  }
  
  return(result)
}

# Need to measure fitness effect relative to actual background then rather than 
# from a baseline, since the fitness effect relative to baseline isn't experienced 
# by the population
# so need to get the moltrait values at a fixed effect's given timepoint, and 
# subtract the fixed effect from it to measure the effect in context of the 
# population's background - can use mean pop phenotype and fitness to measure the
# background

CalcNARPhenotypeEffects <- function(dat, isFixed = T, dat_fixed = dat) {
  dat <- dat %>% filter(modelindex == 2)
  dat_fixed <- dat_fixed %>% filter(modelindex == 2)
  
  # calculate cumulative molecular component values at each step due to only 
  # fixed effects
  # multiply by 2 because diploid
  dat <- dat %>%
    group_by(gen, seed) %>%
    mutate(fixEffectSum_aZ = 2 * sum(dat_fixed[dat_fixed$gen <= cur_group()$gen &
                                                 dat_fixed$mutType == 3 &
                                                 dat_fixed$seed == cur_group()$seed,]$value),
           fixEffectSum_bZ = 2 * sum(dat_fixed[dat_fixed$gen <= cur_group()$gen &
                                                 dat_fixed$mutType == 4 &
                                                 dat_fixed$seed == cur_group()$seed,]$value))
  # Transform to exp scale
  dat$fixEffectSum_aZ <- exp(dat$fixEffectSum_aZ)
  dat$fixEffectSum_bZ <- exp(dat$fixEffectSum_bZ)
  
  dat$rowID <- as.integer(rownames(dat))
  
  # Get phenotypes with the mutation
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, fixEffectSum_aZ, fixEffectSum_bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_popfx <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  if (isFixed) {
    # For fixed comparisons:
    # AA = 1+s; Aa = 1+hs; aa = 1
    # AA = d_popfx$fixEffectSum
    # Aa = d_popfx$fixEffectSum - value
    # aa = d_popfx$fixEffectSum - 2 * value
    
    # Now take away the fixed effect: calculate seperately for alpha/beta
    d_dat_withoutFX_aZ <- dat %>% filter(mutType == 3)
    d_dat_withoutFX_bZ <- dat %>% filter(mutType == 4)
    d_dat_withoutFX_aZ$aZ <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_aZ) - d_dat_withoutFX_aZ$value)
    d_dat_withoutFX_aZ$bZ <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_bZ))
    d_dat_withoutFX_bZ$aZ <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_aZ))
    d_dat_withoutFX_bZ$bZ <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_bZ) - d_dat_withoutFX_bZ$value)
    
    # Homozygous estimation - for dominance calculation
    d_dat_withoutFX_aZ$aZ_aa <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_aZ) - 2 * d_dat_withoutFX_aZ$value)
    d_dat_withoutFX_aZ$bZ_aa <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_bZ))
    d_dat_withoutFX_bZ$aZ_aa <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_aZ))
    d_dat_withoutFX_bZ$bZ_aa <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_bZ) - 2 * d_dat_withoutFX_bZ$value)
    
    d_dat_withoutFX <- rbind(d_dat_withoutFX_aZ, d_dat_withoutFX_bZ)

    write.table(d_dat_withoutFX %>% ungroup() %>% 
                  dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
                "d_grid.csv", sep = ",", col.names = F, row.names = F)
    Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)

    write.table(d_dat_withoutFX %>% ungroup() %>% 
                  dplyr::select(rowID, aZ_aa, bZ_aa, KZ, KXZ), 
                "d_grid.csv", sep = ",", col.names = F, row.names = F)
    aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)

    # Get the effect size by taking away the phenotype missing that fixation
    # Ensure that the tables are aligned by id before we join them
    dat <- dat %>% arrange(rowID)
    Aa <- Aa %>% arrange(id)
    aa <- aa %>% arrange(id)
    
    dat$avFX <- d_popfx$pheno - Aa$pheno
    dat$avFit <- d_popfx$fitness - Aa$fitness
    dat$avFX_AA <- d_popfx$pheno - aa$pheno
    dat$avFit_AA <- d_popfx$fitness - aa$fitness
    dat$wAA <- d_popfx$fitness
    dat$wAa <- Aa$fitness
    dat$waa <- aa$fitness
    return(dat)
  }
  
  # Segregating mutation calculations
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = d_popfx$fixEffectSum + 2 * value
  # Aa = d_popfx$fixEffectSum + value
  # aa = d_popfx$fixEffectSum
  
  # Add on the segregating effect to the fixation effects
  d_dat_withFX_aZ <- dat %>% filter(mutType == 3)
  d_dat_withFX_bZ <- dat %>% filter(mutType == 4)
  
  d_dat_withFX_aZ$aZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + d_dat_withFX_aZ$value)
  d_dat_withFX_aZ$bZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + d_dat_withFX_bZ$value)
  
  # homozygous effect - can measure dominance coefficient
  d_dat_withFX_aZ$aZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + 2 * d_dat_withFX_aZ$value)
  d_dat_withFX_aZ$bZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + 2 * d_dat_withFX_bZ$value)
  
  d_dat_withFX <- rbind(d_dat_withFX_aZ, d_dat_withFX_bZ)
  
  # Get phenotypes with the mutation
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ_AA, bZ_AA, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  AA <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Ensure that the tables are aligned by id before we join them
  dat <- dat %>% arrange(rowID)
  Aa <- Aa %>% arrange(id)
  AA <- AA %>% arrange(id)
  
  # Get effect
  dat$avFX <- Aa$pheno - d_popfx$pheno
  dat$avFit <- Aa$fitness - d_popfx$fitness
  dat$avFX_AA <- AA$pheno - d_popfx$pheno
  dat$avFit_AA <- AA$fitness - d_popfx$fitness
  dat$wAA <- AA$fitness
  dat$wAa <- Aa$fitness
  dat$waa <- d_popfx$fitness
  return(dat)
}

d_fix_nar <- d_fix_adapted %>% ungroup() %>% mutate(r = rownames(.)) %>% 
  filter(modelindex == 2)

# First get matched mol trait data for generations where we have fixations
d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_nar$gen, d_fix_nar$seed, d_fix_nar$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_nar <- inner_join(d_fix_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))
d_fix_nar <- CalcNARPhenotypeEffects(d_fix_nar) %>% filter(gen >= 50000)
d_fix_nar <- CalcDominance(d_fix_nar)

# ranked mutations and average effects
RankFixations <- function(dat, isNAR) {
  index <- as.integer(isNAR) + 1
  
  if (isNAR) {
    d_fix_ranked <- dat %>%
      group_by(seed, modelindex) %>%
      arrange(gen, .by_group = T) %>%
      mutate(rank = row_number()) %>%
      dplyr::select(c(gen, rank, seed, modelindex, mutType, 
                      value, aZ, bZ, phenomean, w, avFX, avFit, avFit_AA, avFX_AA, 
                      wAA, wAa, waa, s, h))
    
    step0_pheno <- d_adapted %>% 
      filter(modelindex == index, gen == 49500, interaction(seed, modelindex) %in%
               interaction(d_fix_ranked$seed, d_fix_ranked$modelindex)) %>%
      mutate(rank = 0, value = NA, avFit = NA) %>%
      dplyr::select(gen, rank, seed, modelindex, value, 
                    phenomean, w, avFit)
  } else {
    d_fix_ranked <- dat %>%
      group_by(seed, modelindex) %>%
      arrange(gen, .by_group = T) %>%
      mutate(rank = row_number()) %>%
      dplyr::select(c(gen, rank, seed, modelindex, mutType, 
                      value, value_AA, aZ, bZ, phenomean, w, avFit, avFit_AA, 
                      wAA, wAa, waa, s, h))
    
    step0_pheno <- d_adapted %>% 
      filter(modelindex == index, gen == 49500, interaction(seed, modelindex) %in%
               interaction(d_fix_ranked$seed, d_fix_ranked$modelindex)) %>%
      mutate(rank = 0, value = NA, avFit = NA) %>%
      dplyr::select(gen, rank, seed, modelindex, value, 
                    phenomean, w, avFit)
  }
  
  return(rbind(d_fix_ranked, step0_pheno))
}

d_fix_ranked <- RankFixations(d_fix_nar, T)


# additive: attach step 0 (phenomean from before the first step in the walk)
d_fix_ranked_add <- RankFixations(d_fix_add, F)

# Summarise
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

# get segregating variants
d_muts_adapted %>% 
  filter(modelindex == 2, Freq < 1, interaction(gen, seed, modelindex) %in%
           interaction(d_fix_ranked$gen, d_fix_ranked$seed, 
                       d_fix_ranked$modelindex),
         Freq > 0) -> d_seg_ranked

d_seg_ranked %>%
  group_by(gen, seed) %>%
  distinct() %>%
  mutate(rank = tail(d_fix_ranked[d_fix_ranked$gen == cur_group()$gen &
                               d_fix_ranked$seed == cur_group()$seed,]$rank, 1)) -> d_seg_ranked

d_muts_adapted %>% 
  filter(modelindex == 1, Freq < 1, interaction(gen, seed, modelindex) %in%
           interaction(d_fix_ranked_add$gen, d_fix_ranked_add$seed, 
                       d_fix_ranked_add$modelindex),
         Freq > 0) -> d_seg_ranked_add

d_seg_ranked_add %>%
  distinct() %>%
  group_by(gen, seed) %>%
  mutate(rank = tail(d_fix_ranked_add[d_fix_ranked_add$gen == cur_group()$gen & 
                               d_fix_ranked_add$seed == cur_group()$seed,]$rank, 1)) -> d_seg_ranked_add

# Calculate fitness etc.
d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_seg_ranked$gen, d_seg_ranked$seed, d_seg_ranked$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_ranked <- inner_join(d_seg_ranked, d_qg_matched_seg, 
                           by = c("gen", "seed", "modelindex"))
d_seg_ranked <- CalcNARPhenotypeEffects(d_seg_ranked, F, d_fix_adapted)
d_seg_ranked <- CalcDominance(d_seg_ranked)


d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_seg_ranked_add$gen, d_seg_ranked_add$seed, 
                       d_seg_ranked_add$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_ranked_add <- inner_join(d_seg_ranked_add, d_qg_matched_seg, 
                               by = c("gen", "seed", "modelindex"))
d_seg_ranked_add <- CalcAddEffects(d_seg_ranked_add, F, d_fix_adapted)
d_seg_ranked_add <- CalcDominance(d_seg_ranked_add)

# contribution to trait seg vs fixed
GetSegFixContributions <- function(seg, fix, isNAR) {
  # weight effects by frequency
  seg %>%
    filter(rank > 0) %>%
    group_by(gen, seed, modelindex, rank) %>%
    summarise(segEffectSum = sum(abs(seg[seg$gen == cur_group()$gen & 
                                   seg$seed == cur_group()$seed,][[ifelse({{ isNAR }}, "avFX_AA", "value_AA")]])*
                                   seg[seg$gen == cur_group()$gen & 
                                         seg$seed == cur_group()$seed,]$Freq),
              weightSumSeg = sum(Freq),
              segWeightedAverage = segEffectSum/weightSumSeg) -> d_segFX

  # Select all fixations before or equal to this one (cumulative sum)
  fix %>%
    filter(rank > 0) %>%
    group_by(gen, seed, modelindex, rank) %>%
    summarise(fixEffectSum = sum(abs(fix[fix$rank <= cur_group()$rank & fix$rank > 0 &
                                           fix$seed == cur_group()$seed,][[ifelse({{ isNAR }}, "avFX_AA", "value_AA")]])),
              fixWeightedAverage = fixEffectSum/rank) -> d_fixFX
  
  return(inner_join(d_fixFX, d_segFX, 
                    by = c("gen", "seed", "modelindex", "rank")) %>%
    group_by(gen, seed, rank) %>%
      # weighted average: rank is the sum of fixation frequencies (1), weightSumFix
    mutate(percFix = (fixEffectSum)/(fixEffectSum + segEffectSum)) %>%
    ungroup(gen, seed) %>%
    summarise(meanPercFix = mean(percFix, na.rm = T),
              CIPercFix = CI(percFix, na.rm = T)))
}

d_segFixRat_add_sum <- GetSegFixContributions(d_seg_ranked_add, d_fix_ranked_add, F)
d_segFixRat_sum <- GetSegFixContributions(d_seg_ranked, d_fix_ranked, T)

# seg effects weighted by frequency
d_seg_ranked$weighteds <- d_seg_ranked$s * d_seg_ranked$Freq
d_seg_ranked_add$weighteds <- d_seg_ranked_add$s * d_seg_ranked_add$Freq

# Adapted and Maladapted
# d_maladapted <- d_qg %>% filter(!isAdapted)
# 
# d_maladapted %>% 
#   group_by(gen, modelindex) %>%
#   summarise(meanPheno = mean(phenomean),
#             CIPheno = CI(phenomean)) -> d_maladapted_sum
# 
# d_fix_maladapted <- d_fix %>% filter(interaction(seed, modelindex) %in% 
#                                     interaction(d_maladapted$seed, d_maladapted$modelindex))
# 
# # If we want to combine adapted and maladapted
# d_fix_maladapted <- d_fix 
# 
# d_fix_maladapted$fixTime <- d_fix_maladapted$gen - d_fix_maladapted$originGen
# 
# d_fix_mal_add <- d_fix_maladapted %>% filter(modelindex == 1, gen >= 50000)
# 
# d_qg_matched_fix <- d_qg %>% 
#   filter(interaction(gen, seed, modelindex) %in% 
#            interaction(d_fix_mal_add$gen, d_fix_mal_add$seed, d_fix_mal_add$modelindex)) %>%
#   dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w, isAdapted) %>% distinct()
# 
# d_fix_mal_add <- inner_join(d_fix_mal_add, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))
# d_fix_mal_add <- CalcAddEffects(d_fix_mal_add)
# 
# d_fix_mal_nar <- d_fix_maladapted %>% filter(modelindex == 2, gen >= 50000)
# 
# # First get matched mol trait data for generations where we have fixations
# d_qg_matched_fix <- d_qg %>% 
#   filter(interaction(gen, seed, modelindex) %in% 
#            interaction(d_fix_mal_nar$gen, d_fix_mal_nar$seed, d_fix_mal_nar$modelindex)) %>%
#   dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w, isAdapted) %>% distinct()
# 
# d_fix_mal_nar <- inner_join(d_fix_mal_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))
# d_fix_mal_nar <- CalcNARPhenotypeEffects(d_fix_mal_nar)
# 
# # ranked mutations and average effects
# d_fix_mal_ranked <- d_fix_mal_nar %>%
#   group_by(seed, modelindex) %>%
#   arrange(gen, .by_group = T) %>%
#   mutate(rank = row_number()) %>%
#   dplyr::select(c(rank, seed, modelindex, mutType, value, aZ, bZ, 
#            phenomean, w, avFX, avFit, isAdapted))
# 
# step0_pheno <- d_qg %>% 
#   filter(modelindex == 2, gen == 49500, interaction(seed, modelindex) %in%
#            interaction(d_fix_mal_ranked$seed, d_fix_mal_ranked$modelindex))
# 
# step0_pheno$rank <- 0
# step0_pheno$value <- NA
# step0_pheno$avFit <- NA
# 
# d_fix_mal_ranked <- rbind(d_fix_mal_ranked, step0_pheno %>% 
#                         dplyr::select(rank, seed, modelindex, value, phenomean, w, avFit, isAdapted))
# 
# 
# # additive: attach step 0 (phenomean from before the first step in the walk)
# d_fix_ranked_add_mal <- d_fix_mal_add %>%
#   group_by(seed, modelindex) %>%
#   arrange(gen, .by_group = T) %>%
#   mutate(rank = row_number()) %>%
#   dplyr::select(c(rank, seed, modelindex, value, phenomean, w, avFit, isAdapted)) %>%
#   ungroup()
# 
# step0_pheno <- d_qg %>% 
#   filter(modelindex == 1, gen == 49500, interaction(seed, modelindex) %in%
#            interaction(d_fix_ranked_add_mal$seed, d_fix_ranked_add_mal$modelindex))
# 
# step0_pheno$rank <- 0
# step0_pheno$value <- NA
# step0_pheno$avFit <- NA
# 
# d_fix_ranked_add_mal <- rbind(d_fix_ranked_add_mal, step0_pheno %>% 
#                             dplyr::select(rank, seed, modelindex, value, phenomean, w, 
#                                    avFit, isAdapted))

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$")
)

# save data frames
d_fix_combined <- rbind(d_fix_add, d_fix_nar)
d_fix_combined <- d_fix_combined %>% dplyr::select(-c(fixGen, constraint, chi, Count))
write_csv(rbind(d_fix_add, d_fix_nar), "d_fix_combined.csv")

# d_fix_mal_combined <- rbind(d_fix_mal_add, d_fix_mal_nar)
# write_csv(d_fix_mal_combined, "d_fix_combined_all.csv")

# d_ranked_combined <- rbind(d_fix_mal_ranked, d_fix_ranked_add_mal)
# write_csv(d_ranked_combined, "d_ranked_combined_all.csv")

d_rank_combined_tbl <- d_ranked_combined %>% filter(rank != 0, modelindex == 2) %>%
  group_by(rank, modelindex, isAdapted, mutType) %>% 
  summarise(n = n()) %>% ungroup() %>% 
  complete(rank, isAdapted, mutType)
d_rank_combined_tbl <- d_rank_combined_tbl[-c(15, 16),]
d_rank_combined_tbl

d_ranked_combined %>% filter(modelindex == 2) %>%
  group_by(seed, isAdapted) %>%
  mutate(aZbZ = aZ/bZ) -> d_ratio

