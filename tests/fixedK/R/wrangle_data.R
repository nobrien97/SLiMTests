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

estimate_mode <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  if (length(x) < 2)
    return(NA)
  d <- density(x)
  d$x[which.max(d$y)]
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

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() %>%
  filter(gen == 49500, isAdapted) %>%
  summarise(n())

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.95, 2.05))) %>%
  ungroup() %>%
  filter(gen == 49500, isAdapted) %>%
  summarise(n())


# 1455 additive models adapted total
# 1181 network models adapted total

d_adapted <- d_qg %>% filter(isAdapted)

d_adapted %>% 
  mutate(model = if_else(modelindex == 1, "Additive", "NAR")) %>%
  group_by(gen, model) %>%
  summarise(meanPheno = mean(phenomean),
            CIPheno = CI(phenomean),
            medPheno = median(phenomean),
            modePheno = estimate_mode(phenomean),
            meanGenomeH = mean(meanH),
            CIGenomeH = CI(meanH)) -> d_adapted_sum

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

# every model had at least 1 fixation
# length(unique(interaction(d_fix$seed, d_fix$modelindex)))

# observed heterozygosity
d_het <- read.table("/mnt/d/SLiMTests/tests/fixedK/calcHet/slim_locusHo.csv", 
                     header = F, sep = ",", 
                     colClasses = c("integer", "factor", "factor", "numeric", "numeric"), 
                                    col.names = c("gen", "seed", "modelindex", "Ho_l1", "Ho_l2"), 
                                    fill = T)
  
d_Ho <- d_het %>% 
  pivot_longer(cols = c(Ho_l1, Ho_l2), names_to = "locus", values_to = "Ho") %>%
  mutate(model = ifelse(modelindex == 1, "Additive", "NAR"))

d_Ho %>%
  group_by(gen, model) %>%
  summarise(meanHo = mean(Ho),
            CIHo = CI(Ho)) -> d_Ho_sum

View(d_Ho %>%
  group_by(model) %>%
  summarise(meanHo = mean(Ho),
            CIHo = CI(Ho)))



  

# is the repeat simulation identical to the original?
d_qg_het <- read.table("/mnt/d/SLiMTests/tests/fixedK/calcHet/slim_qg.csv", 
                       sep = ",", colClasses = c("integer", "factor", "factor", 
                                                 rep("numeric", times = 12)), 
                       col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                     "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                     "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                       fill = T)

d_qg_het %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg_het

d_adapted_het <- d_qg_het %>% filter(isAdapted)


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
    dat$AA_pheno <- dat$fixEffectSum
    dat$Aa_pheno <- dat$fixEffectSum - dat$value
    dat$aa_pheno <- dat$fixEffectSum - 2 * dat$value
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
  
  dat$AA_pheno <- dat$fixEffectSum + 2 * dat$value
  dat$Aa_pheno <- dat$fixEffectSum + dat$value
  dat$aa_pheno <- dat$fixEffectSum
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
    
    dat$AA_pheno <- d_popfx$pheno
    dat$Aa_pheno <- Aa$pheno
    dat$aa_pheno <- aa$pheno
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
  dat$AA_pheno <- AA$pheno
  dat$Aa_pheno <- Aa$pheno
  dat$aa_pheno <- d_popfx$pheno
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
RankFixations <- function(dat, isNAR, dat_burnInFX = dat) {
  index <- as.integer(isNAR) + 1
  
  if (isNAR) {
    # Get fixed effects up to each rank
    dat <- dat %>%
      group_by(gen, seed) %>%
      mutate(fixEffectSum_aZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 3 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value)),
             fixEffectSum_bZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 4 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value)))
    
    d_fix_ranked <- dat %>%
      group_by(seed, modelindex) %>%
      arrange(gen, .by_group = T) %>%
      mutate(rank = row_number()) %>%
      dplyr::select(c(gen, rank, seed, modelindex, mutID, mutType, originGen,
                      value, aZ, bZ, phenomean, w, fixEffectSum_aZ, fixEffectSum_bZ,
                      avFX, avFit, avFit_AA, avFX_AA, AA_pheno, Aa_pheno, aa_pheno,
                      wAA, wAa, waa, s, h))
    
    step0_pheno <- d_adapted %>% 
      filter(modelindex == index, gen == 49500, interaction(seed, modelindex) %in%
               interaction(d_fix_ranked$seed, d_fix_ranked$modelindex)) %>%
      group_by(gen, seed) %>%
      mutate(rank = 0, value = NA, mutID = NA, avFit = NA,
             fixEffectSum_aZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 3 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value)),
             fixEffectSum_bZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 4 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value))) %>%
      dplyr::select(gen, rank, seed, modelindex, mutID, value, 
                    aZ, bZ, phenomean, w, fixEffectSum_aZ, fixEffectSum_bZ, avFit)
  } else {
    d_fix_ranked <- dat %>%
      group_by(seed, modelindex) %>%
      arrange(gen, .by_group = T) %>%
      mutate(rank = row_number()) %>%
      dplyr::select(c(gen, rank, seed, modelindex, mutID, mutType, originGen,
                      fixEffectSum, value, value_AA, aZ, bZ, phenomean, w, 
                      avFit, avFit_AA, AA_pheno, Aa_pheno, aa_pheno, 
                      wAA, wAa, waa, s, h))
    
    step0_pheno <- d_adapted %>% 
      filter(modelindex == index, gen == 49500, interaction(seed, modelindex) %in%
               interaction(d_fix_ranked$seed, d_fix_ranked$modelindex)) %>%
      group_by(gen, seed) %>%
      mutate(rank = 0, value = NA, mutID = NA, avFit = NA,
             fixEffectSum = 2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                     dat_burnInFX$seed == cur_group()$seed,]$value)) %>%
      dplyr::select(gen, rank, seed, modelindex, mutID, value,
                    aZ, bZ, phenomean, w, fixEffectSum, avFit)
  }
  
  return(rbind(d_fix_ranked, step0_pheno))
}

d_fix_ranked <- RankFixations(d_fix_nar, T, d_fix %>% filter(modelindex == 2))

# additive: attach step 0 (phenomean from before the first step in the walk)
d_fix_ranked_add <- RankFixations(d_fix_add, F, d_fix %>% filter(modelindex == 1))

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

# Since there aren't many populations with steps >3, we'll organise the groups
# into 1, 2, >=3
d_fix_ranked %>% 
  mutate(rankFactor = ifelse(rank > 2, "\\geq 3", as.character(rank))) -> d_fix_ranked

d_fix_ranked$rankFactor <- factor(d_fix_ranked$rankFactor, 
                                  levels = c("0", "1", "2", "\\geq 3"))

d_fix_ranked_add %>% 
  mutate(rankFactor = ifelse(rank > 2, "\\geq 3", as.character(rank))) -> d_fix_ranked_add

d_fix_ranked_add$rankFactor <- factor(d_fix_ranked_add$rankFactor, 
                                  levels = c("0", "1", "2", "\\geq 3"))


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




# are overall phenomeans skewed by walks with only 1 fixation?
d_fix_ranked_combined <- rbind(d_fix_ranked, d_fix_ranked_add)
d_fix_ranked_combined$model <- if_else(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")

d_fix_ranked_combined %>%
  group_by(seed, modelindex) %>%
  mutate(manyFixed = any(rank > 1)) %>%
  filter(!manyFixed) -> d_1fixation

d_fix_ranked_combined %>%
  group_by(seed, modelindex) %>%
  mutate(manyFixed = any(rank > 1)) %>%
  filter(manyFixed) -> d_manyfixation

d_adapted_1fix <- d_adapted %>% 
  filter(interaction(seed, modelindex) %in% 
           interaction(d_1fixation$seed, d_1fixation$modelindex))

d_adapted_manyfix <- d_adapted %>% 
  filter(interaction(seed, modelindex) %in% 
           interaction(d_manyfixation$seed, d_manyfixation$modelindex))

d_adapted_manyfixsum <- d_adapted_manyfix %>%
  group_by(gen, modelindex) %>%
  summarise(meanPheno = mean(phenomean),
            CIPheno = CI(phenomean))

# upon plotting these above phenomeans, doesn't make any difference: trait evo
# similar between all of them - must be the phenomeans in the fix_ranked calculation


# Some fixations are deleterious in early steps - why?
# Could be they are adaptive in burn-in environment, reached a high freq
# there and then drifted the rest of the way
# so need to calculate fitness effect of the mutation relative to that environment:
# where they originated, and when the mutation first reached >50% frequency
d_fix_ranked_combined %>% filter(s < 0) -> d_fix_del

# need to find simulations with those muts
d_muts_adapted %>% filter(interaction(seed, modelindex) %in%
                            interaction(d_fix_del$seed, 
                                        d_fix_del$modelindex)) -> d_muts_del

d_muts_del %>% filter(mutID %in% d_fix_del$mutID) %>% distinct() -> d_muts_del

# Calculate fitness effects when they arose, 
# and when they first reached 50% or greater freq - not fixed
d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_muts_del$gen, d_muts_del$seed, d_muts_del$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_del <- inner_join(d_muts_del, d_qg_matched_seg, 
                        by = c("gen", "seed", "modelindex"))

# Calculate phenotypes
d_seg_del_add <- CalcAddEffects(d_seg_del %>% filter(modelindex == 1), isFixed = F,
                                dat_fixed = d_fix_adapted %>% filter(modelindex == 1))
d_seg_del <- CalcNARPhenotypeEffects(d_seg_del %>% filter(modelindex == 2), isFixed = F, 
                        dat_fixed = d_fix_adapted %>% filter(modelindex == 2))

# Fitness calculations aren't correct for those before optimum shift: recalculate
d_seg_del[d_seg_del$gen < 50000,]$wAA <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$AA_pheno, 1, 0.05)
d_seg_del[d_seg_del$gen < 50000,]$wAa <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$Aa_pheno, 1, 0.05)
d_seg_del[d_seg_del$gen < 50000,]$waa <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$aa_pheno, 1, 0.05)
d_seg_del$avFit <- d_seg_del$wAa - d_seg_del$waa
d_seg_del$avFit_AA <- d_seg_del$wAA - d_seg_del$waa
d_seg_del$s <- d_seg_del$wAA - d_seg_del$waa

d_seg_del_add[d_seg_del_add$gen < 50000,]$wAA <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$AA_pheno, 1, 0.05)
d_seg_del_add[d_seg_del_add$gen < 50000,]$wAa <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$Aa_pheno, 1, 0.05)
d_seg_del_add[d_seg_del_add$gen < 50000,]$waa <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$aa_pheno, 1, 0.05)
d_seg_del_add$avFit <- d_seg_del_add$wAa - d_seg_del_add$waa
d_seg_del_add$avFit_AA <- d_seg_del_add$wAA - d_seg_del_add$waa
d_seg_del_add$s <- d_seg_del_add$wAA - d_seg_del_add$waa

d_seg_del <- rbind(d_seg_del, d_seg_del_add)
d_seg_del$model <- ifelse(d_seg_del$modelindex == 2, "NAR", "Additive")
rm(d_seg_del_add)

# Sample a few mutations
eg_muts <- sample(unique(d_seg_del$mutID), 10)

# Plot frequency 
ggplot(d_seg_del %>% filter(mutID %in% eg_muts),
       aes(x = gen, y = Freq, linetype = model, 
           group = as.factor(mutID), colour = as.factor(mutID))) + 
  geom_line(show.legend = F) +
  theme_bw()

# Plot effect over time
ggplot(d_seg_del %>% filter(mutID %in% eg_muts),
       aes(x = gen, y = s, linetype = model, 
           group = as.factor(mutID), colour = as.factor(mutID))) + 
  geom_line(show.legend = F) +
  theme_bw()

# Get the difference in selection coefficient after optimum shift
# as well as mean frequency just prior to optimum shift
d_seg_del %>%
  filter(gen == 49500 | gen == 50000) %>%
  arrange(gen) %>%
  group_by(seed, model, mutID) %>%
  filter(n() > 1) %>%
  summarise(diff_s = s[2] - s[1],
            Freq = Freq[1]) -> d_del_diffs


# d_muts_adapted %>% 
#   filter(interaction(seed, modelindex, mutID) %in% 
#                        interaction(d_fix_del$seed, 
#                                    d_fix_del$modelindex, 
#                                    d_fix_del$mutID)) -> d_muts_del
# 
# d_fix_del$wAA = calcAddFitness(d_fix_del$AA_pheno, 1, 0.05)
# d_fix_del$wAa = calcAddFitness(d_fix_del$Aa_pheno, 1, 0.05)
# d_fix_del$waa = calcAddFitness(d_fix_del$aa_pheno, 1, 0.05)
# 
# d_fix_del <- CalcDominance(d_fix_del)



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

# d_rank_combined_tbl <- d_ranked_combined %>% filter(rank != 0, modelindex == 2) %>%
#   group_by(rank, modelindex, isAdapted, mutType) %>% 
#   summarise(n = n()) %>% ungroup() %>% 
#   complete(rank, isAdapted, mutType)
# d_rank_combined_tbl <- d_rank_combined_tbl[-c(15, 16),]
# d_rank_combined_tbl
# 
# d_fix_ranked_combined %>% filter(modelindex == 2) %>%
#   mutate(aZbZ = fixEffectSum_aZ/fixEffectSum_bZ) %>%
#   ungroup() -> d_ratio

d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(molCompDiff = sum(abs(2 * value_aZ), na.rm = T) - sum(abs(2 * value_bZ), na.rm = T),
         molCompCorr = cor(abs(2 * value_aZ), abs(2 * value_bZ))) %>%
  ungroup() %>%
  dplyr::select(seed, molCompDiff, molCompCorr) %>%
  distinct(seed, .keep_all = T) -> d_molCompDiff

d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(evoBybZ = all(value_aZ == 0, na.rm = T),
         evoByaZ = all(value_bZ == 0, na.rm = T)) %>%
  ungroup() %>%
  distinct(seed, .keep_all = T) %>% 
  summarise(percEvoByaZ = mean(evoByaZ),
            percEvoBybZ = mean(evoBybZ),
            percEvoByBoth = 1 - (percEvoByaZ + percEvoBybZ),
            countEvoByaZ = sum(evoByaZ),
            countEvoBybZ = sum(evoBybZ),
            countEvoByBoth = n() - (countEvoByaZ + countEvoBybZ))
