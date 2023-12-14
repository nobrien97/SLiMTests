# Load packages
packageList <- c("data.table", "dtplyr", "dplyr", "ggplot2", "tibble", "tidyr", "cowplot", "ggridges", 
"ggpmisc", "deSolve", "DescTools", "paletteer", "latex2exp", "readr", "RColorBrewer",
"ggh4x")

lapply(packageList, require, character.only = T)

# Variables
mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$"),
  TeX("$K_Z$"),
  TeX("$K_{XZ}$")
)

# Misc functions
se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

# Adds the parameter combination to a dataframe
AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           nloci = d_combos$nloci[as.numeric(levels(modelindex))[modelindex]],
           tau = d_combos$tau[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]])
}

# Fitness effect calculation functions
## Gaussian fitness function, used for both models
calcAddFitness <- function(phenotypes, optimum, width) {
  dists <- (phenotypes - optimum)^2
  return(exp(-(dists * width)))
}

## Calculate the fitness effects in additive populations
### Fitness/phenotype is calculated relative to "wildtype"
### an individual with no mutations at all OR if there are fixations,
### those count to the wildtype
CalcAddEffects <- function(dat, dat_fixed) {
  # If we are calculating fitness for segregating sites, need to evaluate fitness
  # vs the fixed effect background at the timepoint (i.e. all fixations at timepoint gen)
  # Fixation effect is multiplied by 2 because diploid
  dat <- as.data.table(dat)
  dat_fixed <- as.data.table(dat_fixed)
  
  dat_fixed <- dat_fixed %>%
    group_by(gen, seed, modelindex) %>%
    summarise(fixEffectSum = 2 * sum(value)) %>%
    select(gen, seed, modelindex, fixEffectSum) %>%
    ungroup()
  
  dat <- dat %>% inner_join(dat_fixed, by = c("gen", "seed", "modelindex"))
  
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = fixEffectSum + 2 * value
  # Aa = fixEffectSum + value
  # aa = fixEffectSum
  # Get effect
  Aa <- calcAddFitness(dat$fixEffectSum + dat$value, 2, 0.05)
  AA <- calcAddFitness(dat$fixEffectSum + dat$value * 2, 2, 0.05)
  aa <- calcAddFitness(dat$fixEffectSum, 2, 0.05)
  
  dat %>% 
    mutate(AA_pheno = fixEffectSum + 2 * value,
           Aa_pheno = fixEffectSum + value,
           aa_pheno = fixEffectSum,
           value_AA = value * 2,
)
  dat$avFit <- Aa - aa
  dat$avFit_AA <- AA - aa
  dat$value_AA <- dat$value * 2
  dat$wAA <- AA
  dat$wAa <- Aa
  dat$waa <- aa
  dat$s <- AA - aa
  
  return(as_tibble(dat))
}

## Run the ODELandscaper tool to evaluate phenotype and fitness
## for many individuals at once.
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

## Calculate fitness in network models
CalcNARPhenotypeEffects <- function(dat, dat_fixed) {
  # calculate cumulative molecular component values at each step due to only 
  # fixed effects
  # multiply by 2 because diploid
  dat <- as.data.table(dat)
  dat_fixed <- as.data.table(dat_fixed)
  
  dat_fixed <- dat_fixed %>%
    group_by(gen, seed, modelindex, mutType) %>%
    summarise(fixEffectSum = exp(2 * sum(value))) %>%
    select(gen, seed, modelindex, mutType, fixEffectSum) %>%
    ungroup()
  
  # Pivot fixEffectSums and values 
  dat_fixed <- dat_fixed %>% 
    pivot_wider(names_from = mutType, values_from = fixEffectSum,
                names_glue = "{.value}_{mutType}", values_fill = 1)
  
  # HACK ALERT! HACK ALERT!
  ## use a temporary mutType2 column so we keep mutType when pivot_wider does
  ## its thing
  dat <- dat %>%
    mutate(mutType2 = mutType) %>%
    pivot_wider(names_from = mutType2, values_from = value,
                names_glue = "{.value}_{mutType2}", values_fill = 0)
  
  dat <- dat %>% inner_join(dat_fixed, 
                            by = c("gen", "seed", "modelindex"))
  
  dat$rowID <- as.integer(rownames(dat))
  
  # Get phenotypes without the mutation
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, fixEffectSum_3, fixEffectSum_4, fixEffectSum_5, fixEffectSum_6), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_popfx <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Segregating mutation calculations
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = d_popfx$fixEffectSum + 2 * value
  # Aa = d_popfx$fixEffectSum + value
  # aa = d_popfx$fixEffectSum
  
  # Add on the segregating effect to the fixation effects
  d_dat_withFX <- dat
  d_dat_withFX$aZ <- exp(log(d_dat_withFX$fixEffectSum_3) + d_dat_withFX$value_3)
  d_dat_withFX$bZ <- exp(log(d_dat_withFX$fixEffectSum_4) + d_dat_withFX$value_4)
  d_dat_withFX$KZ <- exp(log(d_dat_withFX$fixEffectSum_5) + d_dat_withFX$value_5)
  d_dat_withFX$KXZ <- exp(log(d_dat_withFX$fixEffectSum_6) + d_dat_withFX$value_6)
  
  # homozygous effect
  d_dat_withFX$aZ_AA <- exp(log(d_dat_withFX$fixEffectSum_3) + 2 * d_dat_withFX$value_3)
  d_dat_withFX$bZ_AA <- exp(log(d_dat_withFX$fixEffectSum_4) + 2 * d_dat_withFX$value_4)
  d_dat_withFX$KZ_AA <- exp(log(d_dat_withFX$fixEffectSum_5) + 2 * d_dat_withFX$value_5)
  d_dat_withFX$KXZ_AA <- exp(log(d_dat_withFX$fixEffectSum_6) + 2 * d_dat_withFX$value_6)
  
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
  dat$s <- dat$wAA - dat$waa
  return(dat)
}

# Calculate pairwise epistasis between combinations of additive mutational effects
# relative to their wildtype background
# a and b are n x 4 dataframes with gen, seed, modelindex for matching to dat_fixed
PairwiseEpistasisAdditive <- function(dat_fixed, a, b) {
  
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)
  
  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex) %>%
    summarise(fixEffectSum = 2 * sum(value)) %>%
    select(gen, seed, modelindex, fixEffectSum) %>%
    ungroup()
  
  # Add on a and b columns
  dat <- dat %>%
    inner_join(., a, by = c("gen", "seed", "modelindex"), 
               relationship = "many-to-many") %>%
    inner_join(., b, by = c("gen", "seed", "modelindex"),
               relationship = "many-to-many")
  
  # Calculate phenotype and fitness effects
  Pw <- dat$fixEffectSum
  Pa <- dat$fixEffectSum + dat$a
  Pb <- dat$fixEffectSum + dat$b
  Pab <- dat$fixEffectSum + dat$a + dat$b
  wa <- calcAddFitness(Pa, 2, 0.05)
  wb <- calcAddFitness(Pb, 2, 0.05)
  wab <- calcAddFitness(Pab, 2, 0.05)
  
  # Epistasis (fitness and trait)
  ew <- log(wab) - log(wa) - log(wb)
  # ep <- Pab - (dat$fixEffectSum + dat$a + dat$b)  
  # e_p = (effect due to ab) - (effect due to a + effect due to b)
  ep <- ( Pab - Pw ) - ( ( Pa - Pw ) + ( Pb - Pw ) ) # should always be zero for additive
  
  # Account for floating point error
  # ep[ep != 0] <- 0
    
  return(tibble(gen = dat$gen,
                seed = dat$seed,
                modelindex = dat$modelindex,
                a = dat$a,
                b = dat$b,
                ew = ew,
                ep = ep))  
}

# Calculate pairwise epistasis between combinations of network mutational effects
# relative to their wildtype background
# a and b are n x 5 dataframes with gen, seed, modelindex, mutType for 
# matching to dat_fixed
PairwiseEpistasisNAR <- function(dat_fixed, a, b) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)
  
  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex, mutType) %>%
    summarise(fixEffectSum = exp(2 * sum(value))) %>%
    select(gen, seed, modelindex, mutType, fixEffectSum) %>%
    ungroup()
  
  dat <- dat %>% 
    pivot_wider(names_from = mutType, values_from = fixEffectSum,
                names_glue = "{.value}_{mutType}", values_fill = 1)
  
  
  # Add on a and b columns
  dat <- dat %>%
    inner_join(., a, by = c("gen", "seed", "modelindex"), 
               relationship = "many-to-many") %>%
    inner_join(., b, by = c("gen", "seed", "modelindex"),
               relationship = "many-to-many") %>%
    rename(mutType_a = mutType.x,
           mutType_b = mutType.y)
  
  abNames <- paste(c(rep("a", times = 4), rep("b", times = 4)), 3:6, sep = "_")
  dat[,abNames] <- 0
  
  # initialize a and b values for the right molecular component
  for (i in unique(dat$mutType_a)) {
    dat[dat$mutType_a == paste0(i), paste0("a_", i)] <- dat[dat$mutType_a == paste0(i), "a"]
    dat[dat$mutType_b == paste0(i), paste0("b_", i)] <- dat[dat$mutType_b == paste0(i), "b"]
  }
  
  # Add on a, b, ab to the base effect
  for (i in unique(dat$mutType_a)) {
    dat[, paste0("value_a_", i)] <- exp(log(dat[,paste0("fixEffectSum_", i)]) + dat[,paste0("a_", i)])
    dat[, paste0("value_b_", i)] <- exp(log(dat[,paste0("fixEffectSum_", i)]) + dat[,paste0("b_", i)])
    dat[, paste0("value_ab_", i)] <- exp(log(dat[,paste0("fixEffectSum_", i)]) + dat[,paste0("a_", i)] + dat[,paste0("b_", i)])
  }
  
  dat$rowID <- as.integer(rownames(dat))
  
  # Run wildtype phenotypes
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, fixEffectSum_3, fixEffectSum_4, fixEffectSum_5, fixEffectSum_6), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_wildtype <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)

  # a phenotypes
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, value_a_3, value_a_4, value_a_5, value_a_6), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_a <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # b phenotypes
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, value_b_3, value_b_4, value_b_5, value_b_6), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_b <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # ab phenotypes
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, value_ab_3, value_ab_4, value_ab_5, value_ab_6), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_ab <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Ensure that the tables are aligned by id before we join them
  dat <- dat %>% arrange(rowID)
  d_a <- d_a %>% arrange(id)
  d_b <- d_b %>% arrange(id)
  d_ab <- d_ab %>% arrange(id)
  d_wildtype <- d_wildtype %>% arrange(id)
  
  # Epistasis (fitness and trait)
  ew <- log(d_ab$fitness) - log(d_a$fitness) - log(d_b$fitness)
  ep <- ( d_ab$pheno - d_wildtype$pheno ) - ( ( d_a$pheno - d_wildtype$pheno ) + ( d_b$pheno - d_wildtype$pheno ) ) 

    return(tibble(gen = dat$gen,
                seed = dat$seed,
                modelindex = dat$modelindex,
                mutType_a = dat$mutType_a,
                mutType_b = dat$mutType_b,
                a = dat$a,
                b = dat$b,
                ew = ew,
                ep = ep))  
}

# Calculates the site frequency spectra for mutations
CalcSFS <- function(dat) {
  dat$FreqBin <- cut(dat$Freq, breaks = 10)
  dat$optPerc <- dat$phenomean - 1
  dat$optPerc <- cut(dat$optPerc, c(-Inf, 0.25, 0.5, 0.75, 0.95, Inf))
  
  dat %>% 
    select(optPerc, seed, modelindex, model, nloci, r, tau, 
           mutID, mutType, value, Freq, FreqBin)
}
