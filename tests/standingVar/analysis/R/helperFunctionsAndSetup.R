# Load packages
packageList <- c("data.table", "dtplyr", "dplyr", "ggplot2", "tibble", "tidyr", "cowplot", "ggridges", 
"ggpmisc", "deSolve", "DescTools", "paletteer", "latex2exp", "readr", "RColorBrewer",
"ggh4x", "progress")

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
  
  if (any(dat$model == "ODE")) {
    dat$fixEffectSum_5 <- 1
    dat$fixEffectSum_6 <- 1
    dat$value_5 <- 0
    dat$value_6 <- 0
    }
    
  
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
  data.table::fwrite(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  data.table::fwrite(d_dat_withFX %>% ungroup() %>% 
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
# muts is a N x 4 dataframe of mutations with gen, seed, modelindex for matching 
# to dat_fixed where N is the number of mutants to test
# will have to do bootstrap methods, since there are a lot of mutations,
# and we're really only interested in an average
# n is the number of iterations to do, m = number of mutations to sample each iteration
PairwiseEpistasisAdditive <- function(dat_fixed, muts, n = 10, m = 100) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)
  
  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex) %>%
    summarise(fixEffectSum = 2 * sum(value)) %>%
    select(gen, seed, modelindex, fixEffectSum) %>%
    ungroup()
  
  # output dataframe: number of generations/seeds/modelindices * iterations
  output_len <- nrow(dat) * n
  out <- tibble(gen = numeric(output_len),
                seed = rep(dat$seed, each = n),
                modelindex = rep(dat$modelindex, each = n),
                meanEP = numeric(output_len),
                meanEW = numeric(output_len),
                seEP = numeric(output_len),
                seEW = numeric(output_len))
                
  # Iterate bootstrap
  i = 1
  while (i <= n) {
    # Sample m mutations from each of muts for a and b: note: chance to sample the
    # same mutation twice, so some epistasis might be dominance: chance is low though,
    # p = 2 * (1 - ((m-1)/m)^2 - 2 * 1/m * ((m-1)/m))
    # for m = 100, p = 0.0002: will probably happen sometimes, but rarely
    a <- muts %>% group_by(gen, seed, modelindex) %>% slice_sample(n = m)
    b <- muts %>% group_by(gen, seed, modelindex) %>% slice_sample(n = m)
    
    # Join a and b and add fixed effects
    result <- a %>% inner_join(., b, by = c("gen", "seed", "modelindex"),
                               relationship = "many-to-many")
    result <- result %>% rename(a = value.x, b = value.y)
    result <- inner_join(result, dat, by = c("gen", "seed", "modelindex"))

    # Calculate phenotype and fitness effects
    Pw <- result$fixEffectSum
    Pa <- result$fixEffectSum + result$a
    Pb <- result$fixEffectSum + result$b
    Pab <- result$fixEffectSum + result$a + result$b
    wa <- calcAddFitness(Pa, 2, 0.05)
    wb <- calcAddFitness(Pb, 2, 0.05)
    wab <- calcAddFitness(Pab, 2, 0.05)
    
    # Epistasis (fitness and trait)
    result$ew <- log(wab) - log(wa) - log(wb)
    # ep <- Pab - (dat$fixEffectSum + dat$a + dat$b)  
    # e_p = (effect due to ab) - (effect due to a + effect due to b)
    result$ep <- ( Pab - Pw ) - ( ( Pa - Pw ) + ( Pb - Pw ) ) # should always be zero for additive
    
    # Account for floating point error
    # ep[ep != 0] <- 0
    
  
    # Calculate mean and se for this iteration
    # put into output vector
    thisIterRange <- ( (i-1) * nrow(dat) + 1 ):( i * nrow(dat) )
    out[thisIterRange,] <- result %>%
      group_by(gen, seed, modelindex) %>%
      summarise(    mean_ep <- mean(ep),
                    mean_ew <- mean(ew),
                    se_ep <- se(ep),
                    se_ew <- se(ew))
    i <- i + 1
  }
  
  return(out)
}

# Calculate pairwise epistasis between combinations of network mutational effects
# relative to their wildtype background
# muts is a N x 4 dataframe of mutations with gen, seed, modelindex for matching 
# to dat_fixed where N is the number of mutants to test
# will have to do bootstrap methods, since there are a lot of mutations,
# and we're really only interested in an average
# n is the number of iterations to do, m = number of mutations to sample each iteration
PairwiseEpistasisNAR <- function(dat_fixed, muts, n = 10, m = 100) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)
  
  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex, mutType) %>%
    summarise(fixEffectSum = exp(2 * sum(value))) %>%
    select(gen, seed, modelindex, mutType, fixEffectSum) %>%
    ungroup()
  
  nMutTypes <- nlevels(muts$mutType) # could be 2 or 4, depending on model
  nPossibleMutTypes <- 4

  # output dataframe: number of generations/seeds/modelindices * iterations
  output_len <- nrow(dat) * n * nMutTypes
  out <- tibble(gen = numeric(output_len),
                seed = rep(dat$seed, each = n * nMutTypes),
                modelindex = rep(dat$modelindex, each = n * nMutTypes),
                mutType_ab = rep(as.character(dat$mutType), each = n * nMutTypes),
                meanEP = numeric(output_len),
                meanEW = numeric(output_len),
                seEP = numeric(output_len),
                seEW = numeric(output_len))
  
  # Pivot wider for easier access to fixed effects for the result vector
  dat <- dat %>% 
    pivot_wider(names_from = mutType, values_from = fixEffectSum,
                names_glue = "{.value}_{mutType}", values_fill = 1)
  
  
  # Iterate bootstrap
  j = 1
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent)", total = n)
  pb$tick(0)
  while (j <= n) {
    # Sample m mutations from each of muts for a and b: note: chance to sample the
    # same mutation twice, so some epistasis might be dominance: chance is low though,
    # p = 2 * (1 - ((m-1)/m)^2 - 2 * 1/m * ((m-1)/m))
    # for m = 100, p = 0.0002: will probably happen sometimes, but rarely
    a <<- muts %>% group_by(gen, seed, modelindex, mutType) %>% slice_sample(n = m/nMutTypes)
    b <<- muts %>% group_by(gen, seed, modelindex, mutType) %>% slice_sample(n = m/nMutTypes)
    
    # Join a and b and add fixed effects
    result <- a %>% inner_join(., b, by = c("gen", "seed", "modelindex"),
                               relationship = "many-to-many")
    result <- result %>% rename(a = value.x, b = value.y,
                                mutType_a = mutType.x, mutType_b = mutType.y)
    result <- inner_join(result, dat, by = c("gen", "seed", "modelindex")) %>%
      select(gen, seed, modelindex, mutType_a, mutType_b, a, b, 
             starts_with("fixEffectSum")) %>%
      mutate(mutType_ab = paste(mutType_a, mutType_b, sep = "_"))
  
  abNames <- paste(c(rep("a", times = 4), rep("b", times = 4)), 3:6, sep = "_")
  result[,abNames] <- 0
  
  # if this is an ODE/not K model, we need to add on fixed effects for KZ and KXZ
  # TODO: handle this better
  if (any(dat_fixed$model == "ODE")) {
    result$fixEffectSum_5 <- 1
    result$fixEffectSum_6 <- 1
  }
  
  # initialize a and b values for the right molecular component
  for (i in (1:nPossibleMutTypes) + 2) {
    result[result$mutType_a == paste0(i), paste0("a_", i)] <- result[result$mutType_a == paste0(i), "a"]
    result[result$mutType_b == paste0(i), paste0("b_", i)] <- result[result$mutType_b == paste0(i), "b"]
  }
  
  # Add on a, b, ab to the base effect
  for (i in (1:nPossibleMutTypes) + 2) {
    result[, paste0("a_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("a_", i)])
    result[, paste0("b_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("b_", i)])
    result[, paste0("ab_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("a_", i)] + result[,paste0("b_", i)])
  }
  result$rowID <- as.integer(rownames(result))
  result$mutGroup <- rep(1:(nrow(result)/m), each = m)
  
  # Split the result into wt, a, b, and ab to reduce non-unique solutions
  d_wildtype <- result %>%
    group_by(gen, seed, modelindex) %>%
    filter(row_number() == 1) %>%
    select(gen, seed, modelindex, rowID, starts_with("fixEffectSum")) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 1)))
    
  d_a <- result %>%
    group_by(gen, seed, modelindex) %>%
    distinct(., mutGroup, .keep_all = T) %>%
    select(gen, seed, modelindex,
           rowID, a_molComp_3, a_molComp_4, a_molComp_5, a_molComp_6) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 2)))

  d_b <- result %>%
    group_by(gen, seed, modelindex) %>%
    #distinct(., mutGroup, .keep_all = T) %>%
    select(gen, seed, modelindex, 
           rowID, b_molComp_3, b_molComp_4, b_molComp_5, b_molComp_6) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 3)))
  
  
  d_ab <- result %>%
    group_by(gen, seed, modelindex) %>%
    select(gen, seed, modelindex, 
           rowID, ab_molComp_3, ab_molComp_4, ab_molComp_5, ab_molComp_6) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 4)))

  # Remove column names so we can join them
  colnames(d_wildtype) <- paste0("v", 1:5)
  colnames(d_a) <- paste0("v", 1:5)
  colnames(d_b) <- paste0("v", 1:5)
  colnames(d_ab) <- paste0("v", 1:5)

  d_landscaper <- rbind(d_wildtype, d_a, d_b, d_ab)
  
  # Run landscaper
  data.table::fwrite(d_landscaper, 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_phenos <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 16, TRUE)
  

  # Ensure that the tables are aligned by id before we join them
  result <- result %>% arrange(rowID)

  # Separate phenos by the rowID value that we assigned earlier and revert to
  # original id
  d_wildtype <- d_phenos %>% filter(id %% 10 == 1) %>% 
    mutate(id = (id - 1) / 10) %>% arrange(id)
  d_a <- d_phenos %>% filter(id %% 10 == 2) %>% 
    mutate(id = (id - 2) / 10) %>% arrange(id)
  d_b <- d_phenos %>% filter(id %% 10 == 3) %>% 
    mutate(id = (id - 3) / 10) %>% arrange(id)
  d_ab <- d_phenos %>% filter(id %% 10 == 4) %>% 
    mutate(id = (id - 4) / 10) %>% arrange(id)
  
  # Fill out results to match ab size
  d_a <- d_a[rep(seq_len(nrow(d_a)), each = m), ]
  d_a$id <- d_ab$id
  
  d_b <- d_b[rep(seq_len(nrow(d_b)), each = m), ]
  d_b$id <- d_ab$id
  
  d_wildtype <- d_wildtype[rep(seq_len(nrow(d_wildtype)), each = m*m), ]
  d_wildtype$id <- d_ab$id
  
  
  # Epistasis (fitness and trait)
  result$ew <- log(d_ab$fitness) - log(d_a$fitness) - log(d_b$fitness)
  result$ep <- ( d_ab$pheno - d_wildtype$pheno ) - ( ( d_a$pheno - d_wildtype$pheno ) + ( d_b$pheno - d_wildtype$pheno ) ) 
  
  thisIterRange <- ( (j-1) * (output_len / (n)) + 1 ):( j * (output_len / (n)) )
  out[thisIterRange,] <- result %>%
    group_by(gen, seed, modelindex, mutType_ab) %>%
    summarise(    mean_ep <- mean(ep),
                  mean_ew <- mean(ew),
                  se_ep <- se(ep),
                  se_ew <- se(ew), .groups = "keep")
  pb$tick(1)
  j <- j + 1
  }
  return(out)
}

# Calculates the site frequency spectra for mutations
CalcSFS <- function(dat) {
  dat$FreqBin <- cut(dat$Freq, breaks = 10)
  dat$optPerc <- dat$phenomean - 1
  dat$optPerc <- cut(dat$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))
  
  dat %>% 
    select(optPerc, seed, modelindex, model, nloci, r, tau, 
           mutID, mutType, value, Freq, FreqBin)
}
