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

CalcPhenotypeEffects <- function(dat, dat_fixed) {
  if (nrow(dat_fixed) == 0) {
    dat_fixed[1,] <- dat[1,]
    dat_fixed %>% mutate(value = 0)
  }
  
  if (dat[1, "model"] == "Add") {
    return(CalcAddEffects(dat, dat_fixed))
  }
  
  CalcNARPhenotypeEffects(dat, dat_fixed)
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
  command <- "~/Tools/odeLandscape/ODELandscaper -i %s -o ./%s -w %f -p %f -t %i"
  #command <- "ODELandscaper -i %s -o ./%s -w %f -p %f -t %i"
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
  # calculate cumulative molecular component values due to fixations,
  # add on the sampled mutation and recalculate phenotype
  # multiply by 2 because diploid
  dat <- as.data.table(dat)
  dat_fixed <- as.data.table(dat_fixed)
  model <- as.character(dat$modelindex)[1]

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
              paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  d_popfx <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 0.05, 2, 4, TRUE)
  
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
              paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 0.05, 2, 4, TRUE)

  data.table::fwrite(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ_AA, bZ_AA, KZ, KXZ), 
              paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  AA <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 0.05, 2, 4, TRUE)

  
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

PairwiseFitnessRank <- function(dat_fixed, muts, A_ids, B_ids) {
  if (muts[1, "model"] == "Add") {
    return(PairwiseFitnessRankAdditive(dat_fixed, muts, A_ids, B_ids))
  }
  
  PairwiseFitnessRankNAR(dat_fixed, muts, A_ids, B_ids)
}

PairwiseEpistasis <- function(dat_fixed, muts, n = 1000, m = 10, 
                              returnAverage = F, weightABByFreq = F) {
  if (dat_fixed[1, "model"] == "Add") {
    return(PairwiseEpistasisAdditive(dat_fixed, muts, n, m, 
                                     returnAverage, weightABByFreq))
  }
  
  PairwiseEpistasisNAR(dat_fixed, muts, n, m, 
                       returnAverage, weightABByFreq)
}


# Calculate pairwise epistasis between combinations of additive mutational effects
# relative to their wildtype background
# muts is a N x 4 dataframe of mutations with gen, seed, modelindex for matching 
# to dat_fixed where N is the number of mutants to test
# will have to do bootstrap methods, since there are a lot of mutations,
# and we're really only interested in an average
# n is the number of iterations to do, m = number of mutations to sample each iteration
PairwiseEpistasisAdditive <- function(dat_fixed, muts, n = 1000, m = 10, 
                                      returnAverage = F, weightABByFreq = F) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)
  
  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex) %>%
    summarise(fixEffectSum = 2 * sum(value)) %>%
    select(gen, seed, modelindex, fixEffectSum) %>%
    ungroup()
  
  # output dataframe: number of generations/seeds/modelindices * iterations
  output_len <- nrow(dat) * n * m
  out <- tibble(gen = numeric(output_len),
                seed = rep(dat$seed, each = n * m),
                modelindex = rep(dat$modelindex, each = n * m),
                wa = numeric(output_len),
                wb = numeric(output_len),
                wab = numeric(output_len),
                Pwt = numeric(output_len),
                Pa = numeric(output_len),
                Pb = numeric(output_len),
                Pab = numeric(output_len),
                ew = numeric(output_len),
                ep = numeric(output_len))
                
  # Iterate bootstrap
  i = 1
  pb <- progress::progress_bar$new(
    format = "[:bar] :current/:total (:percent eta: :eta)", total = n)
  pb$tick(0)
  while (i <= n) {
    # Sample m mutations from each of muts for a and b: note: chance to sample the
    # same mutation twice, so some epistasis might be dominance: chance is low though,
    # p = 2 * (1 - ((m-1)/m)^2 - 2 * 1/m * ((m-1)/m))
    # for m = 100, p = 0.0002: will probably happen sometimes, but rarely
    a <- muts %>% group_by(gen, seed, modelindex) %>% 
      slice_sample(n = m, replace = T,
                   weight_by = case_when(weightABByFreq == T ~ freq, 
                                         weightABByFreq == F ~ rep(1, times = n()))) %>%
      rename(a = value)
    b <- muts %>% group_by(gen, seed, modelindex) %>% 
      slice_sample(n = m, replace = T,
                   weight_by = case_when(weightABByFreq == T ~ freq, 
                                         weightABByFreq == F ~ rep(1, times = n()))) %>%
      rename(b = value)
    
    # Join a and b and add fixed effects
    result <- a %>% select(-freq)
    result$b <- b$b
    result <- inner_join(result, dat, by = c("gen", "seed", "modelindex"))

    # Calculate phenotype and fitness effects
    result$Pwt <- result$fixEffectSum
    result$Pa <- result$fixEffectSum + result$a
    result$Pb <- result$fixEffectSum + result$b
    result$Pab <- result$fixEffectSum + result$a + result$b
    result$wa <- calcAddFitness(result$Pa, 2, 0.05)
    result$wb <- calcAddFitness(result$Pb, 2, 0.05)
    result$wab <- calcAddFitness(result$Pab, 2, 0.05)
    
    # Epistasis (fitness and trait)
    result$ew <- log(result$wab) - log(result$wa) - log(result$wb)
    result$ep <- ( result$Pab - result$Pwt ) - 
      ( ( result$Pa - result$Pwt ) + ( result$Pb - result$Pwt ) ) # should always be zero for additive
        
  
    # Calculate mean and se for this iteration
    # put into output vector
    thisIterRange <- ( (i-1) * (nrow(dat) * m) + 1 ):( i * (nrow(dat) * m) )
    out[thisIterRange,] <- result %>%
      select(gen, seed, modelindex, wa, wb, wab, Pwt, Pa, Pb, Pab, ew, ep)
    pb$tick(1)
    i <- i + 1
  }
    
  if (returnAverage) {
    out <- out %>%
      group_by(gen, seed, modelindex) %>%
      summarise(meanEP = mean(ep),
                meanEW = mean(ew),
                meanPwt = mean(Pwt),
                meanPa = mean(Pa),
                meanPb = mean(Pb),
                meanPab = mean(Pab),
                meanwa = mean(wa),
                meanwb = mean(wb),
                meanwab = mean(wab),
                sdEP = sd(ep),
                sdEW = sd(ew),
                sdPwt = sd(Pwt),
                sdPa = sd(Pa),
                sdPb = sd(Pb),
                sdPab = sd(Pab),
                sdwa = sd(wa),
                sdwb = sd(wb),
                sdwab = sd(wab))
  }
  
  return(out)
}

# Calculate pairwise epistasis between combinations of network mutational effects
# relative to their wildtype background
# muts is a N x 4 dataframe of mutations with gen, seed, modelindex for matching 
# to dat_fixed where N is the number of mutants to test
# will have to do bootstrap methods, since there are a lot of mutations
# calculate average over the 10 sampled mutations
# n is the number of iterations to do, m = number of mutations to sample each iteration
# returnAverage = T will return a dataframe of the average epistasis over the n 
# repetitions and m mutations per repetition. When returnAverage = F, the result is
# not averaged over the n repetitions (i.e. each replicate gets its own average over
# the m mutations)
PairwiseEpistasisNAR <- function(dat_fixed, muts, n = 1000, m = 10, 
                                 returnAverage = F, weightABByFreq = F) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)

  model <- as.character(dat_fixed$modelindex)[1]
  
  fixEffectDat <- dat_fixed %>%
    group_by(gen, seed, modelindex, mutType) %>%
    summarise(fixEffectSum = exp(2 * sum(value))) %>%
    select(gen, seed, modelindex, mutType, fixEffectSum) %>%
    ungroup()
  
  nMutTypes <- length(unique(muts$mutType)) # could be 2 or 4, depending on model
  nPossibleMutTypes <- 4                    # always 4

  # output dataframe: number of generations/seeds/modelindices * iterations
  output_len <- nrow(fixEffectDat) / nMutTypes * n * m
  output_len_each <- ( n * m ) / nMutTypes
  
  out <- tibble(gen = numeric(output_len),
                seed = rep(fixEffectDat$seed, each = output_len_each),
                modelindex = rep(fixEffectDat$modelindex, each = output_len_each),
                mutType_ab = rep(as.character(fixEffectDat$mutType), each = output_len_each),
                wa = numeric(output_len),
                wb = numeric(output_len),
                wab = numeric(output_len),
                Pwt = numeric(output_len),
                Pa = numeric(output_len),
                Pb = numeric(output_len),
                Pab = numeric(output_len),
                ew = numeric(output_len),
                ep = numeric(output_len))
  
  # Pivot wider for easier access to fixed effects for the result vector
  fixEffectDat <- fixEffectDat %>% 
    pivot_wider(names_from = mutType, values_from = fixEffectSum,
                names_glue = "{.value}_{mutType}", values_fill = 1)
  
  
  # Iterate bootstrap
  j = 1
  pb <- progress::progress_bar$new(
   format = "[:bar] :current/:total (:percent eta: :eta)", total = n)
  pb$tick(0)
  while (j <= n) {
    # Sample m mutations from each of muts for a and b: note: chance to sample the
    # same mutation twice, so some epistasis might be dominance: chance is low though,
    # p = 2 * (1 - ((m-1)/m)^2 - 2 * 1/m * ((m-1)/m))
    # for m = 100, p = 0.0002: will probably happen sometimes, but rarely
    a <- muts %>% group_by(gen, seed, modelindex) %>% 
      slice_sample(n = m, replace = T,
                   weight_by = case_when(weightABByFreq == T ~ freq, 
                                         weightABByFreq == F ~ rep(1, times = n())))
    b <- muts %>% group_by(gen, seed, modelindex) %>% 
      slice_sample(n = m, replace = T,
                   weight_by = case_when(weightABByFreq == T ~ freq, 
                                         weightABByFreq == F ~ rep(1, times = n())))

    
    # Join a and b and add fixed effects
    result <- a %>% select(-freq)
    result <- result %>% rename(a = value, mutType_a = mutType)
    result$mutType_b <- b$mutType
    result$b <- b$value
    result <- inner_join(result, fixEffectDat, by = c("gen", "seed", "modelindex")) %>%
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
    
    # Split the result into wt, a, b, and ab to reduce non-unique solutions
    d_wildtype <- result %>%
      group_by(gen, seed, modelindex) %>%
      filter(row_number() == 1) %>%
      select(gen, seed, modelindex, rowID, starts_with("fixEffectSum")) %>%
      ungroup() %>% select(!(gen:modelindex)) %>%
      mutate(rowID = as.numeric(paste0(rowID, 1)))
      
    d_a <- result %>%
      group_by(gen, seed, modelindex) %>%
      select(gen, seed, modelindex,
             rowID, a_molComp_3, a_molComp_4, a_molComp_5, a_molComp_6) %>%
      ungroup() %>% select(!(gen:modelindex)) %>%
      mutate(rowID = as.numeric(paste0(rowID, 2)))
  
    # b is organised differently to a, need to calculate first m for each mutgroup
    # and repeat that for the remaining
    d_b <- result %>%
      group_by(gen, seed, modelindex) %>%
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
                paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
    d_phenos <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 0.05, 2, 4, TRUE)
    
  
    # Ensure that the tables are aligned by id before we join them
    result <- result %>% arrange(rowID)
  
    # Separate phenos by the rowID value that we assigned earlier and revert to
    # original id
    d_wildtype <- d_phenos %>% filter((id - 1) %% 10 == 0) %>% 
      mutate(id = (id - 1) / 10) %>% arrange(id)
    d_a <- d_phenos %>% filter((id - 2) %% 10 == 0) %>% 
      mutate(id = (id - 2) / 10) %>% arrange(id)
    d_b <- d_phenos %>% filter((id - 3) %% 10 == 0) %>% 
      mutate(id = (id - 3) / 10) %>% arrange(id)
    d_ab <- d_phenos %>% filter((id - 4) %% 10 == 0) %>% 
      mutate(id = (id - 4) / 10) %>% arrange(id)
    
    # Repeat wildtype to match the a b pairs
    d_wildtype <- d_wildtype[rep(seq_len(nrow(d_wildtype)), each = m), ]
    d_wildtype$id <- d_ab$id
    
    
    # Epistasis (fitness and trait)
    result$wa <- d_a$fitness
    result$wb <- d_b$fitness
    result$wab <- d_ab$fitness
    result$Pwt <- d_wildtype$pheno
    result$Pa <- d_a$pheno
    result$Pb <- d_b$pheno
    result$Pab <- d_ab$pheno
    
    result$ew <- log(result$wab) - log(result$wa) - log(result$wb)
    result$ep <- ( result$Pab - result$Pwt ) - 
      ( ( result$Pa - result$Pwt ) + ( result$Pb - result$Pwt ) ) 
    
    # Fill output: figure out which range of the output vector to fill
    thisIterRange <- ( (j-1) * (nrow(fixEffectDat) * m) + 1 ):( j * (nrow(fixEffectDat) * m) )
    
    out[thisIterRange,] <- result %>%
      ungroup() %>%
      select(gen, seed, modelindex, mutType_ab, wa, wb, wab, Pwt, Pa, Pb, Pab, ew, ep)

    pb$tick(1)
    j <- j + 1
  }
  
  # if we're averaging over iterations, do that now:
  # identical to grand mean since all replicates have the same sample size
  # Combined SEM via pythagorean theorem 
  if (returnAverage) {
    out <- out %>%
      group_by(gen, seed, modelindex, mutType_ab) %>%
      summarise(meanEP = mean(ep),
                meanEW = mean(ew),
                meanPwt = mean(Pwt),
                meanPa = mean(Pa),
                meanPb = mean(Pb),
                meanPab = mean(Pab),
                meanwa = mean(wa),
                meanwb = mean(wb),
                meanwab = mean(wab),
                sdEP = sd(ep),
                sdEW = sd(ew),
                sdPwt = sd(Pwt),
                sdPa = sd(Pa),
                sdPb = sd(Pb),
                sdPab = sd(Pab),
                sdwa = sd(wa),
                sdwb = sd(wb),
                sdwab = sd(wab))
  }

  return(out)
}

PairwiseFitnessRankAdditive <- function(dat_fixed, muts, A_ids, B_ids) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)
  
  if (nrow(dat_fixed) == 0) {
    dat_fixed <- dat_fixed %>% add_row(muts[1,])
    dat_fixed$value <- 0
  }

  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex) %>%
    summarise(fixEffectSum = 2 * sum(value)) %>%
    select(gen, seed, modelindex, fixEffectSum) %>%
    ungroup()
  
  
  # split mutations into A and B
  A_pos <- match(A_ids, muts$mutID)
  B_pos <- match(B_ids, muts$mutID)
  mutsA <- muts[A_pos[A_pos > 0 & B_pos > 0],]
  mutsB <- muts[B_pos[A_pos > 0 & B_pos > 0],]
  

  # output dataframe: ranking fitness of parental ab/AB
  # parab = parental alleles, parAB = derived alleles
  # solve for fitness of each genotype, which we can then use to 
  # rearrange the LD genotypes according to fitness (so ab lowest fitness, AB highest)
  output_len <- nrow(mutsA)
  out <- tibble(gen = dat$gen, # assumes there is only one gen/seed/modelindex
                seed = dat$seed,
                modelindex = dat$modelindex,
                mutIDA = numeric(output_len),
                mutIDB = numeric(output_len),
                wparab = numeric(output_len),
                wparaB = numeric(output_len),
                wparAb = numeric(output_len),
                wparAB = numeric(output_len)
                )
  
    # Calculate phenotype and fitness effects
    Pab <- dat$fixEffectSum
    PAb <- dat$fixEffectSum + mutsA$value
    PaB <- dat$fixEffectSum + mutsB$value
    PAB <- dat$fixEffectSum + mutsA$value + mutsB$value
    
    out$mutIDA <- mutsA$mutID
    out$mutIDB <- mutsB$mutID
    out$wparab <- calcAddFitness(Pab, 2, 0.05)
    out$wparaB <- calcAddFitness(PaB, 2, 0.05)
    out$wparAb <- calcAddFitness(PAb, 2, 0.05)
    out$wparAB <- calcAddFitness(PAB, 2, 0.05)
    
    # Filter NA results
    out <- out %>% filter(!is.na(mutIDA), !is.na(mutIDB))
    
    return(out)
}

PairwiseFitnessRankNAR <- function(dat_fixed, muts, A_ids, B_ids) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed) %>% distinct()
  muts <- muts %>% distinct()
  
  if (nrow(dat_fixed) == 0) {
    dat_fixed <- dat_fixed %>% add_row(muts[1,])
    dat_fixed$value <- 0
  }
  
  model <- paste0(as.character(dat_fixed$gen)[1], "_", 
                  as.character(dat_fixed$modelindex)[1], "_", 
                  as.character(dat_fixed$seed)[1])
  
  dat <- dat_fixed %>%
    group_by(gen, seed, modelindex, mutType) %>%
    summarise(fixEffectSum = exp(2 * sum(value))) %>%
    select(gen, seed, modelindex, mutType, fixEffectSum) %>%
    ungroup()
  
  nMutTypes <- length(unique(muts$mutType)) # could be 2 or 4, depending on model
  nPossibleMutTypes <- 4                    # always 4
  
  
  # split mutations into A and B
  A_pos <- match(A_ids, muts$mutID, nomatch = 0)
  B_pos <- match(B_ids, muts$mutID, nomatch = 0)
  
  # Remove pairs with at least one 0
  mutsA <- muts[A_pos[A_pos > 0 & B_pos > 0],]
  mutsB <- muts[B_pos[A_pos > 0 & B_pos > 0],]
  
  # output dataframe: ranking fitness of parental ab/AB
  # parab = parental alleles, parAB = derived alleles
  # solve for fitness of each genotype, which we can then use to 
  # rearrange the LD genotypes according to fitness (so ab lowest fitness, AB highest)
  output_len <- nrow(mutsA)
  out <- tibble(gen = dat$gen[1], # assumes there is only one gen/seed/modelindex
                seed = dat$seed[1],
                modelindex = dat$modelindex[1],
                mutType_ab = rep("", times = output_len),
                mutIDA = numeric(output_len),
                mutIDB = numeric(output_len),
                wparab = numeric(output_len),
                wparaB = numeric(output_len),
                wparAb = numeric(output_len),
                wparAB = numeric(output_len)
  )
  
  # Pivot wider for easier access to fixed effects for the result vector
  
  dat <- dat %>% 
    pivot_wider(names_from = mutType, values_from = fixEffectSum,
                names_glue = "{.value}_{mutType}", values_fill = 1)
    
  columns_to_add <- c(
    fixEffectSum_3 = 1,
    fixEffectSum_4 = 1,
    fixEffectSum_5 = 1,
    fixEffectSum_6 = 1
  )
  
  # Add other fixed effects (that might be missing)  
  dat <- dat %>%
    add_column(!!!columns_to_add[!names(columns_to_add) %in% names(.)])

  result <- mutsA
  result <- result %>% rename(a = value, mutType_a = mutType)
  result$mutType_b <- mutsB$mutType
  result$b <- mutsB$value
  result <- inner_join(result, dat, by = c("gen", "seed", "modelindex")) %>%
    select(gen, seed, modelindex, mutType_a, mutType_b, a, b, 
           starts_with("fixEffectSum")) %>%
    mutate(mutType_ab = paste(mutType_a, mutType_b, sep = "_"))
  
  abNames <- paste(c(rep("a", times = 4), rep("b", times = 4)), 3:6, sep = "_")
  result[,abNames] <- 0
  
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
  
  # Split the result into wt, a, b, and ab to reduce non-unique solutions
  d_wildtype <- result %>%
    group_by(gen, seed, modelindex) %>%
    filter(row_number() == 1) %>%
    select(gen, seed, modelindex, rowID, starts_with("fixEffectSum")) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 1)))
  
  d_a <- result %>%
    group_by(gen, seed, modelindex) %>%
    select(gen, seed, modelindex,
           rowID, a_molComp_3, a_molComp_4, a_molComp_5, a_molComp_6) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 2)))
  
  # b is organised differently to a, need to calculate first m for each mutgroup
  # and repeat that for the remaining
  d_b <- result %>%
    group_by(gen, seed, modelindex) %>%
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
                     paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  d_phenos <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 0.05, 2, 4, TRUE)
  
  
  # Ensure that the tables are aligned by id before we join them
  result <- result %>% arrange(rowID)
  
  # Separate phenos by the rowID value that we assigned earlier and revert to
  # original id
  d_wildtype <- d_phenos %>% filter((id - 1) %% 10 == 0) %>% 
    mutate(id = (id - 1) / 10) %>% arrange(id)
  d_a <- d_phenos %>% filter((id - 2) %% 10 == 0) %>% 
    mutate(id = (id - 2) / 10) %>% arrange(id)
  d_b <- d_phenos %>% filter((id - 3) %% 10 == 0) %>% 
    mutate(id = (id - 3) / 10) %>% arrange(id)
  d_ab <- d_phenos %>% filter((id - 4) %% 10 == 0) %>% 
    mutate(id = (id - 4) / 10) %>% arrange(id)
  
  # Repeat wildtype to match the a b pairs
  d_wildtype <- d_wildtype[rep(seq_len(nrow(d_wildtype)), each = nrow(d_a)), ]
  d_wildtype$id <- d_ab$id

  # Save to output
  out$mutIDA <- mutsA$mutID
  out$mutIDB <- mutsB$mutID
  out$mutType_ab <- paste0(mutsA$mutType, "_", mutsB$mutType)
  out$wparab <- d_wildtype$fitness
  out$wparaB <- d_b$fitness
  out$wparAb <- d_a$fitness
  out$wparAB <- d_ab$fitness
  
  out <- out %>% filter(!is.na(mutIDA), !is.na(mutIDB))
  
  return(out)
}

# Calculates the site frequency spectra for mutations
CalcSFS <- function(dat) {
  dat$freqBin <- cut(dat$freq, breaks = 10)
  dat$optPerc <- dat$phenomean - 1
  dat$optPerc <- cut(dat$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))
  
  dat %>% 
    select(optPerc, seed, modelindex, 
           mutID, mutType, value, freqBin)
}
