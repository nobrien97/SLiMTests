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
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]])
}


# Fitness effect calculation functions
## Gaussian fitness function, used for both models
calcAddFitness <- function(phenotypes, optimum, width) {
  dists <- (phenotypes - optimum)^2
  return(exp(-(dists * width)))
}

CalcPhenotypeEffects <- function(dat, dat_fixed, dat_opt) {
  # If there's no fixations, set the fixed value to 0
  # TODO: Need to set to be relative to the baseline reference value (e.g. base = 0.01)
  if (nrow(dat_fixed) == 0) {
    dat_fixed[1,] <- dat[1,]
    dat_fixed %>% mutate(value = 0)
  }
  
  return(CalcNetworkPhenotypeEffects(dat, dat_fixed, dat_opt))
}

## Run the ODELandscaper tool to evaluate phenotype and fitness
## for many individuals at once.
runLandscaper <- function(df_path, output, optimum, motif, threads, useID = FALSE) {
  command <- "~/Tools/odeLandscape/ODELandscaper -i %s -o ./%s -O %s -s %s -t %i"
  #command <- "ODELandscaper -i %s -o ./%s -O %s -s %s -t %i"
  if (useID) {
    command <- paste(command, "-I")
  }
  system(sprintf(command,
                 df_path, output, optimum, motif, threads))
  result <- read_csv(paste0("./", output), col_names = F, col_types = "d")

  result_names <- c("fitness", "trait1", "trait2")

  # Column names depend on the motif 
  switch(motif,
    "NAR"   = { result_names <- c(result_names, "aZ", "bZ", "KZ", "KXZ", "base", "n", "XMult") },
    "PAR"   = { result_names <- c(result_names, "aZ", "bZ", "KZ", "KXZ", "base", "n", "XMult") },
    "FFLC1" = { result_names <- c(result_names, "trait3", "aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult") },
    "FFLI1" = { result_names <- c(result_names, "trait3", "aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult") },
    "FFBH"  = { result_names <- c(result_names, "trait4", "aX", "KZX", "aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult") }
  )

  # Add row id to the names
  if (useID) {
    result_names <- c("id", result_names)
  }

  names(result) <- result_names

  return(result)
}

## Calculate fitness in network models
## TODO: Adjust for multiple motifs, remember that base() needs to be adjusted because the default
## value isn't 1, for the PAR it's 0.01, for others it's 0
CalcNetworkPhenotypeEffects <- function(dat, dat_fixed, dat_opt) {
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
  
  dat <- dat %>% inner_join(dat_fixed, 
                            by = c("gen", "seed", "modelindex"))
  
  dat$rowID <- as.integer(rownames(dat))

  # Write optimum file
  WriteOptimumInputTable(dat_opt, dat)

  # Get phenotypes without the mutation
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, starts_with("fixEffectSum")), 
              paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  d_popfx <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 
                           paste0("d_grid_opt", model, ".csv"), model, 4, TRUE)
  
  # Segregating mutation calculations
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = d_popfx$fixEffectSum + 2 * value
  # Aa = d_popfx$fixEffectSum + value
  # aa = d_popfx$fixEffectSum
  
  # Add on the segregating effect to the fixation effects
  # TODO fix this: programmatically name variables 
  d_dat_withFX <- dat
  nMutTypes <- GetNMutTypes(dat$model[1])

  for (i in seq_len(nMutTypes)) {
    mutType <- as.character(i + 2) # offset by first 2 mut types (not included)
    d_dat_withFX[,paste0("mutValue_", mutType)] <- 
      exp(log(d_dat_withFX[,paste0("fixEffectSum_", mutType)]) + d_dat_withFX[,paste0("value_", mutType)])
    d_dat_withFX[,paste0("mutValueAA_", mutType)] <- 
      exp(log(d_dat_withFX[,paste0("fixEffectSum_", mutType)]) + 2 * d_dat_withFX[,paste0("value_", mutType)])
  }
    
  # Get phenotypes with the mutation
  data.table::fwrite(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, starts_with("mutValue_")), 
              paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 
                      paste0("d_grid_opt", model, ".csv"), model, 4, TRUE)

  data.table::fwrite(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, starts_with("mutValueAA_")), 
              paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  AA <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 
                      paste0("d_grid_opt", model, ".csv"), model, 4, TRUE)
  
  # Ensure that the tables are aligned by id before we join them
  dat <- dat %>% arrange(rowID)
  Aa <- Aa %>% arrange(id)
  AA <- AA %>% arrange(id)
  
  # TODO Get effect
  # Total effect is the sum of differences between the traits
  # But traits are on different scales, so it doesn't make much sense to measure that way
  # Focus on fitness changes and effect size in terms of fitness?
  dat$AA_trait1 <- AA[,"trait1"]
  dat$AA_trait2 <- AA[,"trait2"]
  dat$AA_trait3 <- AA[,"trait3"]
  dat$AA_trait4 <- AA[,"trait4"]
  dat$Aa_trait1 <- Aa[,"trait1"]
  dat$Aa_trait2 <- Aa[,"trait2"]
  dat$Aa_trait3 <- Aa[,"trait3"]
  dat$Aa_trait4 <- Aa[,"trait4"]
  dat$aa_trait1 <- d_popfx[,"trait1"]
  dat$aa_trait2 <- d_popfx[,"trait2"]
  dat$aa_trait3 <- d_popfx[,"trait3"]
  dat$aa_trait4 <- d_popfx[,"trait4"]
  dat$avFit <- Aa$fitness - d_popfx$fitness
  dat$avFit_AA <- AA$fitness - d_popfx$fitness
  dat$wAA <- AA$fitness
  dat$wAa <- Aa$fitness
  dat$waa <- d_popfx$fitness
  dat$s <- dat$wAa - dat$waa
  return(dat)
}

GetModelInvalidTraitIndices <- function(model) {
  switch(model,
    "NAR"   = { return( c(3,4) ) },
    "PAR"   = { return( c(3,4) ) },
    "FFLC1" = { return( c(4) ) },
    "FFLI1" = { return( c(4) ) },
    "FFBH"  = { return( numeric(0) ) }
  )
}

GetNMutTypes <- function(model) {
  switch(model,
    "NAR"   = { return( 7 ) },
    "PAR"   = { return( 7 ) },
    "FFLC1" = { return( 9 ) },
    "FFLI1" = { return( 9 ) },
    "FFBH"  = { return( 11 ) }
  )
}

WriteOptimumInputTable <- function(dat_opt, dat) {
  # Extract the right optima and widths
  dat_opt <- dat_opt %>% filter(modelindex == dat$modelindex[1]) 
  
  exclude_indices <- GetModelInvalidTraitIndices(dat$model[1])
  exclude_indices <- paste0("trait", exclude_indices)

  # Filter out the unneeded traits for non-FFBH models
  dat_opt <- dat_opt %>% filter(!starts_with(exclude_indices))

  # Attach opt and sigma to the dataframe, match by seed and modelindex
  dat <- dat %>% inner_join(dat_opt, by = c("seed", "modelindex"))

  # Now create a table for the optima and sigmas
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, ends_with(c("opt", "sig"))), 
              paste0("d_grid_opt", model, ".csv"), sep = ",", col.names = F, row.names = F)
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
PairwiseEpistasisNetwork <- function(dat_fixed, muts, dat_opt, n = 1000, m = 10, 
                                 returnAverage = F, weightABByFreq = F) {
  # Get fixed effects/wildtype
  dat_fixed <- as.data.table(dat_fixed)

  model <- as.character(dat_fixed$modelindex)[1]
  
  fixEffectDat <- dat_fixed %>%
    group_by(gen, seed, modelindex, model, mutType) %>%
    summarise(fixEffectSum = exp(2 * sum(value))) %>%
    select(gen, seed, modelindex, mutType, fixEffectSum) %>%
    ungroup()
  
  nMutTypes <- GetNMutTypes(dat_fixed$model[1])

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
                wwt = numeric(output_len),
                ew = numeric(output_len))
  
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
      select(gen, seed, modelindex, model, mutType_a, mutType_b, a, b, 
             starts_with("fixEffectSum")) %>%
      mutate(mutType_ab = paste(mutType_a, mutType_b, sep = "_"))
  
    # Column names: mutTypes for a and b mutations
    abNames <- paste(c(rep("a", times = nMutTypes), rep("b", times = nMutTypes)), (1:nMutTypes) + 2, sep = "_")
    result[,abNames] <- 0
        
    # initialize a and b values for the right molecular component
    for (i in (1:nMutTypes) + 2) {
      result[result$mutType_a == paste0(i), paste0("a_", i)] <- result[result$mutType_a == paste0(i), "a"]
      result[result$mutType_b == paste0(i), paste0("b_", i)] <- result[result$mutType_b == paste0(i), "b"]
    }
    
    # Add on a, b, ab to the base effect
    for (i in (1:nMutTypes) + 2) {
      result[, paste0("a_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("a_", i)])
      result[, paste0("b_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("b_", i)])
      result[, paste0("ab_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("a_", i)] + result[,paste0("b_", i)])
    }
    result$rowID <- as.integer(rownames(result))
    
    # Get optimum and write to table
    WriteOptimumInputTable(dat_opt, result)

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
             rowID, starts_with("a_molComp")) %>%
      ungroup() %>% select(!(gen:modelindex)) %>%
      mutate(rowID = as.numeric(paste0(rowID, 2)))
  
    # b is organised differently to a, need to calculate first m for each mutgroup
    # and repeat that for the remaining
    d_b <- result %>%
      group_by(gen, seed, modelindex) %>%
      select(gen, seed, modelindex, 
             rowID, starts_with("b_molComp")) %>%
      ungroup() %>% select(!(gen:modelindex)) %>%
      mutate(rowID = as.numeric(paste0(rowID, 3)))
    
    
    d_ab <- result %>%
      group_by(gen, seed, modelindex) %>%
      select(gen, seed, modelindex, 
             rowID, starts_with("ab_molComp")) %>%
      ungroup() %>% select(!(gen:modelindex)) %>%
      mutate(rowID = as.numeric(paste0(rowID, 4)))
  
    # Remove column names so we can join them
    colnames(d_wildtype) <- paste0("v", 1:ncol(d_wildtype))
    colnames(d_a) <- paste0("v", 1:ncol(d_a))
    colnames(d_b) <- paste0("v", 1:ncol(d_b))
    colnames(d_ab) <- paste0("v", 1:ncol(d_ab))
  
    d_landscaper <- rbind(d_wildtype, d_a, d_b, d_ab)
    
    # Run landscaper
    data.table::fwrite(d_landscaper, 
                paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
    d_phenos <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 
                  paste0("d_grid_opt", model, ".csv"), model, 4, TRUE)
    
  
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
    result$wwt <- d_wildtype$fitness
    
    result$ew <- log(result$wab) - (log(result$wa) + log(result$wb))
    
    # Fill output: figure out which range of the output vector to fill
    thisIterRange <- ( (j-1) * (nrow(fixEffectDat) * m) + 1 ):( j * (nrow(fixEffectDat) * m) )
    
    out[thisIterRange,] <- result %>%
      ungroup() %>%
      select(gen, seed, modelindex, mutType_ab, wa, wb, wab, wwt, ew)

    pb$tick(1)
    j <- j + 1
  }
  
  # if we're averaging over iterations, do that now:
  # identical to grand mean since all replicates have the same sample size
  # Combined SEM via pythagorean theorem 
  if (returnAverage) {
    out <- out %>%
      group_by(gen, seed, modelindex, mutType_ab) %>%
      summarise(meanEW = mean(ew),
                meanwa = mean(wa),
                meanwb = mean(wb),
                meanwab = mean(wab),
                meanwwt = mean(wwt),
                sdEW = sd(ew),
                sdwa = sd(wa),
                sdwb = sd(wb),
                sdwab = sd(wab),
                sdwwt = sd(wwt))
  }

  return(out)
}


PairwiseFitnessRankNetwork <- function(dat_fixed, muts, A_ids, B_ids) {
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
  
  nMutTypes <- GetNMutTypes(muts$model[1])   
  
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
      

  result <- mutsA
  result <- result %>% rename(a = value, mutType_a = mutType)
  result$mutType_b <- mutsB$mutType
  result$b <- mutsB$value
  result <- inner_join(result, dat, by = c("gen", "seed", "modelindex")) %>%
    select(gen, seed, modelindex, model, mutType_a, mutType_b, a, b, 
           starts_with("fixEffectSum")) %>%
    mutate(mutType_ab = paste(mutType_a, mutType_b, sep = "_"))
  
    abNames <- paste(c(rep("a", times = nMutTypes), rep("b", times = nMutTypes)), (1:nMutTypes) + 2, sep = "_")
    result[,abNames] <- 0
  
  # initialize a and b values for the right molecular component
  for (i in (1:nMutTypes) + 2) {
    result[result$mutType_a == paste0(i), paste0("a_", i)] <- result[result$mutType_a == paste0(i), "a"]
    result[result$mutType_b == paste0(i), paste0("b_", i)] <- result[result$mutType_b == paste0(i), "b"]
  }
  
  # Add on a, b, ab to the base effect
  for (i in (1:nMutTypes) + 2) {
    result[, paste0("a_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("a_", i)])
    result[, paste0("b_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("b_", i)])
    result[, paste0("ab_molComp_", i)] <- exp(log(result[,paste0("fixEffectSum_", i)]) + result[,paste0("a_", i)] + result[,paste0("b_", i)])
  }
  result$rowID <- as.integer(rownames(result))

  # Get optimum and write to table
  WriteOptimumInputTable(dat_opt, result)

  
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
           rowID, starts_with("a_molComp")) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 2)))
  
  # b is organised differently to a, need to calculate first m for each mutgroup
  # and repeat that for the remaining
  d_b <- result %>%
    group_by(gen, seed, modelindex) %>%
    select(gen, seed, modelindex, 
           rowID, starts_with("b_molComp")) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 3)))
  
  
  d_ab <- result %>%
    group_by(gen, seed, modelindex) %>%
    select(gen, seed, modelindex, 
           rowID, starts_with("ab_molComp")) %>%
    ungroup() %>% select(!(gen:modelindex)) %>%
    mutate(rowID = as.numeric(paste0(rowID, 4)))
  
  # Remove column names so we can join them
  colnames(d_wildtype) <- paste0("v", 1:ncol(d_wildtype))
  colnames(d_a) <- paste0("v", 1:ncol(d_a))
  colnames(d_b) <- paste0("v", 1:ncol(d_b))
  colnames(d_ab) <- paste0("v", 1:ncol(d_ab))
  
  d_landscaper <- rbind(d_wildtype, d_a, d_b, d_ab)
  

  # Run landscaper
  data.table::fwrite(d_landscaper, 
                paste0("d_grid", model, ".csv"), sep = ",", col.names = F, row.names = F)
  d_phenos <- runLandscaper(paste0("d_grid", model, ".csv"), paste0("data_popfx", model, ".csv"), 
                  paste0("d_grid_opt", model, ".csv"), model, 4, TRUE)
  
  
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
  
  dat %>% 
    select(timePoint, seed, modelindex, 
           mutID, mutType, value, freqBin)
}
