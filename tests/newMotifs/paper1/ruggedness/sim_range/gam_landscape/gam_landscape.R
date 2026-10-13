# Uses GAM to plot the fitness landscape
library(mgcv)

gammamatrix <- function(fitness, traits) {
  # gamma matrix function
  
  # ARGUMENTS: fitness = single-column data frame containing a measure of
  #                 fitness
  #            traits = matrix of traits
  
  # RETURNS:   everything = a list, 
  #                 $results = summary from multiple regression
  #                 $gamma_matrix = gamma matrix
  
  num_traits <- length(traits)
  trait_names <- names(traits)
  dset <- as.data.frame(cbind(fitness, traits))
  
  for (i in 1:num_traits) {
    for (j in i:num_traits) {
      newterm <- traits[ ,i] * traits[ ,j]
      temp_names <- names(traits)
      traits <- cbind(traits, newterm)
      names(traits) <- c(temp_names, paste(names(traits)[i],
                                           names(traits)[j], sep = "_"))
    }
  }
  
  main_effects_L <- paste(names(traits), collapse = "+")
  formL <- as.formula(paste(names(fitness), "~", main_effects_L, sep = ""))
  dset <- as.data.frame(cbind(fitness, traits))
  G <- lm(formL, dset)
  results <- summary(G)
  
  gamma <- matrix(rep(0.0, num_traits * num_traits), nrow = num_traits)
  index <- num_traits + 2
  for (i in 1:num_traits) {
    for (j in i:num_traits) {
      # multiply quadratic terms by 2
      if (i == j) gamma[i, j] <- 2.0 * G$coefficients[index]
      else gamma[i, j] <- G$coefficients[index]
      gamma[j, i] <- gamma[i, j]
      index <- index + 1
    }
  }
  rownames(gamma) <- trait_names
  colnames(gamma) <- trait_names
  everything <- list(results = results, gamma_matrix = gamma)
  return(everything)
}


test <- data.frame(fitness = runif(1000),
                   a = runif(1000),
                   b = runif(1000),
                   c = runif(1000),
                   d = runif(1000))

gam_out <- gam(fitness ~ s(a) + s(b) + s(c) + s(d), data = test)
summary(gam_out)
test_fitted <- fitted(gam_out)

gamma_mat <- gammamatrix(test[,1], as.matrix(test[, -1]))

plot(test_fitted, test$fitness)
