model_names <- c("'NAR'", "'PAR'", "'FFLC1'", 
                 "'FFLI1'", "'FFBH'")



ModelFromIndex <- function(id) {
  motifs <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
  return(motifs[id])
}

ModelFromIndexWithR <- function(id) {
  motifs <- rep(c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), times = 3)
  return(motifs[id])
}

RFromIndex <- function(id) {
  r <- rep(c(1e-10, 1e-5, 1e-1), each = 5)
  return(r[id])
}


se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}


CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi / 180)}

# Calculate evolvability metrics for trait data (Hansen and Houle 2008)
CalcECRATrait <- function(matList, id) {
  require(matrixcalc)
  require(Matrix)
  
  PCAdata <- data.frame(
    ev = numeric(length(matList)),
    cev = numeric(length(matList)),
    res = numeric(length(matList)),
    aut = numeric(length(matList))
  )
  
  PCAdata <- PCAdata %>%
    mutate(timePoint = id$timePoint,
           seed = id$seed,
           modelindex = id$modelindex,
           clus = id$clus,
           isAdapted = id$isAdapted)
  
  Hx <- function(x) {
    1/mean(1/x)
  }
  
  Ix <- function(x) {
    var(x)/mean(x)^2
  }
  
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    idx <- GetMotifTraitRange(as.character(id$model[i]))
    
    # Resize g to the proper dimensions
    g <- g[idx, idx]
    
    # If the matrix isn't positive semi-definite, find the nearest PD
    if (!is.positive.semi.definite(g)) {
      g <- as.matrix(nearPD(g)$mat)
    }
    
    pca <- eigen(g)
    k <- length(pca$values)
    
    # Calculate bTGb evolvability as well
    
    
    PCAdata$ev[i] <- mean(pca$values) #e
    PCAdata$cev[i] <- Hx(pca$values) * (1 + (2*Ix(1/pca$values)) / (k+2) ) #c
    PCAdata$res[i] <- sqrt(mean(pca$values^2)) * (1 - (Ix(pca$values^2) / (4*k+2) ) ) #r
    PCAdata$aut[i] <- (Hx(pca$values) / mean(pca$values)) * (1 + 2 * (Ix(pca$values) + Ix(1/pca$values) - 1 + Hx(pca$values)/mean(pca$values) + 2 * Ix(pca$values) * Ix(1/pca$values)/(k+2))/(k+2)) #a
  }
  
  return(PCAdata)
}


GetMotifTraitRange <- function(model) {
  switch (model,
          "'NAR'"   = { result <- 1:2 },
          "'PAR'"   = { result <- 1:2 },
          "'FFLC1'" = { result <- 1:3 },
          "'FFLI1'" = { result <- 1:3 },
          "'FFBH'"  = { result <- 1:4 },
          "NAR"   = { result <- 1:2 },
          "PAR"   = { result <- 1:2 },
          "FFLC1" = { result <- 1:3 },
          "FFLI1" = { result <- 1:3 },
          "FFBH"  = { result <- 1:4 }
          
  )
  return(result)
} 

# M is list of matrices, b is df of selection vectors
GetCosineSimilarity <- function(matList, bFrame, id) {
  require(matrixcalc)
  require(Matrix)
  
  bFrame <- as.data.frame(bFrame)
  
  result <- data.frame(
    cosSim = numeric(length(matList)),
    bTMb = numeric(length(matList))
  )
  
  result <- result %>%
    mutate(timePoint = id$timePoint,
           seed = id$seed,
           modelindex = id$modelindex,
           clus = id$clus,
           isAdapted = id$isAdapted)
  
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    idx <- GetMotifTraitRange(as.character(id$model[i]))
    
    # Resize g to the proper dimensions
    g <- g[idx, idx]
    b <- unlist(bFrame[i, idx])
    
    # If the matrix isn't positive semi-definite, find the nearest PD
    if (!is.positive.semi.definite(g)) {
      g <- as.matrix(nearPD(g)$mat)
    }
    
    pca <- eigen(g)
    
    # Similarity between leading eigenvector of the matrix and the selection gradient
    result$cosSim[i] <- sum(b * pca$vectors[,1]) / sqrt(sum(b^2) * sum(pca$vectors[,1]^2))
    result$bTMb[i] <- (t(b) %*% g %*% b)
  }
  
  return(result)
}

