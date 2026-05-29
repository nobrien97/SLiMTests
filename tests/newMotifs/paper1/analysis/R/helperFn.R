model_names <- c("'NAR'", "'PAR'", "'FFLC1'", 
                 "'FFLI1'", "'FFBH'")

model_names_noquote <- c("NAR", "PAR", "FFLC1", 
                         "FFLI1", "FFBH")

model_names_labeller <- c("'NAR'" = "NAR", 
                          "'PAR'" = "PAR", 
                          "'FFLC1'" = "FFLC1", 
                          "'FFLI1'" = "FFLI1", 
                          "'FFBH'" = "FFBH")


mutate <- dplyr::mutate
select <- dplyr::select
summarise <- dplyr::summarise
rename <- dplyr::rename

pal <- paletteer_d("nationalparkcolors::Everglades", 5)

# Adds the parameter combination to a dataframe
AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]])
}

select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter

# Cowplot 1.1.3 bug: won't get legend, this fixes
get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

# Convert a row to a matrix
row_to_m <- function(x) {
  # Get the number of traits and covariance terms
  n <- 2
  cov_terms <- 12
  
  if (x$model == "FFLC1" | x$model == "FFLI1") {
    n <- 3
    cov_terms <- c(12, 13, 23)
  }
  
  if (x$model == "FFBH") {
    n <- 4
    cov_terms <- c(12, 13, 14, 23, 24, 34)
  }
  
  # Triangular number for number of covariance terms
  n_cov <- ((n-1) * n) / 2
  
  m <- matrix(NA_real_, nrow = n, ncol = n)
  
  # Variances
  diag(m) <- unlist(x[1,paste0("var_", 1:n)])
  
  # Covariances
  m[lower.tri(m)] <- unlist(x[1,paste0("cov_", cov_terms)])
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  
  return(m)
}

Vrel <- function(l) {
  p <- length(l)
  avg_l <- mean(l)
  
  sum((l - avg_l)^2) / (p * (p-1) * avg_l^2)
}

# Setup correlation matrices
make_matrix <- function(x) {
  # triangular number to get length
  t <- nrow(x) 
  n <- ((-1 + sqrt(1 + 8 * t)) / 2) + 1
  
  cor_mat <- array(1, dim = c(n, n, 3))
  
  # lower CI
  cor_mat[,,1][lower.tri(cor_mat[,,1])] <- x$r_post_ci_lower
  cor_mat[,,1][upper.tri(cor_mat[,,1])] <- t(cor_mat[,,1])[upper.tri(cor_mat[,,1])]
  # mean
  cor_mat[,,2][lower.tri(cor_mat[,,2])] <- x$r_post_mean
  cor_mat[,,2][upper.tri(cor_mat[,,2])] <- t(cor_mat[,,2])[upper.tri(cor_mat[,,2])]
  # upper CI
  cor_mat[,,3][lower.tri(cor_mat[,,3])] <- x$r_post_ci_upper
  cor_mat[,,3][upper.tri(cor_mat[,,3])] <- t(cor_mat[,,3])[upper.tri(cor_mat[,,3])]
  
  return(cor_mat)
}



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




######################
# G matrix functions #
######################
# transform IDs from matrices to a list form
GetMatrixIDs <- function(matList) {
  lapply(matList, function(x) {
    data.frame(timePoint = x$timePoint, 
               seed = x$seed, 
               modelindex = x$modelindex, 
               isAdapted = x$isAdapted)}) -> matList
  
  
  lapply(matList, function(x) {
    split(x, seq(nrow(x)))
  }) -> matList
  # unlist to full form
  matList <- unlist(matList, recursive = F)
  return(matList)
}

# Above but with the Dataset column
GetMatrixIDsWithDataset <- function(matList) {
  lapply(matList, function(x) {
    data.frame(timePoint = x$timePoint, 
               seed = x$seed, 
               modelindex = x$modelindex, 
               isAdapted = x$isAdapted,
               dataset = x$dataset)}) -> matList
  
  
  lapply(matList, function(x) {
    split(x, seq(nrow(x)))
  }) -> matList
  # unlist to full form
  matList <- unlist(matList, recursive = F)
  return(matList)
}


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
  
    # cev and aut can be NaN if Ix(pca$values^2) gives an infinity (really small/0 eigenvalue)
    # in that case, set cev to 0
    if (is.nan(PCAdata$cev[i]))
      PCAdata$cev[i] <- 0.0
      PCAdata$aut[i] <- 0.0
    }
  
  return(PCAdata)
}


GetMotifTraitRange <- function(model) {
  result <- 1:4
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

# similarity between first eigenvectors of two matrices (lists of matrices)
GetCosineSimilarityTwoMats <- function(mat1, mat2, id) {
  #mat1 will be transformed in the original function, need to find
  # e_max for mat2
  
  eig <- lapply(mat2, function(x) return(eigen(x)$vectors[,1]))
  eig <- lapply(eig, function(x) {
    result <- numeric(4)
    result[1:length(x)] <- x
    return(result)
  })
  
  d_eig <- as.data.frame(t(as.data.frame(eig)))

  return(GetCosineSimilarity(mat1, d_eig, id))
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
    

    # If matrix isn't symmetric, force it to be
    if (!is.symmetric.matrix(g)) {
      g[lower.tri(g)] <- t(g)[lower.tri(g)] 
    }
    
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

# https://stackoverflow.com/a/28459434
equal_breaks <- function(n = 3, s = 0.05){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}

# https://stackoverflow.com/a/73416408
process_the_Boruta_data <- function(x, whichShadow=c(TRUE,TRUE,TRUE),
                                    colCode=c('green','yellow','red','blue'),
                                    col=NULL) {
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')
  
  #Removal of -Infs and conversion to a list
  lz <- lapply(1:ncol(x$ImpHistory),
               function(i) x$ImpHistory[is.finite(x$ImpHistory[,i]),i])
  colnames(x$ImpHistory) -> names(lz)
  
  #Selection of shadow meta-attributes
  numShadow <- sum(whichShadow)
  lz[c(rep(TRUE,length(x$finalDecision)),whichShadow)] -> lz
  
  generateCol<-function(x,colCode,col,numShadow){
    #Checking arguments
    if(is.null(col) & length(colCode)!=4)
      stop('colCode should have 4 elements.')
    #Generating col
    if(is.null(col)){
      rep(colCode[4],length(x$finalDecision)+numShadow)->cc
      cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1]
      cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2]
      cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3]
      col=cc
    }
    return(col)
  }
  
  #Generating color vector
  col <- generateCol(x, colCode, col, numShadow)
  
  #Ordering boxes due to attribute median importance
  ii<-order(sapply(lz,stats::median))
  lz[ii] -> lz
  col <- col[ii]
  lz_df <- do.call(rbind.data.frame, lz)
  df <- as.data.frame(t(lz_df))
  names(df) <- names(lz)
  rownames(df) <- NULL
  return(df)
}
