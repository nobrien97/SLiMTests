model_names <- c("'NAR'", "'PAR'", "'FFLC1'", 
                 "'FFLI1'", "'FFBH'")

model_names_noquote <- c("NAR", "PAR", "FFLC1", 
                         "FFLI1", "FFBH")

model_names_labeller <- c("'NAR'" = "NAR", 
                          "'PAR'" = "PAR", 
                          "'FFLC1'" = "FFLC1", 
                          "'FFLI1'" = "FFLI1", 
                          "'FFBH'" = "FFBH")

molComp_names <- list("NAR" = c(
  # NAR and PAR
  "aZ",
  "bZ",
  "KZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),
"FFBH" = c(
  # FFBH
  "aX",
  "KZX",
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
)
)



molComp_labels <- list("NAR" = c(
  # NAR and PAR
  "aZ" = TeX("$\\alpha_Z$", output = "character"),
  "bZ" = TeX("$\\beta_Z$", output = "character"),
  "KZ" = TeX("$K_Z$", output = "character"),
  "KXZ" = TeX("$K_{XZ}$", output = "character"),
  "zZ" = TeX("$\\zeta_Z$", output = "character"), # baseline expression
  "h" = TeX("$h$", output = "character"), # hill coefficient
  "gX" = TeX("$\\gamma_X$", output = "character") # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY" = TeX("$\\alpha_Y$", output = "character"),
  "bY" = TeX("$\\beta_Y$", output = "character"),
  "KY" = TeX("$K_Y$", output = "character"),
  "aZ" = TeX("$\\alpha_Z$", output = "character"),
  "bZ" = TeX("$\\beta_Z$", output = "character"),
  "KXZ" = TeX("$K_{XZ}$", output = "character"),
  "zZ" = TeX("$\\zeta_Z$", output = "character"), # baseline expression
  "h" = TeX("$h$", output = "character"), # hill coefficient
  "gX" = TeX("$\\gamma_X$", output = "character") # X multiplier
),
"FFBH" = c(
  # FFBH
  "aX" = TeX("$\\alpha_X$", output = "character"),
  "KZX" = TeX("$K_{ZX}$", output = "character"),
  "aY" = TeX("$\\alpha_Y$", output = "character"),
  "bY" = TeX("$\\beta_Y$", output = "character"),
  "KY" = TeX("$K_Y$", output = "character"),
  "aZ" = TeX("$\\alpha_Z$", output = "character"),
  "bZ" = TeX("$\\beta_Z$", output = "character"),
  "KXZ" = TeX("$K_{XZ}$", output = "character"),
  "zZ" = TeX("$\\zeta_Z$", output = "character"), # baseline expression
  "h" = TeX("$h$", output = "character"), # hill coefficient
  "gX" = TeX("$\\gamma_X$", output = "character") # X multiplier
)
)

all_molcomp_features <- c(
  "aX" = TeX("$\\alpha_X$", output = "character"),
  "KZX" = TeX("$K_{ZX}$", output = "character"),
  "aY" = TeX("$\\alpha_Y$", output = "character"),
  "bY" = TeX("$\\beta_Y$", output = "character"),
  "KY" = TeX("$K_Y$", output = "character"),
  "aZ" = TeX("$\\alpha_Z$", output = "character"),
  "bZ" = TeX("$\\beta_Z$", output = "character"),
  "KZ" = TeX("$K_{Z}$", output = "character"),
  "KXZ" = TeX("$K_{XZ}$", output = "character"),
  "zZ" = TeX("$\\zeta_Z$", output = "character"), # baseline expression
  "h" = TeX("$h$", output = "character"), # hill coefficient
  "gX" = TeX("$\\gamma_X$", output = "character") # X multiplier
  
)

molComp_names[["PAR"]] <- molComp_names[["NAR"]]
molComp_names[["FFLI1"]] <- molComp_names[["FFLC1"]]

molComp_labels[["PAR"]] <- molComp_labels[["NAR"]]
molComp_labels[["FFLI1"]] <- molComp_labels[["FFLC1"]]

t_mc_combos <- expand.grid(1:11, 1:4)
mc_permodel <- list("NAR" = t_mc_combos[(t_mc_combos$Var1 < 8) & !(t_mc_combos$Var2 %in% c(3,4)),],
                    "PAR" = t_mc_combos[(t_mc_combos$Var1 < 8) & !(t_mc_combos$Var2 %in% c(3,4)),],
                    "FFLC1" = t_mc_combos[(t_mc_combos$Var1 < 10) & !(t_mc_combos$Var2 %in% c(4)),],
                    "FFLI1" = t_mc_combos[(t_mc_combos$Var1 < 10) & !(t_mc_combos$Var2 %in% c(4)),],
                    "FFBH" = t_mc_combos)

h2_colnames <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_aX", "VA_KZX", 
                 "VA_aY", "VA_bY", "VA_KY", "VA_aZ", "VA_bZ", "VA_KZ", "VA_KXZ", 
                 "VA_base", "VA_n", "VA_XMult", "CVA_aX_KZX", "CVA_aX_aY", 
                 "CVA_aX_bY", "CVA_aX_KY", "CVA_aX_aZ", "CVA_aX_bZ", "CVA_aX_KZ", 
                 "CVA_aX_KXZ", "CVA_aX_base", "CVA_aX_n", "CVA_aX_XMult", 
                 "CVA_KZX_aY", "CVA_KZX_bY", "CVA_KZX_KY", "CVA_KZX_aZ", 
                 "CVA_KZX_bZ", "CVA_KZX_KZ", "CVA_KZX_KXZ", "CVA_KZX_base", 
                 "CVA_KZX_n", "CVA_KZX_XMult", "CVA_aY_bY", "CVA_aY_KY", 
                 "CVA_aY_aZ", "CVA_aY_bZ", "CVA_aY_KZ", "CVA_aY_KXZ", 
                 "CVA_aY_base", "CVA_aY_n", "CVA_aY_XMult", "CVA_bY_KY", 
                 "CVA_bY_aZ", "CVA_bY_bZ", "CVA_bY_KZ", "CVA_bY_KXZ", 
                 "CVA_bY_base", "CVA_bY_n", "CVA_bY_XMult", "CVA_KY_aZ", 
                 "CVA_KY_bZ", "CVA_KY_KZ", "CVA_KY_KXZ", "CVA_KY_base", 
                 "CVA_KY_n", "CVA_KY_XMult", "CVA_aZ_bZ", "CVA_aZ_KZ", 
                 "CVA_aZ_KXZ", "CVA_aZ_base", "CVA_aZ_n", "CVA_aZ_XMult", 
                 "CVA_bZ_KZ", "CVA_bZ_KXZ", "CVA_bZ_base", "CVA_bZ_n", 
                 "CVA_bZ_XMult", "CVA_KZ_KXZ", "CVA_KZ_base", "CVA_KZ_n", 
                 "CVA_KZ_XMult", "CVA_KXZ_base", "CVA_KXZ_n", "CVA_KXZ_XMult", 
                 "CVA_base_n", "CVA_base_XMult", "CVA_n_XMult", "h2_aX", "h2_KZX", "h2_aY", "h2_bY", 
                 "h2_KY", "h2_aZ", "h2_bZ", "h2_KZ", "h2_KXZ", "h2_base", "h2_n", 
                 "h2_XMult")





mutate <- dplyr::mutate
select <- dplyr::select
summarise <- dplyr::summarise
rename <- dplyr::rename
filter <- dplyr::filter
desc <- dplyr::desc
arrange <- dplyr::arrange

pal <- paletteer_d("nationalparkcolors::Everglades", 5)

# Adds the parameter combination to a dataframe
AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]])
}


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

generateCol<-function(x,colCode,col,whichShadow=c(T,T,T)){
  lz <- lapply(1:ncol(x$ImpHistory),
               function(i) x$ImpHistory[is.finite(x$ImpHistory[,i]),i])
  colnames(x$ImpHistory) -> names(lz)
  
  #Selection of shadow meta-attributes
  numShadow <- sum(whichShadow)
  lz[c(rep(TRUE,length(x$finalDecision)),whichShadow)] -> lz
  
  
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
  
  # Rearrange by median
  ii<-order(sapply(lz,stats::median))
  
  return(col[ii])
}


RunRandomForestMolCompPerMotif <- function(dataset, seed = NULL, train.test = c(0.7, 0.3)) {
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }
  
  motifs <- levels(dataset$model)
  result <- vector(mode = "list")
  
  for (motif in motifs) {
    set.seed(seed)
    d_rf <- dataset %>%
      filter(model == motif) %>%
      select(isAdapted, dataset, paste0("meanMC", unique(mc_permodel[[motif]]$Var1))) %>%
      rename_at(dplyr::vars(starts_with("meanMC")), function(x) {
        i <- gsub("meanMC", "", x)
        return(molComp_names[[motif]][as.numeric(i)])
      }
      )
    
    adapted_counts <- table(d_rf$isAdapted)
    total_counts <- sum(adapted_counts)
    num_responses <- length(adapted_counts)
    adapted_weights <- total_counts / (num_responses * adapted_counts)
    names(adapted_weights) <- levels(d_rf$isAdapted)
    
    
    idx <- sample(2, nrow(d_rf), replace = T, prob = train.test)
    train <- d_rf[idx == 1,]
    test <- d_rf[idx == 2,]
    
    
    # no balancing
    rf_nobal <- randomForest(formula = isAdapted ~ .,
                             data = train,
                             ntree = 500,
                             proximity = T,
                             importance = T,
                             type = "classification")
    
    print(rf_nobal)
    
    # With balancing (class weights)
    rf_bal <- randomForest(formula = isAdapted ~ .,
                           data = train,
                           strata = train$isAdapted,
                           classwt = adapted_weights,
                           ntree = 500,
                           proximity = T,
                           importance = T,
                           type = "classification")
    
    print(rf_bal)
    
    # Training data
    p_train_bal <- predict(rf_bal, train)
    caret::confusionMatrix(p_train_bal, train$isAdapted)
    
    p_train_nobal <- predict(rf_nobal, train)
    caret::confusionMatrix(p_train_nobal, train$isAdapted)
    
    
    # Test data
    p_test_bal <- predict(rf_bal, test)
    p_test_bal_probs <- predict(rf_bal, test,
                                type = "prob")[,1]
    
    
    
    result[[motif]][["cMat_bal"]] <- caret::confusionMatrix(p_test_bal, test$isAdapted)
    caret::confusionMatrix(p_test_bal, test$isAdapted)
    
    p_test_nobal <- predict(rf_nobal, test)
    p_test_nobal_probs <- predict(rf_nobal, test,
                                  type = "prob")[,1]
    
    result[[motif]][["cMat_nobal"]] <- caret::confusionMatrix(p_test_nobal, test$isAdapted)
    caret::confusionMatrix(p_test_nobal, test$isAdapted)
    
    
    # roc
    d_roc_bal <- roc(response = test$isAdapted,
                     predictor = p_test_bal_probs)
    
    d_roc_nobal <- roc(response = test$isAdapted,
                       predictor = p_test_nobal_probs,
                       levels = rev(levels(test$isAdapted)))
    
    d_rocs <- data.frame(model = c(rep("Weighted", times = length(d_roc_bal$sensitivities)),
                                   rep("Unbalanced", times = length(d_roc_nobal$sensitivities))),
                         sens = c(d_roc_bal$sensitivities, d_roc_nobal$sensitivities),
                         spec = c(d_roc_bal$specificities, d_roc_nobal$specificities))
    
    
    roc_aucs <- c(pROC::auc(d_roc_nobal), pROC::auc(d_roc_bal))
    
    
    ggplot(d_rocs,
           aes(x = 1 - spec, y = sens, colour = model)) +
      geom_line() + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
      annotate("text", x = c(0.75, 0.75), y = c(0.375, 0.25), 
               label = paste("AUC:", round(roc_aucs, digits = 3)),
               colour = pal[1:2]) +
      theme_bw() +
      scale_colour_manual(values = pal) +
      ggtitle(motif) +
      labs(x = "1 - Specificity", y = "Sensitivity", colour = "RF Model") +
      theme(legend.position = "bottom",
            text = element_text(size = 12)) -> plt_roc
    ggsave(paste0("plt_RF_ROC_", motif, ".png"), plt_roc,
           device = png, width = 4, height = 4, bg = "white")
    
    result[[motif]][["d_roc"]] <- d_rocs
    result[[motif]][["pltROC"]] <- plt_roc
    
    # Importance measures
    ## Boruta, permutation importance, sobol MDA
    bor <- Boruta::Boruta(isAdapted ~ ., data = d_rf)
    bor
    plot(bor)
    
    d_bor <- process_the_Boruta_data(bor)
    
    result[[motif]][["bor"]] <- d_bor
    
    shadow_names <- c("shadowMin" = TeX("Shadow Variable (min)", output = "character"),
                      "shadowMean" = TeX("Shadow Variable (mean)", output = "character"),
                      "shadowMax" = TeX("Shadow Variable (max)", output = "character"))
    
    # Sort variables by Boruta median for all other importance plots
    bor_order <- names(sort(unlist(d_bor %>%
                                     summarise_all(median))))
    
    # Axis labels
    motif_labels <- c(molComp_labels[[motif]], 
                      "dataset" = TeX("Trait/selection alignment", output = "character"), 
                      shadow_names)[bor_order]
    
    pal_boruta <- generateCol(bor, colCode=c("#00A000","#EECC00","#DD0000","#00C0EA"),
                              col = NULL)
    
    ggplot(d_bor %>% pivot_longer(everything()) %>%
             mutate(x = fct_reorder(name, value, median)),
           aes(x = x, y = value, fill = x)) +
      geom_boxplot(show.legend = F, linewidth = 0.25) +
      theme_bw() +
      scale_fill_manual(values = pal_boruta) +
      scale_x_discrete(labels = parse(text = motif_labels),
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Boruta Importance") +
      theme(text = element_text(size = 12)) -> plt_boruta_imp
    plt_boruta_imp
    ggsave(paste0("plt_boruta_", motif, ".png"), plt_boruta_imp, 
           device = png, bg = "white",
           width = 12, height = 8)
    
    
    # Permutation
    predictor <- iml::Predictor$new(rf_bal, 
                                    data = test[, 2:(ncol(test))], 
                                    y = test$isAdapted,
                                    type = "prob")
    
    # Need to set the option future globals maxsize
    options(future.globals.maxSize = 3221225472)
    imp <- iml::FeatureImp$new(predictor,
                               loss = "ce",
                               n.repetitions = 100)
    
    result[[motif]][["FeatImp"]] <- imp
    
    
    motif_labels_noshadow <- motif_labels[names(motif_labels) %in% c(molComp_names[[motif]], "dataset")]
    
    ggplot(imp$results %>%
             mutate(feature = factor(feature, levels = names(motif_labels_noshadow)))
           ,
           aes(x = feature, y = importance)) +
      geom_point() +
      geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                    width = 0.2) +
      scale_x_discrete(labels = parse(text=motif_labels_noshadow),
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature/Component", y = "Permutation Importance") +
      theme_bw() +
      theme(text = element_text(size = 12)) -> plt_perm_imp
    plt_perm_imp
    ggsave(paste0("plt_perm_molcomp_", motif,".png"), 
           device = png, width = 9, height = 5, bg = "white")
    
    # Interaction strengths
    ia <- Interaction$new(predictor)
    
    result[[motif]][["pred"]] <- predictor
    result[[motif]][["ia"]] <- ia
    
    
    # Sobol MDA
    rf_sob <- sobolMDA::ranger(isAdapted ~ .,
                               data = train, num.trees = 500, 
                               importance = "sobolMDA")
    sob <- rf_sob$variable.importance
    d_sob <- data.frame(feature = names(sob),
                        sobelMDA = sob)
    
    d_sob$feature <- factor(d_sob$feature, levels = names(motif_labels_noshadow))
    
    ggplot(d_sob,
           aes(x = feature, y = sobelMDA)) +
      geom_point() +
      geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
                   linewidth = 0.5) +
      theme_bw() +
      scale_x_discrete(labels = parse(text = motif_labels_noshadow),
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Sobel MDA") +
      theme(text = element_text(size = 12)) -> plt_sob
    plt_sob
    result[[motif]][["d_sob"]] <- d_sob
    
    layout <- "
AAAA
AAAA
AAAA
BBCC
BBCC
"
    plt_featimp <- plt_boruta_imp +
      plt_perm_imp +
      plt_sob +
      plot_layout(design = layout) +
      plot_annotation(tag_levels = 'A',
                      title = paste("Feature importance for", motif, "motif")) &
      theme(plot.tag = element_text(face = "bold"))
    plt_featimp
    result[[motif]][["plt_featimp"]] <- plt_featimp
    
    ggsave(paste0("plt_featimp_", motif, ".png"), plt_featimp, 
           device = png, width = 12, height = 10, bg = "white")
    
    
    # Accumulated local effects
    ale <- FeatureEffects$new(predictor, grid.size = 10)
    ale$plot()
    
    result[[motif]][["ale"]] <- ale
    
    ale_plots <- vector(mode = "list", length = length(molComp_names[[motif]]) + 1)
    
    ale_labels <- motif_labels[names(ale$results)]
    
    for (i in seq_along(ale_plots)) {
      d_ale <- ale$results[[i]]
      d_ale <- d_ale %>% filter(.class == "Adapted")
      x_label <- ale_labels[d_ale$.feature[1]] 
      
      if (!is.numeric(d_ale$.borders[1])) {
        geom_fn <- geom_lollipop
        scale_fn <- scale_x_discrete
      } else {
        geom_fn <- geom_line
        scale_fn <- scale_x_continuous
        # Remove outliers
        d_ale <- d_ale %>% filter(.borders < 1000)
      }
      ale_plots[[i]] <- ggplot(d_ale,
                               aes(x = .borders, y = .value)) +
        geom_fn() +
        scale_fn() +
        theme_bw() +
        labs(x = parse(text = x_label), y = "ALE of adaptation probability")
    }
    
    alePlot <- plot_grid(plotlist = ale_plots,
                         labels= "AUTO")
    ggsave(paste0("plt_ale_molcomp_", motif, ".png"), alePlot, device = png, bg = "white",
           width = 12, height = 9)
    
    result[[motif]][["alePlot"]] <- alePlot
    
  }
  
  return(result)
}

