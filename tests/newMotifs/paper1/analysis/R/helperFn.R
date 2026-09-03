library(tidyverse)
library(latex2exp)
library(paletteer)
library(ggh4x)
library(Rcpp)
library(ggbeeswarm)
library(emmeans)
library(cowplot)


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
                 "VA_zZ", "VA_h", "VA_gX", "CVA_aX_KZX", "CVA_aX_aY", 
                 "CVA_aX_bY", "CVA_aX_KY", "CVA_aX_aZ", "CVA_aX_bZ", "CVA_aX_KZ", 
                 "CVA_aX_KXZ", "CVA_aX_zZ", "CVA_aX_h", "CVA_aX_gX", 
                 "CVA_KZX_aY", "CVA_KZX_bY", "CVA_KZX_KY", "CVA_KZX_aZ", 
                 "CVA_KZX_bZ", "CVA_KZX_KZ", "CVA_KZX_KXZ", "CVA_KZX_zZ", 
                 "CVA_KZX_h", "CVA_KZX_gX", "CVA_aY_bY", "CVA_aY_KY", 
                 "CVA_aY_aZ", "CVA_aY_bZ", "CVA_aY_KZ", "CVA_aY_KXZ", 
                 "CVA_aY_zZ", "CVA_aY_h", "CVA_aY_gX", "CVA_bY_KY", 
                 "CVA_bY_aZ", "CVA_bY_bZ", "CVA_bY_KZ", "CVA_bY_KXZ", 
                 "CVA_bY_zZ", "CVA_bY_h", "CVA_bY_gX", "CVA_KY_aZ", 
                 "CVA_KY_bZ", "CVA_KY_KZ", "CVA_KY_KXZ", "CVA_KY_zZ", 
                 "CVA_KY_h", "CVA_KY_gX", "CVA_aZ_bZ", "CVA_aZ_KZ", 
                 "CVA_aZ_KXZ", "CVA_aZ_zZ", "CVA_aZ_h", "CVA_aZ_gX", 
                 "CVA_bZ_KZ", "CVA_bZ_KXZ", "CVA_bZ_zZ", "CVA_bZ_h", 
                 "CVA_bZ_gX", "CVA_KZ_KXZ", "CVA_KZ_zZ", "CVA_KZ_h", 
                 "CVA_KZ_gX", "CVA_KXZ_zZ", "CVA_KXZ_h", "CVA_KXZ_gX", 
                 "CVA_zZ_h", "CVA_zZ_gX", "CVA_h_gX", "h2_aX", "h2_KZX", "h2_aY", "h2_bY", 
                 "h2_KY", "h2_aZ", "h2_bZ", "h2_KZ", "h2_KXZ", "h2_zZ", "h2_h", 
                 "h2_gX")





mutate <- dplyr::mutate
select <- dplyr::select
summarise <- dplyr::summarise
rename <- dplyr::rename
filter <- dplyr::filter
desc <- dplyr::desc
arrange <- dplyr::arrange

pal <- paletteer_d("nationalparkcolors::Everglades", 5)
pal_tol <- c("#D2D1E2", "#C1A7C4", "#AD7BA4",
             "#93C7EC", "#6DA1D6", "#3E79BD",
             "#B8D29E", "#7BB47D", "#2C925B",
             "#F4A937", "#EB7121", "#D9070D")


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


ReadH2Data <- function(path, trait_type = "trait", calc_mode = "mkr") {
  require(readr)
  x <- read_csv(path)
  
  if (trait_type == "trait") {
    colnames(x) <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_t1",
                                       "VA_t2", "VA_t3", "VA_t4", "CVA_t1_t2", "CVA_t1_t3",
                                       "CVA_t1_t4", "CVA_t2_t3", "CVA_t2_t4", "CVA_t3_t4",
                                       "h2_t1", "h2_t2", "h2_t3", "h2_t4")
  } else if (trait_type == "mc") {
    colnames(x) <- h2_colnames
  }
  
  x$calcMode <- calc_mode
  
}



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

GetMeanMatrixIDs <- function(matList) {
  lapply(matList, function(x) {
    data.frame(model = x$model,
               dataset = x$dataset,
               timePoint = x$timePoint,
               isAdapted = x$isAdapted)}) -> matList
  
  
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
           dataset = id$dataset,
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

RunRandomForestPerMotif <- function(dataset, seed = NULL, train.test = c(0.7, 0.3), type = "classification") {
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }
  
  motifs <- levels(dataset$model)
  result <- vector(mode = "list")
  
  for (motif in motifs) {
    set.seed(seed)
    d_rf <- dataset %>%
      filter(model == motif) %>%
      select(-model)
    
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
                             type = type)
    
    print(rf_nobal)
    
    # With balancing (class weights)
    rf_bal <- randomForest(formula = isAdapted ~ .,
                           data = train,
                           strata = train$isAdapted,
                           classwt = adapted_weights,
                           ntree = 500,
                           proximity = T,
                           importance = T,
                           type = type)
    
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
    ggsave(paste0("plt_RF_ROC_permod_", motif, ".png"), plt_roc,
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
    # motif_labels <- c(molComp_labels[[motif]], 
    #                   "dataset" = TeX("Trait/selection alignment", output = "character"), 
    #                   shadow_names)[bor_order]
    
    pal_boruta <- generateCol(bor, colCode=c("#00A000","#EECC00","#DD0000","#00C0EA"),
                              col = NULL)
    
    ggplot(d_bor %>% pivot_longer(everything()) %>%
             mutate(x = fct_reorder(name, value, median)),
           aes(x = x, y = value, fill = x)) +
      geom_boxplot(show.legend = F, linewidth = 0.25) +
      theme_bw() +
      scale_fill_manual(values = pal_boruta) +
      scale_x_discrete(#labels = parse(text = motif_labels),
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Boruta Importance") +
      theme(text = element_text(size = 12)) -> plt_boruta_imp
    plt_boruta_imp
    ggsave(paste0("plt_boruta_permod_", motif, ".png"), plt_boruta_imp, 
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
    
    
    #motif_labels_noshadow <- motif_labels[names(motif_labels) %in% c(molComp_names[[motif]], "dataset")]
    
    ggplot(imp$results #%>%
             #mutate(feature = factor(feature, levels = names(motif_labels_noshadow)))
           ,
           aes(x = feature, y = importance)) +
      geom_point() +
      geom_errorbar(aes(ymin = importance.05, ymax = importance.95),
                    width = 0.2) +
      scale_x_discrete(#labels = parse(text=motif_labels_noshadow),
                       guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Permutation Importance") +
      theme_bw() +
      theme(text = element_text(size = 12)) -> plt_perm_imp
    plt_perm_imp
    ggsave(paste0("plt_perm_permod_", motif,".png"), 
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
    
    #d_sob$feature <- factor(d_sob$feature, levels = names(motif_labels_noshadow))
    
    ggplot(d_sob,
           aes(x = feature, y = sobelMDA)) +
      geom_point() +
      geom_segment(aes(xend = feature, y = 0, yend = sobelMDA),
                   linewidth = 0.5) +
      theme_bw() +
      scale_x_discrete(#labels = parse(text = motif_labels_noshadow),
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
    
    ggsave(paste0("plt_featimp_permod_", motif, ".png"), plt_featimp, 
           device = png, width = 12, height = 10, bg = "white")
    
    
    # Accumulated local effects
    ale <- FeatureEffects$new(predictor, grid.size = 10)
    ale$plot()
    
    result[[motif]][["ale"]] <- ale
    
    ale_plots <- vector(mode = "list", length = length(ale$features))
    
    #ale_labels <- motif_labels[names(ale$results)]
    
    for (i in seq_along(ale_plots)) {
      d_ale <- ale$results[[i]]
      d_ale <- d_ale %>% filter(.class == "Adapted")
      x_label <- ale$features[i] 
      
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
        labs(x = x_label, #parse(text = x_label), 
             y = "ALE of adaptation probability")
    }
    
    alePlot <- plot_grid(plotlist = ale_plots,
                         labels= "AUTO")
    ggsave(paste0("plt_ale_permod_", motif, ".png"), alePlot, device = png, bg = "white",
           width = 12, height = 9)
    
    result[[motif]][["alePlot"]] <- alePlot
    
  }
  
  return(result)
}

RunRandomForestPerMotifTimeToAdapt <- function(dataset, seed = NULL, train.test = c(0.7, 0.3), type = "regression") {
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }
  
  motifs <- levels(dataset$model)
  result <- vector(mode = "list")
  
  for (motif in motifs) {
    set.seed(seed)
    d_rf <- dataset %>%
      filter(model == motif) %>%
      select(-model)
    
    idx <- sample(2, nrow(d_rf), replace = T, prob = train.test)
    train <- d_rf[idx == 1,]
    test <- d_rf[idx == 2,]
    
    
    # no balancing
    rf_nobal <- randomForest(formula = timeToAdapt ~ .,
                             data = train,
                             ntree = 500,
                             proximity = T,
                             importance = T,
                             type = type)
    
    print(rf_nobal)
    
    result[[motif]][["rf"]] <- rf_nobal
    
    # Test data
    p_test_nobal <- predict(rf_nobal, test)

    result[[motif]][["prediction_plot"]] <- plot(test$timeToAdapt, p_test_nobal)

    # Importance measures
    ## Boruta, permutation importance, sobol MDA
    bor <- Boruta::Boruta(timeToAdapt ~ ., data = d_rf)
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
    # motif_labels <- c(molComp_labels[[motif]], 
    #                   "dataset" = TeX("Trait/selection alignment", output = "character"), 
    #                   shadow_names)[bor_order]
    
    pal_boruta <- generateCol(bor, colCode=c("#00A000","#EECC00","#DD0000","#00C0EA"),
                              col = NULL)
    
    ggplot(d_bor %>% pivot_longer(everything()) %>%
             mutate(x = fct_reorder(name, value, median)),
           aes(x = x, y = value, fill = x)) +
      geom_boxplot(show.legend = F, linewidth = 0.25) +
      theme_bw() +
      scale_fill_manual(values = pal_boruta) +
      scale_x_discrete(#labels = parse(text = motif_labels),
        guide = guide_axis(n.dodge = 2)) +
      labs(x = "Feature", y = "Boruta Importance") +
      theme(text = element_text(size = 12)) -> plt_boruta_imp
    plt_boruta_imp
    ggsave(paste0("plt_boruta_permod_tta_", motif, ".png"), plt_boruta_imp, 
           device = png, bg = "white",
           width = 12, height = 8)
    
    
    # Permutation
    predictor <- iml::Predictor$new(rf_nobal, 
                                    data = test[which(names(test) != "timeToAdapt")], 
                                    y = test$timeToAdapt)
    
    # Need to set the option future globals maxsize
    options(future.globals.maxSize = 3221225472)

    # Interaction strengths
    ia <- Interaction$new(predictor)
    
    result[[motif]][["pred"]] <- predictor
    result[[motif]][["ia"]] <- ia
    
    
    # Accumulated local effects
    ale <- FeatureEffects$new(predictor, grid.size = 10)
    ale$plot()
    
    result[[motif]][["ale"]] <- ale
    
    ale_plots <- vector(mode = "list", length = length(ale$features))
    
    #ale_labels <- motif_labels[names(ale$results)]
    
    for (i in seq_along(ale_plots)) {
      d_ale <- ale$results[[i]]
      x_label <- ale$features[i] 
      
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
        labs(x = x_label, #parse(text = x_label), 
             y = "ALE of adaptation probability")
    }
    
    alePlot <- plot_grid(plotlist = ale_plots,
                         labels= "AUTO")
    ggsave(paste0("plt_ale_permod_tta_", motif, ".png"), alePlot, device = png, bg = "white",
           width = 12, height = 9)
    
    result[[motif]][["alePlot"]] <- alePlot
    
  }
  
  return(result)
}

RunXGBPerMotif <- function(dataset, seed = NULL, train.test = c(0.7, 0.3)) {
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }
  
  motifs <- levels(dataset$model)
  result <- vector(mode = "list")
  
  for (motif in motifs) {
    set.seed(seed)
    d_rf <- dataset %>%
      filter(model == motif) %>%
      select(-model)
    
    idx <- sample(2, nrow(d_rf), replace = T, prob = train.test)
    train <- d_rf[idx == 1,]
    test <- d_rf[idx == 2,]
    
    # convert to proper format
    train_matrix <- as.data.frame(train %>% mutate(across(where(~is.factor(.)),
                                               ~as.numeric(.))))
    train_matrix <- xgb.DMatrix(data = as.matrix(train_matrix %>% select(-isAdapted)),
                                label = train_matrix$isAdapted) 
    
    test_matrix <- as.data.frame(test %>% mutate(across(where(~is.factor(.)),
                                                         ~as.numeric(.))))
    
    # Tune XGB  hyperparameters
    opt_weight <- sum(train$isAdapted == "Maladapted") / sum(train$isAdapted == "Adapted")
    
    hyper_grid <- expand.grid(
      eta = c(0.3, 0.1, 0.05),
      max_depth = 3,
      min_child_weight = 3,
      scale_pos_weight = c(1, opt_weight, 100, 1000),
      subsample = 0.5,
      colsample_bytree = 0.5,
      gamma = c(0.0, 1, 10),
      lambda = c(0.01, 0.1, 1),
      alpha = c(0.01, 0.1, 1),
      RMSE = NA,
      trees = NA
      )
    
    if (weight == Inf) {
      message(paste("No maladapted populations for motif", motif))
      next
    }
    
    pb <- progress::progress_bar$new(
      format = "  Running [:bar] :percent in :elapsedfull",
      total = nrow(hyper_grid), clear = FALSE, width = 60)
    
    for (i in seq_len(nrow(hyper_grid))) {
      seed <- 416155406
    
      xgb_model <- xgb.cv(data = train_matrix,
                    nrounds = 1000,
                    early_stopping_rounds = 50,
                    nfold = 10,
                    verbose = 0,
                    params = list(
                      objective = "reg:squarederror",
                      scale_pos_weight = hyper_grid$scale_pos_weight[i],
                      eta = hyper_grid$eta[i],
                      max_depth = hyper_grid$max_depth[i],
                      min_child_weight = hyper_grid$min_child_weight[i],
                      subsample = hyper_grid$subsample[i],
                      colsample_bytree = hyper_grid$colsample_bytree[i],
                      gamma = hyper_grid$gamma[i],
                      lambda = hyper_grid$lambda[i],
                      alpha = hyper_grid$alpha[i]
                      )
                    )
    
    
    hyper_grid$RMSE[i] <- min(xgb_model$evaluation_log$test_rmse_mean)
    hyper_grid$trees[i] <- xgb_model$early_stop$best_iteration

    pb$tick()
    }
    
    # Pick the best learning rate to continue with
    hyper_grid %>%
      filter(RMSE > 0) %>%
      arrange(RMSE) %>%
      glimpse()
    
    best_params <- as.list(dplyr::arrange(hyper_grid, RMSE)[1,] %>% 
                             select(-c(trees, RMSE)))
    best_params[["objective"]] <- "reg:squarederror"
    best_params[["scale_pos_weight"]] <- 0.0001
    best_nrounds <- dplyr::arrange(hyper_grid, RMSE)[1, "trees"]
    
    xgb_model_final <- xgb.train(
      params = best_params,
      data = train_matrix,
      nrounds = best_nrounds,
      verbose = 0
    )
    
    # Confusion matrix w/ test set
    pred <- predict(xgb_model_final, as.matrix(test_matrix %>% select(-isAdapted)))
    pred <- as.numeric(pred > 0.5)
    
    xgb_confusion <- careWt::confusionMatrix(factor(pred), factor(test_matrix$isAdapted))
    xgb_confusion
    
    # Save output model
    result[[motif]][["XGB_model"]] <- xgb_model_final
    result[[motif]][["xgb_confusion"]] <- xgb_confusion
    
    # Look at feature importance
    bor <- Boruta::Boruta(isAdapted ~ .,
                   data = d_rf)
    result[[motif]][["boruta"]] <- bor
    
    # Accumulated local effects
    predictor <- iml::Predictor$new(xgb_model_final, 
                                    data = test[which(names(test) != "isAdapted")], 
                                    y = test$isAdapted)
    
    # Need to set the option future globals maxsize
    options(future.globals.maxSize = 3221225472)
    
    # Interaction strengths
    ia <- Interaction$new(predictor)
    
    result[[motif]][["pred"]] <- predictor
    result[[motif]][["ia"]] <- ia
    
    
    # Accumulated local effects
    ale <- FeatureEffects$new(predictor, grid.size = 10)
    ale$plot()
    
    result[[motif]][["ale"]] <- ale
    
    
  }
  return(result)
}

triNumber <- function(n) {
  return (as.integer((n * (n + 1)) / 2))
}

# Read in phenotypic molecular component variance/covariance matrix
getMCVar <- function(x, means) {
  x_id <- x %>% select(c(gen, seed, modelindex, dataset, model, r))
  x <- x %>% select(-c(gen, seed, modelindex, dataset, model, r))
  
  columnsToAdd <- paste0("var_", unique(unlist(molComp_names)))
  result <- x_id
  result[,columnsToAdd] <- NA
  
  for (i in seq_len(nrow(result))) {
    x_i <- as.numeric(unlist(x[i,]))
    x_id_i <- x_id[i,]
    means_i <- as.numeric(unlist(means[i,] %>% select(starts_with("meanMC"))))
    motif <- as.character(unlist(x_id_i %>% select(model)) )
    n <- length(molComp_labels[[motif]])
    covMat <- matrix(numeric(n * n), nrow = n)
    
    # Output is stored as a triangular number, decode
    tri.n <- triNumber(n)
    
    # x is encoded in rowwise order 
    covMat[lower.tri(covMat, T)] <- x_i[1:tri.n]
    covMat <- t(covMat)
#    covMat[lower.tri(covMat, F)] <- t(covMat)[lower.tri(covMat, F)]
    
    relevantColumns <- paste0("var_", molComp_names[[motif]]) 
    
    # divide by means for scaled variance
    result[i, relevantColumns] <- as.list(diag(covMat) / means_i[1:n])
    
  }
  
  return(result)
}

getMCCov <- function(x) {
  x_id <- x %>% select(c(gen, seed, modelindex, dataset, model, r))
  x <- x %>% select(-c(gen, seed, modelindex, dataset, model, r))
  
  result <- vector(mode = "list", length = nrow(x_id))
  pb <- progress::progress_bar$new(
    format = "  Running [:bar] :percent in :elapsedfull eta: :eta",
    total = nrow(x), clear = FALSE, width = 60)
  
  for (i in seq_len(nrow(x))) {
    pb$tick()
    x_i <- as.numeric(unlist(x[i,]))
    x_id_i <- x_id[i,]
    motif <- as.character(unlist(x_id_i %>% select(model)) )
    n <- length(molComp_labels[[motif]])
    covMat <- matrix(numeric(n * n), nrow = n)
    
    # Output is stored as a triangular number, decode
    tri.n <- triNumber(n)
    
    # x is encoded in rowwise order 
    covMat[lower.tri(covMat, T)] <- x_i[1:tri.n]
    covMat <- t(covMat)
    covMat[lower.tri(covMat, F)] <- t(covMat)[lower.tri(covMat, F)]

    if (!matrixcalc::is.positive.definite(covMat)) {
      covMat <- as.matrix(Matrix::nearPD(covMat)$mat)
    }
    covMat <- (covMat + t(covMat)) / 2
    result[[i]] <- covMat
  }
  
  return(result)
}

readMCMatrix <- function(x, motif) {
  if (is.list(x)) {
    x <- unlist(x)
  }
  
  n <- length(molComp_labels[[motif]])
  result <- matrix(numeric(n * n), nrow = n)
  
  # Output is stored as a triangular number, decode
  tri.n <- triNumber(n)
  
  # x is encoded in rowwise order 
  result[lower.tri(result, T)] <- x[1:tri.n]
  result <- t(result)
  result[lower.tri(result, F)] <- t(result)[lower.tri(result, F)]
  return(result)
}

# Finds the most common
EigenTensorExperiment <- function(g_list, ids, n = 10000, seed = 123) {
  pb <- progress::progress_bar$new(
    format = "[:bar] :current/:total (:percent eta: :eta)", total = n)
  
  set.seed(seed)
  
  result <- data.frame(bootstrap = 1:n,
                       id_1 = 1:n,
                       id_2 = 1:n,
                       id_3 = 1:n,
                       mc_1 = character(n),
                       mc_2 = character(n),
                       mc_3 = character(n)
                       )
  
  # Sample a matrix from each treatment, compare them all
  # Dataset - parallel/orthogonal/randomised
  samples <- ids %>%
    slice_sample(n = n, by = "dataset", replace = T)
  samples$label <- as.numeric(samples$label)
  
  # Split so we have n groups
  samples <- split(samples, rep(1:n, times = ceiling(nrow(samples) / n),
                                length.out = nrow(samples)))
  
  # Iterate through samples
  for (i in seq_len(n)) {
    pb$tick()
    labels <- samples[[i]]$label
    x <- g_list[labels]
    et <- EigenTensorDecomposition(x)
    
    # Look for smallest ET, what MCs describe the smallest divergence between
    # populations
    eg <- eigen(et$matrices[,,dim(et$matrices)[3]])
    rownames(eg$vectors) <- rownames(et$matrices[,,1])
    top3 <- sort(abs(eg$vectors[,1]), decreasing = T)[1:3]
    
    
    result[i,] <- c(i, labels,
                    names(top3))
  }
  return(result)
}

cEvolPerTrait <- function(G) {
  n <- nrow(G)

  result <- numeric(n)
  names(result) <- colnames(G)

    for (i in seq_len(n)) {
      Gy <- diag(G)[i]
      Gx <- G[-i, -i]
      Gyx <- matrix(G[i,-i], nrow = 1)
      Gxy <- t(Gyx)
      
      result[i] <- Gy - (Gyx %*% solve(Gx) %*% Gxy)
    }
  
  return(result)
}

autPerTrait <- function(G) {
  n <- nrow(G)
  Ginv <- solve(G)
  result <- 1 / (diag(Ginv) * diag(G))
  
  return(result)
}

ConditionalEvolvabilityExperiment <- function(G_list, id, beta_list) {
  n <- length(G_list)
  
  result <- data.frame(seed = id$seed,
                       modelindex = id$modelindex,
                       dataset = id$dataset)
  
  result[, paste0("cev_", names(all_molcomp_features))] <- NA
  result[, paste0("aut_", names(all_molcomp_features))] <- NA
  
  
  for (i in seq_len(n)) {
    cev <- cEvolPerTrait(G_list[[i]])
    aut <- autPerTrait(G_list[[i]])
    result[i, paste0("cev_", names(cev))] <- cev
    result[i, paste0("aut_", names(cev))] <- aut
  }
  return(result)
}

CalculateMCEffects <- function(dataset, formula_list, seed = 123) {
  out <- list()
  
  pb <- progress::progress_bar$new(
    format = "[:bar] :current/:total (:percent eta: :eta)", total = length(formula_list),
    show_after = 0)
  
  pb$tick(0)
  
  for (formula in formula_list) {
    response <- all.vars(formula)[1]
    predictors <- all.vars(formula)[-1]
    
    bor <- Boruta::Boruta(as.formula(paste0(response, "~",
                                            paste(predictors, collapse = "+"))), 
                          dataset)
    
    # Choose top three most important predictors, use these for the beta regression
    bor_sorted <- colMeans(bor$ImpHistory)
    bor_sorted <- bor_sorted[!str_detect(names(bor_sorted), "shadow")]
    bor_sorted <- bor_sorted[!is.infinite(bor_sorted)]
    bor_sorted <- sort(bor_sorted, decreasing = T)[1:3]
    top_predictors <- names(bor_sorted)
    
    out[[response]][["boruta"]] <- bor
    
    #browser()
    beta_model <- NA
    null_model <- NA
    # Start with beta regression
    beta_model <- betareg::betareg(as.formula(paste(response, "~", 
                                                    paste(top_predictors, collapse = "*"))), dataset)
    null_model <- betareg::betareg(paste(response, "~ 1"), dataset)
    out[[response]][["beta_model"]] <- beta_model
    out[[response]][["beta_model_lrtest"]]  <- lmtest::lrtest(beta_model, null_model)
    
    # tryCatch({
    #   beta_model <- betareg::betareg(as.formula(paste(response, "~", 
    #                                        paste(top_predictors, collapse = "*"))), dataset)
    #   null_model <- betareg::betareg(paste(response, "~ 1"), dataset)
    # },
    # error = function(e) {
    #   print(paste0("Error: ", e))
    #   #beta_model <- NA
    #   #null_model <- NA
    # },
    # warning = function(w) {
    #   #beta_model <- NA
    # },
    # finally = {
    #   #browser()
    #   if (!all(is.na(beta_model))) {
    #     out[[response]][["beta_model"]] <- beta_model
    #     out[[response]][["beta_model_lrtest"]]  <- lmtest::lrtest(beta_model, null_model)
    #   }
    # })
    
    # Shapley value regression to find relative importance of each feature
    sv <- ShapleyValue::shapleyvalue(unlist(dataset[,response]),
                                        as.data.frame(dataset[,predictors]))
    out[[response]][["shapley"]] <- sv
    
    # LMG R2 decomposition
    lmg <- sensitivity::lmg(as.data.frame(dataset[,predictors]),
                            unlist(dataset[,response]))
    out[[response]][["lmg"]] <- lmg
    
    # Random forest
    set.seed(seed)
    
    idx <- sample(2, nrow(dataset), replace = T, prob = c(0.7, 0.3))
    train <- dataset[idx == 1,c(response, predictors)]
    test <- dataset[idx == 2,c(response, predictors)]
    
    rf <- randomForest(formula = formula,
                                 data = train,
                                 ntree = 500,
                                 proximity = T,
                                 importance = T,
                                 type = "regression")
    out[[response]][["rf"]] <- rf
    
    # Test prediction accuracy
    prediction <- predict(rf, test)
    pred_dat <- data.frame(obs = unlist(test[,response]),
                           pred = prediction)
    out[[response]][["rf_rmse"]] <- caret::defaultSummary(pred_dat)
    
    # Variable importance
    predictor <- iml::Predictor$new(rf, 
                                    data = test[, predictors], 
                                    y = test[,response])
    out[[response]][["predictor"]] <- predictor
    
    # Need to set the option future globals maxsize
    options(future.globals.maxSize = 3221225472)
    imp <- iml::FeatureImp$new(predictor,
                               loss = "mae",
                               n.repetitions = 100)
    out[[response]][["imp"]] <- imp
    
    ale <- FeatureEffects$new(predictor, method = "ale", grid.size = 10)
    out[[response]][["ale"]] <- ale
    interact <- Interaction$new(predictor, grid.size = 10)
    out[[response]][["interact"]] <- interact
    pb$tick()
  }
  
  return(out)
}

AdjMatFromCorMat <- function(m, threshold) {
  result <- matrix(0, nrow = nrow(m), ncol = ncol(m))
  
  result[abs(m)>=threshold] <- 1
  
  return(result)
}

# PCA similarity
bootKrzCorFn <- function(x, group = "", PCASim = F) {
  require(evolqg)
  require(dplyr)
  
  fn <- ifelse(PCASim, evolqg::PCAsimilarity, evolqg::KrzCor)
  
  if (group != "") {
    grps <- unique(x[,group])
    nGrps <- length(grps)
    
    
    # output data frame
    res <- data.frame(group1 = character(length(grps)^2),
                      group2 = character(length(grps)^2),
                      krzCor = numeric(length(grps)^2))
    
    # Temporary data frame for filling inner loop
    res_tmp <- data.frame(group1 = character(length(grps)),
                          group2 = character(length(grps)),
                          krzCor = numeric(length(grps)))
    
    for (i in seq_along(grps)) {
      for (j in seq_along(grps)) {
        # Sample matrices in different groups
        if (i != j) {
          g_1 <- slice_sample(x[group == grps[i]], n = 1)
          g_2 <- slice_sample(x[group == grps[j]], n = 1)
        } else 
        {
          g <- slice_sample(x[group == grps[i]], n = 2, replace = F) # Avoid sampling the same matrix
          g_1 <- g[1,]
          g_2 <- g[2,]
        }
        res_tmp$group1[j] <- as.character(g_1[1,group])
        res_tmp$group2[j] <- as.character(g_2[1,group])
        res_tmp$krzCor[j] <- fn(g_1$g[[1]], g_2$g[[1]])
      }
      indices <- (nGrps*(i-1) + 1):(nGrps*i)
      res[indices,] <- res_tmp
    }
    return(res)
  }
  
  # If group is "", sample two random matrices and return that
  g <- slice_sample(x, n = 2)
  g1 <- g[1,]
  g2 <- g[2,]
  return(fn(g1$g[[1]], g2$g[[1]]))
}

# Split by model comparison
GetModelComparison <- function(model1, model2, model_names) {
  result <- character(length(model1))
  # assign by priority (NAR > PAR > FFLC1 > FFLI1 > FFBH)
  model1_index <- match(model1, model_names)
  model2_index <- match(model2, model_names)
  
  if (any(is.na(model1_index)) | any(is.na(model2_index))) {
    return(NA)
  }
  
  # If model2's index is greater than model1,
  # we put model1 first
  model1First_idx <- model2_index > model1_index
  model2First_idx <- model1_index >= model2_index
  result[model2_index > model1_index] <- paste(model1[model1First_idx], "vs", model2[model1First_idx])
  # otherwise model2 comes first, if models are the same these are
  # included here too
  result[model1_index >= model2_index] <- paste(model2[model2First_idx], "vs", model1[model2First_idx])
  return(result)
}


BootSelResponseDiff <- function(x, group, nSkewers = 100) {
  require(evolqg)
  require(dplyr)

  if (group != "") {
    grps <- unique(x[,group])
    nGrps <- length(grps)
    
    
    # output data frame
    res <- data.frame(group1 = character(length(grps)^2),
                      group2 = character(length(grps)^2),
                      SRD = numeric(length(grps)^2))
    
    # Temporary data frame for filling inner loop
    res_tmp <- data.frame(group1 = character(length(grps)),
                          group2 = character(length(grps)),
                          SRD = numeric(length(grps)))
    
    for (i in seq_along(grps)) {
      for (j in seq_along(grps)) {
        # Sample matrices in different groups
        if (i != j) {
          g_1 <- slice_sample(x[group == grps[i]], n = 1)
          g_2 <- slice_sample(x[group == grps[j]], n = 1)
        } else 
        {
          g <- slice_sample(x[group == grps[i]], n = 2, replace = F) # Avoid sampling the same matrix
          g_1 <- g[1,]
          g_2 <- g[2,]
        }
        res_tmp$group1[j] <- as.character(g_1[1,group])
        res_tmp$group2[j] <- as.character(g_2[1,group])
        res_tmp$SRD[j] <- evolqg::SRD(g_1$g[[1]], g_2$g[[1]], iterations = nSkewers)
      }
      indices <- (nGrps*(i-1) + 1):(nGrps*i)
      res[indices,] <- res_tmp
    }
    return(res)
  }
  
  # If group is "", sample two random matrices and return that
  g <- slice_sample(x, n = 2)
  g1 <- g[1,]
  g2 <- g[2,]
  return(evolqg::SRD(g1$g[[1]], g2$g[[1]], iterations = nSkewers))
}


GetTraces <- function(matlist, id) {
  result <- purrr::map(seq_along(matlist), function(i) {
    g <- matlist[[i]]
    
    trace_i <- matrixcalc::matrix.trace(g)
    
    result <- id[i,]
    
    result$trace <- trace_i
    return(result)
  })
  
  return(data.table::rbindlist(result))
}
