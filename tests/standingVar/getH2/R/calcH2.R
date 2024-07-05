# R Script to calculate heritability given some data
library(bWGR)
library(purrr)

# Get command line arguments
## 1: run
## 2: chunk
args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])
chunk <- as.numeric(args[2])
run_chunk <- paste(run, chunk, sep = "_")

# Path to write output
WRITE_PATH_MRR <- paste0("/scratch/ht96/nb9894/standingVar/getH2/out_h2_", run_chunk, "_mrr.csv")
WRITE_PATH_MKR <- paste0("/scratch/ht96/nb9894/standingVar/getH2/out_h2_", run_chunk, "_mkr.csv")

# Load functions for loading relatedness/haplotype matrices
source("~/tests/standingVar/getH2/R/helpFns.R")

# Two methods: kernel and ridge regression based on pedigree or loci 
# Pedigree is really an estimate of breeding values, loci should be more
# accurate: but we will try both

# Load in haplotypes
haplos <- scan(paste0("slim_haplo_sbst_", run, ".csv"), sep = ",")[-(1:3)]
haplos_fix <- scan(paste0("slim_haplo_fix_sbst_", run, ".csv"), sep = ",")[-(1:3)]

haplos <- decompressHap(haplos, 2000, 1024)
names(haplos) <- paste0("q_", seq_len(ncol(haplos)))
haplos <- addFixedHaplos(haplos, haplos_fix)

# Get genotypes from haplos (0, 1, 2 coded)
genos <- hapToGen(haplos)
geno_names <- paste0("i_", seq_len(nrow(genos)))
# Centralise genotypes
X <- CNT(genos)

# Load relatedness matrix and molecular components
relPos <- scan(paste0("slim_relPos_sbst_", run, ".csv"), what = numeric(), sep = ",")
relVals <- scan(paste0("slim_relVals_sbst_", run, ".csv"), what = numeric(), sep = ",")[-(1:3)]

# From the phenos file, extract gen, seed, modelindex and remove them
## Use these later to identify output
run_info <- relPos[1:3]
relPos <- relPos[-(1:3)]

# Additive relatedness matrix (Wright coefficients)
A <- decompressRel(relPos, relVals, 1000)

# Error catching
mkrOrError <- possibly(mkr, otherwise = NA)
mrrOrError <- possibly(mrr, otherwise = NA)

# Get model info from modulo: type of model offset by 3
model_mod <- (run_info[3] %% 3) + 1

switch (model_mod,
  model <- "K",
  model <- "Add",
  model <- "ODE"
)

if (model == "Add") {
  # No molecular components, only phenotype
  relPheno <- scan(paste0("slim_pheno_sbst_", run, ".csv"), sep = ",")[-(1:3)]
  names(relPheno) <- NULL
  
  relPheno_dat <- data.frame(Z   = relPheno)
  # Scale and center
  relPheno_dat_scaled <- scale(relPheno_dat[,1])

  relPheno_mat <- as.matrix(relPheno_dat_scaled[,1])
  mkr_result <- mkrOrError(relPheno_mat, A)
  mrr_result <- mrrOrError(relPheno_mat, X)
  
} else {
  # Molecular components
  relPheno <- scan(paste0("slim_moltrait_sbst_", run, ".csv"), sep = ",")[-(1:3)]
  names(relPheno) <- NULL
  ind_names <- paste0(1:1000)
  
  # Convert to data frame: log mol trait values to transform to normal
  relPheno_dat <- data.frame(Z   = relPheno[seq(1, length(relPheno), by = 5)],
                             aZ  = log(relPheno[seq(2, length(relPheno), by = 5)]),
                             bZ  = log(relPheno[seq(3, length(relPheno), by = 5)]),
                             KZ  = log(relPheno[seq(4, length(relPheno), by = 5)]),
                             KXZ = log(relPheno[seq(5, length(relPheno), by = 5)]))
  
  # Scale and center phenotypes
  relPheno_dat_scaled <- scale(relPheno_dat[,1:5])

  # Run kernel regression w/ eigendecomposition depending on the model
  relPheno_mat <- as.matrix(relPheno_dat_scaled[,1:5])
}

if (model == "ODE") {
  # No K values for the regular models
  mkr_result <- mkrOrError(relPheno_mat[,1:3], A)
  mrr_result <- mrrOrError(relPheno_mat[,1:3], X)
  
} else if (model == "K") {
  # All K values
  mkr_result <- mkrOrError(relPheno_mat[,1:5], A)
  mrr_result <- mrrOrError(relPheno_mat[,1:5], X)
} else if (model != "Add") {
  print(paste("Couldn't find model type in run ", run, "- closing R."))
  q(save = "no")
}

# Print output
if (is.na(mkr_result[1]) & is.na(mrr_result[1])) {
  print(paste("Couldn't solve model in run ", run, "- closing R."))
  q(save = "no")
}

# Write results to separate files
if (!is.na(mkr_result[1])) {
  # Expand results with NA for Additive and ODE cases
  molTraitNames <- c("Z", "aZ", "bZ", "KZ", "KXZ")
  colnames(mkr_result$Vb) <- molTraitNames[1:ncol(mkr_result$Vb)]
  rownames(mkr_result$Vb) <- colnames(mkr_result$Vb)
  
  G <- matrix(NA, nrow = 5, ncol = 5)
  colnames(G) <- molTraitNames
  rownames(G) <- colnames(G)
  
  G[rownames(mkr_result$Vb), colnames(mkr_result$Vb)] <- mkr_result$Vb
  
  h2 <- rep(NA, 5)
  h2[1:length(mkr_result$h2)] <- mkr_result$h2
  
  cov_terms <- combn(molTraitNames, 2)
  cov_terms <- paste(cov_terms[1,], cov_terms[2,], sep = "_") 
  
  #Extract results
  mkr_result_df <- data.frame(gen = run_info[1],
                              seed = as.character(run_info[2]),
                              modelindex = run_info[3]
  )
  # Get genetic variances
  mkr_result_df[1, paste("VA", rownames(G), sep = "_")] <- diag(G)
  
  # get genetic covariances
  mkr_result_df[1, cov_terms] <- G[t(upper.tri(G))]
  
  # heritability
  mkr_result_df[1, paste("h2", molTraitNames, sep = "_")] <- h2
  
  # Output file
  print(paste0("Writing mkr output for model ", run, "..."))
  write.table(mkr_result_df, WRITE_PATH_MKR, 
              sep = ",", row.names = F, col.names = F)
}

if (!is.na(mrr_result[1])) {
  # Expand results with NA for Additive and ODE cases
  molTraitNames <- c("Z", "aZ", "bZ", "KZ", "KXZ")
  colnames(mrr_result$Vb) <- molTraitNames[1:ncol(mrr_result$Vb)]
  rownames(mrr_result$Vb) <- colnames(mrr_result$Vb)
  
  G <- matrix(NA, nrow = 5, ncol = 5)
  colnames(G) <- molTraitNames
  rownames(G) <- colnames(G)
  
  G[rownames(mrr_result$Vb), colnames(mrr_result$Vb)] <- mrr_result$Vb
  
  h2 <- rep(NA, 5)
  h2[1:length(mrr_result$h2)] <- mrr_result$h2
  
  cov_terms <- combn(molTraitNames, 2)
  cov_terms <- paste(cov_terms[1,], cov_terms[2,], sep = "_") 
  
  #Extract results
  mrr_result_df <- data.frame(gen = run_info[1],
                              seed = as.character(run_info[2]),
                              modelindex = run_info[3]
  )
  # Get genetic variances
  mrr_result_df[1, paste("VA", rownames(G), sep = "_")] <- diag(G)
  
  # get genetic covariances
  mrr_result_df[1, cov_terms] <- G[t(upper.tri(G))]
  
  # heritability
  mrr_result_df[1, paste("h2", molTraitNames, sep = "_")] <- h2
  
  # Output file
  print(paste0("Writing mrr output for model ", run, "..."))
  write.table(mrr_result_df, WRITE_PATH_MRR, 
              sep = ",", row.names = F, col.names = F)
}
######################################################################


