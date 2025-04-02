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
WRITE_PATH_MRR <- paste0("/scratch/ht96/nb9894/newMotifs/h2/getH2/out_h2_", run_chunk, "_mrr.csv")
WRITE_PATH_MKR <- paste0("/scratch/ht96/nb9894/newMotifs/h2/getH2/out_h2_", run_chunk, "_mkr.csv")

# Load functions for loading relatedness/haplotype matrices
source("~/tests/standingVar/getH2/R/helpFns.R")
#source("../../../standingVar/getH2/R/helpFns.R")

# Two methods: kernel and ridge regression based on pedigree or loci 
# Pedigree is really an estimate of breeding values, loci should be more
# accurate: but we will try both

# Load in haplotypes
haplos <- scan(paste0("slim_haplo_sbst_", run, ".csv"), sep = ",")[-(1:3)]
haplos_fix <- scan(paste0("slim_haplo_fix_sbst_", run, ".csv"), sep = ",")[-(1:3)]

haplos <- decompressHap(haplos, 2000, 1386)
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

# Get model info from modulo: type of model offset by 5
model_mod <- (run_info[3] %% 5)

if (model_mod == 0) {
  model_mod <- 5
}

switch (model_mod,
  model <- "NAR",
  model <- "PAR",
  model <- "FFLC1",
  model <- "FFLI1",
  model <- "FFBH"
)

# Molecular components
relPheno <- scan(paste0("slim_moltrait_sbst_", run, ".csv"), sep = ",")[-(1:3)]
names(relPheno) <- NULL
ind_names <- paste0(1:1000)


if (model == "NAR" | model == "PAR") {
  # Get molecular component dataframe
  relPheno_dat <- data.frame(aZ    = log(relPheno[seq(1, length(relPheno), by = 7)]),
                             bZ    = log(relPheno[seq(2, length(relPheno), by = 7)]),
                             KZ    = log(relPheno[seq(3, length(relPheno), by = 7)]),
                             KXZ   = log(relPheno[seq(4, length(relPheno), by = 7)]),
                             base  = log(relPheno[seq(5, length(relPheno), by = 7)]),
                             n     = log(relPheno[seq(6, length(relPheno), by = 7)]),
                             XMult = log(relPheno[seq(7, length(relPheno), by = 7)]))
} else if (model == "FFLC1" | model == "FFLI1") {
  # Get molecular component dataframe
  relPheno_dat <- data.frame(aY    = log(relPheno[seq(1, length(relPheno), by = 9)]),
                             bY    = log(relPheno[seq(2, length(relPheno), by = 9)]),
                             KY    = log(relPheno[seq(3, length(relPheno), by = 9)]),
                             aZ    = log(relPheno[seq(4, length(relPheno), by = 9)]),
                             bZ    = log(relPheno[seq(5, length(relPheno), by = 9)]),
                             KXZ   = log(relPheno[seq(6, length(relPheno), by = 9)]),
                             base  = log(relPheno[seq(7, length(relPheno), by = 9)]),
                             n     = log(relPheno[seq(8, length(relPheno), by = 9)]),
                             XMult = log(relPheno[seq(9, length(relPheno), by = 9)]))
} else if (model == "FFBH") {
  relPheno_dat <- data.frame(aX    = log(relPheno[seq(1, length(relPheno), by = 11)]),
                             KZX   = log(relPheno[seq(2, length(relPheno), by = 11)]),
                             aY    = log(relPheno[seq(3, length(relPheno), by = 11)]),
                             bY    = log(relPheno[seq(4, length(relPheno), by = 11)]),
                             KY    = log(relPheno[seq(5, length(relPheno), by = 11)]),
                             aZ    = log(relPheno[seq(6, length(relPheno), by = 11)]),
                             bZ    = log(relPheno[seq(7, length(relPheno), by = 11)]),
                             KXZ   = log(relPheno[seq(8, length(relPheno), by = 11)]),
                             base  = log(relPheno[seq(9, length(relPheno), by = 11)]),
                             n     = log(relPheno[seq(10, length(relPheno), by = 11)]),
                             XMult = log(relPheno[seq(11, length(relPheno), by = 11)]))  
}

# Scale and center
relPheno_dat_scaled <- scale(relPheno_dat)

# Run kernel regression w/ eigendecomposition depending on the model
relPheno_mat <- as.matrix(relPheno_dat_scaled)
mkr_result <- mkrOrError(relPheno_mat, A)
mrr_result <- mrrOrError(relPheno_mat, X)

# Do the same with fitness
relFitness <- scan(paste0("slim_pheno_sbst_", run, ".csv"), sep = ",")[-(1:3)]
# only read first 1000 values, these are fitnesses
relFitness <- relFitness[1:1000]
names(relFitness) <- NULL
relFitness_dat <- data.frame(w    = relFitness)
relFitness_dat_scaled <- scale(relFitness_dat)

# Run kernel regression w/ eigendecomposition depending on the model
relFitness_mat <- as.matrix(relFitness_dat_scaled)
mkr_fitness_result <- mkrOrError(relFitness_mat, A)
mrr_fitness_result <- mrrOrError(relFitness_mat, X)

if (is.na(mkr_result[1]) & is.na(mrr_result[1])) {
  print(paste("Couldn't solve model in run ", run, "- closing R."))
  q(save = "no")
}

# Print output

#Extract results
mkr_result_df <- data.frame(gen = run_info[1],
                            seed = as.character(run_info[2]),
                            modelindex = run_info[3]
)

mrr_result_df <- data.frame(gen = run_info[1],
                            seed = as.character(run_info[2]),
                            modelindex = run_info[3]
)


# fitness variance
if (!is.na(mkr_fitness_result[1])) {
  mkr_result_df[1, "VA_w"] <- mkr_fitness_result$Vb
  mkr_result_df[1, "h2_w"] <- mkr_fitness_result$h2
}

if (!is.na(mrr_fitness_result[1])) {
  mrr_result_df[1, "VA_w"] <- mrr_fitness_result$Vb
  mrr_result_df[1, "h2_w"] <- mrr_fitness_result$h2
}

allMolTraitNames <- c("aX", "KZX", "aY", "bY", "KY", "aZ", "bZ", "KZ", "KXZ", "base", "n", "XMult")

# Write results to separate files
if (!is.na(mkr_result[1])) {
  # Expand results with NA for Additive and ODE cases
  if (model == "NAR" | model == "PAR") {
    molTraitNames <- c("aZ", "bZ", "KZ", "KXZ", "base", "n", "XMult")
    colnames(mkr_result$Vb) <- molTraitNames
  }
  
  if (model == "FFLC1" | model == "FFLI1") {
    molTraitNames <- c("aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult")
    colnames(mkr_result$Vb) <- molTraitNames
  }

  if (model == "FFBH") {
    molTraitNames <- c("aX", "KZX", "aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult")
    colnames(mkr_result$Vb) <- molTraitNames
  }

  rownames(mkr_result$Vb) <- colnames(mkr_result$Vb)
  G <- matrix(NA, nrow = 12, ncol = 12)
  colnames(G) <- allMolTraitNames
  rownames(G) <- colnames(G)
  
  G[rownames(mkr_result$Vb), colnames(mkr_result$Vb)] <- mkr_result$Vb
  
  h2 <- rep(NA, 12)
  names(h2) <- allMolTraitNames
  h2[molTraitNames] <- mkr_result$h2
  
  cov_terms <- combn(allMolTraitNames, 2)
  cov_terms <- paste(cov_terms[1,], cov_terms[2,], sep = "_") 
    
  # Get genetic variances
  mkr_result_df[1, paste("VA", rownames(G), sep = "_")] <- diag(G)
  
  # get genetic covariances
  mkr_result_df[1, cov_terms] <- G[t(upper.tri(G))]
  
  # heritability
  mkr_result_df[1, paste("h2", allMolTraitNames, sep = "_")] <- h2
  

  # Output file
  print(paste0("Writing mkr output for model ", run, "..."))
  write.table(mkr_result_df, WRITE_PATH_MKR, 
              sep = ",", row.names = F, col.names = F)
}

if (!is.na(mrr_result[1])) {
  # Expand results with NA for Additive and ODE cases
  if (model == "NAR" | model == "PAR") {
    molTraitNames <- c("aZ", "bZ", "KZ", "KXZ", "base", "n", "XMult")
    colnames(mrr_result$Vb) <- molTraitNames
  }
  
  if (model == "FFLC1" | model == "FFLI1") {
    molTraitNames <- c("aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult")
    colnames(mrr_result$Vb) <- molTraitNames
  }

  if (model == "FFBH") {
    molTraitNames <- c("aX", "KZX", "aY", "bY", "KY", "aZ", "bZ", "KXZ", "base", "n", "XMult")
    colnames(mrr_result$Vb) <- molTraitNames
  }

  rownames(mrr_result$Vb) <- colnames(mrr_result$Vb)

  
  G <- matrix(NA, nrow = 12, ncol = 12)
  colnames(G) <- allMolTraitNames
  rownames(G) <- colnames(G)
  
  G[rownames(mrr_result$Vb), colnames(mrr_result$Vb)] <- mrr_result$Vb
  
  h2 <- rep(NA, 12)
  names(h2) <- allMolTraitNames
  h2[molTraitNames] <- mrr_result$h2
  
  cov_terms <- combn(allMolTraitNames, 2)
  cov_terms <- paste(cov_terms[1,], cov_terms[2,], sep = "_") 
  
  # Get genetic variances
  mrr_result_df[1, paste("VA", rownames(G), sep = "_")] <- diag(G)
  
  # get genetic covariances
  mrr_result_df[1, cov_terms] <- G[t(upper.tri(G))]
  
  # heritability
  mrr_result_df[1, paste("h2", allMolTraitNames, sep = "_")] <- h2
  
  # Output file
  print(paste0("Writing mrr output for model ", run, "..."))
  write.table(mrr_result_df, WRITE_PATH_MRR, 
              sep = ",", row.names = F, col.names = F)
}
######################################################################


