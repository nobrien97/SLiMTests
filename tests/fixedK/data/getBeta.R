library(sommer)
library(tidyverse)

# load data
setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")

source("/mnt/c/GitHub/SLiMTests/tests/h2/R/helpFns.R")

haplos <- read_csv("/mnt/d/SLiMTests/tests/fixedK/slim_haplo_sbst.csv", col_names = F)[, -(1:3)]
phenos <- unlist(read_csv("/mnt/d/SLiMTests/tests/fixedK/slim_indPheno_sbst.csv", col_names = F))

names(haplos) <- paste0("q_", seq_len(ncol(haplos)))
names(phenos) <- NULL
run_info <- phenos[1:3]
phenos <- phenos[-(1:3)]

# Get genotypes from haplos
genos <- hapToGen(haplos)
geno_names <- paste0("i_", seq_len(nrow(genos)))
pheno_dat <- data.frame(y = phenos,
                        id = as.factor(geno_names))

# Calculate genomic relationship matrices
A <- A.mat(as.matrix(genos)-1)
D <- D.mat(as.matrix(genos)-1)
AA <- A*A

# Add names to dims
dimnames(A) <- list(geno_names, geno_names)
dimnames(D) <- list(geno_names, geno_names)
dimnames(AA) <- list(geno_names, geno_names)

# Add extra IDs for D and AA
pheno_dat$idd <- pheno_dat$id
pheno_dat$idaa <- pheno_dat$id

######################################################################
# This code is adapted from Onigo et al. 2021, thanks to Akio Onogi
# for sharing code samples

# Fit models using A, D, and AA GRMs only to avoid overfitting
Result.A_D_AA <- mmer(y~1,
                      random = ~vsr(id, Gu = A) +
                        vsr(idd, Gu = D) +
                        vsr(idaa, Gu = AA),
                      rcov = ~ units,
                      data = pheno_dat)

Result_GWAS <- GWAS(y~1,
     random = ~vsr(id, Gu = A) +
       vsr(idd, Gu = D) +
       vsr(idaa, Gu = AA),
     rcov = ~ units,
     data = pheno_dat, 
     M = as.matrix(genos)-1,
     gTerm = "u:id")
