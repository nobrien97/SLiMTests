# heritability using lme4GS
library(lme4GS)
library(sommer)
source("helpFns.R")


haplos <- read.csv("~/Desktop/haplos.csv", header = F)
names(haplos) <- paste0("q_", 1:ncol(haplos))
ped <- read.csv("~/Desktop/ped.csv")
genMap <- read.csv("~/Desktop/test.csv")
pheno_means <- read.csv("~/Desktop/mean.csv", header = F)

# Organise data
genos <- hapToGen(haplos)
rel_k <- read.csv("../k_mat.csv", header = F)
rel_k <- as.matrix(rel_k)
colnames(rel_k) <- paste0("i_", 1:nrow(genos))
rownames(rel_k) <- paste0("i_", 1:nrow(genos))
phenos <- unlist(read.csv("~/Desktop/phenos.csv", header=F))
names(phenos) <- NULL
geno_names <- paste0("i_", 1:nrow(genos))

# Set-up phenotype dataframe
pheno_dat <- data.frame(y = phenos, 
                        GCA1 = as.character(ped$mother),
                        GCA2 = as.character(ped$father),
                        SCA = as.character(ped$id))

# Calculate A matrix
A <- A.mat(as.matrix(genos)-1)
dimnames(A) <- list(pheno_dat$GCA1, pheno_dat$GCA2)


# Genomic relationship matrix for parent 1
GCA1 <- unique(pheno_dat$GCA1)
selected <- rownames(A) %in% GCA1
K1 <- A[selected, selected]

# Genomic relationship matrix for parent 2
GCA2 <- unique(pheno_dat$GCA2)
selected <- rownames(A) %in% GCA2
K2 <- A[selected, selected]

# kronecker
K3 <- kronecker(K1, K2, make.dimnames=T)

