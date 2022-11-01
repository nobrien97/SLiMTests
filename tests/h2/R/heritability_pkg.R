# Using heritability package to calc h2
library(heritability)
source("helpFns.R")

# Load data
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

# Calc h2
h2 <- marker_h2(scale(phenos), geno_names, K = rel_k, max.iter = 1000)

