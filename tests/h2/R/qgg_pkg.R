library(qgg)
source("helpFns.R")

# Load data
haplos <- read.csv("~/Desktop/haplos.csv", header = F)
names(haplos) <- paste0("q_", 1:ncol(haplos))
ped <- read.csv("~/Desktop/ped.csv")
genMap <- read.csv("~/Desktop/test.csv")

# Organise data
genos <- hapToGen(haplos)
phenos <- unlist(read.csv("~/Desktop/phenos.csv", header=F))
names(phenos) <- NULL
geno_names <- paste0("i_", 1:nrow(genos))

W <- genos

nMinor <- rowSums(W, na.rm=T)
nAlleles <- rowSums(!is.na(W))
p <- nMinor/(2*nAlleles)
min(p)
max(p)

W <- t(W)
W <- scale(W)
dim(W)

pheno_dat <- data.frame(y = phenos,
                        id = geno_names)
