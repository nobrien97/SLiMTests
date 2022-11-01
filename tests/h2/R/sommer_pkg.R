library(sommer)
source("helpFns.R")


# Load data
haplos <- read.csv("../data/slim_haplo_sbst.csv", header = F)[,-(1:3)][1001:2000,]
names(haplos) <- paste0("q_", 1:ncol(haplos))
ped <- read.csv("../data/slim_pedigree_sbst.csv", header = F)[,-(1:3)][501:1000,]
names(ped) <- c("id", "mother", "father")
genMap <- read.csv("", header = F)[,-(1:3)]

# Organise data
genos <- hapToGen(haplos)
phenos <- unlist(read.csv("../data/slim_sampled_pheno_sbst.csv", header=F)[2,-(1:3)])
names(phenos) <- NULL
geno_names <- paste0("i_", 1:nrow(genos))

pheno_dat <- data.frame(y = phenos,
                        id = geno_names)

# Narrow sense

# calculate relationship matrices using pedigree

A <- A.mat(as.matrix(genos)-1)
D <- D.mat(as.matrix(genos)-1)
AA <- E.mat(as.matrix(genos)-1)
AD <- E.mat(as.matrix(genos)-1, type = "A#D")
DD <- E.mat(as.matrix(genos)-1, type = "D#D")

dimnames(A) <- list(geno_names, geno_names)
dimnames(D) <- list(geno_names, geno_names)
dimnames(AA) <- list(geno_names, geno_names)
dimnames(AD) <- list(geno_names, geno_names)
dimnames(DD) <- list(geno_names, geno_names)

# Add additional columns for dominance and epistasis IDs
pheno_dat$idd <- pheno_dat$id
pheno_dat$idaa <- pheno_dat$id
pheno_dat$idad <- pheno_dat$id
pheno_dat$iddd <- pheno_dat$id
pheno_dat$y <- scale(pheno_dat$y)

ans.ADE <- mmer(y~1,
                random=~vsr(id,Gu=A) + vsr(idd,Gu=D) + vsr(idaa, Gu=AA),
                rcov=~units,
                data=pheno_dat)

summary(ans.ADE)$varcomp
vpredict(ans.ADE, h2 ~ (V1) / ( V1+V2+V3+V4) ) # narrow-sense
vpredict(ans.ADE, h2 ~ (V1+V2+V3) / ( V1+V2+V3+V4) ) # broad-sense

# Using mmec for pedigree information
phenos <- phenos + rnorm(500)

pheno_dat <- data.frame(y = phenos, 
                        GCA1 = as.factor(ped$mother),
                        GCA2 = as.factor(ped$father),
                        SCA = as.factor(ped$id))

pheno_dat$idd <- as.numeric(pheno_dat$SCA)
pheno_dat$ide <- as.numeric(pheno_dat$SCA)


ans.ADE <- mmec(y~1,
                random =~GCA1+GCA2+SCA,
                rcov=~units,
                data=pheno_dat)

suma <- summary(ans.ADE)$varcomp
Vgca <- sum(suma[1:2,1])
Vsca <- suma[3,1]
Ve <- suma[4,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )
