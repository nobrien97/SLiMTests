library(sommer)
source("helpFns.R")


# Load data
haplos <- read.csv("../data/slim_haplo_sbst.csv", header = F)[,-(1:3)][1001:2000,]
names(haplos) <- paste0("q_", seq_len(ncol(haplos)))
ped <- read.csv("../data/slim_pedigree_sbst.csv", header = F)[,-(1:3)][501:1000,]
names(ped) <- c("id", "mother", "father")

# Organise data
genos <- hapToGen(haplos) 
phenos <- unlist(read.csv("../data/slim_sampled_pheno_sbst.csv", header=F)[2,])
names(phenos) <- NULL
run_info <- phenos[1:3]
phenos <- phenos[-(1:3)]


geno_names <- paste0("i_", 1:nrow(genos))

pheno_dat <- data.frame(y = phenos,
                        id = as.factor(geno_names))

# Narrow sense

# calculate relationship matrices using pedigree

#Calculate the genomic relationship matrices##########################
##Geno is a SNP genotype matrix whose entries are coded as 0, 1, and 2
##Geno is a Ni (number of individuals) x Nm (number of markers) matrix
Ni <- nrow(genos)
Nm <- ncol(genos)

Af <- colSums(genos) / (2 * Ni)
Maf <- Af
Maf[Maf > 0.5] <- 1 - Maf[Maf > 0.5]
Use_SNP <- which(Maf >= 0.05)
Nm.use <- length(Use_SNP)

#Additive matrix
Ha <- matrix(0, Ni, Nm.use)
for(i in 1:Nm.use){
  if(i%%100 == 0) cat(i, "\n")
  m <- Use_SNP[i]
  v <- genos[, m]
  paa <- sum(v == 2) / Ni
  pAa <- sum(v == 1) / Ni
  
  Ha[v == 0, i] <- (- pAa - 2 * paa) * (-1)
  Ha[v == 1, i] <- (1 - pAa - 2 * paa) * (-1)
  Ha[v == 2, i] <- (2 - pAa - 2 * paa) * (-1)
}
A <- Ha %*% t(Ha)
Deno.a <- sum(diag(A)) / Ni
A <- A / Deno.a
rownames(A) <- colnames(A) <- ped$id

#Dominance matrix
Hd<-matrix(0, Ni, Nm.use)
for(i in 1:Nm.use){
  m <- Use_SNP[i]
  v <- genos[, m]
  paa <- sum(v == 2) / Ni
  pAa <- sum(v == 1) / Ni
  pAA <- sum(v == 0) / Ni
  w <- pAA + paa - (pAA - paa)^2
  
  Hd[v == 0, i] <- (-1) * 2 * pAa * paa / w
  Hd[v == 1, i] <- 4 * pAA * paa / w
  Hd[v == 2, i] <- (-1) * 2 * pAA * pAa / w
}
D <- Hd %*% t(Hd)
Deno.d <- sum(diag(D)) / Ni
D <- D / Deno.d
rownames(D) <- colnames(D) <- ped$id

#Add 1e-3 to diagonals to avoid singular
diag(A) <- diag(A) + 1e-3
diag(D) <- diag(D) + 1e-3

#Calculate additive by additive (AA)
AA <- A * A
Deno.aa <- sum(diag(AA)) / Ni
AA <- AA / Deno.aa

#Calculate additive by dominance (AD)
AD <- A * D
Deno.ad <- sum(diag(AD)) / Ni
AD <- AD / Deno.ad

#Calculate dominance by dominance (DD)
DD <- D * D
Deno.dd <- sum(diag(DD)) / Ni
DD <- DD / Deno.dd

pheno_dat$ida <- pheno_dat$id
pheno_dat$idd <- pheno_dat$id
pheno_dat$idaa <- pheno_dat$id
pheno_dat$idad <- pheno_dat$id
pheno_dat$iddd <- pheno_dat$id


Result.A_D_AA <- mmer(y~1,
                    random = ~vsr(id, Gu = A) +
                      vsr(idd, Gu = D) +
                      vsr(idaa, Gu = AA),
                      # vsr(idad, Gu = AD) +
                      # vsr(iddd, Gu = DD),
                    rcov = ~ units,
                    data = pheno_dat)

Result.Full <- mmer(y~1,
                      random = ~vsr(id, Gu = A) +
                        vsr(idd, Gu = D) +
                        vsr(idaa, Gu = AA) +
                        vsr(idad, Gu = AD) +
                        vsr(iddd, Gu = DD),
                      rcov = ~ units,
                      data = pheno_dat)

#Extract results
Result.A_D_AA.Extract <- list(VarA = Result.A_D_AA$sigma[[1]],#Variance components
                            VarD = Result.A_D_AA$sigma[[2]],
                            VarAA = Result.A_D_AA$sigma[[3]],
                            # VarAD = Result.A_D_AA$sigma[[4]],
                            # VarDD = Result.A_D_AA$sigma[[5]],
                            VarR = Result.A_D_AA$sigma[[4]],
                            VarSE = Result.A_D_AA$sigmaSE,#SE of variance components
                            H2.A = vpredict(Result.A_D_AA, h2 ~ V1 / (V1+V2+V3+V4)),#Proportions
                            H2.D = vpredict(Result.A_D_AA, h2 ~ V2 / (V1+V2+V3+V4)),
                            H2.AA = vpredict(Result.A_D_AA, h2 ~ V3 / (V1+V2+V3+V4)),
                            # H2.AD = vpredict(Result.A_D_AA, h2 ~ V4 / (V1+V2+V3+V4+V5+V6)),
                            # H2.DD = vpredict(Result.A_D_AA, h2 ~ V5 / (V1+V2+V3+V4+V5+V6)),
                            AIC = Result.A_D_AA$AIC,
                            U.A = Result.A_D_AA$U[[1]][[1]],#Genotypic values 
                            U.D = Result.A_D_AA$U[[2]][[1]],
                            U.AA = Result.A_D_AA$U[[3]][[1]],
                            # U.AD = Result.A_D_AA$U[[4]][[1]],
                            # U.DD = Result.A_D_AA$U[[5]][[1]],
                            PEV.A = diag(Result.A_D_AA$PevU[[1]][[1]]),#PEVs of genotypic values
                            PEV.D = diag(Result.A_D_AA$PevU[[2]][[1]]),
                            PEV.AA = diag(Result.A_D_AA$PevU[[3]][[1]]),
                            # PEV.AD = diag(Result.A_D_AA$PevU[[4]][[1]]),
                            # PEV.DD = diag(Result.A_D_AA$PevU[[5]][[1]]),
                            Beta = Result.A_D_AA$Beta,#Fixed effects
                            VarBeta = Result.A_D_AA$VarBeta,
                            Monitor = Result.A_D_AA$monitor)

# Compare results to calculating GRMs from sommer functions
A2 <- A.mat(as.matrix(genos)-1)
D2 <- D.mat(as.matrix(genos)-1)
AA2 <- E.mat(as.matrix(genos)-1)
AD <- E.mat(as.matrix(genos)-1, type = "A#D")
DD <- E.mat(as.matrix(genos)-1, type = "D#D")

dimnames(A2) <- list(ped$id, ped$id)
dimnames(D2) <- list(ped$id, ped$id)
dimnames(AA2) <- list(ped$id, ped$id)
dimnames(AD) <- list(geno_names, geno_names)
dimnames(DD) <- list(geno_names, geno_names)


Result.A_D_AA2 <- mmer(y~1,
                      random = ~vsr(id, Gu = A2) +
                        vsr(idd, Gu = D2) +
                        vsr(idaa, Gu = AA2),
                      # vsr(idad, Gu = AD) +
                      # vsr(iddd, Gu = DD),
                      rcov = ~ units,
                      data = pheno_dat)

#Extract results
Result.A_D_AA2.Extract <- list(VarA = Result.A_D_AA2$sigma[[1]],#Variance components
                              VarD = Result.A_D_AA2$sigma[[2]],
                              VarAA = Result.A_D_AA2$sigma[[3]],
                              # VarAD = Result.A_D_AA2$sigma[[4]],
                              # VarDD = Result.A_D_AA2$sigma[[5]],
                              VarR = Result.A_D_AA2$sigma[[4]],
                              VarSE = Result.A_D_AA2$sigmaSE,#SE of variance components
                              H2.A = vpredict(Result.A_D_AA2, h2 ~ V1 / (V1+V2+V3+V4)),#Proportions
                              H2.D = vpredict(Result.A_D_AA2, h2 ~ V2 / (V1+V2+V3+V4)),
                              H2.AA = vpredict(Result.A_D_AA2, h2 ~ V3 / (V1+V2+V3+V4)),
                              # H2.AD = vpredict(Result.A_D_AA2, h2 ~ V4 / (V1+V2+V3+V4+V5+V6)),
                              # H2.DD = vpredict(Result.A_D_AA2, h2 ~ V5 / (V1+V2+V3+V4+V5+V6)),
                              AIC = Result.A_D_AA2$AIC,
                              U.A = Result.A_D_AA2$U[[1]][[1]],#Genotypic values 
                              U.D = Result.A_D_AA2$U[[2]][[1]],
                              U.AA = Result.A_D_AA2$U[[3]][[1]],
                              # U.AD = Result.A_D_AA2$U[[4]][[1]],
                              # U.DD = Result.A_D_AA2$U[[5]][[1]],
                              PEV.A = diag(Result.A_D_AA2$PevU[[1]][[1]]),#PEVs of genotypic values
                              PEV.D = diag(Result.A_D_AA2$PevU[[2]][[1]]),
                              PEV.AA = diag(Result.A_D_AA2$PevU[[3]][[1]]),
                              # PEV.AD = diag(Result.A_D_AA2$PevU[[4]][[1]]),
                              # PEV.DD = diag(Result.A_D_AA2$PevU[[5]][[1]]),
                              Beta = Result.A_D_AA2$Beta,#Fixed effects
                              VarBeta = Result.A_D_AA2$VarBeta,
                              Monitor = Result.A_D_AA2$monitor)


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


# Pedigree approach
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
