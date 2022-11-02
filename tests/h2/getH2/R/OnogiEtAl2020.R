#Create the data frame of inputs#####################################
##Pheno is a data frame including IDs of animals, several fixed effects, and phenotypic values
Data <- data.frame(ID1forA = as.factor(Pheno$ID1),
                   ID1forD = as.factor(Pheno$ID1),
                   ID1forAA = as.factor(Pheno$ID1),
                   ID1forAD = as.factor(Pheno$ID1),
                   ID1forDD = as.factor(Pheno$ID1),
                   Sex = as.factor(Pheno$SexClass),#fixed effect
                   Date = as.factor(Pheno$SlaughterDateClass),#fixed effect
                   Farm = as.factor(Pheno$FarmClass),#fixed effect
                   Age = scale(Pheno$AgeDev),#fixed effect (covariate)
                   Hetero = scale(Hetero),#fixed effect (covariate)
                   CW = Pheno$CW,#Phenotypic value
                   REA = Pheno$REA,#Phenotypic value
                   RT = Pheno$RT,#Phenotypic value
                   SFT = Pheno$SFT,#Phenotypic value
                   YI = Pheno$YI,#Phenotypic value
                   BMS = Pheno$BMS)#Phenotypic value

#Calculate the genomic relationship matrices##########################
##Geno is a SNP genotype matrix whose entries are coded as 0, 1, and 2
##Geno is a Ni (number of individuals) x Nm (number of markers) matrix
Ni <- nrow(Geno)
Nm <- ncol(Geno)

Af <- colSums(Geno) / (2 * Ni)
Maf <- Af
Maf[Maf > 0.5] <- 1 - Maf[Maf > 0.5]
Use_SNP <- which(Maf >= 0.05)
Nm.use <- length(Use_SNP)

#Additive matrix
Ha <- matrix(0, Ni, Nm.use)
for(i in 1:Nm.use){
  if(i%%100 == 0) cat(i, "\n")
  m <- Use_SNP[i]
  v <- Geno[, m]
  paa <- sum(v == 2) / Ni
  pAa <- sum(v == 1) / Ni
  
  Ha[v == 0, i] <- (- pAa - 2 * paa) * (-1)
  Ha[v == 1, i] <- (1 - pAa - 2 * paa) * (-1)
  Ha[v == 2, i] <- (2 - pAa - 2 * paa) * (-1)
}
A <- Ha %*% t(Ha)
Deno.a <- sum(diag(A)) / Ni
A <- A / Deno.a
rownames(A) <- colnames(A) <- Pheno$ID1

#Dominance matrix
Hd<-matrix(0, Ni, Nm.use)
for(i in 1:Nm.use){
  m <- Use_SNP[i]
  v <- Geno[, m]
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
rownames(D) <- colnames(D) <- Pheno$ID1

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


#Fit models with all variance components###############################
library(sommer)
#=>I used version 4.09 or 4.12.

#Fit to carcass weight (CW)
Result.Full <- mmer(CW ~ Sex + Date + Farm + Age + Hetero,
                    random = ~vs(ID1forA, Gu = A) +
                      vs(ID1forD, Gu = D) +
                      vs(ID1forAA, Gu = AA) +
                      vs(ID1forAD, Gu = AD) +
                      vs(ID1forDD, Gu = DD),
                    rcov = ~ units,
                    data = Data)

#Extract results
Result.Full.Extract <- list(VarA = Result.Full$sigma[[1]],#Variance components
                            VarD = Result.Full$sigma[[2]],
                            VarAA = Result.Full$sigma[[3]],
                            VarAD = Result.Full$sigma[[4]],
                            VarDD = Result.Full$sigma[[5]],
                            VarR = Result.Full$sigma[[6]],
                            VarSE = Result.Full$sigmaSE,#SE of variance components
                            H2.A = vpredict(Result.Full, h2 ~ V1 / (V1+V2+V3+V4+V5+V6)),#Proportions
                            H2.D = vpredict(Result.Full, h2 ~ V2 / (V1+V2+V3+V4+V5+V6)),
                            H2.AA = vpredict(Result.Full, h2 ~ V3 / (V1+V2+V3+V4+V5+V6)),
                            H2.AD = vpredict(Result.Full, h2 ~ V4 / (V1+V2+V3+V4+V5+V6)),
                            H2.DD = vpredict(Result.Full, h2 ~ V5 / (V1+V2+V3+V4+V5+V6)),
                            AIC = Result.Full$AIC,
                            U.A = Result.Full$U[[1]][[1]],#Genotypic values 
                            U.D = Result.Full$U[[2]][[1]],
                            U.AA = Result.Full$U[[3]][[1]],
                            U.AD = Result.Full$U[[4]][[1]],
                            U.DD = Result.Full$U[[5]][[1]],
                            PEV.A = diag(Result.Full$PevU[[1]][[1]]),#PEVs of genotypic values
                            PEV.D = diag(Result.Full$PevU[[2]][[1]]),
                            PEV.AA = diag(Result.Full$PevU[[3]][[1]]),
                            PEV.AD = diag(Result.Full$PevU[[4]][[1]]),
                            PEV.DD = diag(Result.Full$PevU[[5]][[1]]),
                            Beta = Result.Full$Beta,#Fixed effects
                            VarBeta = Result.Full$VarBeta,
                            Monitor = Result.Full$monitor)


