library(sommer)
source("helpFns.R")

haplos_net <- read.csv("~/Desktop/slim_haplo2930242583378526208_1.csv", header = F)[,-(1:3)]
names(haplos_net) <- paste0("q_", seq_len(ncol(haplos_net)))
ped_net <- read.csv("~/Desktop/slim_pedigree2930242583378526208_1.csv", header = F)[,-(1:3)]
names(ped_net) <- c("id", "mother", "father")

# Organise data
genos_net <- hapToGen(haplos_net) 
phenos_net <- unlist(read.csv("~/Desktop/slim_sampled_pheno2930242583378526208_1.csv", header=F))
names(phenos_net) <- NULL
run_info_net <- phenos_net[1:3]
phenos_net <- phenos_net[-(1:3)]

geno_names <- paste0("i_", 1:nrow(genos_net))

pheno_dat_net <- data.frame(y = phenos_net,
                        id = as.factor(geno_names))

A <- A.mat(as.matrix(genos_net)-1)
D <- D.mat(as.matrix(genos_net)-1)
AA <- A*A

dimnames(A) <- list(geno_names, geno_names)
dimnames(D) <- list(geno_names, geno_names)
dimnames(AA) <- list(geno_names, geno_names)

pheno_dat_net$idd <- pheno_dat_net$id
pheno_dat_net$idaa <- pheno_dat_net$id


Result_net.A_D_AA <- mmer(y~1,
                       random = ~vsr(id, Gu = A) +
                         vsr(idd, Gu = D) +
                         vsr(idaa, Gu = AA),
                       rcov = ~ units,
                       data = pheno_dat_net)

Result_net.A_D_AA.Extract <- data.frame(gen = run_info_net[1],
                                    seed = as.character(run_info_net[2]),
                                    modelindex = run_info_net[3],
                                    VarA = unname(Result_net.A_D_AA$sigma[[1]]), #Variance components
                                    VarD = unname(Result_net.A_D_AA$sigma[[2]]),
                                    VarAA = unname(Result_net.A_D_AA$sigma[[3]]),
                                    VarR = unname(Result_net.A_D_AA$sigma[[4]]),
                                    VarA.SE = sqrt(diag(Result_net.A_D_AA$sigmaSE))[1],#SE of variance components
                                    VarD.SE = sqrt(diag(Result_net.A_D_AA$sigmaSE))[2],
                                    VarAA.SE = sqrt(diag(Result_net.A_D_AA$sigmaSE))[3],
                                    VarR.SE = sqrt(diag(Result_net.A_D_AA$sigmaSE))[4],
                                    H2.A = vpredict(Result_net.A_D_AA, h2 ~ V1 / (V1+V2+V3+V4)),#Proportions
                                    H2.D = vpredict(Result_net.A_D_AA, h2 ~ V2 / (V1+V2+V3+V4)),
                                    H2.AA = vpredict(Result_net.A_D_AA, h2 ~ V3 / (V1+V2+V3+V4)),
                                    AIC = Result_net.A_D_AA$AIC
)
Result_net.A_D_AA.Extract

# Additive
haplos_add <- read.csv("~/Desktop/slim_haplo2873761167128395776_0.csv", header = F)[,-(1:3)]
names(haplos_add) <- paste0("q_", seq_len(ncol(haplos_add)))
ped_add <- read.csv("~/Desktop/slim_pedigree2873761167128395776_0.csv", header = F)[,-(1:3)]
names(ped_add) <- c("id", "mother", "father")

# Organise data
genos_add <- hapToGen(haplos_add) 
phenos_add <- unlist(read.csv("~/Desktop/slim_sampled_pheno2873761167128395776_0.csv", header=F))
names(phenos_add) <- NULL
run_info_add <- phenos_add[1:3]
phenos_add <- phenos_add[-(1:3)]

pheno_dat_add <- data.frame(y = phenos_add,
                            id = as.factor(geno_names))

A <- A.mat(as.matrix(genos_add)-1)
D <- D.mat(as.matrix(genos_add)-1)
AA <- A*A

dimnames(A) <- list(geno_names, geno_names)
dimnames(D) <- list(geno_names, geno_names)
dimnames(AA) <- list(geno_names, geno_names)

pheno_dat_add$idd <- pheno_dat_add$id
pheno_dat_add$idaa <- pheno_dat_add$id


Result_add.A_D_AA <- mmer(y~1,
                          random = ~vsr(id, Gu = A) +
                            vsr(idd, Gu = D) +
                            vsr(idaa, Gu = AA),
                          rcov = ~ units,
                          data = pheno_dat_add)

Result_add.A_D_AA.Extract <- data.frame(gen = run_info_add[1],
                                        seed = as.character(run_info_add[2]),
                                        modelindex = run_info_add[3],
                                        VarA = unname(Result_add.A_D_AA$sigma[[1]]), #Variance components
                                        VarD = unname(Result_add.A_D_AA$sigma[[2]]),
                                        VarAA = unname(Result_add.A_D_AA$sigma[[3]]),
                                        VarR = unname(Result_add.A_D_AA$sigma[[4]]),
                                        VarA.SE = sqrt(diag(Result_add.A_D_AA$sigmaSE))[1],#SE of variance components
                                        VarD.SE = sqrt(diag(Result_add.A_D_AA$sigmaSE))[2],
                                        VarAA.SE = sqrt(diag(Result_add.A_D_AA$sigmaSE))[3],
                                        VarR.SE = sqrt(diag(Result_add.A_D_AA$sigmaSE))[4],
                                        H2.A = vpredict(Result_add.A_D_AA, h2 ~ V1 / (V1+V2+V3+V4)),#Proportions
                                        H2.D = vpredict(Result_add.A_D_AA, h2 ~ V2 / (V1+V2+V3+V4)),
                                        H2.AA = vpredict(Result_add.A_D_AA, h2 ~ V3 / (V1+V2+V3+V4)),
                                        AIC = Result_add.A_D_AA$AIC
)
Result_add.A_D_AA.Extract

# Pedigree approach
pheno_dat <- data.frame(y = phenos_net, 
                        GCA1 = as.factor(ped_net$mother),
                        GCA2 = as.factor(ped_net$father),
                        SCA = as.factor(ped_net$id))

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


pheno_dat <- data.frame(y = phenos_add, 
                        GCA1 = as.factor(ped_add$mother),
                        GCA2 = as.factor(ped_add$father),
                        SCA = as.factor(ped_add$id))

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
