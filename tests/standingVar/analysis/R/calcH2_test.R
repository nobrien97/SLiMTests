# R Script to calculate heritability given some data
library(sommer)
library(MCMCglmm)

# Load function for extracting genotypes
source("helpFns.R")

# Read data and add names
haplos <- scan(paste0("test_haplo.csv"), what = numeric(), sep = ",")[-(1:3)]
relPos <- scan(paste0("test_rel_relPos.csv"), what = numeric(), sep = ",")[-(1:3)]
haplos <- decompressHap(haplos, 2000, 1024)
names(haplos) <- paste0("q_", seq_len(ncol(haplos)))
phenos <- scan(paste0("test_pheno.csv"), sep = ",")
names(phenos) <- NULL
# From the phenos file, extract gen, seed, modelindex and remove them
## Use these later to identify output
run_info <- phenos[1:3]
phenos <- phenos[-(1:3)]


# Get genotypes from haplos
genos <- hapToGen(haplos)
genos <- genos - 1
geno_names <- paste0("i_", seq_len(nrow(genos)))
pheno_dat <- data.frame(y = phenos,
                        id = as.factor(geno_names))

# Calculate genomic relationship matrices
A <- A.mat(genos)
# D <- D.mat(genos)
# AA <- A*A

# Add names to dims
dimnames(A) <- list(geno_names, geno_names)
# dimnames(D) <- list(geno_names, geno_names)
# dimnames(AA) <- list(geno_names, geno_names)

# Add extra IDs for D and AA
# pheno_dat$idd <- pheno_dat$id
# pheno_dat$ide <- pheno_dat$id

######################################################################
# This code is adapted from Onigo et al. 2021, thanks to Akio Onogi
# for sharing code samples

# Fit animal model: additive only
Result.A <- mmer(y~1,
                    random = ~vsr(id, Gu = A),
                    rcov = ~ units,
                    data = pheno_dat)

# Make sure we have some results: try to run a simpler model if we don't
if(!exists("Result.A_D_AA")) {
    print(paste("Unable to generate mmer in model", run, "- trying VA-VD model"))

  Result.A_D_AA <- mmer(y~1,
                    random = ~vsr(id, Gu = A),
                    rcov = ~ units,
                    data = pheno_dat)
}
# Check again
if (!exists("Result.A_D_AA")) {
    print(paste("Unable to generate mmer in model", run, "closing"))
    q(save = "no")
}

#Extract results
Result.A_D_AA.Extract <- data.frame(gen = run_info[1],
                            seed = as.character(run_info[2]),
                            modelindex = run_info[3],
                            VarA = unname(Result.A_D_AA$sigma[[1]]), #Variance components
                            VarD = unname(Result.A_D_AA$sigma[[2]]),
                            VarAA = unname(Result.A_D_AA$sigma[[3]]),
                            VarR = unname(Result.A_D_AA$sigma[[4]]),
                            VarA.SE = sqrt(diag(Result.A_D_AA$sigmaSE))[1],#SE of variance components
                            VarD.SE = sqrt(diag(Result.A_D_AA$sigmaSE))[2],
                            VarAA.SE = sqrt(diag(Result.A_D_AA$sigmaSE))[3],
                            VarR.SE = sqrt(diag(Result.A_D_AA$sigmaSE))[4],
                            H2.A = vpredict(Result.A_D_AA, h2 ~ V1 / (V1+V2+V3+V4)),#Proportions
                            H2.D = vpredict(Result.A_D_AA, h2 ~ V2 / (V1+V2+V3+V4)),
                            H2.AA = vpredict(Result.A_D_AA, h2 ~ V3 / (V1+V2+V3+V4)),
                            AIC = Result.A_D_AA$AIC
                            )
######################################################################

# Measure molecular components
molPhenos <- scan(paste0("test_moltrait.csv"), sep = ",")
names(molPhenos) <- NULL
# From the phenos file, extract gen, seed, modelindex and remove them
## Use these later to identify output
run_info <- molPhenos[1:3]
molPhenos <- molPhenos[-(1:3)]

# Convert to data frame
molPheno_dat <- data.frame(Z   = molPhenos[seq(1, length(molPhenos), by = 5)],
                           aZ  = molPhenos[seq(2, length(molPhenos), by = 5)],
                           bZ  = molPhenos[seq(3, length(molPhenos), by = 5)],
                           KZ  = molPhenos[seq(4, length(molPhenos), by = 5)],
                           KXZ = molPhenos[seq(5, length(molPhenos), by = 5)],
                           id = as.factor(geno_names))

# Add extra IDs for D and AA
# molPheno_dat$idd <- molPheno_dat$id
# molPheno_dat$ide <- molPheno_dat$id

# Only measure additive variance in this one: since we are interested
# only in G, we can measure epistasis on the phenotype separately
# Can't fit all four at once: takes too long
# Singh et al. 2019 (https://doi.org/10.3389/fpls.2019.00394) 
# estimated covariances by fitting multivariate models on two traits at a time

for (i in 1:nrow(covCombos)) {
  molResult.A <- mmer(cbind(aZ, bZ)~1,
                        random=~ vsr(id, Gu=A, Gtc=unsm(2)),
                        rcov=~ units, 
                        data = molPheno_dat)
}

summary(molResult.A)$varcomp

# Get genetic correlation, r_g
covaZbZ <- molResult.A$sigma$`u:id`[upper.tri(molResult.A$sigma$`u:id`)]

rg_aZbZ <- covaZbZ / (sqrt(prod(diag(molResult.A$sigma$`u:id`))))

molResult.A <- mmer(cbind(aZ, bZ)~1,
                    random=~ vsr(id, Gu=A, Gtc=unsm(2)),
                    rcov=~ units, 
                    data = molPheno_dat)

summary(molResult.A)
molResult.A$sigma$`u:id`


rownames(Result.A_D_AA.Extract) <- NULL


# Inverse additive relatedness matrix
A <- A.mat(genos, 0.05)
Ai <- as(solve(A + diag(1e-4,ncol(A),ncol(A))), Class="dgCMatrix")
dimnames(Ai) <- list(geno_names, geno_names)

# prior
# https://devillemereuil.legtux.org/wp-content/uploads/2021/09/tuto_en.pdf
# https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12886
# scale variance by phenotypic variance/number loci
# prior_animal <- list(R = list(V = diag(4), nu = 4),
#                      G = list(G1 = list(V = diag(4), nu = 4)))

# Scale variances
molPheno_dat$aZ_scl <- scale(molPheno_dat$aZ)
molPheno_dat$bZ_scl <- scale(molPheno_dat$bZ)
molPheno_dat$KZ_scl <- scale(molPheno_dat$KZ)
molPheno_dat$KXZ_scl <- scale(molPheno_dat$KXZ)

# # too slow
# prior_test <- list(R = list(V = 1, nu = 0.002),
#                    G = list(G1 = list(V = 1, nu = 0.002)))
# mcmctest <- MCMCglmm(aZ_scl ~ 1,
#                      random = ~ id,
#                      rcov = ~ units,
#                      family = "gaussian",
#                      ginv = list(id = Ai),
#                      data = molPheno_dat, 
#                      prior = prior_test
#                      )

mcmcres <- MCMCglmm(cbind(aZ_scl, bZ_scl, KZ_scl, KXZ_scl) ~ 1,
         random = ~ id,
         rcov = ~ units,
         family = rep("gaussian", times = 4),
         ginv = list(id = Ai),
         data = molPheno_dat, prior = prior_animal)


# Output file
print(paste0("Writing output for model ", run, "..."))
write.table(Result.A_D_AA.Extract, paste0("/scratch/ht96/nb9894/h2_hsfs_nloci/getH2_hsfs_nloci/out_h2_", run_chunk, ".csv"), sep = ",", row.names = F, col.names = F)
