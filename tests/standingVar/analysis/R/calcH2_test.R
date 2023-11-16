# R Script to calculate heritability given some data
library(sommer)

# Get command line arguments
args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])
chunk <- as.numeric(args[2])
run_chunk <- paste(run, chunk, sep = "_")

# Load function for extracting genotypes
source("helpFns.R")

# Read data and add names
haplos <- scan(paste0("test_haplo.csv"), what = integer(), sep = ",")[-(1:3)]
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
pheno_dat$ide <- pheno_dat$id

######################################################################
# This code is adapted from Onigo et al. 2021, thanks to Akio Onogi
# for sharing code samples

# Fit models using A, D, and AA GRMs only to avoid overfitting
Result.A_D_AA <- mmer(y~1,
                    random = ~vsr(id, Gu = A) +
                      vsr(idd, Gu = D) +
                      vsr(ide, Gu = AA),
                    rcov = ~ units,
                    data = pheno_dat)

# Make sure we have some results: try to run a simpler model if we don't
if(!exists("Result.A_D_AA")) {
    print(paste("Unable to generate mmer in model", run, "- trying VA-VD model"))

Result.A_D_AA <- mmer(y~1,
                    random = ~vsr(id, Gu = A) +
                      vsr(idd, Gu = D),
                    rcov = ~ units,
                    data = pheno_dat)

}
# Check again
if(!exists("Result.A_D_AA")) {
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
                            H2.A = predict(Result.A_D_AA, h2 ~ V1 / (V1+V2+V3+V4)),#Proportions
                            H2.D = predict(Result.A_D_AA, h2 ~ V2 / (V1+V2+V3+V4)),
                            H2.AA = predict(Result.A_D_AA, h2 ~ V3 / (V1+V2+V3+V4)),
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
molPheno_dat$idd <- molPheno_dat$id
molPheno_dat$ide <- molPheno_dat$id

Ai <- as(solve(A + diag(1e-4,ncol(A),ncol(A))), Class="dgCMatrix")
Di <- as(solve(D + diag(1e-4,ncol(D),ncol(D))), Class="dgCMatrix")
AAi <-as(solve(A + diag(1e-4,ncol(A),ncol(A))), Class="dgCMatrix")

Result.A_D_AA <- mmec(cbind(aZ, bZ)~1,
                      random=~ vsc(isc(id), Gu=Ai, Gtc=unsm(2)),
                      rcov=~ vsc(units, Gtc=unsm(2)), 
                      data = molPheno_dat)

cov2cor(Result.A_D_AA$sigma$`u:id`)


rownames(Result.A_D_AA.Extract) <- NULL

# Output file
print(paste0("Writing output for model ", run, "..."))
write.table(Result.A_D_AA.Extract, paste0("/scratch/ht96/nb9894/h2_hsfs_nloci/getH2_hsfs_nloci/out_h2_", run_chunk, ".csv"), sep = ",", row.names = F, col.names = F)
