# R Script to calculate heritability given some data
library(sommer)

# Get command line arguments
args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])
chunk <- as.numeric(args[2])
run_chunk <- paste(run, chunk, sep = "_")

# Load function for extracting genotypes
source("~/tests/h2/R/helpFns.R")

# Read data and add names
haplos <- read.csv(paste0("slim_haplo_sbst_", run, ".csv"), header = F)[, -(1:3)]
names(haplos) <- paste0("q_", seq_len(ncol(haplos)))
phenos <- unlist(read.csv(paste0("slim_pheno_sbst_", run, ".csv"), header = F))
names(phenos) <- NULL
# From the phenos file, extract gen, seed, modelindex and remove them
## Use these later to identify output
run_info <- phenos[1:3]
phenos <- phenos[-(1:3)]

# Check to make sure there is variability
if (!max(haplos)) {
   result <- data.frame(gen = run_info[1],
                            seed = as.character(run_info[2]),
                            modelindex = run_info[3],
                            VarA = 0, 
                            VarD = 0,
                            VarAA = 0,
                            VarR = 0,
                            VarA.SE = 0,
                            VarD.SE = 0,
                            VarAA.SE = 0,
                            VarR.SE = 0,
                            H2.A.Estimate = 0,
                            H2.A.SE = 0,
                            H2.D.Estimate = 0,
                            H2.D.SE = 0,
                            H2.AA.Estimate = 0,
                            H2.AA.SE = 0,
                            AIC = 0
                            )
    write.table(result, paste0("/scratch/ht96/nb9894/h2/getH2/out_h2_", run_chunk, ".csv"), sep = ",", row.names = F, col.names = F)
    q(save = "no")
}

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

# Make sure we have some results
if(!exists("Result.A_D_AA")) {
    result <- data.frame(gen = run_info[1],
                            seed = as.character(run_info[2]),
                            modelindex = run_info[3],
                            VarA = 0, 
                            VarD = 0,
                            VarAA = 0,
                            VarR = 0,
                            VarA.SE = 0,
                            VarD.SE = 0,
                            VarAA.SE = 0,
                            VarR.SE = 0,
                            H2.A.Estimate = 0,
                            H2.A.SE = 0,
                            H2.D.Estimate = 0,
                            H2.D.SE = 0,
                            H2.AA.Estimate = 0,
                            H2.AA.SE = 0,
                            AIC = 0
                            )
    write.table(result, paste0("/scratch/ht96/nb9894/h2/getH2/out_h2_", run_chunk, ".csv"), sep = ",", row.names = F, col.names = F)
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

rownames(Result.A_D_AA.Extract) <- NULL

# Output file
write.table(Result.A_D_AA.Extract, paste0("/scratch/ht96/nb9894/h2/getH2/out_h2_", run_chunk, ".csv"), sep = ",", row.names = F, col.names = F)