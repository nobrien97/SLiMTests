library(sommer)
library(readr)
hapToGen <- function(haplos, ploidy = 2L) {
  # Converts a matrix of haplotypes to a genotype matrix
  
  # First make sure we have a proper ploidy for our haplotypes
  stopifnot(nrow(haplos) %% ploidy == 0)
  
  res <- matrix(rep(0, nrow(haplos)%/%ploidy))
  haplo_seq <- seq(from = 1, to = nrow(haplos), by = ploidy)
  res <- haplos[haplo_seq,]
  if (ploidy == 1) {
    return(res)
  }
  for (i in 1:(ploidy-1)) {
    res <- res + haplos[haplo_seq + i,]
  }
  row.names(res) <- NULL
  return(res)
}


haplos <- read_csv("slim_haplo.csv", col_names = F)
names(haplos) <- c("gen", "seed", "modelindex", paste0("q_", 1:(ncol(haplos) - 3)))
ped <- read_csv("slim_pedigree.csv")
phenos <- read_csv("slim_sampled_pheno.csv", col_names = F)
names(phenos)[1:3] <- c("gen", "seed", "modelindex")
phenos$seed <- as.character(phenos$seed)

# For some reason the header is included in the unique() return? Need to remove that
reps <- unique(ped$seed)
reps <- reps[!grepl("^s", reps)]
models <- unique(ped$modelindex)
models <- as.integer(models[!grepl("^m", models)])
gens <- unique(ped$gen)
gens <- as.integer(gens[!grepl("^g", gens)])
nReps <- length(reps)
nModels <- length(models)
nOut <- nReps * nModels * length(gens)
# Set up output vector for the loop
result <- vector(mode = "list", length = nOut)

curResult <- 1

for (i in seq_len(nReps)) {
  for (j in seq_len(nModels)) {
    phenos_smpl <- phenos[phenos$seed == reps[i] & phenos$modelindex == models[j],]
    ped_smpl <- ped[ped$seed == reps[i] & ped$modelindex == models[j],]
    
    for (k in seq_len(nrow(phenos_smpl))) {
      curPhenos <- unlist(phenos_smpl[k, -(1:3)])
      names(curPhenos) <- NULL
      
      if ( var(curPhenos) < 1e-5 ) {
        sprintf("Inadequate variation in phenotypes to predict heritability: Seed = %s, Model = %s, Gen = %i",
                as.character(reps[i]),
                ifelse(models[j], "Network", "Additive"),
                as.integer(gens[k]))
        result[[curResult]] <- rep(NA, 6)
        curResult <- curResult + 1
        next
      }
      curPed <- ped_smpl[ped_smpl$gen == gens[k],(4:6)]
      
      pheno_dat <- data.frame(y = curPhenos, 
                              GCA1 = as.factor(curPed$mother),
                              GCA2 = as.factor(curPed$father),
                              SCA = as.factor(curPed$id))
      
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
      result[[curResult]] <- c(h2, H2, Va, Vd, Vg, Vi)
      curResult <- curResult + 1
    }
    
  }
}