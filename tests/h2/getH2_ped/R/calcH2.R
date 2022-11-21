# R Script to calculate heritability given some data
library(sommer)

# Get command line arguments
args <- commandArgs(trailingOnly = T)
run <- as.numeric(args[1])
chunk <- as.numeric(args[2])
run_chunk <- paste(run, chunk, sep = "_")

# Read data and add names
ped <- read.csv(paste0("slim_ped_sbst_", run, ".csv"), header = F)[, -(1:3)]
names(ped) <- c("id", "mother", "father")
phenos <- unlist(read.csv(paste0("slim_pheno_sbst_", run, ".csv"), header = F))
names(phenos) <- NULL
# From the phenos file, extract gen, seed, modelindex and remove them
## Use these later to identify output
run_info <- phenos[1:3]
phenos <- phenos[-(1:3)]

# Check to make sure there is variability
if (!var(phenos)) {
    sprintf("No variance in run %i, closing R...", run)
    q(save = "no")
}

pheno_dat <- data.frame(y = phenos, 
                        GCA1 = as.factor(ped$mother),
                        GCA2 = as.factor(ped$father),
                        SCA = as.factor(ped$id))

ans.ADE <- mmec(y~1,
                random =~GCA1+GCA2+SCA,
                rcov=~units,
                data=pheno_dat)

# Check that the mmec finished properly
if (!exists("ans.ADE")) {
  q(save = "no")
}

if (length(ans.ADE) != 26) {
  q(save = "no")
}

suma <- summary(ans.ADE)$varcomp
Vgca <- sum(suma[1:2,1])
Vsca <- suma[3,1]
Ve <- suma[4,1]
Va <- 4*Vgca
Vd <- 4*Vsca
Vg <- Va + Vd
H2 <- Vg / (Vg + (Ve))
h2 <- Va / (Vg + (Ve))

#Extract results
result <- data.frame(       gen = run_info[1],
                            seed = as.character(run_info[2]),
                            modelindex = run_info[3],
                            Vgca = Vgca,
                            Vsca = Vsca,
                            Ve = Ve,
                            Va = Va,
                            Vd = Vd,
                            Vg = Vg,
                            H2 = H2,
                            h2 = h2
                            )
######################################################################

rownames(result) <- NULL

# Output file
print(paste("Writing output for run_chunk:", run_chunk))
write.table(result, paste0("/scratch/ht96/nb9894/h2/getH2_ped/out_h2_", run_chunk, ".csv"), sep = ",", row.names = F, col.names = F)