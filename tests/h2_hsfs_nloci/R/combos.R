# Generate cmds.txt
singleRunBashName <- "./h2_hsfs_nlociSR.sh"
seeds <- 1:48
path <- "/mnt/c/GitHub/SLiMTests/tests/h2_hsfs_nloci/R/h2_hsfs_nloci_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", path))
model <- 0:1
nloci <- c(10, 100, 1000)
width <- c(0.05, 0.105, 0.22)
locisigma <- c(0.01, 0.1, 1)

dfCMDs <- expand.grid(model, nloci, width, locisigma)
write.table(dfCMDs, "/mnt/c/GitHub/SLiMTests/tests/h2_hsfs_nloci/R/combos.csv", row.names = F, col.names = F)

cmds <- data.frame(sr = rep(singleRunBashName, times = nrow(dfCMDs)),
                   seed = rep(seeds, times = nrow(dfCMDs)),
                   model = rep(1:nrow(dfCMDs), each = length(seeds)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/h2_hsfs_nloci/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
