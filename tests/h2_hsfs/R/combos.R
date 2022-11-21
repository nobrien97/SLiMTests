# Generate cmds.txt
singleRunBashName <- "./h2_hsfsSR.sh"
seeds <- 1:720
path <- "/mnt/c/GitHub/SLiMTests/tests/h2_hsfs/R/h2_hsfs_seeds.csv"
system(paste0("SeedGenerator -n 720 -t -d ", path))
modelindex <- 0:1

dfCMDs <- expand.grid(seeds, modelindex)

dfCMDs <- cbind(rep(singleRunBashName, nrow(dfCMDs)), dfCMDs)

write.table(dfCMDs, "./cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
