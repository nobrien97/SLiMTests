# Generate cmds.txt
singleRunBashName <- "./h2SR.sh"
seeds <- 1:1008
path <- "/mnt/c/GitHub/SLiMTests/tests/h2/R/h2_seeds.csv"
system(paste0("SeedGenerator -n 1008 -d ", path))
modelindex <- 0:1

dfCMDs <- expand.grid(seeds, modelindex)

dfCMDs <- cbind(rep(singleRunBashName, nrow(dfCMDs)), dfCMDs)

write.table(dfCMDs, "./cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
