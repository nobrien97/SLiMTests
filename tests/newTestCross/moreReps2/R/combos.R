# Generate cmds.txt
singleRunBashName <- "./newTestCrossSR.sh"
seeds <- 1:240
path <- "/mnt/c/GitHub/SLiMTests/tests/newTestCross/moreReps2/R/newTestCross_seeds.csv"
system(paste0("SeedGenerator -n 240 -t -d ", path))
nloci <- c(8, 98, 998)
locisigma <- c(0.1, 1)

dfCMDs <- expand.grid(nloci, locisigma)
write.table(dfCMDs, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/moreReps2/R/combos.csv", row.names = F, col.names = F)

cmds <- data.frame(sr = rep(singleRunBashName, times = nrow(dfCMDs)),
                   seed = rep(seeds, times = nrow(dfCMDs)),
                   model = rep(1:nrow(dfCMDs), each = length(seeds)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/moreReps2/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
