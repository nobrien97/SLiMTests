# Generate cmds.txt
seeds <- 1:48
path <- "/mnt/c/GitHub/SLiMTests/tests/newTestCross/additive/R/newTestCrossfix_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", path))
nloci <- c(10, 100, 1000)
locisigma <- c(0.1, 1)

dfCMDs <- expand.grid(nloci, locisigma)
write.table(dfCMDs, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/additive/R/combos.csv", row.names = F, col.names = F)

cmds <- data.frame(sr = rep(singleRunBashName, times = nrow(dfCMDs)),
                   seed = rep(seeds, times = nrow(dfCMDs)),
                   model = rep(1:nrow(dfCMDs), each = length(seeds)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/additive/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
