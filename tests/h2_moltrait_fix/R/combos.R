# Generate cmds.txt
singleRunBashName <- "./h2_moltrait_fixSR.sh"
seeds <- 1:50
path <- "/mnt/c/GitHub/SLiMTests/tests/h2_moltrait_fix/R/h2_moltrait_fix_seeds.csv"
system(paste0("SeedGenerator -n 50 -t -d ", path))
molTrait <- 0:3
nloci <- c(8, 98, 998)
width <- c(0.05, 0.22)
locisigma <- c(0.1, 1)

dfCMDs <- expand.grid(molTrait, nloci, width, locisigma)
write.table(dfCMDs, "/mnt/c/GitHub/SLiMTests/tests/h2_moltrait_fix/R/combos.csv", row.names = F, col.names = F)

cmds <- data.frame(sr = rep(singleRunBashName, times = nrow(dfCMDs)),
                   seed = rep(seeds, times = nrow(dfCMDs)),
                   model = rep(1:nrow(dfCMDs), each = length(seeds)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/h2_moltrait_fix/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
