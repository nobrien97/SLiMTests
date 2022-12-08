# Generate cmds.txt
singleRunBashName <- "./h2_mol_fix_lvlsSR.sh"
molTrait <- 0:3
seeds <- 1:48

lvls <- c(0.5, 1, 1.5, 2)

dfCMDs <- expand.grid(molTrait, lvls)
write.table(dfCMDs, "/mnt/c/GitHub/SLiMTests/tests/h2_mol_fix_lvls/R/combos.csv", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)


cmds <- data.frame(sr = rep(singleRunBashName, times = nrow(dfCMDs)),
                   seed = rep(seeds, times = nrow(dfCMDs)),
                   model = rep(1:nrow(dfCMDs), each = length(seeds)))


write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/h2_mol_fix_lvls/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
