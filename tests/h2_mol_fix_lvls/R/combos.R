# Generate cmds.txt
singleRunBashName <- "./h2_mol_fix_lvlsSR.sh"
molTrait <- 0:3
seeds <- 1:48

lvls <- c(0.5, 1, 1.5, 2)

dfCMDs <- expand.grid(seeds, molTrait)
dfCMDs <- cbind(rep(singleRunBashName, nrow(dfCMDs)), dfCMDs)

write.table(dfCMDs, "/mnt/c/GitHub/SLiMTests/tests/h2_mol_fix_lvls/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
