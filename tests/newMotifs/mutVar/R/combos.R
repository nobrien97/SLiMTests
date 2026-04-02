# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/mutVar/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./mutVarSR.sh"

models <- 450

seeds <- 48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:models, each = seeds),
                   seed = rep(1:seeds, times = models))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/mutVar/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)


cmds_adjTau <- data.frame(sr = "./mutVar_adjTauSR.sh",
                   model = rep(1:models, each = seeds),
                   seed = rep(1:seeds, times = models))

write.table(cmds_adjTau, "/mnt/c/GitHub/SLiMTests/tests/standingVar/mutVar/PBS/cmds_adjTau.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
