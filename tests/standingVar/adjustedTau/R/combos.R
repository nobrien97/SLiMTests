# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/adjustedTau/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./adjustedTauSR.sh"

models <- 450

seeds <- 48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:models, each = seeds),
                   seed = rep(1:seeds, times = models))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/adjustedTau/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
