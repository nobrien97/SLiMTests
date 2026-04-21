# Create cmds.txt
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/mutVar/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./mutVarSR.sh"

models <- 15

seeds <- 208

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:models, each = seeds),
                   seed = rep(1:seeds, times = models))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/mutVar/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
