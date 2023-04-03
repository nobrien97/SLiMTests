# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/equalK_sweep/equalK_ODE/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./sweepODECalcSR.sh"

# 5233 total adapted models
cmds <- data.frame(sr = singleRunBashName,
                   run = 1:5233)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/equalK_sweep/equalK_ODE/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
