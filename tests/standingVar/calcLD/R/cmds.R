# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/R/"
setwd(path)

# There are 36181 samples to go through, since the 
# model is given to skip argument, it is 0:(36181-1)

# Generate cmds.txt
singleRunBashName <- "./calcLDRSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 0:36180)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/PBS/cmdsR.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

singleRunBashName <- "./calcLD_PRSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 0:36180)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/PBS/cmdsR.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
