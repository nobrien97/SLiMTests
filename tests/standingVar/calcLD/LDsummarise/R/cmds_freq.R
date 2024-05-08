# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/LDsummarise/R/"
setwd(path)

# There are 65249 samples to go through, since the 
# model is given to skip argument, it is 0:(65249-1)

# Generate cmds.txt
singleRunBashName <- "./LDsummariseFreqSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 0:65248)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/LDsummarise/PBS/cmds_freq.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
