# Create cmds.txt and combo file
path <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcLD/R/"
setwd(path)

# Get number of rows in sharedmutfreqs
# wc -l slim_sharedmutfreqs.csv
# 37956

# Generate cmds.txt
singleRunBashName <- "./calcLDSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 0:37955)

write.table(cmds, "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcLD/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

