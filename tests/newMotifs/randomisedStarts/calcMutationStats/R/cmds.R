# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./calcMutationStatsSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:15)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
