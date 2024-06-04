# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/sanityChecks/calcMutationStats/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./calcMutationStatsSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:36)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/sanityChecks/calcMutationStats/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
