# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/epistasisDensity/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./epistasisDensitySR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:15)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/epistasisDensity/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
