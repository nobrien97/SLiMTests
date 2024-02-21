# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./epistasisDensitySR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:450)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
