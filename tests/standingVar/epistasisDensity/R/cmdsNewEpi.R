# Create cmds.txt and combo file
path <- "/mnt/e/Documents/GitHub/SLiMTests/tests/standingVar/epistasisDensity/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./epistasisDensitySRNewEpi.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:450)

write.table(cmds, "/mnt/e/Documents/GitHub/SLiMTests/tests/standingVar/epistasisDensity/PBS/cmdsNewEpi.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
