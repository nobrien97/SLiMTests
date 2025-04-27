# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./epistasisDensityNewEpiSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:450)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/PBS/cmdsNewEpi.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
