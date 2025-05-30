# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/epiPerMolComp/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./epiPerMolCompSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:450)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/epistasisDensity/epiPerMolComp/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
