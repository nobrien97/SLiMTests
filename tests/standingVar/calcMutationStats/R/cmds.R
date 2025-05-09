# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./calcNewFXSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = 1:450)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/PBS/cmds_new_fx.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
