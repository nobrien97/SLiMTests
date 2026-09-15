# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/rhvae/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./rhvaeSR.sh"

# Run for each of 5 models
models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")

cmds <- data.frame(sr = singleRunBashName,
                   run = models)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/rhvae/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
