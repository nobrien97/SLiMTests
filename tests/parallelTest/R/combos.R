# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/parallelTest/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./parallelTestSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/parallelTest/R/parallelTest_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", seed_path))

seeds <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:4, each = length(seeds)),
                   seed = rep(1:length(seeds), times = length(1:4)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/parallelTest/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

singleRunBashName <- "./parallelTestAddSR.sh"
cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:4, each = length(seeds)),
                   seed = rep(1:length(seeds), times = length(1:4)))
write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/parallelTest/PBS/cmds_add.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
