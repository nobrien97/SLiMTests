# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/multODETest/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./multODETestSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/multODETest/R/multODETest_seeds.csv"
system(paste0("SeedGenerator -n 16 -t -d ", seed_path))

library(tidyverse)

seed <- sample(1:.Machine$integer.max, 1)
# sampled seed: 1306481334
set.seed(seed)

# Set model types and nloci for each: ODE has 4 because 2 are reserved for KZ and KXZ: this means
# one locus for aZ and bZ; we match the mutational target size with nloci = 2 for the non-ODE cases
modelType <- c("Add", "ODE", "Mult")
combos <- c(2, 4, 2)

write_delim(data.frame(modelType = modelType, nloci = combos), "combos.csv", col_names = F, quote = "all")

seeds <- 1:16

cmds <- data.frame(sr = singleRunBashName,
                   nloci = rep(1:length(combos), each = length(seeds)),
                   seed = rep(1:length(seeds), times = length(combos)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/fixedK/multODETest/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
