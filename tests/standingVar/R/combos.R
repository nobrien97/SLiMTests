# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/moreReps/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./fixedKSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/moreReps2/R/fixedK_seeds.csv"
system(paste0("SeedGenerator -n 1440 -t -d ", seed_path))

library(tidyverse)

seed <- sample(1:.Machine$integer.max, 1)
# sampled seed: 805432144
set.seed(seed)
combos <- c(2, 4)

write_csv(data.frame(nloci = combos), "combos.csv", col_names = F)

seeds <- 1:1440

cmds <- data.frame(sr = singleRunBashName,
                   nloci = rep(1:length(combos), each = length(seeds)),
                   seed = rep(1:length(seeds), times = length(combos)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/fixedK/moreReps/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
