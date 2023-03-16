# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/equalK_sweep/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./sweepExpSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/equalK_sweep/R/equalK_sweep_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", seed_path))

library(tidyverse)

combos <- read_csv("combos.csv")
combos$nloci <- round(combos$nloci)
write.table(combos %>% select(-1), "combos.csv", row.names = F)

seed <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(combos), each = length(seed)),
                   seed = rep(1:length(seed), times = nrow(combos)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/equalK_sweep/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
