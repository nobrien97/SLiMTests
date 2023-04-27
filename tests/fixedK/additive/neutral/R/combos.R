# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./fixedKSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/R/fixedK_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", seed_path))

library(tidyverse)

seed <- sample(1:.Machine$integer.max, 1)
# sampled seed: 805432144
set.seed(seed)

# In NAR models we have (1,10,100)*4 because even though we aren't using the KXZ and KZ loci,
# they need to have values in the model - this means the aZ and bZ will have 1, 10, 100 each
# In additive models, we don't need to worry about feeding the K loci, so we use 2*nloci instead
# which means the trait will have 2,20,200 loci, which is the same as the trait in the NAR model
combos <- c(1, 10, 100) * 2

write_csv(data.frame(nloci = combos), "combos.csv", col_names = F)

seeds <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   nloci = rep(1:length(combos), each = length(seeds)),
                   seed = rep(1:length(seeds), times = length(combos)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/fixedK/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
