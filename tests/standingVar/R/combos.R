# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./standingVarSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/standingVar_seeds.csv"
system(paste0("SeedGenerator -n 100 -t -d ", seed_path))

library(tidyverse)
library(DoE.wrapper)
library(DiceDesign)

iter <- 0
repeat {
  if (iter >= 100) {
    errorCondition("Unable to find hypercube with max correlation <0.001 in 100 iterations")
    break 
  }
  seed <- sample(1:.Machine$integer.max, 1)
  set.seed(seed)
  lhc <- lhs.design(
    nruns = 256,
    nfactors = 3,
    type = "maximin",
    factor.names = list(
      "nloci" = c(1, 1000),
      "tau"   = c(0.01, 1.5),
      "rwide" = c(0, 0.5)
    ),
    seed = seed
  )
  iter = iter + 1
  maxCor <- max(abs(cor(lhc)[upper.tri(cor(lhc))]))
  if (maxCor < 0.001)
    break
} 

# nloci to integer
lhc$nloci <- ceiling(lhc$nloci)


# final seed: 324263921

# Calculate l2-star discrepancy: values close to 0 indicate good spread, 1 is bad
# measures overall uniformity
discrepancyCriteria(lhc, "L2star")
# 0.006626378

write_csv(lhc, "combos.csv", col_names = F)

seeds <- 1:100

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
