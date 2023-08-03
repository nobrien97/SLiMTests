# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/fixedK/calcHet/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./calcHetSR.sh"

# get list of seeds from fixedK, moreReps, moreReps2
moreReps_seeds <- read.csv("/mnt/c/GitHub/SLiMTests/tests/fixedK/moreReps/R/fixedK_seeds.csv", header = F)
moreReps2_seeds <- read.csv("/mnt/c/GitHub/SLiMTests/tests/fixedK/moreReps2/R/fixedK_seeds.csv", header = F)

seeds <- rbind(moreReps_seeds, moreReps2_seeds)

modelindex <- 1:2
seeds <- data.frame(seed = seeds[rep(seq_len(nrow(seeds)), each = 2), ],
                    modelindex = 1:2)
head(seeds)

seeds$file <- paste0("\"/g/data/ht96/nb9894/fixedK/pop_states/slim_popstate", seeds$seed, "_", seeds$modelindex, ".bin\"")


cmds <- data.frame(sr = singleRunBashName,
                   modelindex = seeds$modelindex,
                   seed = seeds$seed,
                   file = seeds$file)

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/fixedK/calcHet/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
