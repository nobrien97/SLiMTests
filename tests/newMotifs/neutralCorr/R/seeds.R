library(tidyverse)

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/neutralCorr_seeds.csv"
system(paste0("SeedGenerator -n 208 -t -d ", seed_path))

seeds <- read.csv(seed_path, header = F)
hist(seeds$V1)


# Models
models <- 5
seeds <- 208

singleRunBashName <- "./neutralCorrSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:models, each = seeds),
                   seed = rep(1:seeds, times = models))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)


models <- data.frame(model = c("\'NAR\'", "\'PAR\'", "\'FFLC1\'", "\'FFLI1\'", "\'FFBH\'"))

write_delim(models, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/neutralCorr/R/combos.csv", col_names = F)
