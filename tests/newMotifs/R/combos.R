# Create cmds.txt and combo file
path <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./newMotifsSR.sh"

# 48 replicates
seed_path <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/newMotifs_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", seed_path))

library(tidyverse)

models <- c("\'NAR\'", "\'PAR\'", "\'FFLC1\'", "\'FFLI1\'", "\'FFBH\'")
rwide <- c(1e-10, 1e-5, 1e-1)

lhc <- expand.grid(models, rwide)

write_delim(lhc, "combos.csv", col_names = F)

seeds <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
