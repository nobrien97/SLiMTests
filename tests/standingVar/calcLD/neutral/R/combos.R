library(tidyverse)

# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/neutral/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./calcLD_neutralSR.sh"

seeds <- 1:48
# seeds copied from standing var:
file.copy("../../../R/standingVar_seeds.csv", "./standingVar_seeds.csv")

models <- c("\'Add\'", "\'ODE\'", "\'K\'")
rwide <- 10^seq(-10, -1, by = 9)
combos <- expand.grid(models, rwide)
write_delim(combos, "combos.csv", col_names = F)

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(combos), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(combos)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/neutral/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

