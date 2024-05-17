# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/sanityChecks/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./sanityCheckSR.sh"

library(tidyverse)

nloci <- 1024
tau <- c(0.0125, 1.25)
rwide <- c(0, 0.5)
w <- c(0, 0.05, 0.1)
lhc <- expand.grid(nloci, tau, rwide, w)

models <- c("\'Add\'", "\'ODE\'", "\'K\'")

lhc <- lhc %>% slice(rep(1:n(), each = 3))
lhc$model <- rep(models, times = nrow(lhc)/3)
write_delim(lhc, "combos.csv", col_names = F)


seeds <- 1:12

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/sanityChecks/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
