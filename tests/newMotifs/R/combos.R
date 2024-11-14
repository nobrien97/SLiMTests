# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./speedTestSR.sh"

# 19 replicates
seed_path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/speedTest_seeds.csv"
system(paste0("SeedGenerator -n 19 -t -d ", seed_path))

library(tidyverse)

models <- c("\'NAR\'", "\'PAR\'", "\'FFLC1\'", "\'FFLI1\'", "\'FFBH\'")

lhc <- expand.grid(models)

write_delim(lhc, "combos.csv", col_names = F)

seeds <- 1:19

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)


nloci <- 4^(1:5)
tau <- c(0.0125, 0.125, 1.25)
rwide <- 10^seq(-10, -1, by = 1)
lhc <- expand.grid(nloci, tau, rwide)

ggpairs(lhc,
        diag = list(continuous = "barDiag"), 
        columnLabels = c(TeX("Number of loci $(n_{loci})$", output = "character"),
                         TeX("Effect size variance $(\\tau)$", output = "character"), 
                         TeX("Recombination rate $(r)$", output = "character")),
        labeller = "label_parsed"
)
ggsave("factorialSample.png", device = png)

models <- c("\'Add\'", "\'ODE\'", "\'K\'")

lhc <- lhc %>% slice(rep(1:n(), each = 3))
lhc$model <- rep(models, times = 150)
write_delim(lhc, "combos.csv", col_names = F)


seeds <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
