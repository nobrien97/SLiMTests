# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/randomisedStarts/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./randomisedStartsSR.sh"

nloci <- 4^(1:5)
tau <- c(0.0125, 0.125, 1.25)
rwide <- 10^seq(-10, -1, by = 1)
lhc <- expand.grid(nloci, tau, rwide)

models <- c("\'Add\'", "\'ODE\'", "\'K\'")

lhc <- lhc %>% slice(rep(1:n(), each = 3))
lhc$model <- rep(models, times = 150)

seeds <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/randomisedStarts/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
