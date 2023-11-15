# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./standingVarSR.sh"

seed_path <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/standingVar_seeds.csv"
system(paste0("SeedGenerator -n 48 -t -d ", seed_path))

library(tidyverse)
library(DoE.wrapper)
library(DiceDesign)
library(GGally)
library(latex2exp)

iter <- 0
repeat {
  if (iter >= 1000) {
    stop("Unable to find hypercube with max correlation <0.001 in 100 iterations")
    break 
  }
  seed <- sample(1:.Machine$integer.max, 1)
  set.seed(seed)
  lhc <- lhs.design(
    nruns = 144,
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
  if (maxCor < 0.005)
    break
} 

# nloci to integer
lhc$nloci <- ceiling(lhc$nloci)


# final seed: 1605869673
set.seed(1605869673)
lhc <- lhs.design(
  nruns = 144,
  nfactors = 3,
  type = "maximin",
  factor.names = list(
    "nloci" = c(1, 1000),
    "tau"   = c(0.01, 1.5),
    "rwide" = c(0, 0.5)
  ))
lhc$nloci <- ceiling(lhc$nloci)

# plot fit quality
ggpairs(lhc, progress = F,
        lower = list(continuous = wrap("points", size = 0.1)),
        upper = list(continuous = wrap("cor", size = 10)),
        columnLabels = c(TeX("Number of loci $(n_{loci})$", output = "character"), 
                         TeX("Effect size variance ($\\tau$)", output = "character"), 
                         TeX("Recombination rate $(r)$", output = "character")),
        labeller = "label_parsed") +
  theme_classic() +
  theme(text = element_text(size = 16, face = "bold"),
        panel.spacing = unit(1.5, "lines")) -> lhc_pairs
lhc_pairs
ggsave("lhc_pairs.png", lhc_pairs, device = png, width = 8, height = 8)


# Calculate l2-star discrepancy: values close to 0 indicate good spread, 1 is bad
# measures overall uniformity
discrepancyCriteria(lhc, "L2star")
# 0.008073495

models <- c("\'Add\'", "\'ODE\'", "\'K\'")

lhc <- lhc %>% slice(rep(1:n(), each = 3))
lhc$model <- rep(models, times = 144)

write_delim(lhc, "combos.csv", col_names = F)

seeds <- 1:48

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(1:nrow(lhc), each = length(seeds)),
                   seed = rep(1:length(seeds), times = nrow(lhc)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/standingVar/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)


nloci <- 4^(0:4)
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
