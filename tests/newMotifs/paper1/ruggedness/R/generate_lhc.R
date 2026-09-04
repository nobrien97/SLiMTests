library(DoE.wrapper)


# seed <- sample(1:.Machine$integer.max, 1)
# > seed
# [1] 1649063102
seed <- 1649063102

set.seed(seed)

models <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
comps <- c("aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
           "aZ", "bZ", "Hilln", "XMult", "base")
nComps <- length(comps)
NUM_RUNS <- 10000
NUM_BACKGROUNDS <- 10
MAX_COMP_SIZE <- log(3)

# Generate hypercube of parameters
# Hypercube is NUM_RUNs per molecular component, per genetic background
pars <- lhs.design(NUM_RUNS, 1, type = "random")
pars <- pars * MAX_COMP_SIZE
pars <- pars[,1]

# Generate backgrounds
parBackgrounds <- matrix(runif((length(comps)) * NUM_BACKGROUNDS), ncol = (length(comps)))
parBackgrounds <- parBackgrounds * MAX_COMP_SIZE

# Replicate rows for the number of molecular components
parBackgrounds <-  parBackgrounds %x% rep(1, nComps)

# Replicate rows for the number of replicate values
parBackgrounds <-  rep(1, NUM_RUNS) %x% parBackgrounds

# replace diagonal components for each component and background
for (i in seq_len(NUM_RUNS)) {
  # Fill NUM_BACKGROUNDS diagonals
  for (j in seq_len(NUM_BACKGROUNDS)) {
    k <- (i - 1) * nComps + j
    offset_start <- ((i - 1) * nComps * NUM_BACKGROUNDS) + ((j - 1) * nComps) + 1 # number of previous runs
    offset_end <-  offset_start + (nComps - 1)
    diag(parBackgrounds[offset_start:offset_end,]) <- pars[i]
  }
}

colnames(parBackgrounds) <- comps

pars <- as.data.frame(parBackgrounds)
saveRDS(pars, "pars.RDS")


seeds <- sample(1:.Machine$integer.max, NUM_RUNS)
seeds <- rep(seeds, each = NUM_BACKGROUNDS * nComps)
saveRDS(seeds, "seeds.RDS")
