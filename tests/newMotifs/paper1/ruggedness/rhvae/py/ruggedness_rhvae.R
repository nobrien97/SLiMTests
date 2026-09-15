# Run an RHVAE via Julia for a better fitness landscape
library(tidyverse)

# Load data, split into groups
DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/"
d_ruggedness <- data.table::fread(paste0(DATA_PATH, "log3/d_ruggedness_permolcomp.csv"), 
                                       header = F)

setwd(paste0(DATA_PATH, "log3"))

colnames(d_ruggedness) <- c("step", "model", "dataset", "fitness", "startW", 
                                 "endW", "netChangeW", "sumChangeW", "numFitnessHoles", 
                                 "nSteps", "aX", "KZX", "aY", "bY", "KY", "KZ", "KXZ",
                                 "aZ", "bZ", "Hilln", "XMult", "base",
                                 "molComp", "bkg")

# Filter to only the columns we care about
d_ruggedness <- d_ruggedness %>%
  select(2:4, 11:22)

molComp_names <- list("NAR" = c(
  # NAR and PAR
  "aZ",
  "bZ",
  "KZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),
"FFBH" = c(
  # FFBH
  "aX",
  "KZX",
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
)
)
molComp_names[["PAR"]] <- molComp_names[["NAR"]]
molComp_names[["FFLI1"]] <- molComp_names[["FFLC1"]]


d_ruggedness_nar <- d_ruggedness %>% filter(model == "NAR") %>% ungroup() %>%
  select(fitness, molComp_names[["NAR"]])
d_ruggedness_par <- d_ruggedness %>% filter(model == "PAR") %>% ungroup() %>%
  select(fitness, molComp_names[["PAR"]])
d_ruggedness_fflc1 <- d_ruggedness %>% filter(model == "FFLC1") %>% ungroup() %>%
  select(fitness, molComp_names[["FFLC1"]])
d_ruggedness_ffli1 <- d_ruggedness %>% filter(model == "FFLI1") %>% ungroup() %>%
  select(fitness, molComp_names[["FFLI1"]])
d_ruggedness_ffbh <- d_ruggedness %>% filter(model == "FFBH") %>% ungroup() %>%
  select(fitness, molComp_names[["FFBH"]])

write_csv(d_ruggedness_nar, "d_ruggedness_nar.csv")
write_csv(d_ruggedness_par, "d_ruggedness_par.csv")
write_csv(d_ruggedness_fflc1, "d_ruggedness_fflc1.csv")
write_csv(d_ruggedness_ffli1, "d_ruggedness_ffli1.csv")
write_csv(d_ruggedness_ffbh, "d_ruggedness_ffbh.csv")


# Run julia script
system("julia ruggedness_rhvae.jl")