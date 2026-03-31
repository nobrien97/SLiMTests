library(tidyverse)

# Load in data
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/neutralCorr/slim_qg.csv"
read_csv(DATA_PATH, col_Names = c("gen", "seed", "modelindex", "meanH", paste0("phenomean_", 1:4),
                                  paste0("phenovar_", 1:4), paste0("phenocor_", 1:4),))


Mfile = paste(sim.cycle, asString(seed), modelindex, meanH, phenomean, phenovar, phenocor, meanMolTraits, sep=",");
