# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/indTrack/R/"
# Generate cmds.txt
singleRunBashName <- "./indTrackSR.sh"

library(tidyverse)
adapted <- read_csv(paste0(path, "adapted_seeds.csv"))
maladapted <- read_csv(paste0(path, "maladapted_seeds.csv"))
wasadapted <- read_csv(paste0(path, "wasadapted_seeds.csv"))

combos <- rbind(adapted, maladapted, wasadapted)
combos$model <- as.character(combos$model)
combos$seed <- as.character(combos$seed)
write.table(combos, paste0(path, "combos.csv"), row.names = F, col.names = F)


cmds <- data.frame(sr = rep(singleRunBashName, times = nrow(combos)),
                   model = rep(1:nrow(combos)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/indTrack/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
