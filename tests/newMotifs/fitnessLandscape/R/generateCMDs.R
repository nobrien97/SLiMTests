# Generates cmds.txt for nosil per mol comp ruggedness job 
library(tidyverse)

seeds <- readRDS(paste0(DATA_PATH, "seeds.RDS"))
pars <- readRDS(paste0(DATA_PATH, "pars.RDS"))

singleRunBashName = "./ruggednessSR.sh"

cmds <- data.frame(sr = singleRunBashName,
                   index = 1:length(seeds))

write.table(cmds, "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/ruggedness/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
