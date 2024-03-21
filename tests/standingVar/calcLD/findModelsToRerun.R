# Extracts models to rerun from o files showing model replicates that incorrectly
# loaded population files
# The file already is partially processed, removing irrelevant output lines
library(tidyverse)

# Load file
models <- readLines("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/models_to_rerun.txt")
# Remove empty lines
models <- models[models != ""]

# Extract matches
modelindices <- as.integer(str_extract(models, "[0-9]+(?=(, seed))"))
seeds <- as.integer(str_extract(models, "[0-9]+(?=( f))"))

out <- data.frame(modelindex = modelindices,
                  seed = seeds)
# write to file
write_delim(out, "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcLD/rerun_models.csv",
           delim = ",", col_names = F)
