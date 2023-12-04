
repoPath <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/analysis/R/"
dataPath <- "/mnt/d/SLiMTests/tests/standingVar/initial/"

setwd(repoPath)

# Load functions
source("./R/helperFunctionsAndSetup.R")

# Setup data
source("./R/wrangle_data.R")

# Run mutation screen experiment
source("./R/mutationScreenExp.R")

# Load plotting functions
source("./R/figureFunctionsAndSetup.R")

# Create figures
source("./R/figures.R")

# Run stats
source("./R/stats.R")
