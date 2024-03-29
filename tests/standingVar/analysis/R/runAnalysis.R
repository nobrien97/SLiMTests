
repoPath <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/analysis/R/"
dataPath <- "/mnt/d/SLiMTests/tests/standingVar/dec11/"

setwd(repoPath)

# Load functions
source("./helperFunctionsAndSetup.R")

# Setup data
source("./wrangle_data.R")

# Calculate heritability
source("./R/calcH2.R")

# Run mutation screen experiment
source("./R/mutationScreenExp.R")

# Load plotting functions
source("./R/figureFunctionsAndSetup.R")

# Create figures
source("./R/figures.R")

# Run stats
source("./R/stats.R")
