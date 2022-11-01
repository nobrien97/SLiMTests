##############################################################################################################
#  Long test
##############################################################################################################

# Run an experiment with 1000 replicates for a single parameter combination, network and additive

#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.
#  Thanks to David Green at the UQ RCC


# Parse parameters: seed and latin square row number
HOME <- Sys.getenv("HOME")

args <- commandArgs(trailingOnly=TRUE)

testDir <- paste0(HOME,"/tests/h2/")

# seed argument from command line
s <- as.numeric(args[1])
m <- as.numeric(args[2])

# Seeds generated with seedgenerator 
seeds <- read.csv(paste0(testDir, "R/h2_seeds.csv"),header=T)

# m is either 0 or 1
model <- "additive.slim"
if (m) {
  model <- "network.slim"
}

# Check if we are running an additive or network model
  # Run SLiM
  # Sublaunch SLiM with the appropriate values
slim_out <- system(sprintf("$HOME/SLiM/slim -s %s %s/slim/%s",
                              as.character(seeds$Seed[s]),
                              testDir,
                              model), intern=T)


