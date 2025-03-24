# Generates cmds.txt for nosil per mol comp ruggedness job 

singleRunBashName = "./ruggednessSR.sh"

NUM_RUNS <- 10000
NUM_BACKGROUNDS <- 10
nComps <- 12

# Total rows in input
TOTAL_LENGTH <- NUM_RUNS * NUM_BACKGROUNDS * nComps

# 10 replicates per run, each run will return a dataframe with 5000 (5 per model) rows in it
# for 1000 total files to combine
ROWS_PER_RUN <- nComps * NUM_BACKGROUNDS * 10

N_CMDS <- TOTAL_LENGTH / ROWS_PER_RUN


cmds <- data.frame(sr = singleRunBashName,
                   index = 1:N_CMDS)

write.table(cmds, "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/ruggedness/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
