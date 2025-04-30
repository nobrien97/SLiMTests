# File name
filename <- "./getH2SR.sh"

# Run numbers to extract haplotypes
run <- 0:7329

# Chunk numbers to collect output: from 1 to 2, which is 1440 files per chunk

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = 1
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/h2/getH2/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
