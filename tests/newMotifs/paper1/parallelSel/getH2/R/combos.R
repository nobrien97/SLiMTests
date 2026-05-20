# File name
filename <- "./getH2SR.sh"

# Run numbers to extract haplotypes
run <- 0:2079

# Chunk numbers to collect output: from 1 to 2, which is 1440 files per chunk

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = 1
)

write.table(cmd_df, "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/parallelSel/getH2/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
