# File name
filename <- "./getH2SR.sh"

# Run numbers to extract haplotypes
run <- 0:1439

# Chunk numbers to collect output: from 1 to 2, which is 1440 files per chunk

chunk <- 1:2

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 720)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/h2/getH2/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
