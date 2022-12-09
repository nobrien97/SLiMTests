# File name
filename <- "./getH2_newTestCrossSR.sh"

# Run numbers to extract haplotypes
run <- 0:59039

# Chunk numbers to collect output: from 1 to 30, which is 1968 files per chunk

chunk <- 1:30

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 1968)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/getH2_newTestCross/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
