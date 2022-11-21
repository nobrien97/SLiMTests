# File name
filename <- "./getH2_hsfsSR.sh"

# Run numbers to extract haplotypes
run <- 0:29519

# Chunk numbers to collect output: from 1 to 15, which is 1968 files per chunk

chunk <- 1:15

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 1968)
)

write.table(cmd_df, "./cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
