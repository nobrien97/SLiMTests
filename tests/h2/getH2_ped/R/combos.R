# File name
filename <- "./getH2SR.sh"

# Run numbers to extract haplotypes
run <- 0:203615

# Chunk numbers to collect output: from 1 to 101, which is 2016 files per chunk

chunk <- 1:101

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 2016)
)

write.table(cmd_df, "./cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
