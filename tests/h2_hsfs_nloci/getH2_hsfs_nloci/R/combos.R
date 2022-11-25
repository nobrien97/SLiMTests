# File name
filename <- "./getH2_hsfs_nlociSR.sh"

# Run numbers to extract haplotypes
run <- 0:106271

# Chunk numbers to collect output: from 1 to 27, which is 1968 files per chunk

chunk <- 1:54

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 1968)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/h2_hsfs_nloci/getH2_hsfs_nloci/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
