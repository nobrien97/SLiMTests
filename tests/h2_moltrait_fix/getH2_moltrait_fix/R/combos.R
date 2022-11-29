# File name
filename <- "./getH2_moltrait_fixSR.sh"

# Run numbers to extract haplotypes
run <- 0:150551

# Chunk numbers to collect output: from 1 to 27, which is 1968 files per chunk

chunk <- 1:72

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 2091),
    haplos = rep(0:7, each = 18819)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/h2_moltrait_fix/getH2_moltrait_fix/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
