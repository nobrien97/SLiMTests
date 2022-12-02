# File name
filename <- "./getH2_moltrait_fixSR.sh"

# Run numbers to extract haplotypes
run <- 0:98399

# Chunk numbers to collect output: from 1 to 48, which is 2050 files per chunk

chunk <- 1:48

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 2050),
    haplos = rep(0:7, each = 12300)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/h2_moltrait_fix/getH2_moltrait_fix/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
