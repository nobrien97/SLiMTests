# File name
filename <- "./getH2_equalKSR.sh"

#  wc -l slim_haplo.csv -> 118080000/2000 haplotypes = 59040 samples


# Run numbers to extract haplotypes
run <- 0:59039

# Chunk numbers to collect output: from 1 to 30, which is 1968 files per chunk

chunk <- 1:30

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 1968),
    haplos = rep(0:7, each = 7380)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/equalK/getH2/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
