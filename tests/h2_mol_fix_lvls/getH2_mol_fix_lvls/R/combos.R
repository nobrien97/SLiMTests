# File name
filename <- "./getH2_mol_fix_lvlsSR.sh"

# Run numbers to extract haplotypes
run <- 0:5903

# Chunk numbers to collect output: from 1 to 3, which is 1968 files per chunk

chunk <- 1:3

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 1968)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/h2_mol_fix_lvls/getH2_mol_fix_lvls/PBS/cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
