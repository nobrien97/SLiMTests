# File name
filename <- "./getH2SR.sh"

#  wc -l slim_haplo_fix.csv -> 132174 samples / 3147 = 42 chunks 

# When combining outputs, slim_haplo.csv joined both fixations and non-fixations
# Fixations all at the bottom of the file
# So ignore the last 132174 rows of slim_haplo.csv

# Run numbers to extract haplotypes
run <- 0:132173

# Chunk numbers to collect output: from 1 to 42, which is 3147 files per chunk

chunk <- 1:42

cmd_df <- data.frame(
    filename = filename,
    run = run,
    chunk = rep(chunk, each = 3147)
)

write.table(cmd_df, "/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
