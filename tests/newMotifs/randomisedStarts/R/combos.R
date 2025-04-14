# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./newMotifsSR.sh"

# Generate 104 replicates
seed_path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/R/newMotifs_seeds.csv"
system(paste0("SeedGenerator -n 104 -t -d ", seed_path))

# Models are the same as in the newMotifs job, 15 of them
models <- 1:15
seeds <- 1:104

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(models, each = length(seeds)),
                   seed = rep(1:length(seeds), times = length(models)))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
