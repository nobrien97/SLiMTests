# Create cmds.txt and combo file
path <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/orthSel/R/"
setwd(path)
# Generate cmds.txt
singleRunBashName <- "./orthSelSR.sh"

# Models are the same as in the newMotifs job, 15 of them
models <- 15
seeds <- 208

# This job is split in two, we want to avoid double-saving output so we add the job number
njobs <- 2

cmds <- data.frame(sr = singleRunBashName,
                   model = rep(11:models, each = seeds),
                   seed = rep(1:seeds, times = models-10),
                   njob = 1 )#rep(0:(njobs - 1), each = (seeds * models) / njobs))

write.table(cmds, "/mnt/c/GitHub/SLiMTests/tests/newMotifs/paper1/orthSel/PBS/cmds.txt", 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
