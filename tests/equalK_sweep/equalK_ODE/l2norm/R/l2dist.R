library(tidyverse)
library(factoextra)

d_qg <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/checkpoint/d_qg.RDS")

d_qg$id <- paste(d_qg$seed, d_qg$modelindex, sep = "_")
res <- numeric(length(unique(d_qg$id)))

for (i in seq_along(res)) {
  # Get appropriate subset to construct a 140x4 matrix
  d_new <- d_qg[((i-1)*140+1):(i*140),] 
  id <- d_new[1, 18]
  # calculate l2 norm
  res[i] <- norm(as.matrix(d_new %>% select(aZ, bZ, KZ, KXZ)))
  names(res)[i] <- as.character(id)
}

res2 <- as.matrix(res)

# calculate distance matrix - log10 scale because we have those 
# giant distances from KZ
# problem: KZ is contributing a lot to the distance, much more so than the other components
# even though it doesn't contribute very much to the trait value when it's large
distmat <- dist(log10(res2), diag = T)

# convert to dataframe
d_dist <- as.data.frame(as.matrix(distmat))
d_dist <- d_dist %>% rownames_to_column("row") %>% gather("col", "distance", -row)

# Plot
fviz_dist(distmat, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

write.table(solution, paste0("ode_", run, ".csv"), sep = ",", col.names = F, row.names = F)
