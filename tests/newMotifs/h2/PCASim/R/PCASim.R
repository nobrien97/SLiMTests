library(mcreplicate)

# PCA similarity
bootKrzCorFn <- function(x, group = "", PCASim = F) {
require(evolqg)
require(dplyr)

fn <- ifelse(PCASim, evolqg::PCAsimilarity, evolqg::KrzCor)

if (group != "") {
  grps <- unique(x[,group])
  nGrps <- length(grps)
  
  
  # output data frame
  res <- data.frame(group1 = character(length(grps)^2),
                    group2 = character(length(grps)^2),
                    krzCor = numeric(length(grps)^2))
  
  # Temporary data frame for filling inner loop
  res_tmp <- data.frame(group1 = character(length(grps)),
                        group2 = character(length(grps)),
                        krzCor = numeric(length(grps)))
  
  for (i in seq_along(grps)) {
    for (j in seq_along(grps)) {
      # Sample matrices in different groups
      g_1 <- slice_sample(x[group == grps[i]], n = 1)
      g_2 <- slice_sample(x[group == grps[j]], n = 1)
      res_tmp$group1[j] <- as.character(g_1[1,group])
      res_tmp$group2[j] <- as.character(g_2[1,group])
      res_tmp$krzCor[j] <- fn(g_1$g[[1]], g_2$g[[1]])
    }
    indices <- (nGrps*(i-1) + 1):(nGrps*i)
    res[indices,] <- res_tmp
  }
  return(res)
}

# If group is "", sample two random matrices and return that
g1 <- slice_sample(x, n = 1)
g2 <- slice_sample(x, n = 1)
return(fn(g1$g[[1]], g2$g[[1]]))
}

model_names <- c("'NAR'", "'PAR'", "'FFLC1'", "'FFLI1'", "'FFBH'")
tpLevels <- c("Start", "End")

# Read in data
WRITE_PATH <- "/g/data/ht96/nb9894/newMotifs/h2/PCASim/"
DATA_PATH <- "/home/564/nb9894/tests/newMotifs/h2/R/"
krz_in <- readRDS(paste0(DATA_PATH, "pca_tp_in.RDS"))




# Bootstrap in ten parts for RAM reasons
# This is slow: uncomment to run, otherwise read in precalculated data
# Generate seeds
#newseed <- sample(1:.Machine$integer.max, 10)
# 407844323 2049133531  970639651  452738391 1161959903  506461634 2087592727 1740805001
# 1698460455  886904095
newseed <- c(407844323L, 2049133531L, 970639651L, 452738391L, 1161959903L, 
             506461634L, 2087592727L, 1740805001L, 1698460455L, 886904095L)
bootPCASim <- vector(mode = "list", length = 10)
for (i in seq_along(newseed)) {
  # Set seed
  set.seed(newseed[i])
  # Run replicate
  res <- mcreplicate::mc_replicate(1000, bootKrzCorFn(krz_in, "group", T))
  bootPCASim[[i]] <- unnest(as.data.frame(t(res)), cols = everything())
}
# Output list into combined df
bootPCASim2 <- bind_rows(bootPCASim)
bootPCASim <- bootPCASim2 %>%
  separate(group1, c("model1", "r1", "timePoint1"), "\\.",
           extra = "merge") %>%
  separate(group2, c("model2", "r2", "timePoint2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = log10(as.numeric(r1)),
         r2 = log10(as.numeric(r2)),
         model1 = factor(model1, levels = model_names),
         model2 = factor(model2, levels = model_names),
         timePoint1 = factor(timePoint1, levels = tpLevels),
         timePoint2 = factor(timePoint2, levels = tpLevels)) %>%
  rename(PCASim = krzCor)

saveRDS(bootPCASim, paste0(WRITE_PATH, "d_bootPCASim.RDS"))
