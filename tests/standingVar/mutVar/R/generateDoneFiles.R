library(tidyverse)

COMBO_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/"
COMBO_PATH <- "~/tests/standingVar/R/"
FILE_PATH <- "~/tests/standingVar/mutVar/done/"

d_combos <- read.table(paste0(COMBO_PATH, "combos.csv"), header = F,
                       col.names = c("nloci", "tau", "r", "model"))

r_subsample <- c(1e-10, 1e-5, 1e-1)

invalid_combos <- d_combos %>%
    rownames_to_column("id") %>%
    filter(nloci != 1024 | tau != 0.0125 | !(r %in% r_subsample)) %>%
    column_to_rownames("id")

seeds <- 1:48

for (i in 1:nrow(invalid_combos)) {
    modelindex <- rownames(invalid_combos)[i]
    sapply(paste0(FILE_PATH, modelindex, "_", seeds), file.create)
}

# Adjusted tau - we also don't want additive models for these
invalid_combos <- d_combos %>%
    rownames_to_column("id") %>%
    filter(nloci != 1024 | tau != 0.0125 | !(r %in% r_subsample) | model == "Add") %>%
    column_to_rownames("id")

seeds <- 1:48

for (i in 1:nrow(invalid_combos)) {
    modelindex <- rownames(invalid_combos)[i]
    sapply(paste0(FILE_PATH, modelindex, "_", seeds), file.create)
}

