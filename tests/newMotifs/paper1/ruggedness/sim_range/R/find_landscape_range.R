library(tidyverse)
library(latex2exp)

setwd("/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/sim_range")

# Load data
d_qg_tot <- readRDS("/g/data/ht96/nb9894/newMotifs/paper1/d_qg_tot.RDS")


molComp_names <- list("NAR" = c(
  # NAR and PAR
  "aZ",
  "bZ",
  "KZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
),
"FFBH" = c(
  # FFBH
  "aX",
  "KZX",
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "zZ", # baseline expression
  "h", # hill coefficient
  "gX" # X multiplier
)
)

molComp_names[["PAR"]] <- molComp_names[["NAR"]]
molComp_names[["FFLI1"]] <- molComp_names[["FFLC1"]]

model_names_noquote <- c("NAR", "PAR", "FFLC1", 
                         "FFLI1", "FFBH")

all_molcomp_features <- c(
  "aX" = TeX("$\\alpha_X$", output = "character"),
  "KZX" = TeX("$K_{ZX}$", output = "character"),
  "aY" = TeX("$\\alpha_Y$", output = "character"),
  "bY" = TeX("$\\beta_Y$", output = "character"),
  "KY" = TeX("$K_Y$", output = "character"),
  "aZ" = TeX("$\\alpha_Z$", output = "character"),
  "bZ" = TeX("$\\beta_Z$", output = "character"),
  "KZ" = TeX("$K_{Z}$", output = "character"),
  "KXZ" = TeX("$K_{XZ}$", output = "character"),
  "zZ" = TeX("$\\zeta_Z$", output = "character"), # baseline expression
  "h" = TeX("$h$", output = "character"), # hill coefficient
  "gX" = TeX("$\\gamma_X$", output = "character") # X multiplier
)


stat.mode <- function(x) {
  x <- na.omit(x)
  
  if (length(x) < 2) {
    return(NA)
  }
  
  d <- density(log10(x))
  10^(d$x[which.max(d$y)])
}

named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = "."))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}


# Exploring the whole space is very hard with 11 components, so we will construct our 
# fitness landscape along the PC1 and PC2 axes of the real data
# We can save the principal components and randomly generate points within the range of 
# PC1 and PC2 values, translate back to trait space, compute the fitness, then convert
# back to PC space
# This is done on a per simulation basis so we can align our fitness optima with the landscape
# and plot the individual walks

# Rearrange data and filter out very large values
d_filtered_qg <- d_qg_tot %>% filter(isAdapted == T, r == 0.1) %>% 
  select(model, seed, dataset, starts_with("meanMC")) %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = starts_with("meanMC"), names_prefix = "meanMC", names_to = "molComp",
               values_to = "molComp_value") %>%
  mutate(model = factor(model, labels = model_names_noquote),
         molComp = as.integer(molComp)) %>%
  filter(var(molComp_value) > 0) %>%
  group_by(model, seed, dataset) %>%
  mutate(molComp_label = molComp_names[[as.character(model[1])]][molComp]) %>%
  filter(!is.na(molComp_label)) %>%
  select(-molComp) %>%
  group_by(model, seed, dataset, molComp_label) %>%
  mutate(m = 2 * stat.mode(molComp_value),
         min_value = min(molComp_value),
         max_value = (m - min_value + m)) %>%
  ungroup() %>%
  filter(between(molComp_value, min_value, max_value)) %>%
  select(-c(m, min_value, max_value)) %>%
  group_by(model, seed, dataset, id) %>%
  pivot_wider(names_from = molComp_label, 
              values_from = molComp_value) %>%
  select(model, seed, dataset, id, names(all_molcomp_features)) %>%
  ungroup() %>%
  group_by(model, seed, dataset) %>%
  filter(!if_any(all_of(molComp_names[[as.character(cur_group()$model)]]), is.na)) %>%
  ungroup()

# get the pca for each model
d_filtered_split <- d_filtered_qg %>%
  named_group_split(model, seed, dataset)

d_pc_qg <- list()

for (i in seq_along(d_filtered_split)) {
  x <- d_filtered_split[[i]]
  # Group values are in order model, seed, dataset
  grp <- unlist(strsplit(names(d_filtered_split[i]), "\\."))
  result <- prcomp(x[,molComp_names[[grp[1]]]])
  
  # New group is in order model, dataset, seed
  d_pc_qg[[grp[1]]][[grp[[3]]]][[grp[2]]] <- result
}


saveRDS(d_pc_qg, "d_pc_qg.RDS")

