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

# For each model and molcomp, find the minimum and maximum values
d_qg_tot %>%
  pivot_longer(cols = starts_with("meanMC"), names_prefix = "meanMC", names_to = "molComp",
               values_to = "molComp_value") %>%
  mutate(model = factor(model, labels = model_names_noquote),
         molComp = as.integer(molComp)) %>%
  group_by(model, molComp) %>%
  filter(var(molComp_value) > 0) %>% # remove invalid molcomps
  mutate(q2 = quantile(molComp_value, probs = 0.5),
         m = stat.mode(molComp_value),
         min_to_m = m - min(molComp_value)) %>% 
  select(model, molComp, molComp_value, q2, m, min_to_m) -> d_qg_vals

d_qg_vals %>%
  group_by(model, molComp) %>%
  summarise(mean = mean(molComp_value),
            min = min(molComp_value),
            Q1 = quantile(molComp_value, probs = 0.25),
            Q2 = quantile(molComp_value, probs = 0.5),
            Q3 = quantile(molComp_value, probs = 0.75),
            m = stat.mode(molComp_value),
            max = max(molComp_value)
) -> sum_qg_vals


print(sum_qg_vals, n = 43)

# Plot distribution
ggplot(d_qg_vals %>% filter(molComp_value <= (m + min_to_m)), # keep values within half an order of magnitude of the mode
    aes(x = molComp_value)) +
    facet_wrap(model ~ molComp, scales = "free") +
    geom_histogram(bins = 100) + 
    #scale_x_log10() +
    theme_bw() +
    theme(text = element_text(size = 12)) -> plt_dist
ggsave("/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/sim_range/plt_range_dist.png", plt_dist, device = png,
        width = 16, height = 16, dpi = 600)


d_qg_vals %>%
  group_by(model, molComp) %>%
  summarise(molComp_label = molComp_names[[as.character(model[1])]][molComp[1]],
            m = stat.mode(molComp_value),
            min_value = min(molComp_value),
            max_value = m + (m - min_value)) -> d_molcomp_valuerange
nrow(d_molcomp_valuerange)
print(d_molcomp_valuerange, n = 43)


# Exploring the whole space is very hard with 11 components, so we will construct our 
# fitness landscape along the PC1 and PC2 axes of the real data
# We can save the principal components and randomly generate points within the range of 
# PC1 and PC2 values, translate back to trait space, compute the fitness, then convert
# back to PC space

# Rearrange data and filter out very large values
d_filtered_qg <- d_qg_tot %>% filter(isAdapted == T) %>% 
  select(model, starts_with("meanMC")) %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = starts_with("meanMC"), names_prefix = "meanMC", names_to = "molComp",
               values_to = "molComp_value") %>%
  mutate(model = factor(model, labels = model_names_noquote),
         molComp = as.integer(molComp)) %>%
  filter(var(molComp_value) > 0) %>%
  group_by(model) %>%
  mutate(molComp_label = molComp_names[[as.character(model[1])]][molComp]) %>%
  filter(!is.na(molComp_label)) %>%
  select(-molComp) %>%
  group_by(model, molComp_label) %>%
  mutate(m = 2 * stat.mode(molComp_value),
         min_value = min(molComp_value),
         max_value = (m - min_value + m)) %>%
  ungroup() %>%
  filter(between(molComp_value, min_value, max_value)) %>%
  select(-c(m, min_value, max_value)) %>%
  group_by(model, id) %>%
  pivot_wider(names_from = molComp_label, 
              values_from = molComp_value) %>%
  select(model, id, names(all_molcomp_features)) %>%
  ungroup() %>%
  group_by(model) %>%
  filter(!if_any(all_of(molComp_names[[as.character(cur_group()$model)]]), is.na)) %>%
  ungroup()

# get the pca for each model
d_pc_qg <- d_filtered_qg %>%
  group_by(model) %>%
  group_map(~ prcomp(.x[,molComp_names[[as.character(.y$model[1])]]]))
saveRDS(d_pc_qg, "d_pc_qg.RDS")

