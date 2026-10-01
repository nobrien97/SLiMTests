library(tidyverse)

setwd("/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/sim_range")

# Load data
d_qg_tot <- readRDS("/g/data/ht96/nb9894/newMotifs/paper1/d_qg_tot.RDS")


molComp_names <- list("NAR" = c(
  # NAR and PAR
  "aZ",
  "bZ",
  "KZ",
  "KXZ",
  "base", # baseline expression
  "Hilln", # hill coefficient
  "XMult" # X multiplier
),

"FFLC1" = c(
  # FFLC1 and FFLI1
  "aY",
  "bY",
  "KY",
  "aZ",
  "bZ",
  "KXZ",
  "base", # baseline expression
  "Hilln", # hill coefficient
  "XMult" # X multiplier
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
  "base", # baseline expression
  "Hilln", # hill coefficient
  "XMult" # X multiplier
)
)

molComp_names[["PAR"]] <- molComp_names[["NAR"]]
molComp_names[["FFLI1"]] <- molComp_names[["FFLC1"]]

model_names_noquote <- c("NAR", "PAR", "FFLC1", 
                         "FFLI1", "FFBH")


stat.mode <- function(x) {
  d <- density(log10(na.omit(x)))
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


saveRDS(d_molcomp_valuerange, "d_molcomp_maxvals.RDS")