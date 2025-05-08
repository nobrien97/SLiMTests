library(tidyverse)
library(ggh4x)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/mutVar/R")
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/mutVar/"
COMBOS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))


# Read data
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                            col.names = c("replicate", "seed", "modelindex", "variance"))


d_combos <- read.table(paste0(COMBOS_PATH, "combos.csv"), header = F,
                       col.names = c("nloci", "tau", "r", "model"))


d_mutvar <- d_mutvar %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar <- AddCombosToDF(d_mutvar)

d_mutvar_sum <- d_mutvar %>%
  group_by(model, r) %>%
  summarise(meanVar = mean(log10(variance)),
            CIVar = CI(log10(variance)))

# Plot mutational variance
ggplot(d_mutvar, aes(x = model, y = log10(variance), colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_sum, 
             aes(x = model, y = meanVar),
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  labs(x = "Model", 
       y = "Mutational variance (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none")
ggsave("plt_vm.png", device = png, width = 5, height = 9)  

# Plot Vm in adjusted tau runs
d_mutvar_adjtau <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_adjtau.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar_adjtau$normalised <- T

d_mutvar$normalised <- F

d_mutvar_adjtau <- d_mutvar_adjtau %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_adjtau <- AddCombosToDF(d_mutvar_adjtau)

# Join with regular
d_mutvar2 <- full_join(d_mutvar, d_mutvar_adjtau, by = c("replicate", "seed", "modelindex",
                                                         "normalised", "variance", "model",
                                                         "nloci", "tau", "r"))

d_mutvar2 <- d_mutvar2 %>%
  mutate(scaled = if_else(normalised, "Scaled tau", "Unscaled tau"),
         scaled = factor(scaled, levels = c("Unscaled tau", "Scaled tau")))

d_mutvar_sum <- d_mutvar2 %>%
  group_by(model, r, scaled) %>%
  summarise(meanVar = mean(log10(variance)),
            CIVar = CI(log10(variance)))

ggplot(d_mutvar2 %>%
         mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE")), 
       aes(x = interaction(model, scaled), y = log10(variance), colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_sum %>% ungroup() %>%
               mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE")), 
             aes(x = interaction(model, scaled), y = meanVar),
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1)) +
  scale_x_discrete(guide = "axis_nested") +
  labs(x = "Model", 
       y = "Mutational variance (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none")
ggsave("plt_vm_scaled.png", device = png, width = 5, height = 9)  


# Look at adaptive walks
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/adjustedTau/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

# Cowplot 1.1.3 bug: won't get legend, this fixes
get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


d_combos <- read.table("/mnt/c/GitHub/SLiMTests/tests/standingVar/R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg



# Attach additive replicates as well
ADD_DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"

d_qg_add <- data.table::fread(paste0(ADD_DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)


d_qg_add %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg_add

d_qg_add <- AddCombosToDF(d_qg_add) 

r_subsample <- c(1e-10, 1e-5, 1e-1)

d_qg_add <- d_qg_add %>% filter(model == "Add", r %in% r_subsample,
                                nloci == 1024, tau == 0.0125) 

# join
d_adapted <- full_join(d_qg, d_qg_add)

d_adapted_sum <- d_adapted %>% 
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

ggplot(d_adapted_sum,
       aes(x = gen, y = meanPhenomean, colour = model),
       group = as.factor(seed)) +
  facet_grid(log10(r)~.) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean, 
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-"), guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 10),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_adapt_mutScale.png", width = 6, height = 5, device = png)
