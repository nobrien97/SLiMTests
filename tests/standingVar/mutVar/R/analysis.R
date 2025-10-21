library(tidyverse)
library(legendry)
library(ggbeeswarm)
library(paletteer)
library(latex2exp)
library(patchwork)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/mutVar/R")
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/mutVar/"
COMBOS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))


# Read data
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                            col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar_percomp <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_percomp.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", "cov_aZ", "cov_bZ",
                                            "cov_KZ", "cov_KXZ"))


d_combos <- read.table(paste0(COMBOS_PATH, "combos.csv"), header = F,
                       col.names = c("nloci", "tau", "r", "model"))




d_mutvar <- d_mutvar %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar <- AddCombosToDF(d_mutvar)

d_mutvar_sum <- d_mutvar %>% filter(tau == 0.0125) %>%
  group_by(model) %>%
  summarise(meanVar = mean(variance),
            CIVar = CI(variance))

d_mutvar_percomp <- d_mutvar_percomp %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_percomp <- AddCombosToDF(d_mutvar_percomp)

d_mutvar_percomp <- d_mutvar_percomp %>%
  pivot_longer(cols = starts_with("cov_"), names_to = "molComp", values_to = "cov") %>%
  # Remove massive values
  filter(cov < 1e1, cov > -1e1)

# # Remove outliers (Hampel)
# lower_bound <- median(d_mutvar_percomp$cov) - 3 * mad(d_mutvar_percomp$cov, constant = 1)
# upper_bound <- median(d_mutvar_percomp$cov) + 3 * mad(d_mutvar_percomp$cov, constant = 1)
# 
# d_mutvar_percomp <- d_mutvar_percomp %>%
#   filter(cov > lower_bound, cov < upper_bound)


d_mutvar_percomp_sum <- d_mutvar_percomp %>%
  # Average
  ungroup() %>%
  group_by(molComp, model) %>%
  summarise(meanCOV = mean(cov),
            CICOV = CI(cov))

# Plot mutational variance
ggplot(d_mutvar %>% filter(tau == 0.0125), 
       aes(x = model, y = variance, colour = model)) +
#  facet_nested("Recombination rate (log10)" + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_sum,
             aes(x = model, y = meanVar),
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  labs(x = "Model", 
       y = "Mutational variance") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none")
ggsave("plt_vm.png", device = png, width = 5, height = 9)  

# Plot mutational covariance
ggplot(d_mutvar_percomp, 
       aes(x = model, y = cov, colour = molComp)) +
  #facet_nested("Recombination rate (log10)" + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_percomp_sum, 
             aes(x = model, y = meanCOV, group = molComp),
             shape = 3, size = 2, colour = "black", position = position_dodge(0.8), 
             inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, 
                                           direction = -1),
                      labels = c(TeX("$\\alpha_Z$"), 
                                 TeX("$\\beta_Z$"),
                                 TeX("$K_Z$"),
                                 TeX("$K_{XZ}$"))) +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  coord_cartesian(ylim = c(-1e-4, 1e-4)) +
  labs(x = "Model", 
       y = "Mutational covariance",
       colour = "Molecular component") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom")
ggsave("plt_covm.png", device = png, width = 5, height = 9)  


# Plot Vm in adjusted tau runs
d_mutvar_adjtau <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_adjtau.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", "variance"))

d_mutvar_percomp_adjtau <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_percomp_adjtau.csv"), header = F,
                                      col.names = c("replicate", "seed", "modelindex", "cov_aZ", "cov_bZ",
                                                    "cov_KZ", "cov_KXZ"))

d_mutvar_adjtau$normalised <- T
d_mutvar$normalised <- F

d_mutvar_percomp_adjtau$normalised <- T
d_mutvar_percomp$normalised <- F

d_mutvar_adjtau <- d_mutvar_adjtau %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_adjtau <- AddCombosToDF(d_mutvar_adjtau)

d_mutvar_percomp_adjtau <- d_mutvar_percomp_adjtau %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_percomp_adjtau <- AddCombosToDF(d_mutvar_percomp_adjtau)

d_mutvar_percomp_adjtau <- d_mutvar_percomp_adjtau %>%
  pivot_longer(cols = starts_with("cov_"), names_to = "molComp", values_to = "cov") %>%
  # Remove massive values
  filter(cov < 1e1, cov > -1e1)


# Join with regular
d_mutvar2 <- full_join(d_mutvar, d_mutvar_adjtau, by = c("replicate", "seed", "modelindex",
                                                         "normalised", "variance", "model",
                                                         "nloci", "tau", "r"))

d_mutvar2 <- d_mutvar2 %>%
  mutate(scaled = if_else(normalised, "Scaled tau", "Unscaled tau"),
         scaled = factor(scaled, levels = c("Unscaled tau", "Scaled tau")))

d_mutvar_sum <- d_mutvar2 %>% filter(tau == 0.0125) %>%
  group_by(model, scaled) %>%
  summarise(meanVar = mean(variance),
            CIVar = CI(variance))


d_mutvar_percomp2 <- full_join(d_mutvar_percomp, d_mutvar_percomp_adjtau, 
                               by = c("replicate", "seed", "modelindex",
                                           "normalised", "cov", "molComp", "model",
                                           "nloci", "tau", "r"))

# Filter to the tau we tested
d_mutvar_percomp2 <- d_mutvar_percomp2 %>%
  filter(tau == 0.0125)

d_mutvar_percomp2 <- d_mutvar_percomp2 %>%
  mutate(scaled = if_else(normalised, "Scaled tau", "Unscaled tau"),
         scaled = factor(scaled, levels = c("Unscaled tau", "Scaled tau")),
         scaledCOV = scale(cov))

d_mutvar_percomp_sum <- d_mutvar_percomp2 %>%
  group_by(model, molComp, scaled) %>%
  summarise(meanCOV = mean(cov),
            CICOV = CI(cov),
            meanAbsCov = mean(abs(cov)),
            meanScaledCOV = mean(scaledCOV),
            CIUScaledCov = CI(scaledCOV))

ggplot(d_mutvar2 %>% filter(tau == 0.0125) %>%
         mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE")), 
       aes(x = model, y = variance, colour = model)) +
  facet_nested(. ~ scaled, space = "free", scales = "free") +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_sum %>% ungroup() %>%
               mutate(model = fct_recode(model, "Additive" = "Add", "K+" = "K", "K-" = "ODE")), 
             aes(x = model, y = meanVar),
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1)) +
  labs(x = "Model", 
       y = "Mutational variance") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_vm
ggsave("plt_vm_scaled.png", device = png, width = 9, height = 5)  

# Covariance
ggplot(d_mutvar_percomp2 %>% filter(cov != 0.0) %>%
         mutate(model = fct_recode(model, "K+" = "K", "K-" = "ODE")), 
       aes(x = model, y = abs(cov), colour = molComp)) +
  facet_nested(. ~ scaled) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_percomp_sum %>% ungroup() %>% filter(meanAbsCov > 0.0) %>%
               mutate(model = fct_recode(model, "K+" = "K", "K-" = "ODE")),
             aes(x = model, y = meanAbsCov, group = molComp),
             position = position_dodge(0.8),
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 4, 
                                           direction = -1),
                      labels = c(TeX("$\\alpha_Z$"),
                                 TeX("$\\beta_Z$"),
                                 TeX("$K_Z$"),
                                 TeX("$K_{XZ}$"))) +
  #coord_cartesian(ylim = c(-0.00001, 0.00001)) +
  scale_y_log10() +
  labs(x = "Model", 
       y = "Absolute mutational\ncovariance (log10)",
       colour = "Molecular component") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom") -> plt_covm
ggsave("plt_covm_scaled_new.png", device = png, width = 9, height = 5)  


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
d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg_adjTau.csv"), header = F, 
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
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) -> plt_adjtau_pheno
ggsave("plt_adapt_mutScale.png", width = 6, height = 5, device = png)


# Plot as a grid
layout <-
"
AABB
#CC#
"
( ( plt_vm + plt_covm ) / ( plt_adjtau_pheno ) ) -> plt_combined
plt_vm + plt_covm + plt_adjtau_pheno -> plt_combined
plt_combined + plot_layout(design = layout) + 
    plot_annotation(tag_levels = 'A')
ggsave("plt_vm_fig.png", device = png, width = 10, height = 7)
