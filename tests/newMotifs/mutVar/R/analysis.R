library(tidyverse)
library(ggbeeswarm)
library(ggh4x)
library(paletteer)
library(Rcpp)

# Helper functions
AddCombosToDF <- function(df, combos) {
  df %>% ungroup() %>%
    mutate(model = combos$model[as.numeric(levels(modelindex))[modelindex]],
           r = combos$r[as.numeric(levels(modelindex))[modelindex]])
}

se <- function(x, na.rm = F) {
  return(sd(x, na.rm) / sqrt(length(x)))
}

# Mutational variance
DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/mutVar/R/"
COMBOS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/"
ANALYSIS_PATH <- paste0(REPO_PATH, "analysis/")

# Read data: pre-adjusted simulations (nloci = 1024, tau = 0.0125)
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                              col.names = c("replicate", "seed", "modelindex", 
                                            "t1_m", "t2_m", "t3_m", "t4_m", 
                                            "t1_v", "t2_v", "t3_v", "t4_v",
                                            "t1_t2_c", "t1_t3_c", "t1_t4_c",
                                            "t2_t3_c", "t2_t4_c", "t3_t4_c"))

d_mutvar_percomp <- data.table::fread(paste0(DATA_PATH, "slim_mutvar_percomp.csv"), header = F,
                                      col.names = c("replicate", "seed", "modelindex", 
                                                    paste0("t1_m", 1:11),
                                                    paste0("t2_m", 1:11),
                                                    paste0("t3_m", 1:11),
                                                    paste0("t4_m", 1:11)))

# Load combos file
d_combos <- data.table::fread(paste0(COMBOS_PATH, "combos.csv"), header = F,
                              col.names = c("model", "r"))




d_mutvar <- d_mutvar %>%
  mutate(replicate = replicate %/% 2 + 1, # convert replicate id from 1 3 5 to 1 2 3
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar <- AddCombosToDF(d_mutvar, d_combos)

d_mutvar_plt <- d_mutvar %>%
  select(replicate, seed, modelindex, model, r, ends_with("v")) %>%
  pivot_longer(cols = ends_with("v"), names_to = "trait", names_pattern = "t([0-9]*)_v",
               values_to = "variance") %>%
  filter(!is.nan(variance))

d_mutvar_plt2 <- d_mutvar %>%
  select(replicate, seed, modelindex, model, r, ends_with("m")) %>%
  pivot_longer(cols = ends_with("m"), names_to = "trait", names_pattern = "t([0-9]*)_m",
               values_to = "mean")

d_mutvar_plt <- full_join(d_mutvar_plt, d_mutvar_plt2, by = 
                            c("replicate", "seed", "modelindex", "model", "r", "trait"))
rm(d_mutvar_plt2)


d_mutvar_plt <- d_mutvar_plt %>%
  mutate(CVVar = sqrt(variance) / mean)


d_mutvar_sum <- d_mutvar_plt %>%
  filter(variance > 1e-10) %>%
  group_by(model, trait) %>%
  summarise(meanVar = mean(variance),
            meanCVVar = mean(CVVar))

# Plot variance
ggplot(d_mutvar_plt %>% filter(variance > 1e-10),
       aes(x = model, y = CVVar, colour = model)) +
  geom_quasirandom(dodge.width = 0.8) +
  facet_nested("Trait" + trait ~ ., scales = "free") +
  geom_point(data = d_mutvar_sum %>% ungroup(),
             aes(x = model, y = meanCVVar),
             shape = 21, size = 2, fill = "white",
             colour = "black", inherit.aes = F,
             stroke = 1) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5,
                                           direction = -1)) +
  scale_x_discrete(guide = "axis_nested") +
  labs(x = "Model",
       y = "Mutational variance (mean coefficient of variance") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none") -> plt_vm
plt_vm
ggsave("plt_vm_scaled.png", device = png, width = 6, height = 5)

# Split to M matrices
d_m <- d_mutvar %>% select(8:17)
sourceCpp("getCovarianceMatrices.cpp")
extractCovarianceMatrices(d_m) -> m_matrices




d_mutvar_percomp <- d_mutvar_percomp %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar_percomp <- AddCombosToDF(d_mutvar_percomp)

d_mutvar_percomp <- d_mutvar_percomp %>%
  pivot_longer(cols = starts_with("cov_"), names_to = "molComp", values_to = "cov") %>%
  # Clean data
  filter(cov < 1e1, cov > -1e1)




# Covariance
ggplot(d_mutvar_percomp2 %>% filter(cov != 0.0) %>%
         mutate(model = fct_recode(model, "K+" = "K", "K-" = "ODE")),
       aes(x = molComp, y = abs(cov), colour = model)) +
  facet_nested(. ~ scaled, labeller = labeller(scaled = label_parsed)) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_percomp_sum %>% ungroup() %>% filter(meanAbsCov > 0.0) %>%
               mutate(model = fct_recode(model, "K+" = "K", "K-" = "ODE")),
             aes(x = molComp, y = meanAbsCov, group = model),
             position = position_dodge(0.8),
             fill = "white", stroke = 1,
             shape = 21, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 4,
                                           direction = 1)) +
  #coord_cartesian(ylim = c(-0.00001, 0.00001)) +
  scale_x_discrete(labels = c(TeX("$\\alpha_Z$"),
                              TeX("$\\beta_Z$"),
                              TeX("$K_Z$"),
                              TeX("$K_{XZ}$"))) +
  scale_y_log10() +
  labs(x = "Molecular component",
       y = "Absolute mutational\ncovariance (log10)",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom") -> plt_covm
plt_covm
ggsave("plt_covm_scaled_new.png", device = png, width = 9, height = 5)

# Plot adaptive trajectory

# load trait evolution data
d_qg_adjTau <- data.table::fread(paste0(DATA_PATH, "slim_qg_adjTau.csv"), header = F,
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    rep("numeric", times = 12)),
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"),
                          fill = T)

d_qg_adjTau %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg_adjTau

# Attach additive replicates as well
d_qg_add <- AddCombosToDF(d_qg) %>% filter(model == "Add", r %in% r_subsample,
                                nloci == 1024, tau == 0.0125)
# join
d_adapted_adjTau <- full_join(d_qg_adjTau, d_qg_add)

d_adapted_adjTau <- AddCombosToDF(d_adapted_adjTau)

d_adapted_adjTau_sum <- d_adapted_adjTau %>%
  filter(isAdapted, gen >= 49500) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, r) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

ggplot(d_adapted_adjTau_sum,
       aes(x = gen, y = meanPhenomean, colour = model)) +
  facet_grid(log10(r)~.) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_ribbon(aes(ymin = meanPhenomean - sdPhenomean,
                  ymax = meanPhenomean + sdPhenomean, fill = model), colour = NA,
              alpha = 0.2) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)",
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                    labels = c("Additive", "K+", "K-"), guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14),
        panel.spacing = unit(0.75, "lines")) -> plt_adapt_adjTau
plt_adapt_adjTau
ggsave("plt_adapt_mutScale.png", width = 6, height = 5, device = png)

# Combine figures
# Plot as a grid
layout <-
  "
AAC
BBC
"
plt_vm + plt_covm + plt_adapt_adjTau -> plt_combined
plt_combined + patchwork::plot_layout(design = layout) +
  patchwork::plot_annotation(tag_levels = 'A') -> plt_vm_final

ggsave("plt_vm_final.png", plt_vm_final, width = 10, height = 8, device = png, bg = "white")

