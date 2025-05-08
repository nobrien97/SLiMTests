# Plot heritability estimates
library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
library(ape)
library(tidytree)
library(ggtree)
library(phytools)
library(factoextra)
library(report)
library(matrixcalc)
library(Matrix)
library(DMwR2)
library(mcreplicate)
library(betareg)
library(car)
library(emmeans)
library(nlme)
library(MASS)
library(xtable)
library(stargazer)
library(Rcpp)
library(dendextend)

DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/analysis/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))

model_names <- c("'NAR'", "'PAR'", "'FFLC1'", 
                 "'FFLI1'", "'FFBH'")

select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter

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

h2_colnames <- c("gen", "seed", "modelindex", "VA_w", "h2_w", "VA_aX", "VA_KZX", 
                 "VA_aY", "VA_bY", "VA_KY", "VA_aZ", "VA_bZ", "VA_KZ", "VA_KXZ", 
                 "VA_base", "VA_n", "VA_XMult", "CVA_aX_KZX", "CVA_aX_aY", 
                 "CVA_aX_bY", "CVA_aX_KY", "CVA_aX_aZ", "CVA_aX_bZ", "CVA_aX_KZ", 
                 "CVA_aX_KXZ", "CVA_aX_base", "CVA_aX_n", "CVA_aX_XMult", 
                 "CVA_KZX_aY", "CVA_KZX_bY", "CVA_KZX_KY", "CVA_KZX_aZ", 
                 "CVA_KZX_bZ", "CVA_KZX_KZ", "CVA_KZX_KXZ", "CVA_KZX_base", 
                 "CVA_KZX_n", "CVA_KZX_XMult", "CVA_aY_bY", "CVA_aY_KY", 
                 "CVA_aY_aZ", "CVA_aY_bZ", "CVA_aY_KZ", "CVA_aY_KXZ", 
                 "CVA_aY_base", "CVA_aY_n", "CVA_aY_XMult", "CVA_bY_KY", 
                 "CVA_bY_aZ", "CVA_bY_bZ", "CVA_bY_KZ", "CVA_bY_KXZ", 
                 "CVA_bY_base", "CVA_bY_n", "CVA_bY_XMult", "CVA_KY_aZ", 
                 "CVA_KY_bZ", "CVA_KY_KZ", "CVA_KY_KXZ", "CVA_KY_base", 
                 "CVA_KY_n", "CVA_KY_XMult", "CVA_aZ_bZ", "CVA_aZ_KZ", 
                 "CVA_aZ_KXZ", "CVA_aZ_base", "CVA_aZ_n", "CVA_aZ_XMult", 
                 "CVA_bZ_KZ", "CVA_bZ_KXZ", "CVA_bZ_base", "CVA_bZ_n", 
                 "CVA_bZ_XMult", "CVA_KZ_KXZ", "CVA_KZ_base", "CVA_KZ_n", 
                 "CVA_KZ_XMult", "CVA_KXZ_base", "CVA_KXZ_n", "CVA_KXZ_XMult", 
                 "CVA_base_n", "CVA_base_XMult", "CVA_n_XMult", "h2_aX", "h2_KZX", "h2_aY", "h2_bY", 
                 "h2_KY", "h2_aZ", "h2_bZ", "h2_KZ", "h2_KXZ", "h2_base", "h2_n", 
                 "h2_XMult")

d_h2_mrr <- data.table::fread(paste0(DATA_PATH, "out_h2_mrr.csv"), header = F, 
                              sep = ",",
                              col.names = h2_colnames, fill = T)
d_h2_mkr <- data.table::fread(paste0(DATA_PATH, "out_h2_mkr.csv"), header = F, 
                              sep = ",",
                              col.names = h2_colnames, fill = T)

# join
d_h2_mkr$calcMode <- "mkr"
d_h2_mrr$calcMode <- "mrr"

d_h2 <- rbind(d_h2_mkr, d_h2_mrr)

# Remove duplicates
d_h2 %>% distinct() -> d_h2

# Remove rows with invalid estimates
d_h2 <- d_h2 %>%
  filter(VA_w >= 0)


# Add our variables
d_combos <- read_delim('/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/combos.csv', 
                     delim = " ", col_names = F)
names(d_combos) <- c("model", "r")

d_h2 %>% mutate(model = d_combos$model[.$modelindex],
                r = d_combos$r[.$modelindex]) -> d_h2

d_h2 |>
  dplyr::summarise(n = dplyr::n(), .by = c(gen, seed, modelindex, model, r, calcMode)) |>
  dplyr::filter(n > 1L)

d_h2 <- d_h2 %>%
  mutate(model = factor(model, levels = model_names))

ggplot(d_h2 %>% distinct(gen, seed, modelindex, calcMode, .keep_all = T) %>% 
         select(gen, seed, modelindex, model, r, h2_w, calcMode) %>%
         distinct() %>%
         pivot_wider(names_from = calcMode, values_from = h2_w), 
       aes(x = mkr, y = mrr, colour = model)) +
  facet_wrap(.~model) +
  geom_point(shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, direction = -1)) +
  labs(x = TeX("Kernel regression heritability $(h^2)$"), 
       y = TeX("Ridge regression heritability $(h^2)$"),
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_h2 %>% distinct(gen, seed, modelindex, calcMode, .keep_all = T) %>% 
         select(gen, seed, modelindex, model, r, VA_w, calcMode) %>%
         distinct() %>%
         pivot_wider(names_from = calcMode, values_from = VA_w), 
       aes(x = mkr, y = mrr, colour = model)) +
  facet_wrap(.~model) +
  geom_point(shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = TeX("Kernel regression additive variance $(VA)$"), 
       y = TeX("Ridge regression additive variance $(VA)$"),
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 14))


# The two estimates are very different: ridge regression is very biased towards
# low heritability, the kernel method is more evenly spread
# Kernel regression seems most accurate: ridge is almost always at 0 heritability

d_h2 <- d_h2 %>% filter(calcMode == "mkr")

boxplot(d_h2$VA_w)

# Attach quant gen data
d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 29)), 
                          col.names = c("gen", "seed", "modelindex", "meanH",
                                        "trait1_mean", "trait2_mean", "trait3_mean",
                                        "trait4_mean", "trait1_var", "trait2_var", 
                                        "trait3_var", "trait4_var", "dist", 
                                        "dist1", "dist2", "dist3", "dist4", "mean_w",
                                        "var_w", "deltaPheno", "deltaW", 
                                        "meanMC1", "meanMC2", "meanMC3", "meanMC4", 
                                        "meanMC5", "meanMC6", "meanMC7", "meanMC8", 
                                        "meanMC9", "meanMC10", "meanMC11"), 
                          fill = T)

# Add predictors
d_qg <- AddCombosToDF(d_qg) 

# Optimum: fitness > 95%
d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.95)) %>%
  mutate(model = factor(model, levels = model_names)) %>%
  ungroup() -> d_qg

d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, isAdapted) %>% filter(gen >= 49500)


d_h2 <- d_h2 %>%
  distinct(gen, seed, modelindex, calcMode, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_w) %>% distinct()

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex"))


# Counts for each model type:
table(d_h2$model, d_h2$isAdapted)

# Discretise generation
d_h2 <- d_h2 %>%
  mutate(timePoint = if_else(gen == 50000, "Start", "End"),
         timePoint = factor(timePoint, levels = c("Start", "End")))

# summarise
d_h2_sum <- d_h2 %>% 
  group_by(timePoint, model, r, isAdapted) %>%
  dplyr::summarise(meanH2w = mean(h2_w, na.rm = T),
            seH2w = se(h2_w, na.rm = T),
            meanVAw = mean(VA_w, na.rm = T),
            seVAw = se(VA_w, na.rm = T))
d_h2_sum$model <- as.factor(d_h2_sum$model)

# Heritability distribution
ggplot(d_h2 %>% 
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?"),
       aes(x = timePoint, y = h2_w, colour = model)) +
  facet_nested(r_title + log10(r) ~ adapted_title + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = meanH2w, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Narrow-sense heritability $(h^2)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Additive variance
# Small effects as separate figure
ggplot(d_h2 %>% 
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?"),
       aes(x = timePoint, y = log10(VA_w), colour = model)) +
  facet_nested(r_title + log10(r) ~ adapted_title + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_sum %>% ungroup() %>% 
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = timePoint, y = log10(meanVAw), group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Time point", 
       y = TeX("Log additive variance in fitness ($log_{10}(VA)$)"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names,
                      guide = guide_legend(override.aes = list(shape = 16,
                                                               size = 3))) +
  #coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_va_together
plt_va_together

ggsave("plt_va.png", device = png, bg = "white",
       width = 12, height = 6)

# VA per molecular component
d_h2_molcomp <- d_h2 %>%
  pivot_longer(starts_with("VA_"), names_to = "molComp", values_to = "VA") %>%
  dplyr::select(gen, seed, modelindex, model, r, molComp, VA, isAdapted, timePoint)
  
d_h2_molcomp$molComp <- str_replace_all(d_h2_molcomp$molComp, "VA_", "")
d_h2_molcomp <- d_h2_molcomp %>%
  filter(molComp != "w") # fitness isn't a mol comp :)

d_h2_molcomp_sum <- d_h2_molcomp %>% 
  group_by(timePoint, model, r, molComp, isAdapted) %>%
  dplyr::summarise(meanVA = mean(VA, na.rm = T),
                   seVA = se(VA, na.rm = T))


ggplot(d_h2_molcomp %>%
         filter(timePoint == "End") %>%
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?") %>%
         mutate(model = str_replace_all(model, "'", "")),
       aes(x = molComp, y = VA, colour = as.factor(r))) +
  facet_nested("Model" + model ~ adapted_title + isAdapted) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = T) +
  geom_point(data = d_h2_molcomp_sum %>% ungroup() %>% filter(timePoint == "End") %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?") %>%
               mutate(model = str_replace_all(model, "'", "")),
             aes(x = molComp, y = meanVA, group = as.factor(r)), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Molecular component", 
       y = TeX("Additive variance $(VA)$"),
       colour = "Recombination rate (log10)") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3),
                      labels = c("-10", "-5", "-1")) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_bw() +
  theme(text = element_text(size = 10),
        legend.position = "bottom")
ggsave("plt_va_percomp.png", device = png, bg = "white",
       width = 560*6, height = (980*4)/3, units = "px")

# What are the identities of the FFBH models which adapted and had high aZ? 
# What does their trajectory look like compare to the other ones?
# High aZ
View(d_h2_molcomp %>%
  filter(model == "'FFBH'", timePoint == "End", molComp == "aZ",
         VA > 0.5))
# Low aZ
View(d_h2_molcomp %>%
       filter(model == "'FFBH'", timePoint == "End", molComp == "aZ",
              VA < 0.25))

# Try to group these up so we can plot their trajectories
d_qg$simID <- interaction(d_qg$seed, d_qg$modelindex)
d_h2_molcomp$simID <- interaction(d_h2_molcomp$seed, d_h2_molcomp$modelindex)

d_h2_ffbh_aZ <- d_h2_molcomp %>%
  mutate(adaptedWithAZ = "Maladapted") %>%
  filter(model == "'FFBH'", timePoint == "End", molComp == "aZ", isAdapted == T) %>%
  drop_na(VA)
d_h2_ffbh_aZ[(d_h2_ffbh_aZ$VA > 0.5),]$adaptedWithAZ <- "With aZ"
d_h2_ffbh_aZ[(d_h2_ffbh_aZ$VA < 0.25),]$adaptedWithAZ <- "Without aZ"

d_qg_ffbh_aZ <- d_qg %>%
  filter(model == "'FFBH'")

# Now attach these labels to the qg
d_qg_ffbh_aZ$adaptedWithAZ <- "Maladapted"
for (i in unique(d_h2_ffbh_aZ$simID)) {
  d_qg_ffbh_aZ[d_qg_ffbh_aZ$simID == i,]$adaptedWithAZ <-
    d_h2_ffbh_aZ[d_h2_ffbh_aZ$simID == i,]$adaptedWithAZ
}

# Plot adaptive walk for some example models
ggplot(d_qg %>% filter((seed == 3245990319 & modelindex == 15) | # Adapted with aZ
                         (seed == 2954984637 & modelindex == 15) | # Adapted without aZ
                         seed == 1049021301 & modelindex == 15) %>% # Maladapted
       mutate(seed = factor(seed, levels = c(3245990319, 2954984637, 1049021301))) %>%
       filter(gen > 40000) %>%
         mutate(gen = gen - 50000),
       aes(x = gen, y = mean_w, colour = seed, group = seed)) +
  geom_line() +
  labs(x = "Generations post-optimum shift", y = "Mean population fitness", 
       colour = "Simulation outcome") +
  scale_colour_manual(values = c(paletteer_d("nationalparkcolors::Everglades", 2), "#666"), 
                      labels = c("Adapted (high aZ)", 
                                 "Adapted (low aZ)",
                                 "Maladapted")) +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_example_sims.png", device = png, bg = "white",
       width = 9, height = 4.5)

# Now for all of them
# Plot adaptive walk for these models
maladapted_examples <- 
  sample(unique(d_qg_ffbh_aZ[!d_qg_ffbh_aZ$isAdapted,]$simID),
         5)

ggplot(d_qg_ffbh_aZ %>% filter(adaptedWithAZ != "Maladapted" | 
                                 simID %in% maladapted_examples) %>%
         filter(gen > 40000) %>%
         mutate(gen = gen - 50000),
       aes(x = gen, y = mean_w, colour = adaptedWithAZ, group = simID)) +
  geom_line() +
  labs(x = "Generations post-optimum shift", y = "Mean population fitness", 
       colour = "Simulation outcome") +
  scale_colour_manual(values = c("#666", paletteer_d("nationalparkcolors::Everglades", 2))) +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_example_sims.png", device = png, bg = "white",
       width = 9, height = 4.5)

# Are these models adapting because of relaxed constraints? Which traits are optimised?
# load in optima
d_opt <- data.table::fread(paste0(DATA_PATH, "slim_opt.csv"), header = F, 
                          sep = ",", colClasses = c("factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          fill = T)
colnames(d_opt)[1:2] <- c("seed", "modelindex")
d_opt <- AddCombosToDF(d_opt)
d_opt$simID <- interaction(d_opt$seed, d_opt$modelindex)

# Select replicates
adapted_examples <- unique((d_qg_ffbh_aZ %>% filter(isAdapted == T))$simID)
View(d_opt %>% filter(simID %in% maladapted_examples))

# Average
d_qg_ffbh_sum <- d_qg_ffbh_aZ %>% 
filter(model == "FFBH", adaptedWithAZ != "Maladapted" | 
               simID %in% maladapted_examples) %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, adaptedWithAZ) %>%
  summarise(meanFitness = mean(mean_w),
            SEFitness = se(mean_w),
            meanFitnessVar = mean(var_w),
            SEFitnessVar = se(var_w))

ggplot(d_qg_ffbh_sum %>% 
         filter(gen > -10500),
       aes(x = gen, y = meanFitness, colour = adaptedWithAZ)) +
  geom_line() +
  geom_ribbon(aes(ymin = meanFitness - SEFitness, 
                  ymax = meanFitness + SEFitness, fill = adaptedWithAZ), colour = NA,
              alpha = 0.2) +
  labs(x = "Generations post-optimum shift", y = "Mean population fitness", 
       colour = "Simulation outcome") +
  scale_colour_manual(values = c("#666", paletteer_d("nationalparkcolors::Everglades", 2))) +
  scale_fill_manual(values = c("#666", paletteer_d("nationalparkcolors::Everglades", 2)),
                    guide = "none") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) 
ggsave("plt_av_ffbh_sims.png", device = png, bg = "white",
       width = 9, height = 4.5)

# So no real difference between the two groups in terms of the adaptive walks




# Plot relative additive variance? What proportion of total is each contributing?

# Infinitesimal expects zero change in additive variance due to selection
# So see how much it changes between timepoints
# Scale by the total variance as well -> a large effect model will produce a lot
# of variance, so the differences are more likely to be greater
# Also account for drift: estimate Ne via Hossjer et al.

# Total change from start to end
d_h2 %>%
  group_by(model, seed, r, isAdapted) %>%
  filter(n() > 1) %>%
  summarise(totalDeltaVA = sum(diff(VA_w))) -> d_h2_deltaVA

# Change per molecular component
d_h2 %>%
  pivot_longer(starts_with("VA_"), names_to = "molComp", values_to = "VA") %>%
  group_by(model, seed, r, molComp, isAdapted) %>%
  filter(n() > 1) %>%
  summarise(totalDeltaVA = sum(diff(VA))) -> d_h2_deltaVA_molComp
boxplot(d_h2_deltaVA_molComp$totalDeltaVA)

# total distribution
boxplot(d_h2_deltaVA$totalDeltaVA)

d_h2_deltaVA %>%
  group_by(model, r, isAdapted) %>%
  summarise(meanDeltaVA = mean(totalDeltaVA, na.rm = T),
            seDeltaVA = se(totalDeltaVA, na.rm = T)) -> d_h2_deltaVA_sum

d_h2_deltaVA_molComp %>%
  group_by(model, r, molComp, isAdapted) %>%
  summarise(meanDeltaVA = mean(totalDeltaVA, na.rm = T),
            seDeltaVA = se(totalDeltaVA, na.rm = T)) -> d_h2_deltaVA_molComp_sum

ggplot(d_h2_deltaVA %>% 
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?"),
       aes(x = model, y = totalDeltaVA, colour = model)) +
  facet_nested(r_title + log10(r) ~ adapted_title + isAdapted) +
  geom_quasirandom(dodge.width = 0.9) +
  #coord_cartesian(ylim = c(0, 1)) +
  geom_point(data = d_h2_deltaVA_sum %>% 
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = model, y = meanDeltaVA, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Model", 
       y = TeX("Change in additive variance $(\\Delta V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5),
                      guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
ggsave("plt_deltaVA.png", device = png, width = 9, height = 4)

ggplot(d_h2_deltaVA_molComp %>% 
         mutate(r_title = "Recombination rate (log10)",
                adapted_title = "Did the population adapt?"),
       aes(x = molComp, y = totalDeltaVA, colour = model)) +
  facet_nested(r_title + log10(r) ~ adapted_title + isAdapted) +
  geom_quasirandom(dodge.width = 0.9) +
  #coord_cartesian(ylim = c(0, 1)) +
  geom_point(data = d_h2_deltaVA_molComp_sum %>% filter(!is.nan(meanDeltaVA)) %>%
               mutate(r_title = "Recombination rate (log10)",
                      adapted_title = "Did the population adapt?"),
             aes(x = molComp, y = meanDeltaVA, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Molecular component", 
       y = TeX("Change in additive variance $(\\Delta V_A)$"),
       colour = "Model") +
  #scale_x_discrete(labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
ggsave("plt_deltaVA_molComp.png", device = png, width = 9, height = 4)

# did the mean VA change with changing recombination and between adapted/nonadapted?
library(nlme)
summary(gls.dva <- gls(totalDeltaVA ~ model + as.factor(r) + isAdapted, 
                       d_h2_deltaVA, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.dva)
# No difference between models in the change in VA between start and end


d_h2 %>% filter(isAdapted) %>%
  select(!VA_w) %>%  # Remove fitness (since its a different measurement)
  filter(rowSums(is.na(select(., 6:83))) < (83 - 6 + 1)) %>%  # Drop rows with no variance
  group_by(modelindex, timePoint, isAdapted) %>%
  group_split(.) -> split_h2


# Separate into model indices
# each sublist is replicates of a model index
sourceCpp("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R/getCovarianceMatrices.cpp")
lapply(split_h2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices
lapply(split_h2, function(x) {data.frame(timePoint = x$timePoint, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex

# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat <- unlist(cov_matrices, recursive = F)

# get ids from the matrix
cov_matrix_modelindex <- GetMatrixIDs(split_h2)


# Distance between G matrices
# Analysis across all timepoints, timepoint doesn't affect the tree structure
sourceCpp("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R/distanceFunctions.cpp")

dist_matrix <- distanceMatrix(h2_mat)
colnames(dist_matrix) <- paste("Matrix", 1:nrow(dist_matrix))
rownames(dist_matrix) <- colnames(dist_matrix)

# Clustering based on distance
hc <- hclust(as.dist(dist_matrix), method="average")
plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 5 seems to be the best
# elbow plot
fviz_nbclust(dist_matrix, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
png(file = "dendrogram_totaldist.png",
    height = 2250/2, width = 2250, units = "px")
par(lwd=2, mar=c(8,8,8,8))
plot(hc, main = "Power Euclidean distances between molecular G matrices", labels = F,
     cex.main = 2, cex.lab = 2, xlab = "", sub = "", axes = F)
rect.hclust(hc, 5, border = 2)
axis(2, lwd = 2, cex.axis = 2)
dev.off()


clus <- cutree(hc, 5)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)

id <- rbindlist(cov_matrix_modelindex, 
                fill = T)
id$label <- as.character(1:nrow(id))
id$modelindex <- as.factor(id$modelindex)
id <- AddCombosToDF(id)
id$model <- factor(id$model, levels = model_names)
id$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id$clus, id$timePoint, id$model, id$r)
names(dimnames(tab)) <- c("cluster", "timePoint", "model", "r")
tab <- as.data.frame(tab)

# Model describes the clustering
glm.clus <- glm(Freq~model+r,family=poisson(),data=tab)
summary(glm.clus)
report::report(glm.clus)

# Models and recombination rate matter, but not the time point and no interaction
# All models are pretty similar

id %>% ungroup() %>%
  group_by(model, timePoint, r, clus) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  mutate(prop = n/sum(n)) -> cluster_percs

phylo <- full_join(as.phylo(phylo), id, by = "label")

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(size = 2) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH", "")) +
  labs(colour = "Model") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(shape=16, size = 5,
                                                   linetype = 0))) -> tree_full

tree_full
ggsave("plt_tree_gmatrix_full.png", device = png, bg = "white",
       width = 7/2, height = 9/2)


# Evolvability metrics

# First convert to nearest positive definite matrix
h2_pd <- lapply(h2_mat, function(x) {
  if (!is.positive.definite(x)) {as.matrix(nearPD(x)$mat)}
})


d_ecr <- CalcECRA(h2_pd, id)
d_ecr <- AddCombosToDF(d_ecr)

# Refactor model
d_ecr <- d_ecr %>%
  mutate(model = factor(model, levels = model_names))

# Need to calculate cev means separately for the different models
# K- shouldn't mean over cev_KZ and KXZ
d_ecr_sum <- d_ecr %>%
  group_by(model, r) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = log10(cev), colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = log10(cev_mean), group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Model", y = "Mean conditional\nevolvability (log10)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = log10(res), colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = log10(res_mean), group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Model", y = "Mean respondability (log10)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = aut, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = aut_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Model", y = "Mean autonomy",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = model, y = log10(ev), colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = log10(ev_mean), group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5),
                      labels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), 
                      breaks = model_names) +
  labs(x = "Model", y = "Mean evolvability (log10)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "none", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                      plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO", label_size = 12)

plt_evol <- plot_grid(plt_evol,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol.png", device = png, bg = "white",
       width = 10, height = 7)

# Compare models among recombination rates 
# Variance differs between groups, use gls to account for unequal variance
# timepoint doesn't matter
library(nlme)
summary(gls.cev <- gls(cev ~ model + as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.cev)
report(gls.cev)

anova(gls.cev)

# Marginal means
library(emmeans)
library(xtable)
em.cev <- emmeans(gls.cev, ~ model * r)
pairs(em.cev, simple = "model")
pairs(em.cev, simple = "r")
plot(em.cev, comparisons = T)

xtable(em.cev)


# Format to a nice table
library(stargazer)
stargazer(gls.cev)

# Repeat for autonomy, respondability, evolvability
summary(gls.aut <- gls(aut ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))

anova(gls.aut)
plot(gls.aut)
em.aut <- emmeans(gls.aut, ~ model * r)
pairs(em.aut, simple = "model")
pairs(em.aut, simple = "r")
plot(em.aut)

xtable(em.aut)

stargazer(gls.aut)

summary(gls.res <- gls(res ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.res)
anova(gls.res)

em.res <- emmeans(gls.res, ~ model * r)
pairs(em.res, simple = "model")
pairs(em.res, simple = "r")
plot(em.res, comparisons = T)
stargazer(gls.res)

xtable(em.res)


summary(gls.ev <- gls(ev ~ model * as.factor(r), d_ecr, 
                      weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.ev)
anova(gls.ev)
em.ev <- emmeans(gls.ev, ~ model * r, type = "response")
pairs(em.ev, simple = "model")
pairs(em.ev, simple = "r")
plot(em.ev, comparisons = T)

stargazer(gls.ev)
xtable(em.ev)

##############################
# Comparing the shape of G matrices with PCA similarity

# Bootstrap PCA similarity factor test:
# sample two matrices at random within group, get PCA similarity
# compare mean similarity between groups
# data input: dataframe with ids and a column with the matrix

# Each model has different molecular components, so really only makes sense to compare
# within models

krz_in <- id %>%
  mutate(g = h2_pd,
         group = interaction(model, log10(r)))

krz_in_timePoint <- id %>%
  mutate(g = h2_pd,
         group = interaction(log10(r), timePoint))

# Remove null matrices (no nearest matrix found)
krz_in <- krz_in[!sapply(krz_in$g,is.null)];
krz_in_timePoint <- krz_in_timePoint[!sapply(krz_in_timePoint$g,is.null)];

# Save krz_in: run this part on HPC
saveRDS(krz_in, "pca_in.RDS")
saveRDS(krz_in_timePoint, "pca_tp_in.RDS")

# Bootstrap in ten parts for RAM reasons
# This is slow: uncomment to run, otherwise read in precalculated data
# Generate seeds
#newseed <- sample(1:.Machine$integer.max, 10)
# 407844323 2049133531  970639651  452738391 1161959903  506461634 2087592727 1740805001
# 1698460455  886904095
newseed <- c(407844323L, 2049133531L, 970639651L, 452738391L, 1161959903L, 
             506461634L, 2087592727L, 1740805001L, 1698460455L, 886904095L)
bootPCASim <- vector(mode = "list", length = length(newseed))

# Per model inputs
krz_in_timePoint_NAR <- krz_in_timePoint %>% filter(model == "'NAR'")
krz_in_timePoint_PAR <- krz_in_timePoint %>% filter(model == "'PAR'")
krz_in_timePoint_FFLC1 <- krz_in_timePoint %>% filter(model == "'FFLC1'")
krz_in_timePoint_FFLI1 <- krz_in_timePoint %>% filter(model == "'FFLI1'")
krz_in_timePoint_FFBH <- krz_in_timePoint %>% filter(model == "'FFBH'")

  
for (i in seq_along(newseed)) {
  # Set seed
  set.seed(newseed[i])
  # Run replicate but only within models
  res_NAR <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_timePoint_NAR, "group", T))
  res_PAR <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_timePoint_PAR, "group", T))
  res_FFLC1 <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_timePoint_FFLC1, "group", T))
  res_FFLI1 <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_timePoint_FFLI1, "group", T))
  res_FFBH <- mcreplicate::mc_replicate(50, bootKrzCorFn(krz_in_timePoint_FFBH, "group", T))
  
  # To data.frame
  res_NAR <- unnest(as.data.frame(t(res_NAR)), cols = everything()) %>%
    mutate(model = "'NAR'")
  res_PAR <- unnest(as.data.frame(t(res_PAR)), cols = everything()) %>%
    mutate(model = "'PAR'")
  res_FFLC1 <- unnest(as.data.frame(t(res_FFLC1)), cols = everything()) %>%
    mutate(model = "'FFLC1'")
  res_FFLI1 <- unnest(as.data.frame(t(res_FFLI1)), cols = everything()) %>%
    mutate(model = "'FFLI1'")
  res_FFBH <- unnest(as.data.frame(t(res_FFBH)), cols = everything()) %>%
    mutate(model = "'FFBH'")
  
  # Combine to output
  bootPCASim[[i]] <- rbind(res_NAR, res_PAR, res_FFLC1, res_FFLI1, res_FFBH)
}

tpLevels <- c("Start", "End")

# Output list into combined df
bootPCASim2 <- bind_rows(bootPCASim)
bootPCASim <- bootPCASim2 %>%
  separate(group1, c("r1", "timePoint1"), "\\.",
           extra = "merge") %>%
  separate(group2, c("r2", "timePoint2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = as.numeric(r1),
         r2 = as.numeric(r2),
         timePoint1 = factor(timePoint1, levels = tpLevels),
         timePoint2 = factor(timePoint2, levels = tpLevels)) %>%
  rename(PCASim = krzCor)

saveRDS(bootPCASim, paste0(DATA_PATH, "d_bootPCASim_timepoints.RDS"))
bootPCASim <- readRDS(paste0(DATA_PATH, "d_bootPCASim_timepoints.RDS"))

# Plot
# Split by model comparison

GetModelComparison <- function(model1, model2, model_names) {
  result <- character(length(model1))
  # assign by priority (NAR > PAR > FFLC1 > FFLI1 > FFBH)
  model1_index <- match(model1, model_names)
  model2_index <- match(model2, model_names)
  
  if (any(is.na(model1_index)) | any(is.na(model2_index))) {
    return(NA)
  }
  
  # If model2's index is greater than model1,
  # we put model1 first
  model1First_idx <- model2_index > model1_index
  model2First_idx <- model1_index >= model2_index
  result[model2_index > model1_index] <- paste(model1[model1First_idx], "vs", model2[model1First_idx])
  # otherwise model2 comes first, if models are the same these are
  # included here too
  result[model1_index >= model2_index] <- paste(model2[model2First_idx], "vs", model1[model2First_idx])
  return(result)
}

bootPCASim <- bootPCASim %>%
  mutate(#tpCombo = GetModelComparison(timePoint1, timePoint2, tpLevels),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)))

# recomb by modelCombo
bootPCASim_sum <- bootPCASim %>%
  group_by(r1, r2, model) %>%
  summarise(meanPCASim = mean(PCASim),
            ciPCASim = CI(PCASim))


ggplot(bootPCASim_sum, aes(
  x = as.factor(r1), y = as.factor(r2)
)) +
  facet_nested("Model" + model ~ .) + 
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 10), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_r.png", device = png, width = 15, height = 4)

# beta regression
# Distributions
ggplot(bootPCASim, 
       aes(x = rCombo, y = PCASim)) +
  geom_quasirandom(dodge.width = 0.9)

plot(density(bootPCASim$PCASim))

# Definitely different between models, none are normally distributed

# Floating point error: clamp to 1
bootPCASim <- bootPCASim %>%
  mutate(PCASim = raster::clamp(PCASim, 0, 1))

# Run regression: this is slow! Uncomment to run, otherwise load in saved object
# br.pcasim <- betareg(PCASim ~ modelCombo * as.factor(rCombo) | modelCombo * as.factor(rCombo), 
#                      bootPCASim)

# Save output
# saveRDS(br.pcasim, "betareg_pcaSim_big.RDS")
br.pcasim <- readRDS(paste0(DATA_PATH, "betareg_pcaSim_big.RDS"))
summary(br.pcasim)
plot(br.pcasim)

car::Anova(br.pcasim, type = 3, test.statistic = "F")

# Response is on log-odds scale
em.pcasim <- emmeans(br.pcasim, ~modelCombo * as.factor(rCombo))
em.pcasim <- regrid(em.pcasim)
pairs(em.pcasim, simple = "modelCombo")
pairs(em.pcasim, simple = "rCombo")
plot(em.pcasim, comparisons = T)
pwpp(em.pcasim, by = "modelCombo", type = "response")
emmip(br.pcasim,  ~ modelCombo | rCombo)
joint_tests(em.pcasim)
xtable(em.pcasim)

# fractional logit: easier to get confidence intervals in response scale, 
# similar results to betareg
fl.pcasim <- glm(PCASim ~ modelCombo * as.factor(rCombo),
                 data = bootPCASim,
                 family = quasibinomial())
summary(fl.pcasim)
plot(fl.pcasim)
em.fl.pcasim <- emmeans(fl.pcasim, ~modelCombo * as.factor(rCombo))
em.fl.pcasim <- regrid(em.fl.pcasim)
pairs(em.fl.pcasim, simple = "modelCombo")
pairs(em.fl.pcasim, simple = "rCombo")
pwpp(em.fl.pcasim, by = "modelCombo")
emmip(em.fl.pcasim,  ~ modelCombo | rCombo)
xtable(em.fl.pcasim)

