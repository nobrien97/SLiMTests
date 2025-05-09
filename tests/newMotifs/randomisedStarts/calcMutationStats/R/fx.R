# Analysis of effect sizes for each molecular component
library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
library(betareg)
library(emmeans)
library(nlme)
library(MASS)
library(xtable)

DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/"

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

# Use the right mutate/summarise functions
mutate <- dplyr::mutate
summarise <- dplyr::summarise
select <- dplyr::select

d_combos <- read.table("../../../R/combos.csv", header = F,
                       col.names = c("model", "r"))

model_levels <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
model_labels <- c("NAR", "PAR", "FFL-C1", "FFL-I1", "FFBH")

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
d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.95)) %>%
  mutate(model = factor(model, levels = model_levels)) %>%
  ungroup() -> d_qg

d_qg <- d_qg %>% filter(gen >= 49500)

# Fixations
d_fixations <- data.table::fread(paste0(DATA_PATH, "calcMutationStats/d_fixed.csv"),
                                 header = F, sep = ",", colClasses = c("integer",
                                 "factor", "factor", "integer"),
                                 col.names = c("gen", "modelindex", "seed", "nFixed"),
                                 fill = T)


d_fixations <- left_join(d_fixations, d_qg, by = c("gen", "seed", "modelindex"))

d_fixations <- d_fixations %>%
  mutate(model = fct_relevel(model, model_levels))

# Number of fixations should decline with ruggedness 
d_fixations_sum <- d_fixations %>%
  group_by(gen, model, r, isAdapted) %>%
  summarise(meanFixations = mean(nFixed),
            CIFixations = CI(nFixed))


ggplot(d_fixations_sum %>% mutate(gen = (gen - 50000) / 1000), 
       aes(x = interaction(gen, model), y = meanFixations, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r)~
                 "Did the population adapt?" + isAdapted) +
  geom_point() +
  scale_x_discrete(guide = "axis_nested") +
  geom_errorbar(aes(ymin = meanFixations - CIFixations,
                    ymax = meanFixations + CIFixations), position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           5, direction = 1),
                      labels = model_labels) +
  labs(x = TeX("Generations post-optimum shift ($x10^3$) / Model"), 
       y = "Mean number of cumulative fixations",
       colour = "Model") +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines")) -> plt_fixations 
plt_fixations
ggsave("plt_fixations.png", plt_fixations, device = png, 
       width = 12, height = 6)



d_fx <- data.table::fread(paste0(DATA_PATH, "calcMutationStats/d_fx.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    "factor", "character", "numeric"), 
                          col.names = c("gen", "seed", "modelindex", 
                                        "mutType", "mutID", "s"), 
                          fill = T)

# Join with phenotypic data
d_fx <- d_fx %>% distinct()
d_fx <- left_join(d_fx, d_qg, by = c("gen", "seed", "modelindex"))

# plot distribution of fitness effects

ggplot(d_fx %>%
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = s)) +
  facet_nested("Model" + model ~ "Mutation Type" + mutType) +
  geom_density() +
  labs(x = "Selection coefficient (s)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

driftBarrier <- 1 / (2 * 5000)

# More extreme ruggedness = more neutral mutations, more canalisation
# Also fewer deleterious mutations
# So look at proportion of deleterious mutations and those above barrier
# Run on HPC - memory usage
# d_fx_sum <- d_fx %>%
#   group_by(gen, seed, model, r, isAdapted) %>%
#   mutate(nMuts = n(),
#          propNeutral = sum((abs(s) < driftBarrier)) / nMuts,
#          propDel = sum(s < -driftBarrier) / nMuts,
#          propBen = 1 - (propNeutral + propDel)) %>%
#   ungroup() %>%
#   pivot_longer(cols = starts_with("prop"), 
#                names_to = "mutClass", values_to = "prop") %>%
#   group_by(gen, model, r, isAdapted, mutClass) %>%
#   summarise(meanProp = mean(prop),
#             CIProp = CI(prop))
# saveRDS(d_fx_sum, "d_fx_sum.RDS")
d_fx_sum <- readRDS(paste0(DATA_PATH, "calcMutationStats/d_fx_sum.RDS"))

# Create secondary axis
d_fx_sum <- d_fx_sum %>%
  mutate(gen_k = (gen - 50000) / 1000)  # Convert to "k" units

# Assign numeric x-axis positions
gen_levels <- sort(unique(d_fx_sum$gen_k))
gap <- 1  # space between model blocks

d_fx_sum <- d_fx_sum %>%
  group_by(model) %>%
  mutate(gen_index = as.integer(factor(gen_k, levels = gen_levels))) %>%
  ungroup() %>%
  mutate(model_index = as.integer(factor(model, levels = model_levels)),
         x_numeric = (model_index - 1) * (length(gen_levels) + gap) + gen_index)

gen_labels <- d_fx_sum %>%
  group_by(model, gen_k) %>%
  summarise(x = mean(x_numeric), .groups = "drop")

fx_model_labels <- d_fx_sum %>%
  group_by(model, r, isAdapted) %>%
  summarise(x = mean(x_numeric), .groups = "drop")

ggplot(d_fx_sum, 
       aes(x = x_numeric, y = meanProp, fill = mutClass)) +
  facet_nested("Recombination rate (log10)" + log10(r)~
                 "Did the population adapt?" + isAdapted) +
  geom_bar(stat = "identity", width = 1) +
  # geom_errorbar(aes(ymin = meanProp - CIProp,
  #                   ymax = meanProp + CIProp), 
  #               position = position_dodge(0.9)) +
  #scale_x_discrete(guide = "axis_nested") +
  #scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(ylim = c(0, 1), expand = F, clip = "off") +
  scale_fill_manual(values = paletteer_d("ggprism::viridis", 
                                           5, direction = 1)[c(5, 1, 3)],
                    labels = c("Beneficial", "Deleterious", "Neutral")) +
  labs(x = TeX("Generations post-optimum shift ($x10^3$) / Model"), 
       y = "Mean proportion of mutations",
       fill = "Mutation class") +
  scale_x_continuous(breaks = gen_labels$x, labels = gen_labels$gen_k) +
  theme_bw() +
  guides(colour = guide_legend(position = "bottom",
                               override.aes=list(linewidth = 5))) +
  geom_text(data = fx_model_labels %>% filter(r == 0.1), aes(x = x, y = -0.15, 
           label = model), inherit.aes = F) +
  geom_segment(data = fx_model_labels %>% filter(r == 0.1),
               aes(x = x - 5,
                   xend = x + 5,
                   y = -0.085, yend = -0.085),
               inherit.aes = FALSE,
               color = "black", size = 0.4) +
  theme(text = element_text(size = 12),
        panel.spacing = unit(0.75, "lines"),
        legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 28))) -> plt_prop
plt_prop

ggsave("plt_propFX.png", device = png, bg = "white",
       width = 16, height = 10)

# Plot the number of beneficial mutations per molecular component
d_fx_ben <- d_fx %>% filter(s > driftBarrier)
d_fx_del <- d_fx %>% filter(s < driftBarrier)

mutTypes_vec <- c(TeX("Additive", output = "character"),
                  TeX("$\\alpha_Z$", output = "character"),
                  TeX("$\\beta_Z$", output = "character"),
                  TeX("$K_Z$", output = "character"),
                  TeX("$K_{XZ}$", output = "character"))

d_fx_ben <- AddCombosToDF(d_fx_ben) 
d_fx_del <- AddCombosToDF(d_fx_del)

# nloci doesn't matter
d_fx_ben <- d_fx_ben %>%
  group_by(gen, seed, model, r) %>%
  mutate(countBen = n()) %>%
         #mutType = mutTypes_vec[as.numeric(mutType)]) %>%
  ungroup()



d_fx_ben_sum <- d_fx_ben %>%
  group_by(model, r) %>%
  summarise(meanCountBen = mean(countBen),
            CICountBen = CI(countBen),
            meanBen = mean(s),
            CIBen = CI(s))

ggplot(d_fx_ben %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = model, y = countBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_ben_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = model, y = meanCountBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  #scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Number of beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# TODO:
# Proportion of mutations that are beneficial in each model
d_fx_propBen <- d_fx %>%
  group_by(optPerc, seed, model, mutType, nloci, tau, r) %>%
  mutate(isBen = (s > 0)) %>%
  summarise(propBen = sum(isBen)/n()) %>%
  ungroup()

d_fx_propBen_sum <- d_fx_propBen %>%
  group_by(optPerc, model, mutType, r) %>%
  summarise(meanPropBen = mean(propBen),
            CIPropBen = CI(propBen))

ggplot(d_fx_propBen %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = propBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Progress to the optimum" + optPerc) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_propBen_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanPropBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Proportion of beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_propben_muts
plt_propben_muts

# Proportion of beneficial mutations doesn't appear to change over the walk
# because effect sizes are so small -> Geometric model, approaches 50%
# average over opt perc
d_fx_propBen_sum <- d_fx_propBen %>%
  group_by(model, mutType, r) %>%
  summarise(meanPropBen = mean(propBen),
            CIPropBen = CI(propBen))

ggplot(d_fx_propBen %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = propBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_propBen_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanPropBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Proportion of beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_propben_muts_wholewalk
plt_propben_muts_wholewalk
ggsave("plt_propben_muts_wholewalk.png", plt_propben_muts_wholewalk, 
       width = 10, height = 4, device = png)

# Average across molecular components for the average % beneficial 
# across all mutations
d_fx_propBen_sum <- d_fx_propBen %>%
  group_by(model, r) %>%
  summarise(meanPropBen = mean(propBen),
            CIPropBen = CI(propBen))

ggplot(d_fx_propBen %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = model, y = propBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_propBen_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanPropBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  labs(x = "Model", 
       y = "Proportion of beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_propben_muts
plt_propben_muts
ggsave("plt_propben.png", plt_propben_muts, 
       width = 10, height = 4, device = png)

# KZ drags the mean down, recombination rate doesn't seem to have much effect

# beta regression

# Is there a difference between models in % beneficial mutations?
# recombination has no effect, remove
# beta regression

# adjust for 0/1 values - need to inflate (Smithson + Verkuilen 2006)
# Set up combinations
d_fx_propBen <- transform(d_fx_propBen, modelMutType = factor(interaction(model, mutType)))
d_fx_propBen$propBen_adj <- (d_fx_propBen$propBen * (nrow(d_fx_propBen)-1) + 0.5)/nrow(d_fx_propBen)
d_fx_propBen$model <- as.factor(d_fx_propBen$model)

br.benmut <- betareg(propBen_adj ~ mutType, d_fx_propBen)

plot(density(d_fx_propBen$propBen))
summary(br.benmut)
plot(br.benmut)

saveRDS(br.benmut, "betareg_benmut.RDS")
br.benmut <- readRDS("betareg_benmut.RDS")

anova(br.benmut)
em.benmut <- emmeans(br.benmut, ~ modelMutType)


# In additive only one mutation type, so calc mean and CI
mean.benmut.add <- d_fx_propBen %>% filter(model == "Add") %>%
  summarise(meanPropBen = mean(propBen_adj),
            SEPropBen = se(propBen_adj),
            CIPropBen = CI(propBen_adj),
            upperCL = meanPropBen + CIPropBen,
            lowerCL = meanPropBen - CIPropBen)

br.benmut.km <- betareg(propBen_adj ~ mutType, d_fx_propBen %>% filter(model == "ODE"))
br.benmut.kp <- betareg(propBen_adj ~ mutType, d_fx_propBen %>% filter(model == "K"))

em.benmut.km <- emmeans(br.benmut.km, ~ mutType)
em.benmut.kp <- emmeans(br.benmut.kp, ~ mutType)

joint_tests(em.benmut.km)
joint_tests(em.benmut.kp)

xtable(em.benmut.km)
xtable(em.benmut.kp)

pairs(em.benmut, simple = "modelMutType")
plot(em.benmut, comparisons = T)
emmip(em.benmut,  ~ modelMutType)

# Differences between models (apart from KZ) range from <1% to ~5%
# KZ is about 50% change, almost never beneficial

# What about effect size: if K_Z has very large effects, generating that variation
# could be important
d_fx_ben <- d_fx_ben %>%
  group_by(gen, seed, model, r) %>%
  mutate(countBen = n()) %>%
  #mutType = mutTypes_vec[as.numeric(mutType)]) %>%
  ungroup()

d_fx_ben_sum <- d_fx_ben %>%
  mutate(model = factor(model, levels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"))) %>%
  group_by(model, r) %>%
  summarise(meanCountBen = mean(countBen),
            CICountBen = CI(countBen),
            meanBen = mean(s),
            CIBen = CI(s))


ggplot(d_fx_ben_sum %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = model, y = meanBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  #geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_ben_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanBen, group = model),
             size = 1, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanBen - CIBen, 
                    ymax = meanBen + CIBen), width = 0.1,
                position = position_dodge(0.9)) +
  #scale_x_discrete(labels = model_labels) +
  labs(x = "Model", 
       y = "Average fitness effect\nof beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_ben_muts_s
plt_ben_muts_s
ggsave("plt_ben_muts_s.png", plt_ben_muts_s, width = 6, height = 4, device = png)

# Deleterious mutations
d_fx_del <- d_fx_del %>%
  group_by(gen, seed, model, r) %>%
  mutate(countDel = n()) %>%
  #mutType = mutTypes_vec[as.numeric(mutType)]) %>%
  ungroup()

d_fx_del_sum <- d_fx_del %>%
  mutate(model = factor(model, levels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"))) %>%
  group_by(model, r) %>%
  summarise(meanCountDel = mean(countDel),
            CICountDel = CI(countDel),
            meanDel = mean(s),
            CIDel = CI(s))

ggplot(d_fx_del_sum %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = model, y = meanDel, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  #geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_del_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = model, y = meanDel, group = model),
             size = 1, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanDel - CIDel, 
                    ymax = meanDel + CIDel), width = 0.1,
                position = position_dodge(0.9)) +
  #scale_x_discrete(labels = model_labels) +
  labs(x = "Model", 
       y = "Average fitness effect\nof deleterious mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_del_muts_s
plt_del_muts_s
ggsave("plt_del_muts_s.png", plt_ben_muts_s, width = 6, height = 4, device = png)

# Both
d_fx_notNeutral <- d_fx %>% filter(s < driftBarrier | s > driftBarrier)
d_fx_notNeutral <- AddCombosToDF(d_fx_notNeutral) 
d_fx_notNeutral <- d_fx_notNeutral %>%
  group_by(gen, seed, model, r) %>%
  mutate(countMuts = n()) %>%
  #mutType = mutTypes_vec[as.numeric(mutType)]) %>%
  ungroup()
d_fx_notNeutral_sum <- d_fx_notNeutral %>%
  mutate(model = factor(model, levels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"))) %>%
  mutate(mutationClass = if_else(s > driftBarrier, "Beneficial", "Deleterious")) %>%
  group_by(model, r, mutationClass) %>%
  summarise(meanCount = mean(countMuts),
            CICount = CI(countMuts),
            meanS = mean(s),
            CIS = CI(s))

ggplot(d_fx_notNeutral_sum %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)"),
       aes(x = interaction(mutationClass, model), y = meanS, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  scale_x_discrete(guide = "axis_nested") + 
  
  geom_point(size = 1, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanS - CIS, 
                    ymax = meanS + CIS), width = 0.1,
                position = position_dodge(0.9)) +
  labs(x = "Model / mutation type", 
       y = "Average fitness effect",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5), 
                      guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 10),
        legend.position = "bottom") -> plt_muts_s
plt_muts_s
ggsave("plt_muts_s.png", plt_muts_s, width = 8, height = 6, device = png)



ggplot(d_fx_ben %>%
         ungroup() %>%
         filter(log10(r) != -5, 
                optPerc != levels(d_fx_ben$optPerc)[2],
                optPerc != levels(d_fx_ben$optPerc)[3]) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = s, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Progress to the optimum" + optPerc,
               labeller = labeller(optPerc = label_parsed)) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_ben_sum %>% ungroup() %>%
               filter(log10(r) != -5,
                      optPerc != levels(d_fx_ben$optPerc)[2],
                      optPerc != levels(d_fx_ben$optPerc)[3]) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Average fitness effect\nof beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_ben_muts_s_pres
plt_ben_muts_s_pres
ggsave("plt_ben_muts_s_pres.png", plt_ben_muts_s_pres, 
       width = 6, height = 6, device = png)


# GLS - fitness effect of beneficial mutations
# optPerc doesn't appear to matter

summary(gls.s.km <- gls(s ~ mutType, (d_fx_ben %>% filter(model == "ODE")), 
                       weights = varIdent(form = ~ 1 | mutType)))

# No difference between alpha and beta in effect size, so can just calc the mean
# like in the additive model across both alpha and beta
mean.s.km <- d_fx_ben %>% filter(model == "ODE") %>%
  summarise(meanS = mean(s),
            SES = se(s),
            CIS = CI(s),
            upperCL = meanS + CIS,
            lowerCL = meanS - CIS)

summary(gls.s.kp <- gls(s ~ mutType, (d_fx_ben %>% filter(model == "K")), 
                        weights = 
                          varIdent(form = ~ 1 | mutType)))

em.s.kp <- emmeans(gls.s.kp, ~ mutType)

anova(gls.s.kp)
anova(gls.s.km)

em.s.kp

mean.s.add <- d_fx_ben %>% filter(model == "Add") %>%
  summarise(meanS = mean(s),
            SES = se(s),
            CIS = CI(s),
            lowerCL = meanS - CIS,
            upperCL = meanS + CIS)

xtable(mean.s.add, digits = 6)
xtable(mean.s.km, digits = 6)
xtable(em.s.kp, digits = 9)


# Proportion of deleterious mutations
d_fx_propDel <- d_fx %>%
  group_by(optPerc, seed, model, mutType, nloci, tau, r) %>%
  mutate(isDel = (s < 0)) %>%
  summarise(propDel = sum(isDel)/n()) %>%
  ungroup()

d_fx_propDel_sum <- d_fx_propDel %>%
  group_by(optPerc, model, mutType, r) %>%
  summarise(meanPropDel = mean(propDel),
            CIPropDel = CI(propDel))

ggplot(d_fx_propDel %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = propDel, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Progress to the optimum" + optPerc) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_propDel_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanPropDel, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Proportion of deleterious mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_propdel_muts
plt_propdel_muts

# Proportion of beneficial mutations doesn't appear to change over the walk
# because effect sizes are so small -> Geometric model, approaches 50%
# average over opt perc
d_fx_propDel_sum <- d_fx_propDel %>%
  group_by(model, mutType, r) %>%
  summarise(meanPropDel = mean(propDel),
            CIPropDel = CI(propDel))

ggplot(d_fx_propDel %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = propDel, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_propDel_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanPropDel, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Proportion of deleterious mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_propdel_muts_wholewalk
plt_propdel_muts_wholewalk
ggsave("plt_propdel_muts_wholewalk.png", plt_propdel_muts_wholewalk, 
       width = 10, height = 4, device = png)



# SFS: maybe KZ should be rare?

# allele frequency spectrum
d_SFS <- data.table::fread(paste0(DATA_PATH, "calcMutationStats/d_SFS.csv"), header = F, 
                           sep = ",", colClasses = c("integer", "factor", 
                                                     "factor", "factor", 
                                                     rep("numeric", times = 3)), 
                           col.names = c("timePoint", "modelindex", "mutType", "freqBin",
                                         "countFreqBin", "meanValue", "sdValue"), 
                           fill = T)

d_SFS <- AddCombosToDF(d_SFS)

freqBins <- seq(from = 0.1, to = 1, by = 0.1)
d_SFS <- d_SFS %>%
  mutate(freqBin = freqBins[as.numeric(freqBin)]) %>%
  mutate(model = factor(model, levels = c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")))


# change in SFS over the walk
d_deltaSFS <- d_SFS %>%
  group_by(mutType, freqBin, model, r) %>%
  summarise(deltaCount = sum(diff(countFreqBin))) %>%
  ungroup()

d_deltaSFS_mean <- d_deltaSFS %>%
  group_by(mutType, freqBin, model, r) %>%
  summarise(meanDeltaCount = mean(deltaCount),
            CIDeltaCount = CI(deltaCount)) %>%
  ungroup()

ggplot(d_deltaSFS_mean %>% 
         mutate(freqBin = freqBin - 0.1),
       aes(x = freqBin, y = meanDeltaCount, fill = model)) +
  facet_nested(log10(r) ~ model + mutType) +
  geom_col(position = position_nudge(x = 0.05)) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, 
                                         direction = 1),
                    labels = model_labels) +
  labs(x = "Allele frequency", y = "Mean change in number\nof mutations over walk", 
       fill = "Model") +
  #coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))
# low recombination seems to have 

# Average 

d_SFS$mutType <- factor(d_SFS$mutType)

d_SFS <- d_SFS %>%
  mutate(mutType = fct_recode(mutType,
                              "alpha[Z]*'/Additive'" = "3",
                              "beta[Z]" = "4",
                              "K[Z]" = "5",
                              "K[XZ]" = "6"))

ggplot(d_SFS %>% 
         mutate(freqBin = freqBin - 0.1) %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(nloci == 1024, r %in% r_subsample) %>%
         uncount(countFreqBin),
       aes(x = freqBin, y = optPerc, fill = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ model + mutType, 
               labeller = labeller(mutType = label_parsed)) +
  stat_binline(bins = 10, binwidth = 0.1, position = position_nudge(x = 0.05), 
               scale = 0.95) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_y_discrete(labels = c("25%", "50%", "75%", "100%")) +
  labs(x = "Allele frequency", y = "Progress to the optimum", 
       fill = "Model") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12))
  
ggsave("plt_sfs_walk.png", device = png, width = 18, height = 6)

# Higher recombination rates look like they have more rare alleles and fewer
# fixations
# lots of segregating KZ mutations still, and plenty of fixations

# Fit model, hard to see in big grid
# Group bins into rare, intermediate, common/fixed
d_SFS <- d_SFS %>%
  mutate(freqBin_sml = cut(freqBin, breaks = c(0, 0.25, 0.75, 1),
                           labels = c("Rare", "Intermediate", "Common")))

d_SFS_freqClpse <- d_SFS %>%
  group_by(timePoint, mutType, model, r, freqBin_sml) %>%
  summarise(countFreqBin = sum(countFreqBin)) %>%
  ungroup()

ggplot(d_SFS_freqClpse, 
       aes(x = interaction(timePoint, as.factor(freqBin_sml)), y = countFreqBin, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ model + mutType, 
               labeller = labeller(mutType = label_parsed)) +
  scale_x_discrete(guide = "axis_nested") +
  geom_quasirandom() +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 5, 
                                         direction = -1)) + #,
                    #labels = c("Additive", "K+", "K-")) +
  labs(x = "Allele frequency", y = "Count of mutations", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12))
ggsave("plt_sfs_collapsed.png", device = png, width = 18, height = 6)


summary(glm.countSFS <- glm.nb(countFreqBin ~ freqBin_sml *
                                 model, 
                               data = d_SFS_freqClpse))

report(glm.countSFS)
plot(glm.countSFS)

# Look at residuals
library(statmod)
res.countSFS <- qresid(glm.countSFS)
qqnorm(res.countSFS)
qqline(res.countSFS)

performance::check_model(glm.countSFS)
