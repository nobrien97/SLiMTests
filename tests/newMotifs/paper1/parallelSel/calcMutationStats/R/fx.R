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
library(DMwR2)


DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/randomisedStarts/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/"
COMBO_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/"

DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/randomisedStarts/"
R_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/randomisedStarts/calcMutationStats/R/"
COMBO_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/R/"


DATA_PATH <- "/g/data/ht96/nb9894/newMotifs/randomisedStarts/"
 R_PATH <- "~/tests/newMotifs/randomisedStarts/calcMutationStats/R/"
 COMBO_PATH <- "~/tests/newMotifs/R/"
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

d_combos <- read.table(paste0(COMBO_PATH, "combos.csv"), header = F,
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

# Average at end of sim
d_fixations_end_sum <- d_fixations %>% filter(gen == 59000) %>%
  group_by(model, r, isAdapted) %>%
  summarise(meanFixations = mean(nFixed),
            CIFixations = CI(nFixed))

print(xtable(d_fixations_end_sum %>% 
               select(model, isAdapted, r, meanFixations, CIFixations) %>%
               mutate(r = as.integer(log10(r)))),
      include.rownames = F)

# Average across all recombination rates
d_fixations_end_sum <- d_fixations %>% filter(gen == 59000) %>%
  group_by(model, isAdapted) %>%
  summarise(meanFixations = mean(nFixed),
            CIFixations = CI(nFixed))

print(xtable(d_fixations_end_sum %>% 
               select(model, isAdapted, meanFixations, CIFixations)),
      include.rownames = F)




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

# Fixations at end of sim
ggplot(d_fixations_sum %>% filter(gen == 59000), 
       aes(x = model, y = meanFixations, colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~
                 "Adaptation outcome" + isAdapted,
               labeller = labeller(isAdapted = c("FALSE" = "Maladapted", "TRUE" = "Adapted"))) +
  geom_point() +
  geom_errorbar(aes(ymin = meanFixations - CIFixations,
                    ymax = meanFixations + CIFixations), position = position_dodge(0.9),
                width = 0.5) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades",
                                           5, direction = 1),
                      labels = model_labels, guide = "none") +
  scale_x_discrete(labels = model_labels) +
  labs(x = "Model", 
       y = "Mean number of cumulative\nfixations at end of simulation") +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.spacing = unit(0.75, "lines")) -> plt_fixations_end
plt_fixations_end
ggsave("plt_fixations_end.png", plt_fixations_end, device = png, 
       width = 9, height = 4.5)

# Is there any difference between the models in the number of fixations?
d_fixations_sbst <- d_fixations %>% filter(gen == 59000) %>%
  mutate(r = as.factor(r))

d_fixations_sbst %>%
  group_by(model, r, isAdapted) %>%
  summarise(n = n()) %>%
  spread(model, n, fill = 0)

lm_numFix <- glm(nFixed ~ model * r * isAdapted, family = "poisson",
                 data = d_fixations_sbst)

# Check dispersion (we dispersed)
AER::dispersiontest(lm_numFix)

# Do glm instead, remove outliers for model fitting
# Detect outliers: Hampel filter
lofscores <- lofactor(scale(d_fixations_sbst$nFixed), 80)
threshold <- 1.01
outliers <- lofscores > threshold
plot(density(lofscores[lofscores < 1.01]))

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score")
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)
boxplot(d_fixations_sbst[!outliers,]$nFixed)

d_fixations_glm <- d_fixations_sbst[!outliers,]

# Model not significant, removed 
glm_numFix <- glm.nb(nFixed ~ r * isAdapted, 
                 data = d_fixations_glm)

summary(glm_numFix)
gratia::appraise(glm_numFix)
plot(glm_numFix)

report::report(glm_numFix)

anova(glm_numFix)

# Marginal means
em_numFix <- emmeans(glm_numFix, ~ r * isAdapted)
logemm.src <- regrid(emmeans(glm_numFix, ~ r * isAdapted), transform = "log")
confint(logemm.src, type = "response")
pairs(logemm.src, simple = "r", type = "response")
plot(em_numFix, comparisons = T)

xtable(em_numFix)


########
# effect sizes
########

d_fx <- data.table::fread(paste0(DATA_PATH, "calcMutationStats/d_fx.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    "factor", "character", "numeric"), 
                          col.names = c("gen", "seed", "modelindex", 
                                        "mutType", "mutID", "s"), 
                          fill = T)

# Join with phenotypic data
d_fx <- d_fx %>% distinct()
d_fx <- left_join(d_fx, d_qg, by = c("gen", "seed", "modelindex"))

# Drop any NAs and Infs
d_fx <- d_fx %>% drop_na(s)

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
#   drop_na(s) %>%
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
# 
# saveRDS(d_fx_sum, paste0(DATA_PATH, "calcMutationStats/d_fx_sum.RDS"))
d_fx_sum <- readRDS(paste0(DATA_PATH, "calcMutationStats/d_fx_sum.RDS"))

# Table of propNeutral, del, ben across all time points
d_fx_props <- d_fx %>%
  drop_na(s) %>%
  group_by(seed, model, r, isAdapted) %>%
  summarise(nMuts = n(),
         propNeutral = sum((abs(s) < driftBarrier)) / nMuts,
         propDel = sum(s < -driftBarrier) / nMuts,
         propBen = 1 - (propNeutral + propDel)) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("prop"),
               names_to = "mutClass", values_to = "prop")


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

gen_labels <- d_fx_sum %>% filter(gen_k != 0) %>%
  group_by(model, gen_k) %>%
  summarise(x = mean(x_numeric), .groups = "drop")

fx_model_labels <- d_fx_sum %>% filter(gen_k != 0) %>%
  group_by(model, r, isAdapted) %>%
  summarise(x = mean(x_numeric), .groups = "drop")

ggplot(d_fx_sum %>% filter(gen_k > 0), 
       aes(x = x_numeric, y = meanProp, fill = mutClass)) +
  facet_nested("Recombination rate (log10)" + log10(r)~
                 "Adaptation outcome" + isAdapted,
               labeller = labeller(isAdapted = c("FALSE" = "Maladapted", "TRUE" = "Adapted"))) +
  geom_bar(stat = "identity", width = 1) +
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
  geom_text(data = fx_model_labels %>% filter(r == 0.1), aes(x = x, y = -0.15 * 2.5, 
           label = model), inherit.aes = F) +
  geom_segment(data = fx_model_labels %>% filter(r == 0.1),
               aes(x = x - 5,
                   xend = x + 5,
                   y = -0.085 * 2.5, yend = -0.085 * 2.5),
               inherit.aes = FALSE,
               color = "black", size = 0.4) +
  theme(text = element_text(size = 14),
        panel.spacing = unit(0.75, "lines"),
        legend.position = "bottom",
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = 28))) -> plt_prop
plt_prop

ggsave("plt_propFX.png", device = png, bg = "white",
       width = 10, height = 5)

# Table of proportions
print(xtable(d_fx_sum %>% filter(gen_index == 10) %>%
               mutate(r = as.integer(log10(r)),
                      isAdapted = if_else(isAdapted, "Adapted", "Maladapted"),
                      mutClass = case_match(
                        mutClass,
                        "propBen" ~ "Beneficial",
                        "propDel" ~ "Deleterious",
                        "propNeutral" ~ "Neutral"
                      )) %>%
               select(model, isAdapted, r, mutClass, meanProp, CIProp),
             digits = 3),
      include.rownames = F)
