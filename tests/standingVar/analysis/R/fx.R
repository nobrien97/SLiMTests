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

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/analysis/R")
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
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

# Use the right mutate/summarise functions
mutate <- dplyr::mutate
summarise <- dplyr::summarise
select <- dplyr::select

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

# Add predictors
d_qg <- AddCombosToDF(d_qg) 

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg

d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, optPerc, isAdapted) %>% filter(gen >= 49500)


d_fx <- data.table::fread(paste0(DATA_PATH, "mutationStats/d_fx.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor",
                                                    "factor", "character", "numeric"), 
                          col.names = c("gen", "seed", "modelindex", 
                                        "mutType", "mutID", "s"), 
                          fill = T)

# Join with phenotypic data
d_fx <- d_fx %>% distinct()
d_fx <- left_join(d_fx, d_qg_optPerc, by = c("gen", "seed", "modelindex"))
d_fx <- AddCombosToDF(d_fx)

r_subsample <- c(1e-10, 1e-5, 1e-1)

# Filter to groups we're looking at
d_fx <- d_fx %>%
  filter(isAdapted, tau == 0.0125, r %in% r_subsample)

# plot distribution of fitness effects

ggplot(d_fx %>% filter(mutType != "5") %>%
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = s)) +
  facet_nested("Model" + model ~ "Mutation Type" + mutType) +
  geom_density() +
  labs(x = "Selection coefficient (s)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")


# Plot the number of beneficial mutations per molecular component
d_fx_ben <- d_fx %>% filter(s > 0)

mutTypes_vec <- c(TeX("$\\alpha_Z$/Additive", output = "character"),
                  TeX("$\\beta_Z$", output = "character"),
                  TeX("$K_Z$", output = "character"),
                  TeX("$K_{XZ}$", output = "character"))

# nloci doesn't matter
d_fx_ben <- d_fx_ben %>%
  group_by(optPerc, seed, model, mutType, nloci, tau, r) %>%
  mutate(countBen = n()) %>%
         #mutType = mutTypes_vec[as.numeric(mutType)]) %>%
  ungroup()

d_fx_ben_sum <- d_fx_ben %>%
  group_by(optPerc, model, mutType, r) %>%
  summarise(meanCountBen = mean(countBen),
            CICountBen = CI(countBen),
            meanBen = mean(s),
            CIBen = CI(s))

ggplot(d_fx_ben %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = countBen, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Progress to the optimum" + optPerc) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_ben_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = mutType, y = meanCountBen, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_x_discrete(labels = parse(text = mutTypes_vec)) +
  labs(x = "Molecular component", 
       y = "Number of beneficial mutations",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_ben_muts
plt_ben_muts
ggsave("plt_ben_muts.png", plt_ben_muts, width = 10, height = 4, device = png)

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
em.benmut <- emmeans(br.benmut, ~ mutType)

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

xtable(em.benmut.km)
xtable(em.benmut.kp)

pairs(em.benmut, simple = "modelMutType")
plot(em.benmut, comparisons = T)
emmip(em.benmut,  ~ modelMutType)

# Differences between models (apart from KZ) range from <1% to ~5%
# KZ is about 50% change, almost never beneficial

# What about effect size: if K_Z has very large effects, generating that variation
# could be important
d_fx_ben_sum <- d_fx_ben %>%
  group_by(optPerc, model, mutType, r) %>%
  summarise(meanCountBen = mean(countBen),
            CICountBen = CI(countBen),
            meanBen = mean(s),
            CIBen = CI(s))

levels(d_fx_ben$optPerc) <- c(
  TeX("$\\leq 25\\%$"),
  TeX("$\\leq 50\\%$"),
  TeX("$\\leq 75\\%$"),
  TeX("$\\leq 100\\%$")
)

levels(d_fx_ben_sum$optPerc) <- c(
  TeX("$\\leq 25\\%$"),
  TeX("$\\leq 50\\%$"),
  TeX("$\\leq 75\\%$"),
  TeX("$\\leq 100\\%$")
)


ggplot(d_fx_ben %>% 
         ungroup() %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = mutType, y = s, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Progress to the optimum" + optPerc,
               labeller = labeller(optPerc = label_parsed)) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, varwidth = T, na.rm = F) +
  geom_point(data = d_fx_ben_sum %>% ungroup() %>%
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
        legend.position = "bottom") -> plt_ben_muts_s
plt_ben_muts_s
ggsave("plt_ben_muts_s.png", plt_ben_muts_s, width = 10, height = 4, device = png)


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


# SFS: maybe KZ should be rare?

# allele frequency spectrum
d_SFS <- data.table::fread(paste0(DATA_PATH, "mutationStats/d_SFS.csv"), header = F, 
                           sep = ",", colClasses = c("factor", "factor", "numeric", 
                                                     "factor",
                                                     rep("numeric", times = 3)), 
                           col.names = c("optPerc", "modelindex", "mutType", "freqBin",
                                         "countFreqBin", "meanValue", "sdValue"), 
                           fill = T)

d_SFS <- AddCombosToDF(d_SFS)

freqBins <- seq(from = 0.1, to = 1, by = 0.1)
d_SFS <- d_SFS %>%
  mutate(freqBin = freqBins[as.numeric(freqBin)])

r_subsample <- c(1e-10, 1e-5, 1e-1)

# Filter to groups we're looking at
d_SFS <- d_SFS %>%
  filter(tau == 0.0125, r %in% r_subsample)

# change in SFS over the walk
d_deltaSFS <- d_SFS %>%
  group_by(mutType, freqBin, model, nloci, r) %>%
  summarise(deltaCount = sum(diff(countFreqBin))) %>%
  ungroup()

d_deltaSFS_mean <- d_deltaSFS %>%
  group_by(mutType, freqBin, model, r) %>%
  summarise(meanDeltaCount = mean(deltaCount),
            CIDeltaCount = CI(deltaCount)) %>%
  ungroup()

ggplot(d_deltaSFS_mean %>% 
         mutate(freqBin = freqBin - 0.1) %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE"),
                mutType = as.factor(mutType),
                mutType = fct_recode(mutType, "$\\alpha_Z$" = "3",
                                     "$\\beta_Z$" = "4",
                                     "K_Z" = "5",
                                     "K_{XZ}" = "6")),
       aes(x = freqBin, y = meanDeltaCount, fill = model)) +
  facet_nested(log10(r) ~ model + mutType) +
  geom_col(position = position_nudge(x = 0.05)) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1),
                    labels = c("Additive", "K+", "K-")) +
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
  group_by(optPerc, mutType, model, nloci, tau, r, freqBin_sml) %>%
  summarise(countFreqBin = sum(countFreqBin)) %>%
  ungroup()

ggplot(d_SFS_freqClpse %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")), 
       aes(x = freqBin_sml, y = countFreqBin, colour = model,
           shape = optPerc)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ model + mutType, 
               labeller = labeller(mutType = label_parsed)) +
  geom_quasirandom() +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                         direction = -1),
                    labels = c("Additive", "K+", "K-")) +
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
