library(tidyverse)
library(gghighlight)
library(paletteer)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/sanityChecks/R")

source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

d_combos <- read.table("combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "width", "model"))

AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           nloci = d_combos$nloci[as.numeric(levels(modelindex))[modelindex]],
           tau = d_combos$tau[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]],
           width = d_combos$width[as.numeric(levels(modelindex))[modelindex]])
}


DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/sanityChecks/"

d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  filter(gen >= 49500) %>% distinct() %>%
 # mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg


d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg <- AddCombosToDF(d_qg)

d_qg_sum <- d_qg %>% 
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r, width) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

# phenotype: looks reasonable - 
# really very little difference between 0.5 and 0.1 recombination, same with 
# 1e-10 and 0
# neutral model behaves as expected for 1024 loci, little change over time
# stronger selection makes adaptation faster
ggplot(d_qg %>% filter(tau == 0.0125),
       aes(x = gen, y = phenomean, colour = model, 
           group = interaction(modelindex, as.factor(seed)))) +
  facet_grid(r~width) +
  geom_line() +
  gghighlight(seed == sample(unique(d_qg$seed), 1), calculate_per_facet = T,
              unhighlighted_params = list(colour = NULL, alpha = 0.1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 0.0125") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Selection strength", 
                                         breaks = NULL, labels = NULL)) +
  #coord_cartesian(ylim = c(0, 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Generations post-optimum shift", y = "Population mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))

