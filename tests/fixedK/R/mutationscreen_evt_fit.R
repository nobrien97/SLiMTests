library(tidyverse)
library(paletteer)
library(cowplot)
library(ggridges)
source("mutationScreenExp.R")

# Set seed for reproducibility - 6673185 
seed <- sample(0:.Machine$integer.max, 1)
#set.seed(seed)
set.seed(6673185)

test_add_ben <- test_add %>% filter(s > 0)
test_ben <- test %>% filter(s > 0)

fit_add <- fevd(test_add_ben$s, method = "Lmoments")
plot(fit_add)
# bootstrap to find CIs
addci <- ci(fit_add, R = 1000, type = "parameter")

summary(fit_add)

fit_nar <- fevd(test_ben$s, method = "Lmoments")
plot(fit_nar)
NARci <- ci(fit_nar, R = 1000, type = "parameter")

summary(fit_nar)

attr(addci, "class") <- NULL
attr(NARci, "class") <- NULL

shape_plt <- data.frame(effectType = c(TeX("$Additive$", output = "character"), 
                                       TeX("NAR", output = "character")),
                        estimate = c(addci[3,2], NARci[3,2]),
                        minCI = c(addci[3,1], NARci[3,1]),
                        maxCI = c(addci[3,3], NARci[3,3]))

# Figure S1: bar plot of shape parameters
ggplot(shape_plt, aes(x = effectType, y = estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = minCI, ymax = maxCI), width = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  scale_x_discrete(labels = parse(text = shape_plt$effectType)) +
  labs(x = "Allelic effect type", y = TeX("Estimated shape parameter ($\\xi$)")) +
  theme_bw() +
  theme(text = element_text(size = 16))
