# Note: data loaded using moreReps_figures.R
library(extRemes)
library(logspline)

setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("wrangle_data.R")

# Try to fit GEV
x <- revd(10000, loc = 0, scale = 1, shape = 0)
hist(x,prob=T,xlab="Random Variables from Gumbel (location = 0,scale = 1, shape =0)",
     main="Gumbel Distribution",ylab="f(x)",font=2,family="serif",font.lab=2,cex.lab=1.5)
plot(logspline(x),add=T)

# Set seed for reproducibility - 6673185 
seed <- sample(0:.Machine$integer.max, 1)
#set.seed(seed)
set.seed(6673185)

# Fitness effects
fit_add <- fevd(d_fix_add$avFit, method = "Lmoments")
plot(fit_add)
# bootstrap to find CIs
addci <- ci(fit_add, R = 1000, type = "parameter")

fit_aZ <- fevd(d_fix_aZ$avFit, method = "Lmoments")
plot(fit_aZ)
aZci <- ci(fit_aZ, R = 1000, type = "parameter")

fit_bZ <- fevd(d_fix_bZ$avFit, method = "Lmoments")
plot(fit_bZ)
bZci <- ci(fit_bZ, R = 1000, type = "parameter")

attr(addci, "class") <- NULL
attr(aZci, "class") <- NULL
attr(bZci, "class") <- NULL

shape_plt <- data.frame(effectType = c(TeX("$Additive$", output = "character"), 
                                       TeX("$\\alpha_Z$", output = "character"), 
                                       TeX("$\\beta_Z$", output = "character")),
                        estimate = c(addci[3,2], aZci[3,2], bZci[3,2]),
                        minCI = c(addci[3,1], aZci[3,1], bZci[3,1]),
                        maxCI = c(addci[3,3], aZci[3,3], bZci[3,3]))

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
ggsave("GEVshape_parameter_fitness.png", device = png)

# Phenotype effects
fit_add <- fevd(d_fix_add$value, method = "Lmoments")
plot(fit_add)
# bootstrap to find CIs
ci(fit_add, R = 1000, type = "parameter")

fit_aZ <- fevd(d_fix_aZ$avFX, method = "Lmoments")
plot(fit_aZ)
ci(fit_aZ, R = 1000, type = "parameter")

fit_bZ <- fevd(d_fix_bZ$avFX, method = "Lmoments")
plot(fit_bZ)
ci(fit_bZ, R = 1000, type = "parameter")


# might need to plot entire range of mutations to show that fixations are 
# the extreme values and the overall distribution is in the gumbel domain
d_muts_adapting <- d_muts %>% filter(gen >= 50000) %>% 
  distinct(seed, modelindex, mutID, .keep_all = T)

ggplot(d_muts_adapting %>% filter(modelindex == 1),
       aes(x = value)) +
  geom_density() +
  theme_bw() +
  theme(text = element_text(size = 16))

# This is the molecular component effect, will need to recalculate avFX
ggplot(d_muts_adapting %>% filter(modelindex == 2, mutType == 3),
       aes(x = value)) +
  geom_density() +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_muts_adapting %>% filter(modelindex == 2, mutType == 4),
       aes(x = value)) +
  geom_density() +
  theme_bw() +
  theme(text = element_text(size = 16))
