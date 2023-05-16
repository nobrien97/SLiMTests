# Note: data loaded using moreReps_figures.R
library(extRemes)
library(logspline)

# Try to fit GEV
x <- revd(10000, loc = 0, scale = 1, shape = 0)
hist(x,prob=T,xlab="Random Variables from Gumbel (location = 0,scale = 1, shape =0)",
     main="Gumbel Distribution",ylab="f(x)",font=2,family="serif",font.lab=2,cex.lab=1.5)
plot(logspline(x),add=T)

# Fitness effects
fit_add <- fevd(d_fix_add$avFit, method = "Lmoments")
plot(fit_add)
# bootstrap to find CIs
ci(fit_add, R = 1000, type = "parameter")

fit_aZ <- fevd(d_fix_aZ$avFit, method = "Lmoments")
plot(fit_aZ)
ci(fit_aZ, R = 1000, type = "parameter")

fit_bZ <- fevd(d_fix_bZ$avFit, method = "Lmoments")
plot(fit_bZ)
ci(fit_bZ, R = 1000, type = "parameter")

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
