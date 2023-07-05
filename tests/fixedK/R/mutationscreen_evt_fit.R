library(tidyverse)
library(paletteer)
library(cowplot)
library(ggridges)
library(extRemes)
library(fitdistrplus)
source("mutationScreenExp.R")

# Set seed for reproducibility - 6673185 
seed <- sample(0:.Machine$integer.max, 1)
#set.seed(seed)
set.seed(6673185)

N_BOOT <- 1000
results <- data.frame(model = character(N_BOOT*2),
                      replicate = numeric(N_BOOT*2),
                      location = numeric(N_BOOT*2),
                      scale = numeric(N_BOOT*2),
                      shape = numeric(N_BOOT*2)) 

# bootstrap to fit EVD and find CIs
for (i in seq(from = 1, to = nrow(results), by = 2)) {
  mutExp_add_sub <- mutExp_add %>% ungroup() %>% sample_n(1000) %>% as.data.frame()
  mutExp_sub <- mutExp %>% ungroup() %>% sample_n(1000) %>% as.data.frame()
  fit_add <- fevd(s, mutExp_add_sub %>% dplyr::select(s), method = "Lmoments")
  fit_nar <- fevd(s, mutExp_sub %>% dplyr::select(s), method = "Lmoments")
  
  results[i,] <- c("Additive", i, fit_add$results)
  results[i+1,] <- c("NAR", i+1, fit_nar$results)
}

results %>% group_by(model) %>%
  mutate(location = as.numeric(location),
         scale = as.numeric(scale),
         shape = as.numeric(shape)) %>%
  summarise(meanLoc = mean(location),
            CILoc = CI(location),
            meanScale = mean(scale),
            CIScale = CI(scale),
            meanShape = mean(shape),
            CIShape = CI(shape)) -> evd_fit

plot(devd(seq(-1, 1, length.out = 1000), loc = evd_fit$meanLoc[1], 
          scale = evd_fit$meanScale[1], 
          shape = evd_fit$meanShape[1], threshold = 0,
          type = "GEV"), type = "l")

plot(devd(seq(-1, 1, length.out = 1000), loc = evd_fit$meanLoc[2], 
          scale = evd_fit$meanScale[2], 
          shape = evd_fit$meanShape[2], threshold = 0,
          type = "GEV"), type = "l")

ggplot(results,
       aes(x = as.numeric(scale), fill = model)) +
  geom_density(alpha = 0.4)

palette_main <- paletteer_d("ggsci::nrc_npg", 2)
palette_main
# plot mean EVDs
ggplot(mutExp_combined, aes(x = s, colour = as.factor(model))) +
  xlim(-1, 1) +
  geom_density(linetype = "dashed", alpha = 0.4, show.legend = F) +
  geom_function(fun = devd, linewidth = 1.2,
                colour = palette_main[1], n = 1000,
                args = list(loc = evd_fit$meanLoc[1], 
                            scale = evd_fit$meanScale[1], 
                            shape = evd_fit$meanShape[1])) +
  geom_function(fun = devd, linewidth = 1.2,
                colour = palette_main[2], n = 1000,
                args = list(loc = evd_fit$meanLoc[2], 
                            scale = evd_fit$meanScale[2], 
                            shape = evd_fit$meanShape[2])) +
  labs(x = "Fitness effect (s)", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_evd_fit

# Hack together a legend
plt_legend <- ggplot(data.frame(x = 1,2, 
                                y = 1,2,
                                model = c("Additive","NAR")), 
                     aes(x = x, y = y, colour = model)) +
  geom_line() +
  scale_colour_manual(values = palette_main) +
  labs(colour = "Model") + 
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
plt_legend

leg <- get_legend(plt_legend)

plot_grid(plt_evd_fit, leg, nrow = 2, rel_heights = c(1, 0.1))
ggsave("plt_evd_fit.png", device = png, bg = "white")

# Fit to beneficial mutations only

# Set seed for reproducibility - 6673185 
seed <- sample(0:.Machine$integer.max, 1)
#set.seed(seed)
set.seed(6673185)

N_BOOT <- 1000
results <- data.frame(model = character(N_BOOT*2),
                      replicate = numeric(N_BOOT*2),
                      location = numeric(N_BOOT*2),
                      scale = numeric(N_BOOT*2),
                      shape = numeric(N_BOOT*2)) 

mutExp_add_ben <- mutExp_add %>% filter(s > 0)
mutExp_ben <- mutExp %>% filter(s > 0)


# bootstrap to fit EVD and find CIs
for (i in seq(from = 1, to = nrow(results), by = 2)) {
  mutExp_add_sub <- mutExp_add_ben %>% ungroup() %>% sample_n(1000) %>% as.data.frame()
  mutExp_sub <- mutExp_ben %>% ungroup() %>% sample_n(1000) %>% as.data.frame()
  fit_add <- fevd(s, mutExp_add_sub %>% dplyr::select(s), method = "Lmoments")
  fit_nar <- fevd(s, mutExp_sub %>% dplyr::select(s), method = "Lmoments")
  
  results[i,] <- c("Additive", i, fit_add$results)
  results[i+1,] <- c("NAR", i+1, fit_nar$results)
}

results %>% group_by(model) %>%
  mutate(location = as.numeric(location),
         scale = as.numeric(scale),
         shape = as.numeric(shape)) %>%
  summarise(meanLoc = mean(location),
            CILoc = CI(location),
            meanScale = mean(scale),
            CIScale = CI(scale),
            meanShape = mean(shape),
            CIShape = CI(shape)) -> evd_fit

# plot mean EVDs
ggplot(mutExp_combined %>% filter(s > 0), aes(x = s, colour = as.factor(model))) +
  xlim(0, 0.25) +
  geom_density(linetype = "dashed", alpha = 0.4, show.legend = F) +
  geom_function(fun = devd, linewidth = 1.2,
                colour = palette_main[1], n = 1000,
                args = list(loc = evd_fit$meanLoc[1], 
                            scale = evd_fit$meanScale[1], 
                            shape = evd_fit$meanShape[1])) +
  geom_function(fun = devd, linewidth = 1.2,
                colour = palette_main[2], n = 1000,
                args = list(loc = evd_fit$meanLoc[2], 
                            scale = evd_fit$meanScale[2], 
                            shape = evd_fit$meanShape[2])) +
  labs(x = "Fitness effect (s)", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_evd_fit_ben

# Hack together a legend
plt_legend <- ggplot(data.frame(x = 1,2, 
                                y = 1,2,
                                model = c("Additive","NAR")), 
                     aes(x = x, y = y, colour = model)) +
  geom_line() +
  scale_colour_manual(values = palette_main) +
  labs(colour = "Model") + 
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
plt_legend

leg <- get_legend(plt_legend)

plot_grid(plt_evd_fit_ben, leg, nrow = 2, rel_heights = c(1, 0.1))
ggsave("plt_evd_ben_fit.png", device = png, bg = "white")


# Fit exponential to beneficial mutations
fit_nar <- fitdist((mutExp %>% filter(s > 0))$s, 
                   distr = "exp", method = "mle")
fit_add <- fitdist((mutExp_add %>% filter(s > 0))$s, 
                   distr = "exp", method = "mle")
summary(fit_nar)
summary(fit_add)
plot(fit_nar)
plot(fit_add)

# This code from Lebeuf-Taylor et al. 2019
## Testing the shape of the distribution of beneficial mutation ####
## Use both the synonymous and non-synonymous together because their beneficial distributions are not significantly different.

# sample fitness effects and return a sorted vector
sampleFitnessEffects <- function(s, n) {
  # Measure fitness relative to the smallest beneficial selection coefficient
  # as recommended by Beisel et al. 2007
  X = s-min(s)
  X = X[X!=0]
  X = sample(X, n)
  X = sort(X)
  return(X)
}

X = sampleFitnessEffects(mutExp_ben$s, 1000)

# additive
X = sampleFitnessEffects(mutExp_add_ben$s, 1000)


LL.GPD <- function(par, d){ # equations from Beisel et al 2007
  # X is the adjusted selection coefficients, i.e. the data
  X = d
  # tau is a parameter of the GPD
  tau = par[1]
  # kappa is a parameter of the GDP
  kappa = par[2]
  # n_1 is the number of observations
  n_1 = length(X)
  LL.a = -(n_1)*log(tau)
  LL.b = NULL
  # for kappa>0
  if(kappa>0) {
    LL.b = -(kappa+1)/kappa * sum (log(1+(kappa*X)/tau))
  }
  # for kappa<0
  if(kappa<0) {
    if(sum(X>-tau/kappa)){LL.b = -1000000}
    else{LL.b = -(kappa+1)/kappa * sum (log(1+(kappa*X)/tau))}
  }
  # for kappa=0
  if(kappa==0) {
    LL.b = -(1/tau)*sum(X)
  }
  LL = LL.a + LL.b
  return(-LL)
}

## Special case kappa=0 (exponential distribution)
LL.Exp <- function(par, d){ # equations from Beisel et al 2007
  # X is the adjusted selection coefficients, i.e. the data
  X = d
  # tau is a parameter of the GPD
  tau = par
  # kappa is a parameter of the GDP
  kappa = 0
  
  # n_1 is the number of observations
  n_1 = length(X)
  LL.a = -(n_1)*log(tau)
  LL.b = -(1/tau)*sum(X)
  LL = LL.a + LL.b
  return(-LL)
}

#1# Optimization ####
# install.packages("GenSA")
require(GenSA)
require(boot)

bootBeisel <- function(d, nullLR) {
  start.par = c(.1, 0) # starting parameter values for the optimization
  ## 1) Find tau, When kappa = 0
  Exp.opt = GenSA(par = start.par[1], fn = LL.Exp, lower = 0.000001, upper = 100, d = d)
  Exp.opt.tau = Exp.opt$par
  LL.Exp.opt = -Exp.opt$value
  
  ## 2) Find best tau and kappa
  GPD.opt = GenSA(par = start.par, fn = LL.GPD, lower = c(0.000001, -100), upper = c(100, 100), d = d)
  GPD.opt.tau = GPD.opt$par[1]
  GPD.opt.kappa = GPD.opt$par[2]
  LL.GPD.opt = -GPD.opt$value
  
  ## Likelihood ratio ####
  LRT = LL.Exp.opt/LL.GPD.opt
  
  P.value = sum(LRT>nullLR)/length(d)
  return(c(GPD.opt.tau, GPD.opt.kappa, LRT, P.value))
  
}

## Find a null distribution of likelihood ratio ####
B = 10000
LR.null.dist = vector(length = B)
for(b in 1:B){
  X = rexp(n = length(X), rate = 1/Exp.opt.tau)
  Exp.null.opt = optim(par = start.par[1], fn = LL.Exp, method = "Brent", lower = 0.000001, upper = 100, d = X)
  LL.Exp.null.opt = -Exp.null.opt$value
  ## 2) Find best tau and kappa
  GPD.null.opt = optim(par = c(Exp.opt.tau, 0), fn = LL.GPD, method = "L-BFGS-B", lower = c(0.000001, -100), upper = c(100, 100), d = X)
  # print(b)
  # GPD.null.opt = GenSA(par = start.par, fn = LL.GPD, lower = c(0.000001, -100), upper = c(100, 100))
  LL.GPD.null.opt = -GPD.null.opt$value
  ## Likelihood ratio
  LR.null.dist[b] = LL.Exp.null.opt/LL.GPD.null.opt
}

bootBeisel_nar <- data.frame(tau = numeric(B),
                             kappa = numeric(B),
                             LRT = numeric(B),
                             p.value = numeric(B))

for (b in 1:B) {
  d <- sampleFitnessEffects(mutExp_ben$s, 1000)
  bootBeisel_nar[b, ] <- bootBeisel(d, LR.null.dist)
}
write.csv(bootBeisel_nar, "bootBeisel_nar.csv", row.names = F)

bootBeisel_add <- data.frame(tau = numeric(B),
                             kappa = numeric(B),
                             LRT = numeric(B),
                             p.value = numeric(B))

library(future)
library(doParallel)
library(foreach)

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run in parallel
bootBeisel_add <- foreach(b=1:B, .combine = "rbind") %dopar% {
  require(GenSA)
  d <- sampleFitnessEffects(mutExp_add_ben$s, 1000)
  bootBeisel_add[b, ] <- bootBeisel(d, LR.null.dist)
}
stopCluster(cl)

bootBeisel_add <- as.data.frame(bootBeisel_add)
colnames(bootBeisel_add) <- c("tau", "kappa", "LRT", "p.value")

write.csv(bootBeisel_add, "bootBeisel_add.csv", row.names = F)



# calculate mean and CI of bootstrap
mean(bootBeisel_nar$kappa)
CI(bootBeisel_nar$kappa)
mean(bootBeisel_nar$p.value)
CI(bootBeisel_nar$p.value)
mean(bootBeisel_nar$LRT)
CI(bootBeisel_nar$LRT)

mean(bootBeisel_add$kappa)
CI(bootBeisel_add$kappa)
mean(bootBeisel_add$p.value)
CI(bootBeisel_add$p.value)
mean(bootBeisel_add$LRT)
CI(bootBeisel_add$LRT)
## Plot data and best fit distribution ####
install.packages("fExtremes")
require(fExtremes)

## Results plot ####
hist(X, freq = F, main = 'All beneficial mutations', breaks = 15, ylim = c(0, 35), 
     xlab = "s", col = c("mediumpurple1"))
y = seq(from=0.001, to=0.15, by=0.001)
lines(y, dgpd(x = y, xi = GPD.opt.kappa, mu = 0, beta = GPD.opt.tau), col='black', lwd=2)
#lines(y, dexp(x = y, rate = 1/Exp.opt.tau), col='black', lwd=2)
legend(0.04, 25, legend = c( 
  expression(paste("Generalized Pareto (", kappa, " is fit)"))), col = c("black"), lwd = 2, bty = "n")
text(0.08, 15, labels = expression(paste("Weibull (", kappa, " = -0.370, ", tau, " = 0.0330)")))



# likelihood ratio mutExp: is the fit different if we try to fit a Gumbel?
fit_add_gumbel <- fevd(s, mutExp_add_ben %>% select(s), type = "Gumbel")
fit_gumbel <- fevd(s, mutExp_ben %>% select(s), type = "Gumbel")

lr.mutExp(fit_add_gumbel, fit_add)
lr.mutExp(fit_gumbel, fit_nar)


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
ggsave("fig_s1_evdfit.png")
