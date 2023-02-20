library(tidyverse)
library(mgcv)
library(ggeffects)
library(tidymv)
library(latex2exp)

# First need to combine additive/network data with a model type parameter
path_add <- "/mnt/d/SLiMTests/tests/newTestCross/additive/"
path_net <- "/mnt/d/SLiMTests/tests/newTestCross/getH2_newTestCross/data/"

combos_add <- read_delim("/mnt/c/GitHub/SLiMTests/tests/newTestCross/additive/R/combos.csv", delim = " ", col_names = F)
names(combos_add) <- c("nloci", "locisigma")
d_qg_add %>% mutate(nloci = combos_add$nloci[.$modelindex],
                sigma = combos_add$locisigma[.$modelindex]) -> d_qg_add


combos_net <- read_delim("/mnt/c/GitHub/SLiMTests/tests/newTestCross/R/combos.csv", delim = " ", col_names = F)
names(combos_net) <- c("model", "nloci", "locisigma")

d_muts_add <- read_csv(paste0(path_add, "d_combined_after.csv"))
d_muts_add$model <- "Additive"

d_muts_add %>% mutate(nloci = combos_add$nloci[.$modelindex],
                      sigma = combos_add$locisigma[.$modelindex]) -> d_muts_add


d_muts_net <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/d_combined_after.csv")
d_muts_net$model <- "NAR"


d_muts_net %>% 
  mutate(fixedEffect = combos_net$model[.$modelindex],
         nloci = combos_net$nloci[.$modelindex] + 2,
         sigma = combos_net$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_muts_net


d_muts <- bind_rows(d_muts_add, d_muts_net)
rm(d_muts_add, d_muts_net)

d_muts <- d_muts %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex),
         model = as_factor(model),
         nloci = as_factor(nloci),
         sigma = as_factor(sigma),
         fixedEffect = as_factor(fixedEffect))



# Using GAMs to fit out phenotype-time model
# u_z ~ gen + nloci + locisigma + c
d_qg_sbst <- d_muts %>% filter(fixedEffect == -1 | is.na(fixedEffect),
                             gen >= 49500,
                             phenomean < 10,
                             abs(estR) < 5,
                             aZ < 100 | is.na(aZ), bZ < 100 | is.na(bZ), 
                             KZ < 100 | is.na(KZ), KXZ < 100 | is.na(KXZ))  %>%
  distinct(gen, seed, modelindex, .keep_all = T)
  


mod_gam <- gam(phenomean ~ s(gen, k = 12) + model * nloci * sigma, 
               data = d_qg_sbst, family = scat, method = "REML")
summary(mod_gam)
plot(mod_gam)
gam.check(mod_gam)
concurvity(mod_gam)

# Test against a linear model
mod_lm <- gam(phenomean ~ gen + model * nloci * sigma, data = d_qg_sbst, method = "REML")
summary(mod_lm)


# GAM fits better
AIC(mod_gam); AIC(mod_lm)
summary(mod_gam)$sp.criterion; summary(mod_lm)$sp.criterion

plot(ggpredict(mod_gam), facets = T)

pred_gam <- predict_gam(mod_gam)
saveRDS(pred_gam, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/paper_figures/R/pred_gam.RDS")
saveRDS(mod_gam, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/paper_figures/R/mod_gam.RDS")


# Mutations
d_muts_add <- read_csv(paste0(path_add, "d_combined_after.csv"))
d_muts_add$model <- "Additive"

d_muts_add %>% mutate(nloci = combos_add$nloci[.$modelindex],
                    sigma = combos_add$locisigma[.$modelindex]) -> d_muts_add


d_muts_net <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/d_combined_after.csv")
d_muts_net$model <- "NAR"


d_muts_net %>% 
  mutate(fixedEffect = combos_net$model[.$modelindex],
                    nloci = combos_net$nloci[.$modelindex] + 2,
                    sigma = combos_net$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_muts_net


d_muts <- bind_rows(d_muts_add, d_muts_net)
rm(d_muts_add, d_muts_net)

d_muts <- d_muts %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex),
         model = as_factor(model),
         nloci = as_factor(nloci),
         sigma = as_factor(sigma),
         fixedEffect = as_factor(fixedEffect))

saveRDS(d_muts, "/mnt/d/SLiMTests/tests/newTestCross/addNetCombined/d_com_add+net_after.RDS")
d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/addNetCombined/d_com_add+net_after.RDS")
d_com %>% filter(phenomean < 10, abs(estR) < 5 | is.na(estR),
                 !(is.na(AIC) & gen >= 50000),
                 aZ < 100 | is.na(aZ), bZ < 100 | is.na(bZ), 
                 KZ < 100 | is.na(KZ), KXZ < 100 | is.na(KXZ)) -> d_com

# Filter out crazy variances
d_com %>% filter((VarA >= 0 & VarA < 10) | is.na(VarA),
                 (VarD >= 0 & VarD < 10) | is.na(VarD),
                 (VarAA >= 0 & VarAA < 10) | is.na(VarAA),
                 (VarR >= 0) | is.na(VarR)) -> d_com

# Using GAMs to fit out phenotype-time model
# u_z ~ gen + nloci + locisigma + c
d_qg_sbst <- d_com %>% filter(gen >= 49500)  %>%
  mutate(nloci = as.numeric(nloci),
         sigma = as.numeric(nloci)) %>%
  distinct(gen, seed, modelindex, .keep_all = T)



mod_gam <- gam(phenomean ~ s(gen, k = 20) + nloci * sigma * model, 
               data = d_qg_sbst, family = scat, method = "REML")
summary(mod_gam)
plot(mod_gam, all.terms = T, pages = 1, shade = T, seWithMean = T, shift = coef(mod_gam)[1])
gam.check(mod_gam)
concurvity(mod_gam)
concurvity(mod_gam, full = F)

# Test against a linear model
mod_lm <- gam(phenomean ~ gen + model * nloci * sigma, data = d_qg_sbst, method = "REML")
summary(mod_lm)


# GAM fits better
AIC(mod_gam); AIC(mod_lm)
summary(mod_gam)$sp.criterion; summary(mod_lm)$sp.criterion

plot(ggpredict(mod_gam), facets = T)

pred_gam <- predict_gam(mod_gam)
saveRDS(pred_gam, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/paper_figures/R/pred_gam.RDS")
saveRDS(mod_gam, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/paper_figures/R/mod_gam.RDS")

# GAM on allele frequency spectrum at gen 1950
# Load in data from paper_figures.R
d_freqs <- readRDS("d_freqs.RDS")

# normalised count ~ freq + nloci + locisigma + model + c
d_freqs_sbst <- d_freqs %>% filter(gen == 51950)

mod_gam <- gam(freqBinCount ~ s(Freq, k = 4) + model * nloci * sigma, 
               data = d_freqs_sbst, family = betar(), method = "REML")
summary(mod_gam)
plot(mod_gam)
gam.check(mod_gam)
concurvity(mod_gam)

# Test against a linear model
mod_lm <- gam(freqBinCount ~ Freq * model * nloci * sigma, data = d_freqs_sbst, method = "REML")
summary(mod_lm)

# GAM fits better
AIC(mod_gam); AIC(mod_lm)
summary(mod_gam)$sp.criterion; summary(mod_lm)$sp.criterion

plot(ggpredict(mod_gam), facets = T)
pred_gam <- predict_gam(mod_gam)
saveRDS(pred_gam, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/paper_figures/R/pred_gam_freq.RDS")
saveRDS(mod_gam, "/mnt/c/GitHub/SLiMTests/tests/newTestCross/paper_figures/R/mod_gam_freq.RDS")

# Various stats
## Comparing average additive variance
d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/addNetCombined/d_com_add+net_after.RDS")

d_com %>% filter(phenomean < 10, abs(estR) < 5 | is.na(estR),
                 !(is.na(AIC) & gen >= 50000),
                 aZ < 100 | is.na(aZ), bZ < 100 | is.na(bZ), 
                 KZ < 100 | is.na(KZ), KXZ < 100 | is.na(KXZ)) -> d_com
d_com %>% filter(phenomean < 10, abs(estR) < 5 | is.na(estR)) -> d_com

# Average h2 
d_com %>% filter(!is.na(AIC)) %>%
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  group_by(model) %>%
  summarise(meanH2A = mean(H2.A.Estimate),
            ci95H2A = qnorm(0.975) * se(H2.A.Estimate),
            meanH2D = mean(H2.D.Estimate),
            ci95H2D = qnorm(0.975) * se(H2.D.Estimate),
            meanH2AA = mean(H2.AA.Estimate),
            ci95H2AA = qnorm(0.975) * se(H2.AA.Estimate)) -> d_com_h2_table

stargazer(as.data.frame(d_com_h2_table), summary = F, rownames = F)


d_com %>% filter(!is.na(AIC)) %>%
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  group_by(model, nloci, sigma) %>%
  summarise(meanH2A = mean(H2.A.Estimate),
            ci95H2A = qnorm(0.975) * se(H2.A.Estimate),
            meanH2D = mean(H2.D.Estimate),
            ci95H2D = qnorm(0.975) * se(H2.D.Estimate),
            meanH2AA = mean(H2.AA.Estimate),
            ci95H2AA = qnorm(0.975) * se(H2.AA.Estimate)) -> d_com_h2_table

stargazer(as.data.frame(d_com_h2_table), summary = F, rownames = F)

# Mean deviation between observed and expected phenotypic response
d_com %>% filter(abs(estR) < 5) %>%
  filter(!(is.na(AIC) & gen >= 50000)) %>%
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  group_by(gen, model) %>%
  mutate(expPheno = lag(phenomean, default = 0) + lag(estR, default = 0)) %>%
  mutate(devEstR = sqrt((expPheno - phenomean)^2)) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(meanDev = mean(devEstR),
            ci95Dev = qnorm(0.975) * se(devEstR)) -> d_com_resp_dev_table

stargazer(as.data.frame(d_com_resp_dev_table), summary = F, rownames = F)


d_com %>% filter(abs(estR) < 5) %>%
  filter(!(is.na(AIC) & gen >= 50000)) %>%
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  group_by(gen, model, nloci, sigma) %>%
  mutate(expPheno = lag(phenomean, default = 0) + lag(estR, default = 0)) %>%
  mutate(devEstR = sqrt((expPheno - phenomean)^2)) %>% 
  ungroup() %>%
  group_by(model, nloci, sigma) %>%
  summarise(meanDev = mean(devEstR),
            ci95Dev = qnorm(0.975)*se(devEstR)) -> d_com_resp_dev_table

stargazer(as.data.frame(d_com_resp_dev_table), summary = F, rownames = F)
