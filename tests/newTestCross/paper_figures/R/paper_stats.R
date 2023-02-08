library(tidyverse)
library(mgcv)
library(ggeffects)
library(tidymv)
library(latex2exp)

# First need to combine additive/network data with a model type parameter
path_add <- "/mnt/d/SLiMTests/tests/newTestCross/additive/"
path_net <- "/mnt/d/SLiMTests/tests/newTestCross/getH2_newTestCross/data/"

d_qg_add <- read_csv(paste0(path_add, "slim_qg.csv"), col_names = F)
names(d_qg_add) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")
d_qg_add$model <- "Additive"

combos_add <- read_delim("/mnt/c/GitHub/SLiMTests/tests/newTestCross/additive/R/combos.csv", delim = " ", col_names = F)
names(combos_add) <- c("nloci", "locisigma")
d_qg_add %>% mutate(nloci = combos_add$nloci[.$modelindex],
                sigma = combos_add$locisigma[.$modelindex]) -> d_qg_add


d_qg_net <- read_csv(paste0(path_net, "slim_qg_total.csv"), col_names = F)
names(d_qg_net) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
d_qg_net$model <- "NAR"

combos_net <- read_delim("/mnt/c/GitHub/SLiMTests/tests/newTestCross/R/combos.csv", delim = " ", col_names = F)
names(combos_net) <- c("model", "nloci", "locisigma")
d_qg_net %>% mutate(fixedEffect = combos_net$model[.$modelindex],
                nloci = combos_net$nloci[.$modelindex],
                sigma = combos_net$locisigma[.$modelindex]) -> d_qg_net

# Recode fixedEffect to mol trait name
d_qg_net %>% mutate(fixedEffect = recode_factor(fixedEffect, `-1`="None", `0`="\u03B1", `1`="\u03B2", `2`="KZ", `3`="KXZ")) -> d_h2

d_qg_net$nloci <- d_qg_net$nloci + 2

d_qg <- bind_rows(d_qg_add, d_qg_net)
rm(d_qg_add, d_qg_net)

d_qg <- d_qg %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex),
         model = as_factor(model),
         nloci = as_factor(nloci),
         sigma = as_factor(sigma),
         fixedEffect = as_factor(fixedEffect))

# Using GAMs to fit out phenotype-time model
# u_z ~ gen + nloci + locisigma + c
d_qg_sbst <- d_qg %>% filter(fixedEffect == -1 | is.na(fixedEffect),
                             gen >= 50000,
                             phenomean < 10)

mod_gam <- gam(phenomean ~ s(gen, k = 12) + model * nloci * sigma, 
               data = d_qg_sbst, family = scat(link = "identity"), method = "REML")
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