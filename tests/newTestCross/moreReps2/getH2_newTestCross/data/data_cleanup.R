library(tidyverse)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}
# Running on HPC for RAM reasons
path <- "/g/data/ht96/nb9894/newTestCross/"
setwd(path)

# Read in data

# combos
d_combos_net_fixed <- read_delim("~/tests/newTestCross/R/combos.csv", delim = " ", col_names = F)
names(d_combos_net_fixed) <- c("model", "nloci", "locisigma")

d_combos_net <- read_delim("~/tests/newTestCross/moreReps2/R/combos.csv", delim = " ", col_names = F)
names(d_combos_net) <- c("nloci", "locisigma")


d_qg_old <- read_csv(paste0(path, "slim_qg.csv"), col_names = F)
d_qg_new <- read_csv(paste0(path, "moreReps/slim_qg.csv"), col_names = F)
names(d_qg_old) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
names(d_qg_new) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
# filter off any of the fixed molecular trait ones, since the newest dataset lacks those
d_qg_old %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_qg_old


d_qg_new %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_qg_new


d_qg_newest <- read_csv(paste0(path, "moreReps2/slim_qg.csv"), col_names = F)
names(d_qg_newest) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
                  
d_qg_newest %>% 
  mutate(fixedEffect = -1,
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) -> d_qg_newest


d_h2_old <- read_csv(paste0(path, "getH2_newTestCross/out_h2.csv"), col_names = F)
d_h2_new <- read_csv(paste0(path, "moreReps/getH2_newTestCross/out_h2_new.csv"), col_names = F)
names(d_h2_old) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR",
                  "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                  "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                  "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")
names(d_h2_new) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR",
                  "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                  "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                  "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")

d_h2_old %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_h2_old


d_h2_new %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_h2_new


d_h2_newest <- read_csv(paste0(path, "moreReps2/getH2_newTestCross/out_h2.csv"), col_names = F)
names(d_h2_newest) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR",
                  "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                  "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                  "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")
d_h2_newest %>% 
  mutate(fixedEffect = -1,
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) -> d_h2_newest


d_muts_old <- read_csv(paste0(path, "slim_muts.csv"), col_names = F)
names(d_muts_old) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")    
d_muts_old %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_muts_old

d_muts_new <- read_csv(paste0(path, "moreReps/slim_muts.csv"), col_names = F)  
names(d_muts_new) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")                                                                                                                   

d_muts_new %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_muts_new


# Newest doesn't fix any molecular traits and has a different combos.csv
d_muts_newest <- read_csv(paste0(path, "moreReps2/slim_muts.csv"), col_names = F)  
names(d_muts_newest) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")                                                                                                                   


d_muts_newest %>% 
  mutate(fixedEffect = -1,
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) -> d_muts_newest



d_qg <- rbind(d_qg_old, d_qg_new, d_qg_newest)
d_h2 <- rbind(d_h2_old, d_h2_new, d_h2_newest)
d_muts <- rbind(d_muts_old, d_muts_new, d_muts_newest)

rm(d_qg_old, d_qg_new, d_qg_newest, d_h2_old, d_h2_new, d_h2_newest, d_muts_new, d_muts_old, d_muts_newest)

head(d_muts)
head(d_qg)
head(d_h2)

# Factorise seed and modelindex
d_muts <- d_muts %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))
d_qg <- d_qg %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))
d_h2 <- d_h2 %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))



d_phenos_old <- read_csv(paste0(path, "slim_sampled_pheno.csv"), col_names = F)
d_phenos_old %>% pivot_longer(cols = 4:ncol(d_phenos_old), names_to = NULL, values_to = "phenotype") -> d_phenos_old
colnames(d_phenos_old) <- c("gen", "seed", "modelindex", "phenotype")

d_phenos_old %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_phenos_old


d_phenos_new <- read_csv(paste0(path, "moreReps/slim_sampled_pheno.csv"), col_names = F)
d_phenos_new %>% pivot_longer(cols = 4:ncol(d_phenos_new), names_to = NULL, values_to = "phenotype") -> d_phenos_new
colnames(d_phenos_new) <- c("gen", "seed", "modelindex", "phenotype")

d_phenos_new %>% 
  mutate(fixedEffect = d_combos_net_fixed$model[.$modelindex],
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) %>%
  filter(fixedEffect == -1) -> d_phenos_new


d_phenos_newest <- read_csv(paste0(path, "moreReps2/slim_sampled_pheno.csv"), col_names = F)
d_phenos_newest %>% pivot_longer(cols = 4:ncol(d_phenos_newest), names_to = NULL, values_to = "phenotype") -> d_phenos_newest
colnames(d_phenos_newest) <- c("gen", "seed", "modelindex", "phenotype")
d_phenos_newest %>% 
  mutate(fixedEffect = -1,
         nloci = d_combos_net_fixed$nloci[.$modelindex] + 2,
         sigma = d_combos_net_fixed$locisigma[.$modelindex]) -> d_phenos_newest


d_phenos <- rbind(d_phenos_old, d_phenos_new, d_phenos_newest)
rm(d_phenos_new, d_phenos_old, d_phenos_newest)

d_phenos <- d_phenos %>%
  mutate(modelindex = as_factor(modelindex),
         seed = as_factor(seed))

# Calculate selection gradients for R = B * VA
# width here is w^2, so need to convert our 1/2w^2 to this format
width_condensed <- 0.05
width_sqrd <- 1/(width_condensed * 2)
opt <- 2
TIME_BETWEEN_SAMPLES <- 50 # 50 generations between each sample

calcBeta <- function(S, mu, theta) {
  return(-S * (mu - theta))
}

# Use Morrissey and Goudie 2022 to predict beta

d_phenos <- d_phenos %>%
  group_by(gen, seed, nloci, sigma) %>%
  summarise(S = 1/(width_sqrd + var(phenotype)),
         beta = calcBeta(S, mean(phenotype), opt))

d_h2 <- full_join(d_h2, d_phenos, by = c("gen", "seed", "nloci", "sigma"))
d_h2$estR <- (d_h2$VarA * d_h2$beta) * TIME_BETWEEN_SAMPLES

d_h2 <- d_h2 %>%
  filter(!is.na(AIC)) %>%
  distinct(.keep_all = T)

if (!dir.exists("./checkpoint"))
{
  dir.create("checkpoint")
}

saveRDS(d_h2, "./checkpoint/d_h2.RDS")
saveRDS(d_muts, "./checkpoint/d_muts.RDS")
saveRDS(d_qg, "./checkpoint/d_qg.RDS")

d_h2 <- readRDS("./checkpoint/d_h2.RDS")
d_muts <- readRDS("./checkpoint/d_muts.RDS")
d_qg <- readRDS("./checkpoint/d_qg.RDS")

# Combine the data frames, remove the unneeded variables
d_h2 <- d_h2 %>% select(!c(fixedEffect, modelindex))
d_muts <- d_muts %>% select(!c(fixedEffect, modelindex))
d_qg <- d_qg %>% select(!c(fixedEffect, modelindex))

d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "nloci", "sigma"))
rm(d_muts, d_qg)
d_combined <- full_join(d_combined, d_h2, by = c("gen", "seed", "nloci", "sigma"))
rm(d_h2)
# Deviation of observed R to estimated R
d_combined$devEstR <- d_combined$deltaPheno - d_combined$estR

saveRDS(d_combined, paste0(path, "moreReps2/d_combined.RDS"))

d_combined <- readRDS("./checkpoint/d_combined.RDS")



# Write only the data after burn-in
d_combined_after <- d_combined %>% filter(gen >= 49500)
nrow(d_combined_after)
nrow(d_combined)

rm(d_combined)

# d_combined_after <- read_csv(paste0(path, "moreReps2/getH2_newTestCross/data/d_combined_after.csv"))

saveRDS(d_combined_after, paste0(path, "moreReps2/d_combined_after.RDS"))
d_combined_after <- readRDS(paste0(path, "moreReps2/d_combined_after.RDS"))

# Combine with additive
path_add <- paste0(path, "additive/")

combos_add <- read_delim(paste0(path_add, "combos.csv"), delim = " ", col_names = F)
names(combos_add) <- c("nloci", "locisigma")


d_com_add <- read_csv(paste0(path_add, "d_combined_after.csv"))

# read in more reps data
d_com_add2 <- readRDS(paste0(path_add, "moreReps/d_combined_after.RDS"))
d_com_add <- rbind(d_com_add, d_com_add2)
rm(d_com_add2)

d_com_add$model <- "Additive"
d_combined_after$model <- "NAR"

d_com_add$modelindex <- as.numeric(d_com_add$modelindex)

d_com_add %>% mutate(fixedEffect = -1,
                    nloci = combos_add$nloci[modelindex],
                    sigma = combos_add$locisigma[modelindex]) -> d_com_add

d_com_add %>% select(-c(fixedEffect, modelindex)) -> d_com_add

d_com_add <- d_com_add %>%
  mutate(seed = as_factor(seed),
         model = as_factor(model),
         nloci = as_factor(nloci),
         sigma = as_factor(sigma))

d_combined_after <- d_combined_after %>%
  mutate(seed = as_factor(seed),
         model = as_factor(model),
         nloci = as_factor(nloci),
         sigma = as_factor(sigma))


d_com <- bind_rows(d_com_add, d_combined_after)
rm(d_combined_after, d_com_add)

saveRDS(d_com, paste0(path, "d_com_add+net_after.RDS"))

# Have a look at the data
d_com <- readRDS("d_com_add+net_after.RDS")

d_com %>%
  group_by(model, nloci, sigma) %>%
  drop_na(estR) %>%
  summarise(meanR = mean(estR),
          seR = se(estR),
          medianR = median(estR),
          IQRR = IQR(estR),
          minR = min(estR),
          maxR = max(estR)) -> d_range_estR

# Drop the outliers using Mahalanobis distance - 
# https://www.r-bloggers.com/2017/12/combined-outlier-detection-with-dplyr-and-ruler/
# library(robustbase)

# alpha <- .001
# cutoff <- (qchisq(p = 1 - alpha, df = 9))

# d_com_filtered <- d_com %>%
#   drop_na(meanH, value, Freq, phenomean, w, VarA, VarD, VarAA, estR) %>%
#   filter(H2.A.Estimate <= 1, H2.D.Estimate <= 1, H2.AA.Estimate <= 1)

# dist_dat <- d_com_filtered %>% 
#               select(model, nloci, sigma, phenomean, estR)


# mcddat <- by(dist_dat, list(dist_dat$model, dist_dat$nloci, dist_dat$sigma), 
#           function(x)
#           {
#             y <- x %>% select(!c(model, nloci, sigma))
#             centerCov <- covMcd(y)
#             y_len <- nrow(y)
#             c(
#               model = rep(x$model, y_len),
#               nloci = rep(x$nloci, y_len),
#               sigma = rep(x$sigma, y_len),
#               md = mahalanobis(y, centerCov$center, centerCov$cov)
#             )
#           })

# df_md <- do.call(rbind, mcddat)

# d_com_filtered_md <- inner_join(d_com_filtered, df_md, by = c("model", "nloci", "sigma"))

# d_com_filtered_md %>% filter(md < cutoff) -> d_com_filtered

# # glue back on the gen == 49500 entries
# d_com %>% mutate(md = 0) -> d_com
# d_com_filtered <- rbind(d_com %>% filter(gen == 49500), d_com_filtered)

# saveRDS(d_com_filtered, "d_com_add+net_after_prefiltered.RDS")

# ggplot(d_com_filtered, aes(x = nloci, y = estR, colour = model)) +
#   facet_grid(.~sigma) +
#   geom_boxplot() -> plt_dist

# ggsave("estR_MCD_dist.png", plt_dist, width = 10, height = 10)


# IQR 1.5

d_com_filtered <- d_com %>%
  drop_na(meanH, value, Freq, phenomean, w, VarA, VarD, VarAA, estR) %>%
  filter(H2.A.Estimate <= 1, H2.D.Estimate <= 1, H2.AA.Estimate <= 1)

d_com_filtered %>%
  group_by(seed, model, sigma, nloci) %>%
  mutate(IQRestR = IQR(estR),
         IQRpheno = IQR(phenomean),
         outlier_estR_upper = quantile(estR, probs=c( .75), na.rm = FALSE)+1.5*IQRestR,
         outlier_estR_lower = quantile(estR, probs=c( .25), na.rm = FALSE)-1.5*IQRestR,
         outlier_pheno_upper = quantile(phenomean, probs=c( .75), na.rm = FALSE)+1.5*IQRpheno,
         outlier_pheno_lower = quantile(phenomean, probs=c( .25), na.rm = FALSE)-1.5*IQRpheno
        ) %>%
  filter(any(outlier_estR_lower <= estR & outlier_estR_upper > estR),
         any(outlier_pheno_lower <= phenomean & outlier_pheno_upper > phenomean)) -> d_com_filtered

ggplot(d_com_filtered, aes(x = nloci, y = estR, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("IQR_estR_dist.png", plt_dist, width = 10, height = 10)

ggplot(d_com_filtered %>% mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("IQR_pheno_dist.png", plt_dist, width = 10, height = 10)

# Combine with end of burn-in and remove IQR values
d_com_filtered <- rbind(d_com %>% filter(gen == 49500), d_com_filtered %>% select(!c(49:54)))

# filter again to get rid of wacky phenomeans in gen == 49500
d_com_filtered %>%
  filter(H2.A.Estimate <= 1, H2.D.Estimate <= 1, H2.AA.Estimate <= 1,
         phenomean < 5) -> d_com_filtered

saveRDS(d_com_filtered, "d_com_add+net_after_prefiltered.RDS")



d_com %>%
  drop_na(estR) %>%
  filter(between(estR,
          quantile(estR, 0.025), quantile(estR, 0.975))) -> d_com_top95

d_com_top95 %>%
  group_by(model, nloci, sigma) %>%
  drop_na(estR) %>%
  summarise(meanR = mean(estR),
          seR = se(estR),
          medianR = median(estR),
          IQRR = IQR(estR),
          minR = min(estR),
          maxR = max(estR)) -> d_range_estR

saveRDS(d_com_top95, "d_com_add+net_after_prefiltered.RDS")


ggplot(d_com_top95, aes(x = nloci, y = estR, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("estR_dist.png", plt_dist, width = 10, height = 10)


# Too much removed
# d_com %>% filter(phenomean < 10, abs(estR) < 5 | is.na(estR),
#                  !(is.na(AIC) & gen >= 50000),
#                  aZ < 100 | is.na(aZ), bZ < 100 | is.na(bZ), 
#                  KZ < 100 | is.na(KZ), KXZ < 100 | is.na(KXZ)) -> d_com

# # Filter out crazy variances
# d_com %>% filter((VarA >= 0 & VarA < 10) | is.na(VarA),
#                  (VarD >= 0 & VarD < 10) | is.na(VarD),
#                  (VarAA >= 0 & VarAA < 10) | is.na(VarAA),
#                  (VarR >= 0) | is.na(VarR)) -> d_com

saveRDS(d_com, paste0(path, "d_com_add+net_after_prefiltered.RDS"))




# Sample some populations and track them, don't care about h2, need all the data. 
# So we have to reprocess...

# NAR
d_muts <- readRDS("./checkpoint/d_muts.RDS")
d_qg <- readRDS("./checkpoint/d_qg.RDS")

# Combine the data frames
d_com <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex", "fixedEffect", "nloci", "sigma"))
rm(d_muts, d_qg)

#d_com <- d_combined %>% filter(gen >= 49500)
#rm(d_combined)

# Remove unneeded variables
d_com <- d_com %>% select(!c(fixedEffect, modelindex))
d_com$model <- "NAR"

# Additive
path_add <- paste0(path, "additive/")
# Read in data
d_qg <- read_csv(paste0(path_add, "slim_qg.csv"), col_names = F)
d_muts <- read_csv(paste0(path_add, "slim_muts.csv"), col_names = F)                                        
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")
names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")
                                                                                                                   
# Factorise seed and modelindex
d_muts <- d_muts %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))
d_qg <- d_qg %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))

d_combined_add <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex"))
rm(d_muts, d_qg)

combos <- read_delim(paste0(path_add, "combos.csv"), delim = " ", col_names = F)
names(combos) <- c("nloci", "locisigma")

d_combined_add %>% mutate(
                            nloci = combos$nloci[.$modelindex],
                            sigma = combos$locisigma[.$modelindex]) -> d_combined_add
d_combined_add <- d_combined_add %>% select(!modelindex) #%>% filter(gen >= 49500)

d_combined_add <- d_combined_add %>% mutate(aZ = NA, bZ = NA, KZ = NA, KXZ = NA, model = "Additive")

d_com_both <- rbind(d_combined_add, d_com)

# Remove outlier replicates - if they have an outlier at any point we need to remove the whole replicate 
d_com_both %>%
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  group_by(seed, model, nloci, sigma) %>%
  mutate(IQRpheno = IQR(phenomean),
         outlier_pheno_upper = quantile(phenomean, probs=c( .75), na.rm = FALSE)+1.5*IQRpheno,
         outlier_pheno_lower = quantile(phenomean, probs=c( .25), na.rm = FALSE)-1.5*IQRpheno,
         outlier = phenomean < outlier_pheno_lower | phenomean > outlier_pheno_upper
        ) %>%
  filter(!outlier) -> d_com_filtered


# d_com_filtered %>%
#  mutate(nloci = as_factor(nloci),
        # sigma = as_factor(sigma)) -> d_com_filtered

d_com_filtered$nloci <- factor(d_com_filtered$nloci, levels = c(10, 100, 1000))
d_com_filtered$sigma <- factor(d_com_filtered$sigma, levels = c(0.1, 1))


ggplot(d_com_filtered %>% mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("IQR_pheno_dist.png", plt_dist, width = 10, height = 10)

saveRDS(d_com_filtered, "d_com_prefiltered.RDS")


# Filter data into group of those who made it to the optimum: the groups where 
# any of the phenotypes means are within (1.9, 2.1) after generation 51800

d_com_adapted <- d_com_filtered %>%
  group_by(seed, model, nloci, sigma) %>%
  filter(any(gen > 51800)) %>%
  filter(any(gen >= 51800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup()

# Filter data into group of those who didn't make it: the groups where 
# any of the phenotypes means aren't within (1.9, 2.1) at any time
d_com_maladapted <- d_com_filtered %>% 
  group_by(seed, model, nloci, sigma) %>%
  filter(any(gen > 51800)) %>%
  filter(all(phenomean < 1.9 | phenomean > 2.1)) %>%
  ungroup()

# Filter data into group of those who made it and then kept moving: the groups where 
# any of the phenotypes means are within (1.9, 2.1) sometime before generation 51800
# also check that the last value isn't in (1.9, 2.1)
d_com_wasadapted <- anti_join(
  d_com_filtered %>% group_by(seed, model, nloci, sigma) %>%
  filter(any(gen > 51800)), 
  rbind(d_com_adapted, d_com_maladapted), by = c("seed", "model", "nloci", "sigma"))

d_com_wasadapted <- d_com_filtered %>%
  group_by(seed, model, nloci, sigma) %>%
  filter(any(gen > 51800)) %>%
  filter(all( (gen < 51800 & (phenomean >= 1.9 & phenomean <= 2.1)) | 
    (gen >= 51800 & (phenomean < 1.9 | phenomean > 2.1)))) %>%
  ungroup()



saveRDS(d_com_adapted, "d_com_prefiltered_adapted.RDS")
saveRDS(d_com_maladapted, "d_com_prefiltered_maladapted.RDS")
saveRDS(d_com_wasadapted, "d_com_prefiltered_wasadapted.RDS")

# Randomly sample three populations from each group so we can track them
seed <- sample(1:.Machine$integer.max, 1)
# Sampled seed: 1849510312
set.seed(seed)
set.seed(1849510312)
d_com_adapted_eg <- d_com_adapted %>% 
  group_by(model, nloci, sigma) %>%
  filter(seed %in% sample(unique(seed), min(length(unique(seed)), 3)))

d_com_maladapted_eg <- d_com_maladapted %>% 
  group_by(model, nloci, sigma) %>%
  filter(seed %in% sample(unique(seed), min(length(unique(seed)), 3)))

d_com_wasadapted_eg <- d_com_wasadapted %>% 
  group_by(model, nloci, sigma) %>%
  filter(seed %in% sample(unique(seed), min(length(unique(seed)), 3)))


saveRDS(d_com_adapted_eg, "d_com_prefiltered_adapted_eg.RDS")
saveRDS(d_com_maladapted_eg, "d_com_prefiltered_maladapted_eg.RDS")
saveRDS(d_com_wasadapted_eg, "d_com_prefiltered_wasadapted_eg.RDS")
