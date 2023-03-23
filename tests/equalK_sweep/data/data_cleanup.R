library(tidyverse)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}
# Running on HPC for RAM reasons
path <- "/g/data/ht96/nb9894/equalK_sweep/"
setwd(path)

# Read in data

# combos
d_combos <- read_delim("~/tests/equalK_sweep/R/combos.csv", delim = " ", col_names = F)
names(d_combos) <- c("nloci", "locisigma")


d_qg <- read_csv(paste0(path, "slim_qg.csv"), col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
# filter off any of the fixed molecular trait ones, since the newest dataset lacks those
d_qg %>% 
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_qg

d_muts <- read_csv(paste0(path, "slim_muts.csv"), col_names = F)
names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")    
d_muts %>% 
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_muts


head(d_muts)
head(d_qg)

# Factorise seed and modelindex
d_muts <- d_muts %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))
d_qg <- d_qg %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))

d_phenos <- read_csv(paste0(path, "slim_indPheno.csv"), col_names = F)
colnames(d_phenos) <- c("gen", "seed", "modelindex", "ind", "phenotype", "aZ", "bZ", "KZ", "KXZ")

d_phenos %>% 
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_phenos

d_phenos <- d_phenos %>%
  mutate(modelindex = as_factor(modelindex),
         seed = as_factor(seed))

d_phenos %>%
  group_by(gen, seed, modelindex) %>%
  mutate(ind = 1:n()) %>%
  ungroup() %>%
  mutate(ind = as_factor(ind)) -> d_indPheno



if (!dir.exists("./checkpoint"))
{
  dir.create("checkpoint")
}

saveRDS(d_muts, "./checkpoint/d_muts.RDS")
saveRDS(d_qg, "./checkpoint/d_qg.RDS")
saveRDS(d_indPheno, "./checkpoint/d_indPheno.RDS")

d_muts <- readRDS("./checkpoint/d_muts.RDS")
d_qg <- readRDS("./checkpoint/d_qg.RDS")


d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex", "nloci", "sigma"))
rm(d_muts, d_qg)

saveRDS(d_combined, paste0(path, "checkpoint/d_combined.RDS"))

d_combined <- readRDS("./checkpoint/d_combined.RDS")



# Write only the data after burn-in
d_combined_after <- d_combined %>% filter(gen >= 49500)
nrow(d_combined_after)
nrow(d_combined)

rm(d_combined)

# d_combined_after <- read_csv(paste0(path, "moreReps2/getH2_equalK_sweep/data/d_combined_after.csv"))

saveRDS(d_combined_after, paste0(path, "d_combined_after.RDS"))
d_combined_after <- readRDS(paste0(path, "d_combined_after.RDS"))


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

d_com_filtered <- d_combined_after %>%
  drop_na(meanH, value, Freq, phenomean, w)

d_com_filtered %>%
  group_by(seed, modelindex) %>%
  mutate(IQRpheno = IQR(phenomean),
         outlier_pheno_upper = quantile(phenomean, probs=c( .75), na.rm = FALSE)+1.5*IQRpheno,
         outlier_pheno_lower = quantile(phenomean, probs=c( .25), na.rm = FALSE)-1.5*IQRpheno
        ) %>%
  filter(any(outlier_pheno_lower <= phenomean & outlier_pheno_upper > phenomean)) -> d_com_filtered

ggplot(d_com_filtered %>% mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean, colour = sigma)) +
  geom_boxplot() -> plt_dist

ggsave("IQR_pheno_dist.png", plt_dist, width = 10, height = 10)


saveRDS(d_com_filtered, "d_com_after_prefiltered.RDS")
d_com_filtered <- readRDS("d_com_after_prefiltered.RDS")

# Split into adapted/maladapted subgroups

d_com_adapted <- d_com_filtered %>%
  group_by(seed, modelindex) %>%
  mutate(adapted = between(phenomean, 1.9, 2.1)) %>%
  filter(!(gen == 51950 & adapted == F)) %>%
  filter(any(gen == 51950)) %>%
  ungroup()

# Check they are adapted
ggplot(d_com_adapted %>% 
  distinct(gen, seed, modelindex, nloci, sigma, .keep_all = T) %>%
  filter(gen == 51950) %>%
  mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("adapted_pheno_dist.png", plt_dist, width = 10, height = 10)


# Filter data into group of those who didn't make it: the groups where 
# any of the phenotypes means aren't within (1.9, 2.1) at any time
d_com_maladapted <- d_com_filtered %>% 
  group_by(seed, modelindex) %>%
  mutate(adapted = between(phenomean, 1.9, 2.1)) %>%
  filter(any(gen > 51800)) %>%
  filter(all(adapted == F)) %>%
  ungroup()

# Check they are all maladapted
ggplot(d_com_maladapted %>% 
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("maladapted_pheno_dist.png", plt_dist, width = 10, height = 10)

# Filter data into group of those who made it and then kept moving: the groups where 
# any of the phenotypes means are within (1.9, 2.1) sometime before generation 51800
# also check that the last value isn't in (1.9, 2.1)
d_com_wasadapted <- d_com_filtered %>%
  group_by(seed, modelindex) %>%
  mutate(adapted = between(phenomean, 1.9, 2.1)) %>%
  filter(any(gen > 51800)) %>%
  filter(all( (gen < 51800 & adapted | gen < 51800 & !adapted | 
    (gen >= 51800 & !adapted)))) %>%
  ungroup()



saveRDS(d_com_adapted, "d_comh2_prefiltered_adapted.RDS")
saveRDS(d_com_maladapted, "d_comh2_prefiltered_maladapted.RDS")
saveRDS(d_com_wasadapted, "d_comh2_prefiltered_wasadapted.RDS")

# Randomly sample three populations from each group so we can track them
seed <- sample(1:.Machine$integer.max, 1)
# Sampled seed: 
set.seed(seed)
set.seed()

d_com_adapted_eg <- d_com_adapted %>% 
  group_by(modelindex) %>%
  filter(seed %in% sample(unique(seed), min(length(unique(seed)), 1)))

d_com_maladapted_eg <- d_com_maladapted %>% 
  group_by(modelindex) %>%
  filter(seed %in% sample(unique(seed), min(length(unique(seed)), 1)))

d_com_wasadapted_eg <- d_com_wasadapted %>% 
  group_by(modelindex) %>%
  filter(seed %in% sample(unique(seed), min(length(unique(seed)), 1)))


saveRDS(d_com_adapted_eg, "d_com_prefiltered_adapted_eg.RDS")
saveRDS(d_com_maladapted_eg, "d_com_prefiltered_maladapted_eg.RDS")
saveRDS(d_com_wasadapted_eg, "d_com_prefiltered_wasadapted_eg.RDS")

d_com_eg <- read_table("~/tests/indTrack/R/combos.csv", col_names = F)
names(d_com_eg) <- c("model", "nloci", "sigma", "seed")
