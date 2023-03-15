library(tidyverse)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}
# Running on HPC for RAM reasons
path <- "/g/data/ht96/nb9894/newTestCross/equalK/"
setwd(path)

# Read in data

# combos
d_combos <- read_delim("~/tests/newTestCross/equalK/R/combos.csv", delim = " ", col_names = F)
names(d_combos) <- c("nloci", "locisigma")


d_qg <- read_csv(paste0(path, "slim_qg.csv"), col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
# filter off any of the fixed molecular trait ones, since the newest dataset lacks those
d_qg %>% 
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_qg

d_h2 <- read_csv(paste0(path, "getH2/out_h2.csv"), col_names = F)
names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR",
                  "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                  "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                  "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")

d_h2 %>% 
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_h2

d_muts <- read_csv(paste0(path, "slim_muts.csv"), col_names = F)
names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")    
d_muts %>% 
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_muts

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


# Phenotype information
d_phenos <- read_csv(paste0(path, "slim_sampled_pheno.csv"), col_names = F)
d_phenos %>% pivot_longer(cols = 4:ncol(d_phenos), names_to = NULL, values_to = "phenotype") -> d_phenos
colnames(d_phenos) <- c("gen", "seed", "modelindex", "phenotype")

d_phenos %>%
  mutate(nloci = d_combos$nloci[.$modelindex],
         sigma = d_combos$locisigma[.$modelindex]) -> d_phenos

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

saveRDS(d_phenos, "./checkpoint/d_phenos.RDS")
saveRDS(d_h2, "./checkpoint/d_h2.RDS")
saveRDS(d_muts, "./checkpoint/d_muts.RDS")
saveRDS(d_qg, "./checkpoint/d_qg.RDS")

d_h2 <- readRDS("./checkpoint/d_h2.RDS")
d_muts <- readRDS("./checkpoint/d_muts.RDS")
d_qg <- readRDS("./checkpoint/d_qg.RDS")

# Combine the data frames
d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex", "nloci", "sigma"))
rm(d_muts, d_qg)
d_combined <- full_join(d_combined, d_h2, by = c("gen", "seed", "modelindex", "nloci", "sigma"))
rm(d_h2)
# Deviation of observed R to estimated R
d_combined$devEstR <- d_combined$deltaPheno - d_combined$estR

saveRDS(d_combined, "./checkpoint/d_combined.RDS")

d_combined <- readRDS("./checkpoint/d_combined.RDS")



# Write only the data after burn-in
d_combined_after <- d_combined %>% filter(gen >= 49500)
nrow(d_combined_after)
nrow(d_combined)

rm(d_combined)

# d_combined_after <- read_csv(paste0(path, "moreReps2/getH2_newTestCross/data/d_combined_after.csv"))

saveRDS(d_combined_after, paste0(path, "d_combined_after.RDS"))
d_combined_after <- readRDS(paste0(path, "d_combined_after.RDS"))


# Have a look at the data

d_combined_after %>%
  group_by(nloci, sigma) %>%
  drop_na(estR) %>%
  summarise(meanR = mean(estR),
          seR = se(estR),
          medianR = median(estR),
          IQRR = IQR(estR),
          minR = min(estR),
          maxR = max(estR)) -> d_range_estR

# IQR 1.5 - drop outliers

d_com_filtered <- d_combined %>%
  drop_na(meanH, value, Freq, phenomean, w, VarA, VarD, VarAA, estR) %>%
  filter(H2.A.Estimate <= 1, H2.D.Estimate <= 1, H2.AA.Estimate <= 1)

# Remove outlier replicates - if they have an outlier at any point we need to remove the whole replicate 
d_com_filtered %>%
  group_by(seed, nloci, sigma) %>%
  mutate(IQRpheno = IQR(phenomean),
         outlier_pheno_upper = quantile(phenomean, probs=c( .75), na.rm = FALSE)+1.5*IQRpheno,
         outlier_pheno_lower = quantile(phenomean, probs=c( .25), na.rm = FALSE)-1.5*IQRpheno,
         outlier = phenomean < outlier_pheno_lower | phenomean > outlier_pheno_upper
        ) %>%
  filter(!outlier) -> d_com_filtered


ggplot(d_com_filtered %>% mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = estR)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("IQR_estR_dist.png", plt_dist, width = 10, height = 10)

ggplot(d_com_filtered %>% mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("IQR_pheno_dist.png", plt_dist, width = 10, height = 10)


saveRDS(d_com_filtered, "d_com_prefiltered.RDS")
d_com_filtered <- readRDS("d_com_prefiltered.RDS")

# Split into adapted/maladapted subgroups

d_com_adapted <- d_com_filtered %>%
  group_by(seed, model, nloci, sigma) %>%
  mutate(adapted = between(phenomean, 1.9, 2.1)) %>%
  filter(!(gen == 51950 & adapted == F)) %>%
  filter(any(gen == 51950)) %>%
  ungroup()

# Check they are adapted
ggplot(d_com_adapted %>% 
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  filter(gen == 51950) %>%
  mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("adapted_pheno_dist.png", plt_dist, width = 10, height = 10)


# Filter data into group of those who didn't make it: the groups where 
# any of the phenotypes means aren't within (1.9, 2.1) at any time
d_com_maladapted <- d_com_filtered %>% 
  group_by(seed, model, nloci, sigma) %>%
  mutate(adapted = between(phenomean, 1.9, 2.1)) %>%
  filter(any(gen > 51800)) %>%
  filter(all(adapted == F)) %>%
  ungroup()

# Check they are all maladapted
ggplot(d_com_maladapted %>% 
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("maladapted_pheno_dist.png", plt_dist, width = 10, height = 10)

# Filter data into group of those who made it and then kept moving: the groups where 
# any of the phenotypes means are within (1.9, 2.1) sometime before generation 51800
# also check that the last value isn't in (1.9, 2.1)
d_com_wasadapted <- d_com_filtered %>%
  group_by(seed, model, nloci, sigma) %>%
  mutate(adapted = between(phenomean, 1.9, 2.1)) %>%
  filter(any(gen > 51800)) %>%
  filter(all( (gen < 51800 & adapted | gen < 51800 & !adapted | 
    (gen >= 51800 & !adapted)))) %>%
  ungroup()



saveRDS(d_com_adapted, "d_comh2_prefiltered_adapted.RDS")
saveRDS(d_com_maladapted, "d_comh2_prefiltered_maladapted.RDS")
saveRDS(d_com_wasadapted, "d_comh2_prefiltered_wasadapted.RDS")



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

# Check they are adapted
ggplot(d_com_adapted %>% 
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  filter(gen > 51800) %>%
  mutate(nloci = as_factor(nloci), sigma = as_factor(sigma)), 
  aes(x = nloci, y = phenomean, colour = model)) +
  facet_grid(.~sigma) +
  geom_boxplot() -> plt_dist

ggsave("adapted_pheno_dist.png", plt_dist, width = 10, height = 10)


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