
library(tidyverse)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}


path <- "/g/data/ht96/nb9894/newTestCross/additive/moreReps/"
setwd(path)

# Read in data
d_qg <- read_csv(paste0(path, "slim_qg.csv"), col_names = F)
d_h2 <- read_csv(paste0(path, "getH2/out_h2.csv"), col_names = F)
d_muts <- read_csv(paste0(path, "slim_muts.csv"), col_names = F)                                        
names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
                 "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                 "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                 "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")
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
d_h2 <- d_h2 %>%
  mutate(seed = as_factor(seed),
         modelindex = as_factor(modelindex))



# Calculate selection gradients for R = B * VA
# width here is w^2, so need to convert our 1/2w^2 to this format
width_condensed <- 0.05
width_sqrd <- 1/(width_condensed * 2)
opt <- 2
TIME_BETWEEN_SAMPLES <- 50 # 50 generations between each sample

d_phenos <- read_csv(paste0(path, "slim_sampled_pheno.csv"), col_names = F)
d_phenos <- d_phenos %>% pivot_longer(cols = 4:ncol(d_phenos), names_to = NULL, values_to = "phenotype")

colnames(d_phenos) <- c("gen", "seed", "modelindex", "phenotype")

d_phenos <- d_phenos %>%
  mutate(modelindex = as_factor(modelindex),
         seed = as_factor(seed))


calcBeta <- function(S, mu, theta) {
  return(-S * (mu - theta))
}

# Use Morrissey and Goudie 2022 to predict beta

d_phenos <- d_phenos %>%
  group_by(gen, seed, modelindex) %>%
  summarise(S = 1/(width_sqrd + var(phenotype)),
            beta = calcBeta(S, mean(phenotype), opt))

d_h2 <- full_join(d_h2, d_phenos, by = c("gen", "seed", "modelindex"))
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


# Combine the data frames
d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex"))
rm(d_muts, d_qg)
d_combined <- full_join(d_combined, d_h2, by = c("gen", "seed", "modelindex"))
rm(d_h2)
# Deviation of observed R to estimated R
d_combined$devEstR <- d_combined$deltaPheno - d_combined$estR

d_combined_after <- d_combined %>% filter(gen >= 49500)
saveRDS(d_combined_after, paste0(path, "d_combined_after.RDS"))
