
library(tidyverse)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}


path <- "/mnt/d/SLiMTests/tests/newTestCross/"

# Read in data
d_qg_old <- read_csv(paste0(path, "getH2_newTestCross/data/slim_qg.csv"), col_names = F)
d_qg_new <- read_csv(paste0(path, "moreReps/getH2_newTestCross/data/slim_qg_new.csv"), col_names = F)
d_h2_old <- read_csv(paste0(path, "getH2_newTestCross/data/out_h2.csv"), col_names = F)
d_h2_new <- read_csv(paste0(path, "moreReps/getH2_newTestCross/data/out_h2_new.csv"), col_names = F)

d_muts_old <- read_csv(paste0(path, "getH2_newTestCross/data/slim_muts.csv"), col_names = F)
d_muts_new <- read_csv(paste0(path, "moreReps/getH2_newTestCross/data/slim_muts.csv"), col_names = F)  

d_qg <- rbind(d_qg_old, d_qg_new)
d_h2 <- rbind(d_h2_old, d_h2_new)
d_muts <- rbind(d_muts_old, d_muts_new)

rm(d_qg_old, d_qg_new, d_h2_old, d_h2_new, d_muts_new, d_muts_old)

names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR",
                  "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                  "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                  "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")
names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")                                                                                                                   
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


# Calculate selection gradients for R = B * VA
# width here is w^2, so need to convert our 1/2w^2 to this format
width_condensed <- 0.05
width_sqrd <- 1/(width_condensed * 2)
opt <- 2
TIME_BETWEEN_SAMPLES <- 50 # 50 generations between each sample

d_phenos_new <- read_csv(paste0(path, "moreReps/getH2_newTestCross/data/slim_sampled_pheno.csv"), col_names = F)
d_phenos_old <- read_csv(paste0(path, "getH2_newTestCross/data/slim_sampled_pheno.csv"), col_names = F)
d_phenos <- rbind(d_phenos_old, d_phenos_new)
d_phenos <- d_phenos %>% pivot_longer(cols = 4:ncol(d_phenos), names_to = NULL, values_to = "phenotype")
rm(d_phenos_new, d_phenos_old)

colnames(d_phenos) <- c("gen", "seed", "modelindex", "phenotype")

# d_phenos <- distinct(d_phenos, gen, seed, modelindex, .keep_all = T)

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
d_h2$estR <- (d_h2$H2.A.Estimate * d_h2$beta) * TIME_BETWEEN_SAMPLES

d_h2 <- d_h2 %>%
  filter(!is.na(AIC)) %>%
  distinct(.keep_all = T)

saveRDS(d_h2, "./checkpoint/d_h2.RDS")
saveRDS(d_muts, "./checkpoint/d_muts.RDS")
saveRDS(d_qg, "./checkpoint/d_qg.RDS")

d_h2 <- readRDS("./checkpoint/d_h2.RDS")
d_muts <- readRDS("./checkpoint/d_muts.RDS")
d_qg <- readRDS("./checkpoint/d_qg.RDS")

# Combine the data frames
d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex"))
rm(d_muts, d_qg)
d_combined <- full_join(d_combined, d_h2, by = c("gen", "seed", "modelindex"))
rm(d_h2)

saveRDS(d_combined, "./checkpoint/d_combined.RDS")

d_combined <- readRDS("./checkpoint/d_combined.RDS")

# Deviation of observed R to estimated R
d_combined$devEstR <- d_combined$deltaPheno - d_combined$estR


# View the combined data and write to file
head(d_combined) 
colnames(d_combined)
d_combined[1,] %>% as_tibble() %>% print(width=Inf)
d_combined[800000,] %>% as_tibble() %>% print(width=Inf)
data.table::fwrite(d_combined, paste0(path, "moreReps/getH2_newTestCross/data/d_combined_full.csv"))

# Write only the data after burn-in
d_combined_after <- d_combined %>% filter(gen >= 50000)                                                                       
nrow(d_combined_after)
nrow(d_combined_full)
data.table::fwrite(d_combined_after, paste0(path, "moreReps/getH2_newTestCross/data/d_combined_after.csv"))

d_combined_after <- read_csv(paste0(path, "moreReps/getH2_newTestCross/data/d_combined_after.csv"))

combos <- read_delim("../../R/combos.csv", delim = " ", col_names = F)
names(combos) <- c("model", "nloci", "locisigma")

d_combined_after %>% mutate(fixedEffect = combos$model[.$modelindex],
                nloci = combos$nloci[.$modelindex] + 2,
                sigma = combos$locisigma[.$modelindex]) -> d_combined_after


# Write summary data for mean frequencies, effect sizes, and estimated responses to selection
labs <- c("[0,0.25]", "(0.25, 0.5]", "(0.5,0.75]", "(0.75,1]")
d_combined_after %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = labs, include.lowest = T, ordered_result = T),
                            h2D = cut(H2.D.Estimate, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = labs, include.lowest = T, ordered_result = T),
                            h2AA = cut(H2.AA.Estimate, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = labs, include.lowest = T, ordered_result = T),
                            fixTime = fixGen - originGen) -> d_com_end_ctrl


d_com_end_ctrl %>% mutate(nloci = combos$nloci[.$modelindex] + 2,
                          sigma = combos$locisigma[.$modelindex],
                          molTrait = recode_factor(mutType, `1`="Neutral", `2`="Del", `3`="\u03B1", `4`="\u03B2", `5`="KZ", `6`="KXZ")) -> d_com_end_ctrl


d_com_ctrl_h2A_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, molTrait, h2A) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  )

d_com_ctrl_h2D_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, molTrait, h2D) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  )

d_com_ctrl_h2AA_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, molTrait, h2AA) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  )

# include stuff for molecular trait

# d_com_end_ctrl %>% pivot_longer(c(aZ, bZ, KZ, KXZ), names_to = "molTraitVal_name", 
#                             values_to = "molTraitVal_value") -> d_com_end_ctrl_molTraitVal

d_com_ctrl_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, molTrait, h2A, h2D, h2AA) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  )

# d_com_ctrl_molTraitVal_sum <- d_com_end_ctrl_molTraitVal %>%
#   select(-c(fixGen, fixTime)) %>%
#   drop_na() %>%
#   group_by(nloci, sigma, molTraitVal_name, h2A, h2D, h2AA) %>%
#   summarise(meanMolTraitVal = mean(molTraitVal_value),
#             seMolTraitVal = se(molTraitVal_value)
#   )
# data.table::fwrite(d_com_ctrl_molTraitVal_sum, paste0(path, "moreReps/getH2_newTestCross/data/d_com_sum_molTraitVal.csv"))

# Response to selection
path <- "/mnt/d/SLiMTests/tests/newTestCross/"

d_com_resp_h2A <- d_com_end_ctrl %>%
  mutate(h2A = fct_collapse(h2A,
               `≤0.5` = c("[0,0.25]", "(0.25, 0.5]"),
               `>0.5` = c("(0.5,0.75]", "(0.75,1]"))) %>%
  select(-c(fixGen, fixTime)) %>%
  filter(phenomean < 10, fixedEffect == -1) %>%
  drop_na() %>%
  group_by(gen, nloci, h2A) %>%
  summarise(meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  ) %>%
  ungroup() %>%
  group_by(nloci, h2A) %>%
  mutate(expPheno = lag(meanPheno) + meanEstResp,
         errPheno = cumsum(meanDeltaEst))

data.table::fwrite(d_com_resp_h2A, paste0(path, "moreReps/getH2_newTestCross/data/d_com_resp_h2A.csv"))

d_com_resp_h2D <- d_com_end_ctrl %>%
  mutate(h2D = fct_collapse(h2D,
                            `≤0.5` = c("[0,0.25]", "(0.25, 0.5]"),
                            `>0.5` = c("(0.5,0.75]", "(0.75,1]"))) %>%
  select(-c(fixGen, fixTime)) %>%
  filter(phenomean < 10, fixedEffect == -1) %>%
  drop_na() %>%
  group_by(gen, nloci, h2D) %>%
  summarise(meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  ) %>%
  ungroup() %>%
  group_by(nloci, h2D) %>%
  mutate(expPheno = lag(meanPheno) + meanEstResp,
         errPheno = cumsum(meanDeltaEst))

data.table::fwrite(d_com_resp_h2D, paste0(path, "moreReps/getH2_newTestCross/data/d_com_resp_h2D.csv"))

d_com_resp_h2AA <- d_com_end_ctrl %>%
  mutate(h2AA = fct_collapse(h2AA,
                            `≤0.5` = c("[0,0.25]", "(0.25, 0.5]"),
                            `>0.5` = c("(0.5,0.75]", "(0.75,1]"))) %>%
  select(-c(fixGen, fixTime)) %>%
  filter(phenomean < 10, fixedEffect == -1) %>%
  drop_na() %>%
  group_by(gen, nloci, h2AA) %>%
  summarise(meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanEstResp = mean(estR),
            seEstResp = se(estR),
            meanDelta = mean(deltaPheno),
            seDelta = se(deltaPheno),
            meanDeltaEst = mean(devEstR),
            seDeltaEst = se(devEstR)
  ) %>%
  ungroup() %>%
  group_by(nloci, h2AA) %>%
  mutate(expPheno = lag(meanPheno) + meanEstResp,
         errPheno = cumsum(meanDeltaEst))

data.table::fwrite(d_com_resp_h2AA, paste0(path, "moreReps/getH2_newTestCross/data/d_com_resp_h2AA.csv"))



data.table::fwrite(d_com_ctrl_h2A_sum, paste0(path, "moreReps/getH2_newTestCross/data/d_com_sum_h2A.csv"))
data.table::fwrite(d_com_ctrl_h2D_sum, paste0(path, "moreReps/getH2_newTestCross/data/d_com_sum_h2D.csv"))
data.table::fwrite(d_com_ctrl_h2AA_sum, paste0(path, "moreReps/getH2_newTestCross/data/d_com_sum_h2AA.csv"))
data.table::fwrite(d_com_ctrl_sum, paste0(path, "moreReps/getH2_newTestCross/data/d_com_sum.csv"))

  
