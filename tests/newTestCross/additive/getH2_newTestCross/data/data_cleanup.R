
library(tidyverse)

path <- "/mnt/d/SLiMTests/tests/newTestCross/additive/"

# Read in data
d_qg <- read_csv(paste0(path, "slim_qg.csv"), col_names = F)
d_h2 <- read_csv(paste0(path, "getH2_newTestCross/data/out_h2.csv"), col_names = F)
d_muts <- read_csv(paste0(path, "slim_muts.csv"), col_names = F)                                        
names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
                 "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                 "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                 "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")
names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", 
                   "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")
                                                                                                                   
head(d_muts)
head(d_qg)
head(d_h2)

# Factorise seed and modelindex
d_muts$seed <- as.factor(d_muts$seed)
d_muts$modelindex <- as.factor(d_muts$modelindex)
d_qg$modelindex <- as.factor(d_qg$modelindex)
d_qg$seed <- as.factor(d_qg$seed)
d_h2$seed <- as.factor(d_h2$seed)
d_h2$modelindex <- as.factor(d_h2$modelindex)

# Combine the data frames
d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex"))
d_combined_full <- full_join(d_combined, d_h2, by = c("gen", "seed", "modelindex"))

# View the combined data and write to file
head(d_combined_full) 
colnames(d_combined_full)
d_combined_full[1,] %>% as_tibble() %>% print(width=Inf)
d_combined_full[800000,] %>% as_tibble() %>% print(width=Inf)
data.table::fwrite(d_combined_full, "d_combined_full.csv")

# Write only the data after burn-in
d_combined_after <- d_combined_full %>% filter(gen >= 50000)                                                                       
nrow(d_combined_after)
nrow(d_combined_full)
data.table::fwrite(d_combined_after, "d_combined_after.csv")


# Write summary data for mean frequencies

d_com_end %>% filter(fixedEffect == -1) %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                          h2D = cut(H2.D.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                          h2AA = cut(H2.AA.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                          fixTime = fixGen - originGen) -> d_com_end_ctrl

d_com_end_ctrl %>% mutate(nloci = combos$nloci[.$modelindex],
                          sigma = combos$locisigma[.$modelindex],
                          molTrait = recode_factor(mutType, `1`="Neutral", `2`="Del", `3`="\u03B1", `4`="\u03B2", `5`="KZ", `6`="KXZ")) -> d_com_end_ctrl


d_com_ctrl_h2A_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, molTrait, h2A) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value)
  )

d_com_ctrl_h2D_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, molTrait, h2D) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value)
  )

d_com_ctrl_h2AA_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, molTrait, h2AA) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value)
  )

# include stuff for molecular trait

d_com_end_ctrl %>% pivot_longer(c(aZ, bZ, KZ, KXZ), names_to = "molTraitVal_name", 
                            values_to = "molTraitVal_value") -> d_com_end_ctrl_molTraitVal

d_com_ctrl_sum <- d_com_end_ctrl %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, molTrait, h2A, h2D, h2AA) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq),
            meanFX = mean(value), 
            seFX = se(value)
  )

d_com_ctrl_molTraitVal_sum <- d_com_end_ctrl_molTraitVal %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, molTraitVal_name, h2A, h2D, h2AA) %>%
  summarise(meanMolTraitVal = mean(molTraitVal_value),
            seMolTraitVal = se(molTraitVal_value)
  )
data.table::fwrite(d_com_ctrl_molTraitVal_sum, "d_com_sum_molTraitVal.csv")



data.table::fwrite(d_com_ctrl_h2A_sum, "d_com_sum_h2A.csv")
data.table::fwrite(d_com_ctrl_h2D_sum, "d_com_sum_h2D.csv")
data.table::fwrite(d_com_ctrl_h2AA_sum, "d_com_sum_h2AA.csv")
data.table::fwrite(d_com_ctrl_sum, "d_com_sum.csv")


