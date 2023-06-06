source("wrangle_data.R")

# see how segregating variation affects adaptation
singlestep_mdls_nar <- d_fix_ranked %>% 
  group_by(seed) %>% 
  filter(rank != 0) %>% 
  filter(n() == 1)

singlestep_mdls_add <- d_fix_ranked_add %>% 
  group_by(seed) %>% 
  filter(rank != 0) %>% 
  filter(n() == 1)


View(d_fix_adapted %>% filter(modelindex == 2, seed == 2206280217, gen > 49500))
View(d_fix_adapted %>% filter(modelindex == 1, seed == 1448101263, gen > 49500))

View(d_muts %>% filter(modelindex == 2, seed == 2206280217, 
                       Freq < 1, gen > 49500, mutID == 4735407) %>% distinct())
View(d_muts %>% filter(modelindex == 1, seed == 1448101263, 
                       Freq < 1, gen > 49500) %>% distinct())


ggplot(d_adapted %>% filter(modelindex == 2, seed == 2206280217, gen > 49000) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = phenomean)) +
  geom_line() +
  theme_bw() +
  labs(x = "Generations post-optimum shift", y = "Phenotype mean") +
  theme(text = element_text(size = 16))


ggplot(d_adapted %>% filter(modelindex == 1, seed == 1448101263, gen > 49000) %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = phenomean)) +
  geom_line() +
  theme_bw() +
  labs(x = "Generations post-optimum shift", y = "Phenotype mean") +
  theme(text = element_text(size = 16))

ggplot(d_com_adapted %>% filter(modelindex == 2, seed == 2270695859, 
                                gen >= 49000, mutID == 4607309) %>% 
         distinct() %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = Freq)) +
  geom_line() +
  theme_bw() +
  labs(x = "Generations post-optimum shift", y = "Allele frequency") +
  theme(text = element_text(size = 16))

d_seg_ranked %>% filter(seed == 2270695859, 
                        gen >= 49000, mutID == 4607309) %>% select(avFit, avFit_AA, h)

ggplot(d_com_adapted %>% filter(modelindex == 1, seed == 1448101263, 
                                gen >= 49000, (mutID == 4624809)) %>% 
         distinct() %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = Freq)) +
  geom_line() +
  theme_bw() +
  scale_colour_paletteer_d("ggsci::alternating_igv", guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Allele frequency", colour = "Mutation") +
  theme(text = element_text(size = 16))

d_seg_ranked %>% filter(modelindex == 2, seed == 2206280217, 
                        gen >= 49000, mutID == 4735407)


# see what's going on with fitness effect estimation for NAR models
# 
write.table(d_seg_ranked %>% filter(seed == 2016338465, gen == 59700, mutType == 3) %>% ungroup() %>% 
              dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ.csv", sep = ",", col.names = F, row.names = F)
d_dat_aZ_absence <- d_seg_ranked %>% filter(seed == 2016338465, gen == 59700, mutType == 3)
d_dat_aZ_absence$aZ <- exp(log(d_dat_aZ_absence$aZ) - d_dat_aZ_absence$value)

write.table(d_dat_aZ_absence %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ_absence.csv", sep = ",", col.names = F, row.names = F)

d_popfx_aZ_absence <- runLandscaper("d_grid_aZ_absence.csv", "data_popfx_aZ_absence.csv", 0.05, 2, 8)
d_popfx_aZ <- runLandscaper("d_grid_aZ.csv", "data_popfx_aZ.csv", 0.05, 2, 8)
d_popfx_aZ$fitness - d_popfx_aZ_absence$fitness


d_seg_ranked %>%
  filter(rank > 0, gen == 59700, seed == 2016338465) %>%
  mutate(curEffect = abs(.data[["avFX"]])) %>%
  summarise(segEffectSum = sum(curEffect))

d_fix_nar %>%
  filter(gen > 49500, seed == 2016338465) %>% 
  group_by(mutType, seed, gen) %>%
  mutate(fixedSum = sum(.data[["avFX"]]))

d_muts_adapted %>%
  filter(Freq == 1, seed == 2016338465) %>%
  distinct(mutID, .keep_all = T)

View(d_muts_adapted %>% filter(Freq == 1) %>% distinct())

# Verify -7.452538 big effect for fixation in seed 489849073
test <- d_fix_adapted %>% ungroup() %>% mutate(r = rownames(.)) %>% filter(modelindex == 2, seed == 119202918)

# First get matched mol trait data for generations where we have fixations
d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(test$gen, test$seed, test$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

test <- inner_join(test, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

test <- test %>%
  group_by(gen, seed) %>%
  mutate(fixEffectSum_aZ = 2 * sum(test[test$gen <= cur_group()$gen &
                                          test$mutType == 3 &
                                          test$seed == cur_group()$seed,]$value),
         fixEffectSum_bZ = 2 * sum(test[test$gen <= cur_group()$gen &
                                          test$mutType == 4 &
                                          test$seed == cur_group()$seed,]$value))
# Transform to exp scale
test$fixEffectSum_aZ <- exp(test$fixEffectSum_aZ)
test$fixEffectSum_bZ <- exp(test$fixEffectSum_bZ)

# Get phenotypes with the mutation
write.table(test %>% ungroup() %>%
              dplyr::select(fixEffectSum_aZ, fixEffectSum_bZ, KZ, KXZ), 
            "d_grid.csv", sep = ",", col.names = F, row.names = F)
d_popfx <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8)


# For fixed comparisons:
# AA = 1+s; Aa = 1+hs; aa = 1
# AA = d_popfx$fixEffectSum
# Aa = d_popfx$fixEffectSum - value
# aa = d_popfx$fixEffectSum - 2 * value

test_withoutFX_aZ <- test %>% filter(mutType == 3)
test_withoutFX_bZ <- test %>% filter(mutType == 4)
test_withoutFX_aZ$aZ <- exp(log(test_withoutFX_aZ$fixEffectSum_aZ) - test_withoutFX_aZ$value)
test_withoutFX_aZ$bZ <- exp(log(test_withoutFX_aZ$fixEffectSum_bZ))
test_withoutFX_bZ$aZ <- exp(log(test_withoutFX_bZ$fixEffectSum_aZ))
test_withoutFX_bZ$bZ <- exp(log(test_withoutFX_bZ$fixEffectSum_bZ) - test_withoutFX_bZ$value)

# Homozygous estimation - for dominance calculation
test_withoutFX_aZ$aZ_aa <- exp(log(test_withoutFX_aZ$fixEffectSum_aZ) - 2 * test_withoutFX_aZ$value)
test_withoutFX_aZ$bZ_aa <- exp(log(test_withoutFX_aZ$fixEffectSum_bZ))
test_withoutFX_bZ$aZ_aa <- exp(log(test_withoutFX_bZ$fixEffectSum_aZ))
test_withoutFX_bZ$bZ_aa <- exp(log(test_withoutFX_bZ$fixEffectSum_bZ) - 2 * test_withoutFX_bZ$value)

test_withoutFX <- rbind(test_withoutFX_aZ, test_withoutFX_bZ)

write.table(test_withoutFX %>% ungroup() %>% 
              dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid.csv", sep = ",", col.names = F, row.names = F)
Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8)

write.table(test_withoutFX %>% ungroup() %>% 
              dplyr::select(aZ_aa, bZ_aa, KZ, KXZ), 
            "d_grid.csv", sep = ",", col.names = F, row.names = F)
aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8)

# Get the effect size by taking away the phenotype missing that fixation
test$avFX <- d_popfx$pheno - Aa$pheno
test$avFit <- d_popfx$fitness - Aa$fitness
test$avFX_AA <- d_popfx$pheno - aa$pheno
test$avFit_AA <- d_popfx$fitness - aa$fitness
test$wAA <- d_popfx$fitness
test$wAa <- Aa$fitness
test$waa <- aa$fitness



View(d_com_adapted %>% filter(seed == 489849073, gen == 53950) %>% distinct())
write.table(d_seg_ranked %>% filter(seed == 2016338465, gen == 59700, mutType == 3) %>% ungroup() %>% 
              dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ.csv", sep = ",", col.names = F, row.names = F)
d_dat_aZ_absence <- d_seg_ranked %>% filter(seed == 2016338465, gen == 59700, mutType == 3)
d_dat_aZ_absence$aZ <- exp(log(d_dat_aZ_absence$aZ) - d_dat_aZ_absence$value)

write.table(d_dat_aZ_absence %>% ungroup() %>% dplyr::select(aZ, bZ, KZ, KXZ), 
            "d_grid_aZ_absence.csv", sep = ",", col.names = F, row.names = F)

d_popfx_aZ_absence <- runLandscaper("d_grid_aZ_absence.csv", "data_popfx_aZ_absence.csv", 0.05, 2, 8)
d_popfx_aZ <- runLandscaper("d_grid_aZ.csv", "data_popfx_aZ.csv", 0.05, 2, 8)
d_popfx_aZ$fitness - d_popfx_aZ_absence$fitness
