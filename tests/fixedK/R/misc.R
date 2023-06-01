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

ggplot(d_com_adapted %>% filter(modelindex == 2, seed == 2206280217, 
                                gen >= 49000, mutID == 4735407) %>% 
         distinct() %>%
         mutate(gen = gen - 50000), 
       aes(x = gen, y = Freq)) +
  geom_line() +
  theme_bw() +
  labs(x = "Generations post-optimum shift", y = "Allele frequency") +
  theme(text = element_text(size = 16))

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
