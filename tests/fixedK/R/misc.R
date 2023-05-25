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


ggplot(d_adapted %>% filter(modelindex == 2, seed == 2206280217, gen > 49000), 
       aes(x = gen, y = phenomean)) +
  geom_line() +
  theme_bw()

ggplot(d_adapted %>% filter(modelindex == 1, seed == 1448101263, gen > 49000), 
       aes(x = gen, y = phenomean)) +
  geom_line() +
  theme_bw()
