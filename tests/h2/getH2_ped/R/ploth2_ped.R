# Plot heritability estimates
library(tidyverse)

se <- function(x) {
  return(sd(x)/sqrt(length(x)))
}

d_h2 <- read_csv("../data/out_h2_ped.csv", col_names = F)

# Header:
names(d_h2) <- c("gen", "seed", "modelindex", "Vgca", "Vsca", 
                 "Ve", "Va", "Vd", "Vg", "H2", "h2")

# Remove duplicates
d_h2 %>% distinct() -> d_h2

# Remove NAs
d_h2 %>% drop_na() -> d_h2

# Recode modelindex to additive/network
d_h2 %>% mutate(model = recode_factor(modelindex, `0`="Additive", `1`="Network")) -> d_h2

# Remove weird values: huge variance
d_h2 %>% filter(Vgca >= 0, Vgca < 10,
                Vsca >= 0, Vsca < 10) -> d_h2

d_h2_sum <- d_h2 %>% 
  group_by(gen, model) %>% 
  summarise(meanH2 = mean(h2),
            seH2 = se(h2),
            meanVA = mean(Va), 
            seVA = se(Va),
            meanVD = mean(Vd),
            seVD = se(Vd),
            meanVE = mean(Ve),
            seVE = se(Ve)
            )

d_h2_sum$gen <- d_h2_sum$gen - 100000

# Distribution
ggplot(d_h2 %>% filter(model == "Additive"), aes(x = h2, fill = model)) +
  geom_density(alpha = 0.4)

# Mean over time
ggplot(d_h2_sum, aes(x = gen, y = meanH2, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2 - seH2), ymax = (meanH2 + seH2), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average narrow-sense heritability", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=20))

# What about other variance components?
ggplot(d_h2_sum, aes(x = gen, y = meanVD, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanVD - seVD), ymax = (meanVD + seVD), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average dominance variance", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=20))


ggplot(d_h2_sum, aes(x = gen, y = meanVA, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanVA - seVA), ymax = (meanVA + seVA), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average additive variance", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=20))


# stacked percent bar chart
ggplot(d_h2_sum %>% pivot_longer(
  cols = c(meanVA, meanVD, meanVE),
  names_to = "varComp", values_to = "prop"), 
  aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Environmental")) +
  facet_grid(.~model) +
  geom_bar(position="fill", stat="identity", width = 1) +
  labs(x = "Generations after optimum shift", y = "Proportion of total phenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 20))


# Plot trait data over time
d_qg <- read_csv("../../data/slim_qg.csv", col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")

# Only plot from data points that are in d_h2
# Match by gen, seed, modelindex
d_h2 <- d_h2 %>% unite("id", gen:modelindex, remove = F)
d_qg <- d_qg %>% unite("id", gen:modelindex, remove = F)

d_qg <- d_qg[d_qg$id %in% d_h2$id,]
d_qg %>% mutate(model = recode_factor(modelindex, `0`="Additive", `1`="Network")) -> d_qg

d_qg <- d_qg %>% filter(phenomean < 10)

d_qg_sum <- d_qg %>% group_by(gen, model) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

d_qg_sum$gen <- d_qg_sum$gen - 100000

ggplot(d_qg_sum, aes(gen, meanPheno, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Model") +
  theme_bw() +
  theme(text = element_text(size=20))
