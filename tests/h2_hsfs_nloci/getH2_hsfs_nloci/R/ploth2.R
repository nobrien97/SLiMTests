# Plot heritability estimates
library(tidyverse)

se <- function(x) {
  return(sd(x)/sqrt(length(x)))
}

d_h2 <- read_csv("../data/new/out_h2.csv", col_names = F)

names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
                 "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                 "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                 "H2.AA.SE", "AIC")

# Remove duplicates
d_h2 %>% distinct() -> d_h2

# Remove NAs
d_h2 %>% drop_na() -> d_h2

# Recode modelindex to additive/network
d_h2 %>% mutate(model = recode_factor(modelindex, `0`="Additive", `1`="Network")) -> d_h2

# Remove weird values: negative variance, huge variance
d_h2 %>% filter(VarA >= 0, VarA < 10,
                VarD >= 0, VarAA >= 0,
                VarR > 0, sum(VarA, VarD, VarAA, VarR) > 0) -> d_h2

d_h2$VarA_scl <- scale(d_h2$VarA)
d_h2$VarD_scl <- scale(d_h2$VarD)
d_h2$VarAA_scl <- scale(d_h2$VarAA)

d_h2_sum <- d_h2 %>%
  group_by(gen, model) %>% 
  summarise(meanH2A = mean(H2.A.Estimate),
            seH2A = se(H2.A.Estimate),
            meanH2D = mean(H2.D.Estimate), 
            seH2D = se(H2.D.Estimate),
            meanH2AA = mean(H2.AA.Estimate),
            seH2AA = se(H2.AA.Estimate),
            meanVarA = mean(VarA),
            seVarA = se(VarA),
            meanVarD = mean(VarD),
            seVarD = se(VarD),
            meanVarAA = mean(VarAA),
            seVarAA = se(VarAA)
            )

d_h2_sum$gen <- d_h2_sum$gen - 50000

# Distribution
ggplot(d_h2, aes(x = H2.A.Estimate, fill = model)) +
  geom_density(alpha = 0.4)

# Mean over time
ggplot(d_h2_sum, aes(x = gen, y = meanH2A, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2A - seH2A), ymax = (meanH2A + seH2A), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average narrow-sense heritability", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=20))

# What about other variance components?
ggplot(d_h2_sum, aes(x = gen, y = meanH2D, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2D - seH2D), ymax = (meanH2D + seH2D), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average dominance variance proportion", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=20))


ggplot(d_h2_sum, aes(x = gen, y = meanH2AA, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2AA - seH2AA), ymax = (meanH2AA + seH2AA), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average AxA variance proportion", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=20))


# stacked percent bar chart
ggplot(d_h2_sum %>% pivot_longer(
  cols = c(meanH2A, meanH2D, meanH2AA),
  names_to = "varComp", values_to = "prop"), 
  aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  facet_grid(.~model) +
  geom_bar(position="fill", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total phenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 20))

ggplot(d_h2_sum %>% pivot_longer(
  cols = c(meanVarA, meanVarD, meanVarAA),
  names_to = "varComp", values_to = "var"),
  aes(x = gen, y = var)) +
  scale_color_viridis_d() +
  facet_grid(varComp~model) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Variance") +
  theme_bw() +
  theme(text = element_text(size = 20))


# Plot trait data over time
d_qg <- read_csv("../data/new/slim_qg.csv", col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")
d_com <- inner_join(d_h2, d_qg, by = c("gen", "seed", "modelindex"))
d_com <- d_com %>% filter(phenomean < 10)

d_com_sum <- d_com %>% group_by(gen, model) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

d_com_sum$gen <- d_com_sum$gen - 50000

ggplot(d_com_sum, aes(gen, meanPheno, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Model") +
  theme_bw() +
  theme(text = element_text(size=20))

# plots of variance components per phenotype range


