library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R")

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

d_h2_mkr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr.csv"), header = F,
                          col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                        "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                        "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                        "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                        "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                        "h2_KXZ"))

d_h2_mrr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mrr_pt1 <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr_pt1.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mkr_pt1 <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr_pt1.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

# Add on optPerc
d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg <- d_qg %>% select(gen, seed, modelindex, optPerc) %>% filter(gen >= 49500)


# Combine
d_h2_mkr <- rbind(d_h2_mkr, d_h2_mkr_pt1) %>% distinct()
d_h2_mrr <- rbind(d_h2_mrr, d_h2_mrr_pt1) %>% distinct()

d_h2_mkr$method <- "mkr"
d_h2_mrr$method <- "mrr"

d_h2 <- rbind(d_h2_mkr, d_h2_mrr)

# Clean data
d_h2 <- d_h2 %>%
  distinct(gen, seed, modelindex, method, .keep_all = T) %>%
  mutate(modelindex = as.factor(modelindex),
         seed = as.factor(seed)) %>%
  drop_na(VA_Z)

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg, by = c("gen", "seed", "modelindex"))

d_h2 <- AddCombosToDF(d_h2)

# Distribution, how different are the estimates
ggplot(d_h2 %>% 
         select(gen, seed, modelindex, optPerc, h2_Z, method) %>% drop_na() %>%
         distinct() %>%
         pivot_wider(names_from = method, values_from = h2_Z), 
       aes(x = mkr, y = mrr)) +
  geom_point(shape = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = TeX("Kernel regression heritability $(h^2)$"), 
       y = TeX("Ridge regression heritability $(h^2)$")) +
  theme_bw() +
  theme(text = element_text(size = 14))

ggplot(d_h2 %>%
         distinct(), 
       aes(x = method, y = h2_Z)) +
  geom_boxplot() +
  labs(x = TeX("Heritability estimation method"), 
       y = TeX("Narrow-sense heritability $(h^2)$")) +
  theme_bw() +
  theme(text = element_text(size = 14))

# The two estimates are very different: ridge regression is very biased towards
# high heritability, the kernel method is biased towards lower heritability
# Can we trust either? Maybe a parent-offspring regression is the safest bet
# Kernel regression seems most accurate: ridge is almost always at 0 heritability


boxplot(d_h2$VA_Z)
# Remove outliers, summarise
d_h2_sum <- d_h2 %>% 
  group_by(optPerc, model, tau, r, method) %>%
  filter(VA_Z < 50) %>% 
  summarise(meanH2Z = mean(h2_Z, na.rm = T),
            seH2Z = se(h2_Z, na.rm = T),
            meanVAZ = mean(VA_Z, na.rm = T),
            seVAZ = se(VA_Z, na.rm = T))

# Number of loci doesn't seem to affect it too much, average across
ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr"),
       aes(x = optPerc, y = h2_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(shape = 1, dodge.width = 1) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr"), 
             aes(x = optPerc, y = meanH2Z, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(1)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Narrow-sense heritability $(h^2)$"),
       colour = "Model") +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Additive variance
# Again, nloci not too important
ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr"),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(shape = 1, dodge.width = 1) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr"), 
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(1)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Infinitesimal expects zero change in additive variance due to selection
# So see how much it changes between timepoints
d_h2 %>%
  group_by(model, tau, r, nloci, method) %>%
  summarise(totalDeltaVA = )


ggplot(d_h2 %>% 
         group_by(optPerc, model, tau, r, method) %>% 
         summarise(meanVAZ = mean(VA_Z, na.rm = T),
                   seVAZ = se(VA_Z, na.rm = T)) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr"),
       aes(x = optPerc, y = meanVAZ, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  geom_errorbar(aes(ymin = meanVAZ - seVAZ, ymax = meanVAZ + seVAZ),
                width = 0.3) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_colour_paletteer_d("nationalparkcolors::Badlands",
                           labels = c("K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

