library(tidyverse)
library(gghighlight)
library(paletteer)
library(data.table)
library(latex2exp)
library(ggridges)
library(ggh4x)
library(ggbeeswarm)
library(cowplot)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/sanityChecks/R")

source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

d_combos <- read.table("combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "width", "model"))

AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]],
           nloci = d_combos$nloci[as.numeric(levels(modelindex))[modelindex]],
           tau = d_combos$tau[as.numeric(levels(modelindex))[modelindex]],
           r = d_combos$r[as.numeric(levels(modelindex))[modelindex]],
           width = d_combos$width[as.numeric(levels(modelindex))[modelindex]])
}


DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/sanityChecks/"

d_qg <- data.table::fread(paste0(DATA_PATH, "slim_qg.csv"), header = F, 
                          sep = ",", colClasses = c("integer", "factor", "factor", 
                                                    rep("numeric", times = 12)), 
                          col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                        "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                        "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                          fill = T)

d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  filter(gen >= 49500) %>% distinct() %>%
 # mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg


d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg <- AddCombosToDF(d_qg)

d_qg_sum <- d_qg %>% 
  mutate(gen = gen - 50000) %>%
  group_by(gen, model, nloci, tau, r, width) %>%
  summarise(meanPhenomean = mean(phenomean),
            SEPhenomean = se(phenomean),
            sdPhenomean = sd(phenomean),
            meanPhenovar = mean(phenovar),
            sdPhenovar = sd(phenovar))

# phenotype: looks reasonable - 
# really very little difference between 0.5 and 0.1 recombination, same with 
# 1e-10 and 0
# neutral model behaves as expected for 1024 loci, little change over time
# stronger selection makes adaptation faster
ggplot(d_qg %>% filter(tau == 0.0125),
       aes(x = gen, y = phenomean, colour = model, 
           group = interaction(modelindex, as.factor(seed)))) +
  facet_grid(r~width) +
  geom_line() +
  gghighlight(seed == sample(unique(d_qg$seed), 1), calculate_per_facet = T,
              unhighlighted_params = list(colour = NULL, alpha = 0.1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 0.0125") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Selection strength", 
                                         breaks = NULL, labels = NULL)) +
  #coord_cartesian(ylim = c(0, 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Generations post-optimum shift", y = "Population mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))
ggsave("plt_neu_phenomean.png", device = png, width = 9, height = 4)

ggplot(d_qg %>% filter(tau == 1.25),
       aes(x = gen, y = phenomean, colour = model, 
           group = interaction(modelindex, as.factor(seed)))) +
  facet_grid(r~width) +
  geom_line() +
  gghighlight(seed == sample(unique(d_qg$seed), 1), calculate_per_facet = T,
              unhighlighted_params = list(colour = NULL, alpha = 0.1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 1.25") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Selection strength", 
                                         breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Generations post-optimum shift", y = "Population mean phenotype", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))

ggplot(d_qg %>% filter(tau == 0.0125),
       aes(x = gen, y = phenovar, colour = model, 
           group = interaction(modelindex, as.factor(seed)))) +
  facet_grid(r~width) +
  geom_line() +
  gghighlight(seed == sample(unique(d_qg$seed), 1), calculate_per_facet = T,
              unhighlighted_params = list(colour = NULL, alpha = 0.1)) +
  ggtitle("Tau = 0.0125") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Selection strength", 
                                         breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Generations post-optimum shift", y = "Population phenotypic variance", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))
ggsave("plt_neu_phenovar.png", device = png, width = 9, height = 4)

ggplot(d_qg %>% filter(tau == 1.25),
       aes(x = gen, y = phenovar, colour = model, 
           group = interaction(modelindex, as.factor(seed)))) +
  facet_grid(r~width) +
  geom_line() +
  gghighlight(seed == sample(unique(d_qg$seed), 1), calculate_per_facet = T,
              unhighlighted_params = list(colour = NULL, alpha = 0.1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  ggtitle("Tau = 1.25") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(labels = scales::comma, 
                     sec.axis = sec_axis(~ ., name = "Selection strength", 
                                         breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Generations post-optimum shift", y = "Population phenotypic variance", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))

# epistasis figures
d_epi_means <- read.table(paste0(DATA_PATH, "d_epi_mean.csv"), header = F, sep = ",",
                          col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                        "meanEW", "sdEW", "count"))

d_epi_means_sbst <- d_epi_means %>% 
  distinct()

d_epi_means_plt <- AddCombosToDF(d_epi_means_sbst %>% 
                                   mutate(modelindex = as.factor(modelindex)))

# Box plots: fitness epistasis means nothing in the neutral case
ggplot(d_epi_means_plt, aes(x = model, y = meanEW)) +
  facet_grid(.~optPerc) +
  geom_boxplot() +
  labs(x = "Model", y = "Average fitness epistasis") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14))

# Zoom in, ignore big outlier
ggplot(d_epi_means_plt %>% filter(meanEP < 3) %>%
         mutate(r_title = "Recombination rate",
                width_title = "Selection strength"), 
       aes(x = model, y = meanEP, colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_boxplot() +
  labs(x = "Model", y = "Average trait epistasis", colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom")
ggsave("plt_neu_ep.png", device = png)

# frequency adjusted
d_epi_freq_means <- read.table(paste0(DATA_PATH, "d_epi_freqweight_mean.csv"), header = F, sep = ",",
                               col.names = c("optPerc", "modelindex", "meanEP", "sdEP",
                                             "meanEW", "sdEW", "count"))

# Get pairwise differences between in mean epistasis between models 
d_epi_freq_means_sbst <- d_epi_freq_means %>% 
  distinct()

d_epi_freq_means_plt <- AddCombosToDF(d_epi_freq_means_sbst %>% 
                                        mutate(modelindex = as.factor(modelindex)))


d_epi_freq_means_plt$model <- as.factor(d_epi_freq_means_plt$model)
d_epi_freq_means_plt$tau <- factor(d_epi_freq_means_plt$tau, 
                                   levels = c("0.0125", "0.125", "1.25"))
d_epi_freq_means_plt$width_title <- "Selection strength"
d_epi_freq_means_plt$tau_title <- "Mutational effect variance"
d_epi_freq_means_plt$r_title <- "Recombination rate (log10)"

ggplot(d_epi_freq_means_plt %>% filter(meanEP < 3), 
       aes(x = as.factor(tau), 
           y = meanEP, colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_boxplot(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanEP - sdEP/sqrt(count), 
                    ymax = meanEP + sdEP/sqrt(count)), position = position_dodge(0.9)) +
  labs(x = "Mutational effect variance", y = "Average trait epistasis", colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Density
d_epi_dens <- fread(paste0(DATA_PATH, "d_epi_density.csv"), header = F, sep = ",",
                    col.names = c("optPerc", "modelindex", "mutType_ab", 
                                  "wa_x", "wa_y", "wb_x", "wb_y", "wab_x",
                                  "wab_y", "Pa_x", "Pa_y", "Pb_x", "Pb_y",
                                  "Pab_x", "Pab_y",
                                  "ew_x", "ew_y", "ep_x", "ep_y" 
                    ))

d_epi_dens <- d_epi_dens %>%
  distinct() 

d_epi_dens <- AddCombosToDF(d_epi_dens %>% 
                              mutate(modelindex = as.factor(modelindex)))


ggplot(d_epi_dens %>% filter(model == "K"), 
       aes(x = ep_x, y = as.factor(r), height = ep_y)) +
  facet_grid(optPerc~mutType_ab) +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  geom_density_ridges(stat = "identity", scale = 0.95, alpha = 0.4)

# LD: LD sign is meaningless for neutral models
d_ld <- read.table(paste0(DATA_PATH, "out_LD.csv"), header = F,
                   col.names = c("gen", "seed", "modelindex", 
                                 "meanD", "sdD", "nD",
                                 "nDP", "nDN", "nDHalf", 
                                 paste0("n", 1:20)))

d_ld_freq <- data.table::fread(paste0(DATA_PATH, "out_LDf.csv"), header = F,
                               col.names = c("gen", "seed", "modelindex", "freqBin",
                                             "meanD", "sdD", "nD",
                                             "nDP", "nDN", "nDHalf",
                                             paste0("n", 1:20)))

d_ld <- d_ld %>% 
  mutate(seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_ld <- left_join(d_ld, d_qg, by = c("gen", "seed", "modelindex"))

d_ld <- d_ld %>%
  group_by(optPerc, model, width, tau, r) %>%
  mutate(propDP = nDP / nD,
         propDN = nDN / nD)

d_ld_sum <- d_ld %>%
  group_by(optPerc, model, width, tau, r) %>%
  summarise_at(vars(-seed,-gen,-modelindex), list(mean = mean, sd = sd), na.rm = T)

# Outliers: histogram of all estimates
bins <- seq(-0.25, 0.25, length.out = 21)

d_ld_dist_hist <- d_ld %>% select(gen, seed, optPerc, model, width, tau, r, 10:29) %>%
  pivot_longer(cols = matches("n[0-9]"), names_to = "col", values_to = "count") %>%
  group_by(gen, seed, optPerc, model, width, tau, r) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 r_title = "Recombination rate (log10)",
                                 width_title = "Selection strength",
                                 tau_title = "Mutational effect size variance") %>%
         filter(tau == 1.25), 
       aes(x = col, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_boxplot(position = position_identity()) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
    #position = position_dodge(width = 0.9)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  #scale_x_continuous(labels = bins[7:15]) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom")

ggplot(d_ld_dist_hist %>% mutate(col = bins[as.numeric(str_extract(col, "[0-9]+"))],
                                 r_title = "Recombination rate",
                                 width_title = "Selection strength",
                                 tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")), 
       aes(x = col, y = prop, colour = model, group = interaction(col, model))) +
  facet_nested(r_title + r ~ width_title + width + model) +
  geom_boxplot(position = position_identity(), outlier.shape = 1,
               outlier.alpha = 0.2) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(group = model, colour = model)
    #position = position_dodge(width = 0.9)
  ) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-"), guide = "none") +
  labs(x = "D", y = "Proportion of estimates", colour = "Model") +
  #scale_x_continuous(labels = bins[7:15]) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom")
ggsave("plt_neu_ld.png", device = png, width = 14, height = 5)



# Proportion of D estimates > 0.05 or < 0.05: i.e. there is non-zero LD
ggplot(d_ld %>%
         ungroup() %>%
         complete(optPerc, model, width, tau, r) %>%
         group_by(optPerc, model, width, tau, r) %>%
         mutate(propDP = nDP / nD,
                propDN = nDN / nD) %>%
         group_by(optPerc, model, width, r) %>%
         summarise_at(vars(-seed,-gen,-modelindex), 
                      list(mean = mean, sd = se), na.rm = T) %>%
         mutate(r_title = "Recombination rate (log10)",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance"),
       aes(x = optPerc, y = propDP_mean,
           colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = propDP_mean - propDP_sd, ymax = propDP_mean + propDP_sd),
                position = position_dodge(0.9, preserve = 'single')) +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Progress to optimum", y = "Proportion of positive D estimates", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))

# Negative LD
ggplot(d_ld %>%
         ungroup() %>%
         complete(optPerc, model, width, tau, r) %>%
         group_by(optPerc, model, width, tau, r) %>%
         mutate(propDP = nDP / nD,
                propDN = nDN / nD) %>%
         group_by(optPerc, model, width, r) %>%
         summarise_at(vars(-seed,-gen,-modelindex), 
                      list(mean = mean, sd = se), na.rm = T) %>%
         mutate(r_title = "Recombination rate (log10)",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance"),
       aes(x = optPerc, y = propDN_mean,
           colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = propDN_mean - propDN_sd, ymax = propDN_mean + propDN_sd),
                position = position_dodge(0.9, preserve = 'single')) +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~ ., name = "Recombination rate (log10)", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = "Progress to optimum", y = "Proportion of negative D estimates", 
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))


# H2

# Cowplot 1.1.3 bug: won't get legend, this fixes
get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

d_h2_mkr <- data.table::fread(paste0(DATA_PATH, "out_h2_mkr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mrr <- data.table::fread(paste0(DATA_PATH, "out_h2_mrr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_a_b", "CVA_Z_KZ", "CVA_a_KZ", "CVA_b_KZ",
                                            "CVA_Z_KXZ", "CVA_a_KXZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mkr$method <- "mkr"
d_h2_mrr$method <- "mrr"

d_h2 <- rbind(d_h2_mkr, d_h2_mrr)

d_h2 <- d_h2 %>%
  distinct(gen, seed, modelindex, method, .keep_all = T) %>%
  mutate(modelindex = as.factor(modelindex),
         seed = as.factor(seed)) %>%
  drop_na(VA_Z) %>% distinct()

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg, by = c("gen", "seed", "modelindex"))

d_h2 <- AddCombosToDF(d_h2)

d_h2$VA_Z <- as.numeric(d_h2$VA_Z)

table(d_h2$model)

# Distribution, how different are the estimates
ggplot(d_h2 %>% 
         select(gen, seed, modelindex, optPerc, h2_Z, method, model) %>% drop_na() %>%
         distinct() %>%
         pivot_wider(names_from = method, values_from = h2_Z), 
       aes(x = mkr, y = mrr, colour = model)) +
  geom_point(shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  labs(x = TeX("Kernel regression heritability $(h^2)$"), 
       y = TeX("Ridge regression heritability $(h^2)$"),
       colour = "Model") +
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

# In this case the models are more similar, can use both?

boxplot(d_h2$VA_Z)

d_h2_all <- d_h2

# Outlier detection
library(DMwR2)

lofscores <- lofactor(d_h2$VA_Z, 10)
threshold <- 1.6
outliers <- lofscores > threshold

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score", ylim = c(0, 10))
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)

# filter out outliers
boxplot(d_h2[!outliers,]$VA_Z)

d_h2 <- d_h2[!outliers,]

d_h2_sum <- d_h2 %>%
  group_by(optPerc, model, tau, r, width) %>%
  summarise(meanH2Z = mean(h2_Z, na.rm = T),
            seH2Z = se(h2_Z, na.rm = T),
            meanVAZ = mean(VA_Z, na.rm = T),
            seVAZ = se(VA_Z, na.rm = T))

ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance") %>%
         filter(tau == 0.0125),
       aes(x = optPerc, y = h2_Z, colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)",
                      width_title = "Selection strength",
                      tau_title = "Mutational effect size variance") %>%
               filter(tau == 0.0125),
             aes(x = optPerc, y = meanH2Z, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Narrow-sense heritability $(h^2)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Additive variance
# Again, nloci not important
ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance") %>%
         filter(tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      width_title = "Selection strength",
                      tau_title = "Mutational effect size variance") %>%
               filter(tau == 0.0125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_sml
ggsave("plt_neu_va.png", plt_add_va_sml, device = png, width = 9, height = 4)

ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance") %>%
         filter(tau == 1.25),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + r ~ width_title + width) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate",
                      width_title = "Selection strength",
                      tau_title = "Mutational effect size variance") %>%
               filter(tau == 1.25),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 2.5)) +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_lrg

leg <- get_legend(plt_add_va_lrg)

plt_add_va <- plot_grid(plt_add_va_sml + theme(legend.position = "none"),
                        plt_add_va_lrg + theme(legend.position = "none"),
                        ncol = 1, labels = "AUTO")

plt_add_va <- plot_grid(plt_add_va,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_add_va
ggsave("plt_va.png", device = png, bg = "white",
       width = 560*4, height = 980*4, units = "px")

# compare VA to VI
d_epi_means_plt <- d_epi_means_plt %>% 
  mutate(VI = sdEP^2,
         optPerc = as.factor(optPerc))
  

d_VA_VI <- inner_join(d_h2 %>% 
                        group_by(modelindex, model, 
                                 nloci, tau, r, width) %>%
                        summarise(VA = mean(VA_Z)), 
                      d_epi_means_plt %>%
                        group_by(modelindex, model, 
                               nloci, tau, r, width) %>%
                        summarise(VI = mean(VI)),  # mean across time points
                      by = c("modelindex", 
                             "model", "nloci", "tau", "r", "width"))

# Proportions of VA vs VI
ggplot(d_VA_VI %>% 
         rowwise() %>%
         mutate(propVA = VA / ( VI),
                propVI = VI / ( VA)) %>%
         mutate(r_title = "Recombination rate",
                width_title = "Selection strength",
                tau_title = "Mutational effect size variance"),
       aes(x = as.factor(r), y = propVI, colour = model)) +
  facet_nested(tau_title + tau ~ width_title + width, scales = "free") +
  geom_point(position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  #coord_cartesian(ylim = c(0, 0.0001)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
