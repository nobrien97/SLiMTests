library(tidyverse)
library(data.table)
library(latex2exp)
library(paletteer)
library(ggridges)
library(ggh4x)
library(cowplot)
library(ggbeeswarm)
setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/getH2/R")

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

d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, optPerc) %>% filter(gen >= 49500)


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
  drop_na(VA_Z) %>% distinct()

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex"))

d_h2 <- AddCombosToDF(d_h2)

# Counts for each model type: K+ harder to estimate than the other two
table(d_h2$model)

# We have many recombination rates: choose a few
r_subsample <- c(1e-10, 1e-5, 1e-1)

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

# The two estimates are very different: ridge regression is very biased towards
# high heritability, the kernel method is biased towards lower heritability
# Can we trust either? Maybe a parent-offspring regression is the safest bet
# Kernel regression seems most accurate: ridge is almost always at 0 heritability

d_h2 <- d_h2 %>% filter(method == "mkr")

boxplot(d_h2$VA_Z)
d_h2_all <- d_h2

# Detect outliers: Hampel filter
# variance depends on the tau group mainly - check outliers within groups
library(DMwR2)
# lower_bound <- list(3)
# upper_bound <- list(3)
# tau_lvls <- unique(d_h2$tau)
# outlier_ind <- list(3)
# 
# for (i in 1:3) {
#   lower_bound[[i]] <- median(d_h2[d_h2$tau == tau_lvls[i],]$VA_Z) - 3 * mad(d_h2[d_h2$tau == tau_lvls[i],]$VA_Z, constant = 1)
#   upper_bound[[i]] <- median(d_h2[d_h2$tau == tau_lvls[i],]$VA_Z) + 3 * mad(d_h2[d_h2$tau == tau_lvls[i],]$VA_Z, constant = 1)
#   outlier_ind[[i]] <- which(d_h2[d_h2$tau == tau_lvls[i],]$VA_Z < lower_bound[[i]] | d_h2[d_h2$tau == tau_lvls[i],]$VA_Z > upper_bound[[i]])
# }
# 
# # This is quite a lot of data to remove: ~12% are outliers?
# length(unlist(outlier_ind))
# 
# # Remove outliers
# d_h2 <- d_h2 %>%
#   filter((tau == 0.0125 & VA_Z >= lower_bound[[1]] & VA_Z <= upper_bound[[1]])|
#          (tau == 0.125  & VA_Z >= lower_bound[[2]] & VA_Z <= upper_bound[[2]])|
#          (tau == 1.25   & VA_Z >= lower_bound[[3]] & VA_Z <= upper_bound[[3]]))


lofscores <- lofactor(scale(d_h2$VA_Z), 10)
threshold <- 1.6
outliers <- lofscores > threshold

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score")
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)

# filter out outliers
d_h2 <- d_h2[!outliers,]

boxplot(d_h2[!outliers,]$VA_Z)

# summarise
d_h2_sum <- d_h2 %>%
  filter(r %in% r_subsample) %>%
  group_by(optPerc, model, tau, r, method) %>%
  summarise(meanH2Z = mean(h2_Z, na.rm = T),
            seH2Z = se(h2_Z, na.rm = T),
            meanVAZ = mean(VA_Z, na.rm = T),
            seVAZ = se(VA_Z, na.rm = T))

# Number of loci doesn't seem to affect it too much, average across
ggplot(d_h2 %>%
         filter(method == "mkr", r %in% r_subsample) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = optPerc, y = h2_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_sum %>% ungroup() %>%
               filter(method == "mkr", r %in% r_subsample) %>% 
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
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
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
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

ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 0.7)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_med

ggplot(d_h2 %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 1.25),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 1.25),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  coord_cartesian(ylim = c(0, 2.5)) +
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
          plt_add_va_med + theme(legend.position = "none"),
          plt_add_va_lrg + theme(legend.position = "none"),
          ncol = 1, labels = "AUTO")

plt_add_va <- plot_grid(plt_add_va,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_add_va
ggsave("plt_va.png", device = png, bg = "white",
       width = 560*4, height = 980*4, units = "px")


# Paixao and Barton 2016 (regarding a polygenic trait with each QTL under 
# negligible selection): "Drift will disperse allele frequencies, decreasing the
# additive variance by a factor (1 - 1/2Ne) per generation"

# Hossjer et al 2016: Ne estimate from VA:
# 1 - (1 - 1/2Ne)^tau = (VA_t - VA_{t+tau})/VA_t

# Infinitesimal expects zero change in additive variance due to selection
# So see how much it changes between timepoints
# Scale by the total variance as well -> a large effect model will produce a lot
# of variance, so the differences are more likely to be greater
# Also account for drift: estimate Ne via Hossjer et al.

d_h2 %>%
  group_by(model, seed, tau, r, nloci, method) %>%
  filter(n() > 1) %>%
  summarise(totalDeltaVA = sum(diff(VA_Z))/sum(VA_Z)) -> d_h2_deltaVA


# total distribution
boxplot(d_h2_deltaVA$totalDeltaVA)

d_h2_deltaVA %>%
  group_by(model, tau, r, method) %>%
  summarise(meanDeltaVA = mean(totalDeltaVA, na.rm = T),
            seDeltaVA = se(totalDeltaVA, na.rm = T)) -> d_h2_deltaVA_sum


# nloci doesn't matter again
ggplot(d_h2_deltaVA %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = model, y = totalDeltaVA, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.9) +
  #coord_cartesian(ylim = c(0, 1)) +
  geom_point(data = d_h2_deltaVA_sum %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = model, y = meanDeltaVA, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Model", 
       y = TeX("Change in additive variance $(\\Delta V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 
                                           3, direction = -1),
                      guide = "none",
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

ggsave("plt_deltaVA.png", device = png, width = 9, height = 4)

# Correlation of genetic variance to time to adaptation
d_pheno_va <- left_join(d_h2, 
                        d_qg %>% select(gen, seed, modelindex, optPerc,
                                            deltaPheno, deltaw) %>% 
                          filter(gen >= 49500),
                        by = c("gen", "seed", "modelindex", "optPerc"))

d_pheno_va <- d_pheno_va %>%
  filter(as.numeric(optPerc) == 4, VA_Z < 50) %>%
  group_by(seed, model, tau, r, nloci, method) %>%
  mutate(timeToAdaptation = gen - 50000)

d_pheno_va_cor <- d_pheno_va %>%
  group_by(model, tau, r, method) %>%
  summarise(corVAw = cor(timeToAdaptation, VA_Z))

# linear model: time to adaptation vs VA
lm_VA <- lm(timeToAdaptation ~ VA_Z * model * tau + r, data = d_pheno_va)
summary(lm_VA)
plot(timeToAdaptation ~ VA_Z * model * tau + r, data = d_pheno_va)

# plot
ggplot(d_pheno_va_cor %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample),
       aes(x = as.factor(tau), y = corVAw, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_point(position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #coord_cartesian(ylim = c(0, 1)) +
  # geom_point(data = d_h2_deltaVA_sum %>%
  #              mutate(r_title = "Recombination rate (log10)",
  #                     nloci_title = "Number of loci",
  #                     tau_title = "Mutational effect size variance") %>%
  #              filter(method == "mkr", r %in% r_subsample),
  #            aes(x = as.factor(tau), y = meanDeltaVA, group = model), colour = "black",
  #            shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Mutational effect size variance", 
       y = "Correlation between additive\nvariance and fitness",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")
