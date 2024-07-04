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
d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg

d_qg$optPerc <- d_qg$phenomean - 1
d_qg$optPerc <- cut(d_qg$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

d_qg_optPerc <- d_qg %>% select(gen, seed, modelindex, optPerc, isAdapted) %>% filter(gen >= 49500)

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
table(d_h2$model, d_h2$isAdapted)

# We have many recombination rates: choose a few
r_subsample <- c(1e-10, 1e-5, 1e-1)

# Distribution, how different are the estimates
ggplot(d_h2 %>% filter(isAdapted) %>%
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
boxplot(d_h2[!outliers,]$VA_Z)

# filter out outliers
d_h2 <- d_h2[!outliers,]
write.csv(d_h2, "d_h2_outliersremoved.csv")
d_h2 <- read.csv("d_h2_outliersremoved.csv")
d_h2 <- d_h2[,-1]

# summarise
d_h2_sum <- d_h2 %>%
  filter(r %in% r_subsample) %>%
  group_by(optPerc, model, tau, r, method, isAdapted) %>%
  summarise(meanH2Z = mean(h2_Z, na.rm = T),
            seH2Z = se(h2_Z, na.rm = T),
            meanVAZ = mean(VA_Z, na.rm = T),
            seVAZ = se(VA_Z, na.rm = T))

# Number of loci doesn't seem to affect it too much, average across
ggplot(d_h2 %>% filter(isAdapted) %>%
         filter(method == "mkr", r %in% r_subsample) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = optPerc, y = h2_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_sum %>% ungroup() %>%
               filter(method == "mkr", r %in% r_subsample, isAdapted) %>% 
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
ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
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

ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
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

ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 1.25),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
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

# Only small effects
library(lemon)
# Additive variance
# Again, nloci not important
ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         mutate(model = fct_recode(model, "Additive" = "Add", 
                                   "K+" = "K",
                                   "K-" = "ODE")) %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested_wrap(model ~ r_title + log10(r), scales = "free") +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               mutate(model = fct_recode(model, "Additive" = "Add", 
                                         "K+" = "K",
                                         "K-" = "ODE")) %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = optPerc, y = meanVAZ, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.8)) +
  #coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14), panel.spacing.y = unit(1, "line"),
        legend.position = "bottom") -> plt_add_va_sml_scl
plt_add_va_sml_scl
ggsave("VA_sml_scl.png", plt_add_va_sml_scl, width = 8, height = 8, device = png)



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
  group_by(model, seed, tau, r, nloci, isAdapted, method) %>%
  filter(n() > 1) %>%
  summarise(totalDeltaVA = sum(diff(VA_Z))/sum(VA_Z)) -> d_h2_deltaVA


# total distribution
boxplot(d_h2_deltaVA$totalDeltaVA)

d_h2_deltaVA %>%
  group_by(model, tau, r, method, isAdapted) %>%
  summarise(meanDeltaVA = mean(totalDeltaVA, na.rm = T),
            seDeltaVA = se(totalDeltaVA, na.rm = T)) -> d_h2_deltaVA_sum


# nloci doesn't matter again
ggplot(d_h2_deltaVA %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = model, y = totalDeltaVA, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.9) +
  #coord_cartesian(ylim = c(0, 1)) +
  geom_point(data = d_h2_deltaVA_sum %>% filter(isAdapted) %>%
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
d_h2$seed <- as.factor(d_h2$seed)
d_h2$modelindex <- as.factor(d_h2$modelindex)
d_h2$optPerc <- as.factor(d_h2$optPerc)

d_pheno_va <- left_join(d_h2, 
                        d_qg %>% select(gen, seed, modelindex, optPerc, isAdapted,
                                            deltaPheno, deltaw) %>% 
                          filter(gen >= 49500),
                        by = c("gen", "seed", "modelindex", "optPerc", "isAdapted"))

d_pheno_va <- d_pheno_va %>%
  filter(as.numeric(optPerc) == 4, VA_Z < 50, isAdapted) %>%
  group_by(seed, model, tau, r, nloci, method) %>%
  mutate(timeToAdaptation = gen - 50000)

d_pheno_va_cor <- d_pheno_va %>%
  group_by(model, tau, r, method) %>%
  summarise(corVAw = cor(1/timeToAdaptation, VA_Z))

# linear model: time to adaptation vs VA
lm_VA <- lm(timeToAdaptation ~ log(VA_Z) * model * tau + r, data = d_pheno_va)
summary(lm_VA)
plot(timeToAdaptation ~ log(VA_Z), data = d_pheno_va)# * model * tau + r, data = d_pheno_va)
emmeans::emmeans(lm_VA, specs = pairwise ~ VA_Z|model|tau)

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
       y = "Correlation between additive\nvariance and time to adaptation",
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1),
                      labels = c("Additive", "K+", "K-")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# Molecular G: extract to array of matrices
# small cluster of non-adapted pops in cluster 1, otherwise isAdapted
# doesn't predict cluster, removing

d_h2 %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  mutate(tmpCVA = CVA_a_b,
         CVA_a_b = CVA_a_KZ,
         CVA_a_KZ = if_else(model == "ODE", NA, tmpCVA)) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2


# Separate into model indices
# each sublist is replicates of a model index
Rcpp::sourceCpp("./getCovarianceMatrices.cpp")
lapply(split_h2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices
lapply(split_h2, function(x) {list(optPerc = x$optPerc[1], seed = x$seed[1], modelindex = x$modelindex[1], isAdapted = x$isAdapted[1])}) -> cov_matrix_modelindex


# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

# So eigentensor analysis: sample random combinations of seeds to get a distribution
# of eigenvectors telling us the models which have the largest difference in variation
# then projection to find the important components

# First sample a matrix from each group
test_indices <- sapply(cov_matrices, function(x) {
  sample(1:length(x), 1)
  }, simplify = F)

test <- sapply(1:length(cov_matrices), function(x) {
  cov_matrices[[x]][[test_indices[[x]]]]
}, simplify = F)


# Repeat with all matrices
h2_mat <- unlist(cov_matrices, recursive = F)

# get ids
lapply(split_h2, function(x) {
  data.frame(optPerc = x$optPerc, 
       seed = x$seed, 
       modelindex = x$modelindex, 
       isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_full

# Split data frames of replicates to individual lists of dataframes with 1 row
lapply(cov_matrix_modelindex_full, function(x) {
  split(x, seq(nrow(x)))
}) -> cov_matrix_modelindex_full
# unlist to full form
cov_matrix_modelindex_full <- unlist(cov_matrix_modelindex_full, recursive = F)

# PCAS <- PCAsimilarity(test[1:100])

# Distance between G matrices
library(ape)
frobenius_distance <- function(A, B) {
  return(sqrt(sum((A - B)^2)))
}

compute_distance_matrix <- function(list_of_matrices) {
  n <- length(list_of_matrices)
  distance_matrix <- matrix(0, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      distance <- shapes::distcov(list_of_matrices[[i]], list_of_matrices[[j]], "Power")
      distance_matrix[i, j] <- distance
      distance_matrix[j, i] <- distance
    }
  }
  
  rownames(distance_matrix) <- colnames(distance_matrix) <- paste("Matrix", 1:n)
  return(distance_matrix)
}

# microbenchmark::microbenchmark(times = 30, 
#                                compute_distance_matrix(h2_mat[1:100]),
#                                distanceMatrix(h2_mat[1:100]))
# 
# # ensure results are the same - distance should be 0 (with some precision error)
# shapes::distcov(compute_distance_matrix(h2_mat[1:100]), 
#                 distanceMatrix(h2_mat[1:100]), "Power")

library(tidytree)
library(ggtree)
library(phytools)

Rcpp::sourceCpp("./distanceFunctions.cpp")
dist_matrix <- distanceMatrix(h2_mat)
colnames(dist_matrix) <- paste("Matrix", 1:nrow(dist_matrix))
rownames(dist_matrix) <- colnames(dist_matrix)

#fviz_dist(as.dist(dist_matrix), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

hc <- hclust(as.dist(dist_matrix), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 2 seems to be the best
library(factoextra)
# elbow plot
fviz_nbclust(dist_matrix, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 2, border = 2)

# gap stat
#gap_stat <- cluster::clusGap(dist_matrix, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
#fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

# silhouette
#fviz_nbclust(dist_matrix, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot")

clus <- cutree(hc, 2)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id <- rbindlist(cov_matrix_modelindex_full, fill = T)
id$label <- as.character(1:nrow(id))
id$modelindex <- as.factor(id$modelindex)
id <- AddCombosToDF(id)
id$nloci_group <- "[4, 64)"
id$nloci_group[id$nloci >= 64 & id$nloci < 1024] <- "[64, 256]"
id$nloci_group[id$nloci == 1024] <- "[1024]"
id$nloci_group <- factor(id$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id[idx,"clus"] <- i
}




# with id, check how frequent genetic architectures are with the clusters
tab <- table(id$clus, id$nloci_group, id$r, id$model, id$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

model <- glm(Freq~cluster*nloci,family=poisson(),data=tab)
summary(model)

id %>% ungroup() %>%
  group_by(r, model, clus) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  mutate(prop = n/sum(n)) -> cluster_percs_r

id %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  mutate(prop = n/sum(n)) -> cluster_percs_nloci

phylo <- full_join(as.phylo(phylo), id, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_full.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id$clus)) {
  if(length(id$clus[id$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id$label[id$clus == i]), 
                                          fill = clus_palette[i], alpha = 0.2,
                                          type = "encircle", to.bottom = T) 
  
  # Find the position for the annotation
  cluster_tips <- tree_r$data %>% filter(clus == i)
  annotation_x <- min(cluster_tips$x) - 0.15
  annotation_y <- max(cluster_tips$y) + 0.05
  
  # Add text annotation
  for (j in 1:length(cluster_labels)) {
    tree_r <- tree_r + annotate("text", x = annotation_x, y = annotation_y, label = cluster_labels[j],
                                        color = clus_palette[i], hjust = 0, vjust = 1 + j*1.5, size = 3)
  }
}
tree_r


ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(nloci_group)), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Number of loci") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id$clus)) {
  if(length(id$clus[id$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
       x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id$label[id$clus == i]), 
                                          fill = clus_palette[i], alpha = 0.2,
                                          type = "encircle", to.bottom = T) 
  
  # Find the position for the annotation
  cluster_tips <- tree_nloci$data %>% filter(clus == i)
  annotation_x <- min(cluster_tips$x) - 0.15
  annotation_y <- max(cluster_tips$y) + 0.05
  
  # Add text annotation
  for (j in 1:length(cluster_labels)) {
    tree_nloci <- tree_nloci + annotate("text", x = annotation_x, y = annotation_y, label = cluster_labels[j],
                       color = clus_palette[i], hjust = 0, vjust = 1 + j*1.5, size = 3)
  }
}
tree_nloci

plt_trees <- plot_grid(tree_r,
                        tree_nloci,
                        ncol = 2, labels = "AUTO")

plt_trees
ggsave("plt_tree_gmatrix_full_2clus.png", device = png, bg = "white",
       width = 12, height = 6)

# Full one very similar to subset, run analysis on sample


angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

# Want to find the similarities in structure between models in these clusters
# PCA of each covariance matrix to find proportion of variance explained in PC1 PC2,
# contributions of each to PC1 PC2, total variance in each matrix,
# major/minor axis
CovMatrixPCA <- function(matList, id) {
  traitMap = list(
    "KXZ_KZ" = c(1,2),
    "KXZ_Z" = c(1,3),
    "KXZ_a" = c(1,4),
    "KXZ_b" = c(1,5),
    "KZ_Z" = c(2,3),
    "KZ_a" = c(2,4),
    "KZ_b" = c(2,5),
    "Z_a" = c(3,4),
    "Z_b" = c(3,5),
    "a_b" = c(4,5)
  )
  
  singleTraitsMap <- c(
    "KXZ",
    "KZ",
    "Z",
    "a",
    "b"
  )
  
  PCAdata <- data.frame(
    trait1 = numeric(length(matList) * length(traitMap)),
    trait2 = numeric(length(matList) * length(traitMap)),
    totalVariation = numeric(length(matList) * length(traitMap)),
    pc1_prop = numeric(length(matList) * length(traitMap)),
    pc2_prop = numeric(length(matList) * length(traitMap)),
    pc1_contrib_t1 = numeric(length(matList) * length(traitMap)),
    pc1_contrib_t2 = numeric(length(matList) * length(traitMap)),
    pc2_contrib_t1 = numeric(length(matList) * length(traitMap)),
    pc2_contrib_t2 = numeric(length(matList) * length(traitMap)),
    pc_majorlength_t1 = numeric(length(matList) * length(traitMap)),
    pc_majorlength_t2 = numeric(length(matList) * length(traitMap)),
    pc_minorlength_t1 = numeric(length(matList) * length(traitMap)),
    pc_minorlength_t2 = numeric(length(matList) * length(traitMap)),
    pc_majorangle = numeric(length(matList) * length(traitMap))
  )
  
  # Set trait columns: keep as number for indexing
  t1 <- unlist(lapply(traitMap, function(x) {x[1]}))
  t2 <- unlist(lapply(traitMap, function(x) {x[2]}))
  PCAdata$trait1 <- rep(t1, times = length(matList))
  PCAdata$trait2 <- rep(t2, times = length(matList))
  
  # Iterate through matrices for PCA
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    pca <- eigen(g)
    
    # Rows to fill for all traits
    i_range <- ( (i-1)*length(traitMap) + 1 ):(i * length(traitMap))
    
    totalVar <- sum(pca$values)
    PCAdata[i_range,]$totalVariation <- totalVar
    PCAdata[i_range,4:5] <- matrix(rep((pca$values/totalVar)[1:2], 
                                       each = length(i_range)), ncol = 2)
    
    pc_sqr <- pca$vectors^2
    
    pc_contrib <- sweep(pc_sqr, 2, colSums(pc_sqr), FUN="/") * 100
    
    # Contributions
    ## PC1
    PCAdata[i_range,6] <- c(pc_contrib[t1,1])
    PCAdata[i_range,7] <- c(pc_contrib[t2,1])
    
    ## PC2
    PCAdata[i_range,8] <- c(pc_contrib[t1,2])
    PCAdata[i_range,9] <- c(pc_contrib[t2,2])
    
    # Major length
    PCAdata[i_range,10] <- pca$vector[t1,1] * sqrt(pca$values[1])
    PCAdata[i_range,11] <- pca$vector[t2,1] * sqrt(pca$values[1])
    
    # Minor length
    PCAdata[i_range,12] <- pca$vector[t1,2] * sqrt(pca$values[2])
    PCAdata[i_range,13] <- pca$vector[t2,2] * sqrt(pca$values[2])
    
    # Angle of rotation (relative to trait 1, second axis orthogonal)
    PCAdata[i_range,14] <- diag(atan2((pca$values[1] - diag(g)[t1]), g[t1,t2]))
  }
  
  PCAdata <- PCAdata %>%
    mutate(optPerc = rep(id$optPerc, each = length(traitMap)),
           seed = rep(id$seed, each = length(traitMap)),
           modelindex = rep(id$modelindex, each = length(traitMap)),
           clus = rep(id$clus, each = length(traitMap)))
  
  return(PCAdata)
}

test_pca <- eigen(h2_mat[[1]])
test_angle <- angle(test_pca$vectors[,1], test_pca$vectors[,2])

# 95% confidence ellipse axis length
test_major_len <- ( qnorm(0.975) * test_pca$values[1] )
test_minor_len <- ( qnorm(0.975) * test_pca$values[2] )

test_dplot_ellipse <- data.frame(vert_x = cos(test_angle * pi/180) * test_major_len,
                            vert_y = sin(test_angle * pi/180) * test_major_len,
                            covert_x = cos((test_angle - 90) * pi/180) * test_minor_len,
                            covert_y = sin((test_angle - 90) * pi/180) * test_minor_len,
                            mean_t1 = 0,
                            mean_t2 = 0,
                            theta = test_angle,
                            major_len = test_major_len,
                            minor_len = test_minor_len)


ggplot(test_dplot_ellipse, aes(x = mean_t1, y = mean_t2)) +
  geom_point(data = test_scores,
             aes(x = PC1, y = PC2)) + 
  geom_ellipse(aes(x0 = mean_t1, y0 = mean_t2, 
                   a = major_len, b = minor_len, angle = theta * pi/180)) +
  geom_segment(aes(xend = (mean_t1 + vert_x), yend = (mean_t2 + vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - vert_x), yend = (mean_t2 - vert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 + covert_x), yend = (mean_t2 + covert_y)),
               linetype = "dashed") +
  geom_segment(aes(xend = (mean_t1 - covert_x), yend = (mean_t2 - covert_y)),
               linetype = "dashed") +
  theme_bw() +
  labs(x = "Trait 1", y = "Trait 2") +
  coord_fixed() # important! ensures gmax and g2 appear orthogonal. different aspect ratios distort the angles between gmax and g2 


covpca_test <- CovMatrixPCA(h2_mat, id)
covpca_test <- AddCombosToDF(covpca_test)

# Summarise
covpca_sum <- covpca_test %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))

# Calculate ellipse x and y coordinates after rotation
ellipse_points <- function(x_center, y_center, major_axis, minor_axis, angle, n_points = 100) {
  theta <- seq(0, 2 * pi, length.out = n_points)
  x <- major_axis * cos(theta)
  y <- minor_axis * sin(theta)
  rotated_x <- x * cos(angle) - y * sin(angle)
  rotated_y <- x * sin(angle) + y * cos(angle)
  data.frame(x = rotated_x + x_center, y = rotated_y + y_center)
}

# normal direction at point t on an ellipse with major axis length a pointing
# in direction u1, and minor axis length b in direction u2
ellipse_normal <- function(t, a, b, u1, u2) {
  return(a * sin(t) * u2 + b * cos(t) * u1)
}

# Length of the normal of an ellipse at point t with major/minor lengths a and b
ellipse_normal_len <- function(t, a, b) {
  return(sqrt(a^2 * (sin(t)^2) + b^2 * (cos(t)^2)))
}

# Get unit vector directions of major/minor axes
getAxisDirections <- function(angle) {
  u1 <- c(cos(angle), sin(angle))
  u2 <- c(-sin(angle), cos(angle))
  return(list(u1 = u1, u2 = u2))
}

# Solve O(t) for the outer boundary
outer_boundary <- function(t, a, b, u1, u2, h) {
  P_t <- a * cos(t) * u1 + b * sin(t) * u2
  normal <- ellipse_normal(t, a, b, u1, u2)
  normal_len <- ellipse_normal_len(t, a, b)
  O_t <- P_t + (h / normal_len) * normal
  return(O_t)
}

# Solve I(t) for the inner boundary
inner_boundary <- function(t, a, b, u1, u2, h) {
  P_t <- a * cos(t) * u1 + b * sin(t) * u2
  normal <- ellipse_normal(t, a, b, u1, u2)
  normal_len <- ellipse_normal_len(t, a, b)
  I_t <- P_t - (h / normal_len) * normal
  return(I_t)
}

ribbon_points <- function(major_axis_mean, minor_axis_mean, major_axis_se, minor_axis_se, angle_mean, n_points = 100) {
  theta <- seq(0, 2 * pi, length.out = n_points)
  
  u <- getAxisDirections(angle_mean)
  u1 <- u$u1
  u2 <- u$u2
  
  outer_points <- t(sapply(theta, outer_boundary, a = major_axis_mean, 
                         b = minor_axis_mean, u1 = u1, u2 = u2,
                         h = major_axis_se))
  inner_points <- t(sapply(theta, inner_boundary, a = major_axis_mean, 
                         b = minor_axis_mean, u1 = u1, u2 = u2,
                         h = major_axis_se))
  
  upper_df <- data.frame(x = outer_points[,1], y = outer_points[,2])
  lower_df <- data.frame(x = inner_points[,1], y = inner_points[,2])
  
  ribbon_df <- rbind(upper_df, lower_df[nrow(lower_df):1, ])
  return(ribbon_df)
}

# Major and minor lines
axis_lines <- function(major_axis, minor_axis, angle) {
  
  vert_x = cos(angle) * major_axis
  vert_y = sin(angle) * major_axis
  covert_x = cos((angle - (90 * pi/180))) * minor_axis
  covert_y = sin((angle - (90 * pi/180))) * minor_axis

  return(data.frame(maj_xend = vert_x, maj_yend = vert_y,
                              min_xend = covert_x, min_yend = covert_y))

}

# Run ellipse function on data
## Mean ellipse
covpca_dplot <- covpca_sum %>% ungroup() %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  mutate(
  mean_ellipse = pmap(list(0,0,pc_majorlength_t1_mean,pc_minorlength_t1_mean,pc_majorangle_mean), ellipse_points),
  ribbon_ellipse = pmap(list(pc_majorlength_t1_mean, pc_minorlength_t1_mean, pc_majorlength_t1_se/2, pc_minorlength_t1_se/2, pc_majorangle_mean), ribbon_points)) %>%
  pivot_longer(cols = ends_with("ellipse"), names_to = "ellipse_type", values_to = "data") %>%
  unnest(data)

# Average over nloci, doesn't matter; only 1 tau group being looked at
covpca_dplot_axes <- covpca_sum %>% ungroup() %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  mutate(mean_axes = pmap(list(pc_majorlength_t1_mean, 
                               pc_minorlength_t1_mean, pc_majorangle_mean), 
                          axis_lines)) %>%
  unnest(mean_axes)

# Iterate over all 10 combinations of traits
traitCombos <- list(
  c(1,2),
  c(1,3),
  c(1,4),
  c(1,5),
  c(2,3),
  c(2,4),
  c(2,5),
  c(3,4),
  c(3,5),
  c(4,5)
)

traitLabels <- c(TeX("$K_{XZ}$", output = "expression"), 
                 TeX("$K_Z$", output = "expression"),
                 TeX("$Z$", output = "expression"),
                 TeX("$\\alpha_Z$", output = "expression"),
                 TeX("$\\beta_Z$", output = "expression")
                 )
# Iterate over 2 time points (<25% optimum, >=75%)
traitCombos <- rep(traitCombos, each = 2)
res_plt <- vector("list", length(traitCombos) * 2)

for (i in seq_along(traitCombos)) {
  combo <- traitCombos[[i]]
  timepoint <- ifelse(i %% 2 != 0, 1, 4)
  # Subset data
sbst_plt <- covpca_dplot %>% filter(trait1 == combo[1], trait2 == combo[2], 
                        as.numeric(optPerc) == timepoint) %>%
  mutate(r_title = "Recombination rate (log10)",
         model_title = "Model (cluster)",
         model = fct_recode(model, "K+" = "K", "K-" = "ODE"))

axes_sbst_plt <- covpca_dplot_axes %>% filter(trait1 == combo[1], trait2 == combo[2], 
                                         as.numeric(optPerc) == timepoint) %>%
  mutate(r_title = "Recombination rate (log10)",
         model_title = "Model (cluster)",
         model = fct_recode(model, "K+" = "K", "K-" = "ODE"))


# A, B, C for different clusters, shapes of the clusters
# overlaid ellipses for optPerc

# plot
ggplot(sbst_plt,
       aes(x = x, y = y, colour = as.factor(clus))) +
  facet_nested(r_title + log10(r) ~ model_title + model + as.factor(clus)) +
  geom_polygon(data = filter(sbst_plt, ellipse_type == "ribbon_ellipse"),
              aes(fill = as.factor(clus)), linetype = "dashed", colour = NA, alpha = 0.3) +
  geom_polygon(data = filter(sbst_plt, ellipse_type == "mean_ellipse"),
               fill = NA) + 
  geom_segment(data = axes_sbst_plt, x = 0, y = 0, 
               aes(xend = (maj_xend), yend = (maj_yend)),
               linetype = "dashed") +
  geom_segment(data = axes_sbst_plt, x = 0, y = 0, 
               aes(xend = (-maj_xend), yend = (-maj_yend)),
               linetype = "dashed") +
  geom_segment(data = axes_sbst_plt, x = 0, y = 0, 
               aes(xend = (min_xend), yend = (min_yend)),
               linetype = "dashed") +
  geom_segment(data = axes_sbst_plt, x = 0, y = 0, 
               aes(xend = (-min_xend), yend = (-min_yend)),
               linetype = "dashed") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3,
                                           direction = -1),
                      labels = c("Cluster 1", "Cluster 2", "Cluster 3")) +
  scale_fill_manual(values = paletteer_d("nationalparkcolors::Everglades", 3,
                                           direction = -1),
                      guide = "none") +
  #lims(x = c(-0.06, 0.06), y = c(-0.06, 0.06)) +
  coord_equal() +
  labs(title = paste0("Progress to the optimum: ", ifelse(timepoint == 1, "<25%", ">=75%")),
       x = traitLabels[combo[1]], y = traitLabels[combo[2]], colour = "Cluster") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(1, "lines")) -> plt

res_plt[[i]] <- plt
filename <- paste0("plt_gmat_", traitMap[combo[1]], "_", traitMap[combo[2]], 
                   "_", timepoint, ".png")
#ggsave(filename, plt, device = png, width = 8, height = 8)
}


# Angles for each pair of traits, do a comparison plot, look at similarity of
# angle with clusters, nloci, r, tau etc.
# compare proportion contributing to PC1 and PC2: similar within clusters?
# what are the most constrained populations? Is K+ less constrained?


# Trait contributions
# text: percent the trait contributed to pc1 and pc2

traitContribs <- function(matList, id) {
  
  nTraits <- 5
  
  PCAdata <- data.frame(
    trait = numeric(length(matList) * nTraits * 2),
    pc = numeric(length(matList) * nTraits * 2),
    totalVar = numeric(length(matList) * nTraits * 2),
    pc_prop = numeric(length(matList) * nTraits * 2),
    contrib = numeric(length(matList) * nTraits * 2),
    dir = numeric(length(matList) * nTraits * 2)
  )
  
  PCAdata$trait <- rep(1:nTraits, times = 2 * length(matList))
  PCAdata$pc <- rep(1:2, each = nTraits)
  
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    pca <- eigen(g)
    
    # Rows to fill for all traits
    i_range <- ( (i-1)*nTraits * 2 + 1 ):(i * nTraits * 2)
    
    totalVar <- sum(pca$values)
    PCAdata[i_range,]$totalVar <- totalVar
    PCAdata[i_range,4] <- rep((pca$values/totalVar)[1:2], 
                              each = nTraits)
    
    pc_sqr <- pca$vectors^2
    
    pc_contrib <- sweep(pc_sqr, 2, colSums(pc_sqr), FUN="/") * 100
    
    # Contributions
    PCAdata[i_range,5] <- c(pc_contrib[,1:2])
    
    # Eigenvector direction
    PCAdata[i_range,6] <- c(pca$vectors[,1:2])
  }
  
  PCAdata <- PCAdata %>%
    mutate(optPerc = rep(id$optPerc, each = nTraits*2),
           seed = rep(id$seed, each = nTraits*2),
           modelindex = rep(id$modelindex, each = nTraits*2),
           clus = rep(id$clus, each = nTraits*2))
  
  return(PCAdata)
}

# Get trait contributions
covpca_traitprops <- traitContribs(h2_mat, id)
covpca_traitprops <- AddCombosToDF(covpca_traitprops)


# Summarise
covpca_traitprops %>%
  select(-nloci, -seed, -tau) %>%
  group_by(optPerc, model, clus, trait, pc, r) %>%
  summarise_if(is.numeric, list(mean = mean, CI = CI, sd = sd)) -> covpca_traitprops_sum

covpca_traitprops_sum <- covpca_traitprops_sum %>%
  mutate(trait = traitMap[trait])

ggplot(covpca_traitprops_sum %>% filter(as.numeric(optPerc) == 4) %>%
         mutate(r_title = "Recombination rate (log10)",
                model_title = "Model (cluster)",
                model = fct_recode(model, "K+" = "K", "K-" = "ODE")),
       aes(x = trait, y = contrib_mean, colour = model, shape = as.factor(pc))) +
  facet_nested(r_title + log10(r) ~ model_title + model + as.factor(clus)) +
  geom_point(size = 1.5, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = contrib_mean - contrib_CI, ymax = contrib_mean + contrib_CI),
                width = 0.3, position = position_dodge(0.3)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3,
                                           direction = -1)[2:3],
                      labels = c("K+", "K-"), guide = "none") +
  scale_shape_manual(values = c(0, 15)) +
  labs(x = "Molecular Trait", y = "Contribution to PC", 
       shape = "Principal Component") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))

# Contributions are more or less the same across clusters, models, r
# Z is the largest contributor to PC1
# PC2 is mainly a and b, in K+ KZ is in there as well
# KXZ does not contribute much to additive variation

# Eigenvectors
ggplot(covpca_traitprops_sum %>% filter(as.numeric(optPerc) == 4) %>%
         mutate(r_title = "Recombination rate (log10)",
                model_title = "Model (cluster)",
                model = fct_recode(model, "K+" = "K", "K-" = "ODE")),
       aes(x = trait, y = dir_mean, colour = model, shape = as.factor(pc))) +
  facet_nested(r_title + log10(r) ~ model_title + model + as.factor(clus)) +
  geom_point(size = 1.5, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = dir_mean - dir_CI, ymax = dir_mean + dir_CI),
                width = 0.3, position = position_dodge(0.3)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3,
                                           direction = -1)[2:3],
                      labels = c("K+", "K-"), guide = "none") +
  scale_shape_manual(values = c(0, 15)) +
  labs(x = "Molecular Trait", y = "Eigenvector", 
       shape = "Principal Component") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 14))
