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

# Make sure we use the right dplyr functions
summarise <- dplyr::summarise
mutate <- dplyr::mutate

d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/"
source("/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/helperFunctionsAndSetup.R")

d_h2_mkr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                            "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mrr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr.csv"), header = F,
                              col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                            "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                            "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                            "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
                                            "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                            "h2_KXZ"))

d_h2_mrr_pt1 <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr_pt1.csv"), header = F,
                                  col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                                "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                                "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                                "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
                                                "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                                "h2_KXZ"))

d_h2_mkr_pt1 <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr_pt1.csv"), header = F,
                                  col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                                "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                                "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                                "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
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
  dplyr::mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
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
  dplyr::mutate(modelindex = as.factor(modelindex),
         seed = as.factor(seed)) %>%
  drop_na(VA_Z) %>% distinct()

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex"))

d_h2 <- AddCombosToDF(d_h2)

# Counts for each model type:
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
  dplyr::summarise(meanH2Z = mean(h2_Z, na.rm = T),
            seH2Z = se(h2_Z, na.rm = T),
            meanVAZ = mean(VA_Z, na.rm = T),
            seVAZ = se(VA_Z, na.rm = T))
d_h2_sum$model <- as.factor(d_h2_sum$model)

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
# Small effects as separate figure
ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
       aes(x = optPerc, y = VA_Z, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
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
        legend.position = "bottom")
ggsave("plt_va_sml.png", device = png, bg = "white",
       width = 560*4, height = (980*4)/3, units = "px")


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

# did the mean VA change with changing recombination?
library(nlme)
summary(gls.dva <- gls(totalDeltaVA ~ model * as.factor(r), 
                       d_h2_deltaVA %>% filter(isAdapted, r %in% r_subsample), 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.dva)

# Average across recombination rates
d_h2_deltaVA %>% filter(isAdapted) %>%
  group_by(model, tau, method, isAdapted) %>%
  summarise(meanDeltaVA = mean(totalDeltaVA, na.rm = T),
            seDeltaVA = se(totalDeltaVA, na.rm = T)) -> d_h2_deltaVA_sum


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

# Scaled h2 estimates
d_h2_scaled <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr_scaled.csv"), header = F,
                                 col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                               "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                               "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                               "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
                                               "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                               "h2_KXZ"))

# d_h2_mrr_scaled <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr_scaled.csv"), header = F,
#                                      col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
#                                                    "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
#                                                    "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
#                                                    "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
#                                                    "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
#                                                    "h2_KXZ"))

d_h2_scaled$method <- "mkr"

# Clean data
d_h2_scaled <- d_h2_scaled %>%
  distinct(gen, seed, modelindex, method, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  drop_na(VA_Z) %>% distinct()

# inner join optPerc
d_h2_scaled <- left_join(d_h2_scaled, d_qg_optPerc, by = c("gen", "seed", "modelindex"))
d_h2_scaled <- AddCombosToDF(d_h2_scaled)

d_h2_scaled %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 1) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc1

d_h2_scaled %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 2) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc2

d_h2_scaled %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 3) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc3

d_h2_scaled %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 4) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc4



# Separate into model indices
# each sublist is replicates of a model index
Rcpp::sourceCpp("./getCovarianceMatrices.cpp")
lapply(split_h2_optPerc1, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op1
lapply(split_h2_optPerc2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op2
lapply(split_h2_optPerc3, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op3
lapply(split_h2_optPerc4, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op4

lapply(split_h2_optPerc1, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op1
lapply(split_h2_optPerc2, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op2
lapply(split_h2_optPerc3, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op3
lapply(split_h2_optPerc4, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op4


# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

# So eigentensor analysis: sample random combinations of seeds to get a distribution
# of eigenvectors telling us the models which have the largest difference in variation
# then projection to find the important components

# Repeat with all matrices
h2_mat_op1 <- unlist(cov_matrices_op1, recursive = F)
h2_mat_op2 <- unlist(cov_matrices_op2, recursive = F)
h2_mat_op3 <- unlist(cov_matrices_op3, recursive = F)
h2_mat_op4 <- unlist(cov_matrices_op4, recursive = F)

# get ids
GetMatrixIDs <- function(matList) {
  lapply(matList, function(x) {
    data.frame(optPerc = x$optPerc, 
               seed = x$seed, 
               modelindex = x$modelindex, 
               isAdapted = x$isAdapted)}) -> matList
  
  
  lapply(matList, function(x) {
    split(x, seq(nrow(x)))
  }) -> matList
  # unlist to full form
  matList <- unlist(matList, recursive = F)
  return(matList)
}

#cov_matrix_modelindex_full <- GetMatrixIDs(split_h2)

cov_matrix_modelindex_op1 <- GetMatrixIDs(split_h2_optPerc1)
cov_matrix_modelindex_op2 <- GetMatrixIDs(split_h2_optPerc2)
cov_matrix_modelindex_op3 <- GetMatrixIDs(split_h2_optPerc3)
cov_matrix_modelindex_op4 <- GetMatrixIDs(split_h2_optPerc4)


# Distance between G matrices
library(ape)
library(tidytree)
library(ggtree)
library(phytools)

# dist matrix for optperc 1
Rcpp::sourceCpp("./distanceFunctions.cpp")
dist_matrix_op1 <- distanceMatrix(h2_mat_op1)
colnames(dist_matrix_op1) <- paste("Matrix", 1:nrow(dist_matrix_op1))
rownames(dist_matrix_op1) <- colnames(dist_matrix_op1)

#fviz_dist(as.dist(dist_matrix), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

hc <- hclust(as.dist(dist_matrix_op1), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 3 seems to be the best
library(factoextra)
# elbow plot
fviz_nbclust(dist_matrix_op1, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 3, border = 2:3)

clus <- cutree(hc, 3)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id_op1 <- rbindlist(cov_matrix_modelindex_op1, fill = T)
id_op1$label <- as.character(1:nrow(id_op1))
id_op1$modelindex <- as.factor(id_op1$modelindex)
id_op1 <- AddCombosToDF(id_op1)
id_op1$nloci_group <- "[4, 64)"
id_op1$nloci_group[id_op1$nloci >= 64 & id_op1$nloci < 1024] <- "[64, 256]"
id_op1$nloci_group[id_op1$nloci == 1024] <- "[1024]"
id_op1$nloci_group <- factor(id_op1$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id_op1$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id_op1[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id_op1$clus, id_op1$nloci_group, id_op1$r, id_op1$model, id_op1$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

model <- glm(Freq~cluster*nloci,family=poisson(),data=tab)
summary(model)

id_op1 %>% ungroup() %>%
  group_by(r, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_r

id_op1 %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_nloci

id_op1 %>% ungroup() %>%
  group_by(optPerc, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_optPerc


phylo <- full_join(as.phylo(phylo), id_op1, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op1_scaled.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id_op1$clus)) {
  if(length(id_op1$clus[id_op1$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id_op1$label[id_op1$clus == i]), 
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
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id_op1$clus)) {
  if(length(id_op1$clus[id_op1$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
       x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id_op1$label[id_op1$clus == i]), 
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
ggsave("plt_tree_gmatrix_optperc1_scaled.png", device = png, bg = "white",
       width = 12, height = 6)

# Progress to optimum <25%: very distinct clustering by model type
# early in the walk variance-covariance is different between the models
# outlier K+ clusters have lower recombination maybe




# Repeat for opt perc 4
dist_matrix_op4 <- distanceMatrix(h2_mat_op4)
colnames(dist_matrix_op4) <- paste("Matrix", 1:nrow(dist_matrix_op4))
rownames(dist_matrix_op4) <- colnames(dist_matrix_op4)

#fviz_dist(as.dist(dist_matrix), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

hc <- hclust(as.dist(dist_matrix_op4), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 3 seems to be the best
# elbow plot
fviz_nbclust(dist_matrix_op4, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 3, border = 2:3)

clus <- cutree(hc, 3)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id_op4 <- rbindlist(cov_matrix_modelindex_op4, fill = T)
id_op4$label <- as.character(1:nrow(id_op4))
id_op4$modelindex <- as.factor(id_op4$modelindex)
id_op4 <- AddCombosToDF(id_op4)
id_op4$nloci_group <- "[4, 64)"
id_op4$nloci_group[id_op4$nloci >= 64 & id_op4$nloci < 1024] <- "[64, 256]"
id_op4$nloci_group[id_op4$nloci == 1024] <- "[1024]"
id_op4$nloci_group <- factor(id_op4$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id_op4$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id_op4[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id_op4$clus, id_op4$nloci_group, id_op4$r, id_op4$model, id_op4$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

model <- glm(Freq~cluster+r,family=poisson(),data=tab)
summary(model)

id_op4 %>% ungroup() %>%
  group_by(r, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_r

id_op4 %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_nloci

phylo <- full_join(as.phylo(phylo), id_op4, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op4_scaled.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id_op4$clus)) {
  if(length(id_op4$clus[id_op4$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id_op4$label[id_op4$clus == i]), 
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
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id_op4$clus)) {
  if(length(id_op4$clus[id_op4$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id_op4$label[id_op4$clus == i]), 
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
ggsave("plt_tree_gmatrix_optperc4_scaled.png", device = png, bg = "white",
       width = 12, height = 6)

# By the end of the walk the G matrices are more similar: there is still clustering
# between models, but the distance between them is quite small
# the outliers are still in the K+ model, mainly this time driven by high
# recombination rates


angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}
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

covpca_op1 <- CovMatrixPCA(h2_mat_op1, id_op1)
covpca_op1 <- AddCombosToDF(covpca_op1)

covpca_op4 <- CovMatrixPCA(h2_mat_op4, id_op4)
covpca_op4 <- AddCombosToDF(covpca_op4)


# Summarise
covpca_op1_sum <- covpca_op1 %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))

# Run ellipse function on data
## Mean ellipse
covpca_dplot_op1 <- covpca_op1_sum %>% ungroup() %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  mutate(
  mean_ellipse = pmap(list(0,0,pc_majorlength_t1_mean,pc_minorlength_t1_mean,
                           pc_majorangle_mean), ellipse_points),
  ribbon_ellipse = pmap(list(pc_majorlength_t1_mean, pc_minorlength_t1_mean, 
                             pc_majorlength_t1_se/2, pc_minorlength_t1_se/2, pc_majorangle_mean), ribbon_points)) %>%
  pivot_longer(cols = ends_with("ellipse"), names_to = "ellipse_type", values_to = "data") %>%
  unnest(data)

# Average over nloci, doesn't matter; only 1 tau group being looked at
covpca_dplot_op1_axes <- covpca_op1_sum %>% ungroup() %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  mutate(mean_axes = pmap(list(pc_majorlength_t1_mean, 
                               pc_minorlength_t1_mean, pc_majorangle_mean), 
                          axis_lines)) %>%
  unnest(mean_axes)


# Summarise
covpca_op4_sum <- covpca_op4 %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))

# Run ellipse function on data
## Mean ellipse
covpca_dplot_op4 <- covpca_op4_sum %>% ungroup() %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  mutate(
    mean_ellipse = pmap(list(0,0,pc_majorlength_t1_mean,pc_minorlength_t1_mean,pc_majorangle_mean), ellipse_points),
    ribbon_ellipse = pmap(list(pc_majorlength_t1_mean, pc_minorlength_t1_mean, pc_majorlength_t1_se/2, pc_minorlength_t1_se/2, pc_majorangle_mean), ribbon_points)) %>%
  pivot_longer(cols = ends_with("ellipse"), names_to = "ellipse_type", values_to = "data") %>%
  unnest(data)

# Average over nloci, doesn't matter; only 1 tau group being looked at
covpca_dplot_op4_axes <- covpca_op4_sum %>% ungroup() %>%
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

singleTraitsMap <- c(
  "KXZ",
  "KZ",
  "Z",
  "a",
  "b"
)

for (i in seq_along(traitCombos)) {
  combo <- traitCombos[[i]]
  timepoint <- ifelse(i %% 2 != 0, 1, 4)
  
  if (i %% 2 != 0) {
    covpca_dplot <- covpca_dplot_op1
    covpca_dplot_axes <- covpca_dplot_op1_axes
  } else {
    covpca_dplot <- covpca_dplot_op4
    covpca_dplot_axes <- covpca_dplot_op4_axes
  }

  # Subset data
sbst_plt <- covpca_dplot %>% filter(trait1 == combo[1], trait2 == combo[2]) %>%
  mutate(r_title = "Recombination rate (log10)",
         model_title = "Model (cluster)",
         model = fct_recode(model, "K+" = "K", "K-" = "ODE"))

axes_sbst_plt <- covpca_dplot_axes %>% filter(trait1 == combo[1], trait2 == combo[2]) %>%
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
  #lims(x = c(-0.25, 0.25), y = c(-0.25, 0.25)) +
  coord_equal() +
  labs(title = paste0("Progress to the optimum: ", ifelse(timepoint == 1, "<25%", ">=75%")),
       x = traitLabels[combo[1]], y = traitLabels[combo[2]], colour = "Cluster") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 12),
        panel.spacing = unit(1, "lines")) -> plt

res_plt[[i]] <- plt
filename <- paste0("plt_gmat_", singleTraitsMap[combo[1]], "_", singleTraitsMap[combo[2]], 
                   "_", timepoint, "_scaled.png")
ggsave(filename, plt, device = png, width = 8, height = 8)
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
  dplyr::summarise_if(is.numeric, list(mean = mean, CI = CI, sd = sd)) -> covpca_traitprops_sum

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

# Plot additive variances of each model




# Distance matrix: response distances (Hansen and Houle 2008)
dist_matrix_d_op1 <- distanceMatrix(h2_mat_op1, metric = "response-distance")
colnames(dist_matrix_d_op1) <- paste("Matrix", 1:nrow(dist_matrix_d_op1))
rownames(dist_matrix_d_op1) <- colnames(dist_matrix_d_op1)
diag(dist_matrix_d_op1) <- 0

hc <- hclust(as.dist(dist_matrix_d_op1), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 6 seems to be the best
library(factoextra)
# elbow plot
fviz_nbclust(dist_matrix_d_op1, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 6, border = 2:6)

clus <- cutree(hc, 6)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id_op1 <- rbindlist(cov_matrix_modelindex_op1, fill = T)
id_op1$label <- as.character(1:nrow(id_op1))
id_op1$modelindex <- as.factor(id_op1$modelindex)
id_op1 <- AddCombosToDF(id_op1)
id_op1$nloci_group <- "[4, 64)"
id_op1$nloci_group[id_op1$nloci >= 64 & id_op1$nloci < 1024] <- "[64, 256]"
id_op1$nloci_group[id_op1$nloci == 1024] <- "[1024]"
id_op1$nloci_group <- factor(id_op1$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id_op1$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id_op1[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id_op1$clus, id_op1$nloci_group, id_op1$r, id_op1$model, id_op1$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

model <- glm(Freq~cluster,family=poisson(),data=tab)
summary(model)

id_op1 %>% ungroup() %>%
  group_by(r, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_r

id_op1 %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_nloci

phylo <- full_join(as.phylo(phylo), id_op1, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op1_d.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id_op1$clus)) {
  if(length(id_op1$clus[id_op1$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id_op1$label[id_op1$clus == i]), 
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
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id_op1$clus)) {
  if(length(id_op1$clus[id_op1$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id_op1$label[id_op1$clus == i]), 
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
ggsave("plt_tree_gmatrix_optperc1_d.png", device = png, bg = "white",
       width = 12, height = 6)

# Repeat for opt perc 4
dist_matrix_d_op4 <- distanceMatrix(h2_mat_op4, metric = "response-distance")
colnames(dist_matrix_d_op4) <- paste("Matrix", 1:nrow(dist_matrix_d_op4))
rownames(dist_matrix_d_op4) <- colnames(dist_matrix_d_op4)

#fviz_dist(as.dist(dist_matrix), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

hc <- hclust(as.dist(dist_matrix_d_op4), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 3 seems to be the best
# elbow plot
fviz_nbclust(dist_matrix_d_op4, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 6, border = 2:6)

clus <- cutree(hc, 6)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id_op4 <- rbindlist(cov_matrix_modelindex_op4, fill = T)
id_op4$label <- as.character(1:nrow(id_op4))
id_op4$modelindex <- as.factor(id_op4$modelindex)
id_op4 <- AddCombosToDF(id_op4)
id_op4$nloci_group <- "[4, 64)"
id_op4$nloci_group[id_op4$nloci >= 64 & id_op4$nloci < 1024] <- "[64, 256]"
id_op4$nloci_group[id_op4$nloci == 1024] <- "[1024]"
id_op4$nloci_group <- factor(id_op4$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id_op4$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id_op4[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id_op4$clus, id_op4$nloci_group, id_op4$r, id_op4$model, id_op4$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

model <- glm(Freq~cluster+nloci+r,family=poisson(),data=tab)
summary(model)

id_op4 %>% ungroup() %>%
  group_by(r, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_r

id_op4 %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_nloci

phylo <- full_join(as.phylo(phylo), id_op4, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op4_d.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id_op4$clus)) {
  if(length(id_op4$clus[id_op4$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id_op4$label[id_op4$clus == i]), 
                                  fill = clus_palette[i], alpha = 0.2,
                                  type = "encircle", to.bottom = T) 
  
  # Find the position for the annotation
  cluster_tips <- tree_r$data %>% filter(clus == i)
  annotation_x <- min(cluster_tips$x) - 0.5
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
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id_op4$clus)) {
  if(length(id_op4$clus[id_op4$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id_op4$label[id_op4$clus == i]), 
                                          fill = clus_palette[i], alpha = 0.2,
                                          type = "encircle", to.bottom = T) 
  
  # Find the position for the annotation
  cluster_tips <- tree_nloci$data %>% filter(clus == i)
  annotation_x <- min(cluster_tips$x) - 0.5
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
ggsave("plt_tree_gmatrix_optperc4_d.png", device = png, bg = "white",
       width = 12, height = 6)

# Evolvability, respondability, conditional evolvability, autonomy (Hansen Houle 2008)
# Code from Puentes et al. 2016
CalcECRA <- function(matList, id) {
  require(matrixcalc)
  require(Matrix)
  
  PCAdata <- data.frame(
    ev = numeric(length(matList)),
    cev = numeric(length(matList)),
    res = numeric(length(matList)),
    aut = numeric(length(matList)),
    cev_a = numeric(length(matList)),
    cev_b = numeric(length(matList)),
    cev_Z = numeric(length(matList)),
    cev_KXZ = numeric(length(matList)),
    cev_KZ = numeric(length(matList))
  )
  
  PCAdata <- PCAdata %>%
    mutate(optPerc = id$optPerc,
           seed = id$seed,
           modelindex = id$modelindex,
           clus = id$clus)
  
  Hx <- function(x) {
    1/mean(1/x)
  }
  
  Ix <- function(x) {
    var(x)/mean(x)^2
  }
  
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    
    # If its ODE, remove K parameters
    if (all(g[1:2,] == 0)) {
      g <- g[3:5, 3:5]
    }
    # If the matrix isn't positive semi-definite, find the nearest PD
    if (!is.positive.semi.definite(g)) {
      g <- as.matrix(nearPD(g)$mat)
    }
    
    pca <- eigen(g)
    k <- length(pca$values)
    
    PCAdata$ev[i] <- mean(pca$values) #e
    PCAdata$cev[i] <- Hx(pca$values) * (1 + (2*Ix(1/pca$values)) / (k+2) ) #c
    PCAdata$res[i] <- sqrt(mean(pca$values^2)) * (1 - (Ix(pca$values^2) / (4*k+2) ) ) #r
    PCAdata$aut[i] <- (Hx(pca$values) / mean(pca$values)) * (1 + 2 * (Ix(pca$values) + Ix(1/pca$values) - 1 + Hx(pca$values)/mean(pca$values) + 2 * Ix(pca$values) * Ix(1/pca$values)/(k+2))/(k+2)) #a
    
    # What if we look at conditional evolvability of alpha and beta alone?
    # Maybe certain values of KZ and KXZ are more conducive to producing 
    # phenotype changes via alpha and beta mutations?
    # Rearrange
    if (nrow(g) > 3) {
      g <- g[c(3:5, 1:2), c(3:5, 1:2)]
    }
    
    # conditional of trait y on the other traits x
    PCAdata$cev_z[i] <- g[1,1] - g[1,-1] %*% solve(g[-1, -1]) %*% g[-1,1]
    PCAdata$cev_a[i] <- g[2,2] - g[2,-2] %*% solve(g[-2, -2]) %*% g[-2,2]
    PCAdata$cev_b[i] <- g[3,3] - g[3,-3] %*% solve(g[-3, -3]) %*% g[-3,3]
    
    if (nrow(g) > 3) {
      PCAdata$cev_KXZ[i] <- g[4,4] - g[4,-4] %*% solve(g[-4, -4]) %*% g[-4,4]
      PCAdata$cev_KZ[i] <- g[5,5] - g[5,-5] %*% solve(g[-5, -5]) %*% g[-5,5]
    }
    
    }
  
  return(PCAdata)
}

id <- rbindlist(c(cov_matrix_modelindex_op1,
                  cov_matrix_modelindex_op2,
                  cov_matrix_modelindex_op3,
                  cov_matrix_modelindex_op4), 
                fill = T)
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

d_ecr <- CalcECRA(c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4), 
                   id)
d_ecr <- AddCombosToDF(d_ecr)

# Outliers
library(DMwR2)
lofscores <- lofactor(d_ecr$cev, 20)
threshold <- 1.5
outliers <- lofscores > threshold

plot(density(lofscores[lofscores < 4]))

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score")
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)
boxplot(d_ecr[!outliers,]$cev)

# filter out outliers
d_ecr <- d_ecr[!outliers,]

d_ecr_sum <- d_ecr %>%
  group_by(optPerc, model, r) %>%
  dplyr::summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev, colour = model)) +
  facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean conditional evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = res, colour = model)) +
  facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = res_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean respondability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = aut, colour = model)) +
  facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = aut_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean autonomy",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = ev, colour = model)) +
  facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = ev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                        plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                        ncol = 2, labels = "AUTO")

plt_evol <- plot_grid(plt_evol,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol.png", device = png, bg = "white",
       width = 10, height = 7)

# What if we look at conditional evolvability of alpha and beta alone?
# Maybe certain values of KZ and KXZ are more conducive to producing 
# phenotype changes via alpha and beta mutations?
ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev_a, colour = model)) +
  facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_a_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Mean conditional evolvability on $\\alpha_Z$"),
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev_a

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev_b, colour = model)) +
  facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_b_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Mean conditional evolvability on $\\beta_Z$"),
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev_b

leg <- get_legend(plt_cev_a)

plt_cev_ab <- plot_grid(plt_cev_a + theme(legend.position = "none"),
                      plt_cev_b + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO")

plt_cev_ab <- plot_grid(plt_cev_ab,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_cev_ab
ggsave("plt_cev_ab.png", device = png, bg = "white",
       width = 9, height = 4)


# Generate random matrices for null: Morrissey et al. 2019, Walter 2023
# Null differences in genetic effects among treatments (among models)
# sampling breeding values using rbv from a null G-matrix (the average G matrix)
# Resample G matrices across models
d_h2_scaled %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 1) %>%
  mutate(tmpCVA = CVA_a_b,
         CVA_a_b = CVA_a_KZ,
         CVA_a_KZ = if_else(model == "ODE", NA, tmpCVA)) %>%
  group_by(seed, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_seed_optPerc1

d_h2_scaled %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 4) %>%
  mutate(tmpCVA = CVA_a_b,
         CVA_a_b = CVA_a_KZ,
         CVA_a_KZ = if_else(model == "ODE", NA, tmpCVA)) %>%
  group_by(seed, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_seed_optPerc4

lapply(split_h2_seed_optPerc1, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_seed_op1
lapply(split_h2_seed_optPerc4, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_seed_op4

lapply(split_h2_seed_optPerc1, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_seed_op1
lapply(split_h2_seed_optPerc4, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_seed_op4

# Get average across models for each seed
lapply(cov_matrices_seed_op1, function(x) {
  array(unlist(x), dim = c(5, 5, length(x)))
}) -> cov_matrices_seed_op1

h2_avg_op1 <- lapply(cov_matrices_seed_op1, function(x) {
  apply(x, 1:2, mean)
})

lapply(cov_matrices_seed_op4, function(x) {
  array(unlist(x), dim = c(5, 5, length(x)))
}) -> cov_matrices_seed_op4

h2_avg_op4 <- lapply(cov_matrices_seed_op4, function(x) {
  apply(x, 1:2, mean)
})

# Bootstrap krzanowski correlation/subspace test:
# sample two matrices at random, get correlation
# compare to sampling two within r
# data input: dataframe with ids and a column with the matrix

# Nearest positive definite matrix
h2_pd <- lapply(c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4), function(x) {
  if (!is.positive.definite(x)) {as.matrix(nearPD(x)$mat)}
})


krz_in <- id %>%
  mutate(g = h2_pd, #c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4),
         group = interaction(model, r))


# Remove null matrices (no nearest matrix found)
krz_in <- krz_in[!sapply(krz_in$g,is.null)];

bootKrzCor <- function(x, group) {
  require(evolqg)
  require(dplyr)
  grps <- unique(x[,group])
  nGrps <- length(grps)
  # output data frame
  res <- data.frame(group1 = character(length(grps)^2),
                    group2 = character(length(grps)^2),
                    krzCor = numeric(length(grps)^2))
  
  # Temporary data frame for filling inner loop
  res_tmp <- data.frame(group1 = character(length(grps)),
                        group2 = character(length(grps)),
                        krzCor = numeric(length(grps)))
  
  for (i in seq_along(grps)) {
    for (j in seq_along(grps)) {
      # Sample matrices in different groups
      g_1 <- slice_sample(x[group == grps[i]], n = 1)
      g_2 <- slice_sample(x[group == grps[j]], n = 1)
      res_tmp$group1[j] <- as.character(g_1[1,group])
      res_tmp$group2[j] <- as.character(g_2[1,group])
      res_tmp$krzCor[j] <- KrzCor(g_1$g[[1]], g_2$g[[1]])
    }
      indices <- (nGrps*(i-1) + 1):(nGrps*i)
      res[indices,] <- res_tmp
  }
  return(res)
}

library(mcreplicate)
bootKrzCor <- mc_replicate(1000, bootKrzCor(krz_in, "group"))
bootKrzCor <- unnest(as.data.frame(t(bootKrzCor)), cols = everything())

bootKrzCor <- bootKrzCor %>%
  separate(group1, c("model1", "r1"), "\\.",
                       extra = "merge") %>%
  separate(group2, c("model2", "r2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = log10(as.numeric(r1)),
         r2 = log10(as.numeric(r2)),
         model1 = factor(model1, levels = c("ODE", "K")),
         model2 = factor(model2, levels = c("ODE", "K")))

bootKrzCor_sum <- bootKrzCor %>%
  group_by(model1, r1, model2, r2) %>%
  dplyr::summarise(meanKrzCor = mean(krzCor),
            ciKrzCor = CI(krzCor))

ggplot(bootKrzCor_sum, aes(
  x = model1, y = model2
)) +
  facet_nested("Recombination rate 2 (log10))" + 
                 r2 ~ "Recombination rate 1 (log10))" + r1) +
  geom_tile(aes(fill = meanKrzCor)) +
  theme_bw() +
  geom_jitter(data = bootKrzCor, mapping = aes(fill = krzCor),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("K-", "K+")) +
  scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Model 1", y = "Model 2", fill = "Krzanowski Correlation") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))

ggsave("krzcor_r_model.png", device = png, width = 7, height = 5)
  
# Look at just recombination rate - is there any effect?
bootKrzCor_sum <- bootKrzCor %>%
  group_by(r1, r2) %>%
  dplyr::summarise(meanKrzCor = mean(krzCor),
                   ciKrzCor = CI(krzCor))

ggplot(bootKrzCor_sum, aes(
  x = as.factor(r1), y = as.factor(r2)
)) +
  geom_tile(aes(fill = meanKrzCor)) +
  theme_bw() +
  geom_jitter(data = bootKrzCor, mapping = aes(fill = krzCor),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #scale_x_discrete(labels = c("K-", "K+")) +
  #scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", fill = "Krzanowski Correlation") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("krzcor_r.png", device = png, width = 7, height = 5)

# Effect seems to be that decreasing recombination in both models
# makes the major axes of variation less similar, but this is a small effect
# Overall pretty close to 50% mean correlation (expected value)
# Increasing recombination in both models reduces constraint slightly?

# Now just model
bootKrzCor_sum <- bootKrzCor %>%
  group_by(model1, model2) %>%
  dplyr::summarise(meanKrzCor = mean(krzCor),
                   ciKrzCor = CI(krzCor))

ggplot(bootKrzCor_sum, aes(
  x = model1, y = model2
)) +
  geom_tile(aes(fill = meanKrzCor)) +
  theme_bw() +
  geom_jitter(data = bootKrzCor, mapping = aes(fill = krzCor),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("K-", "K+")) +
  scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Model 1", y = "Model 2", fill = "Krzanowski Correlation") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("krzcor_model.png", device = png, width = 7, height = 5)

# This is the strong effect: K- almost always have the same axes
# K+ often have the same axes
# comparing the two are most dissimilar

# ANOVA
bootKrzCor <- bootKrzCor %>%
  mutate(modelCombo = ifelse(model1 != model2, "mix",
                             # paste(as.character(model1), 
                             #       as.character(model2), 
                             #      sep = "_"), 
                             as.character(model1)),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)))

# R combo only explains 2% of the variance, removing
lm.krz <- lm(krzCor ~ modelCombo, data = bootKrzCor)
plot(lm.krz)
summary(lm.krz)
aov.krz <- aov(krzCor ~ as.factor(modelCombo), data = bootKrzCor)
summary(aov.krz)
aov.krz$coefficients
TukeyHSD(aov.krz)
plot(TukeyHSD(aov.krz))

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = krzCor ~ as.factor(modelCombo), data = bootKrzCor)
# 
# $`as.factor(modelCombo)`
#               diff        lwr        upr p adj
# mix-K   -0.1337903 -0.1374235 -0.1301571     0
# ODE-K    0.3851058  0.3809106  0.3893010     0
# ODE-mix  0.5188961  0.5152629  0.5225292     0


# The shared subspace between two matrices is increased by ~0.5 
# when those two matrices are K- instead of one of them being K+
# When both matrices are K- the shared subspace is increased by ~0.4
# relative to both being K+
# When comparing shared subspace when both matrices are K+ to when there
# is one of each, correlation decreases by ~0.13

# Overall, K- subspace is extremely similar, K+ is very similar, and
# there is a difference in the major axis of variation between models
# i.e. the K+ models do use the KZ and KXZ components to adapt, and they
# tend to use them in a similar way (same relative amounts to alpha/beta)

# So what are the PC contributions in the models then?
covpca <- CovMatrixPCA(#h2_pd,
                       c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4), 
                       id)
covpca <- AddCombosToDF(covpca)

covpca_sum <- covpca %>%
  group_by(optPerc, r, model, clus, trait1, trait2) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))


# Rerun analysis on matrices without Z
d_h2_noZ_mrr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mrr_noZ_scaled.csv"), header = F,
                                  col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                                "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                                "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                                "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
                                                "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                                "h2_KXZ"))

d_h2_noZ_mkr <- data.table::fread(paste0(DATA_PATH, "getH2/out_h2_mkr_noZ_scaled.csv"), header = F,
                                  col.names = c("gen", "seed", "modelindex", "VA_Z", "VA_a",
                                                "VA_b", "VA_KZ", "VA_KXZ", "CVA_Z_a", "CVA_Z_b",
                                                "CVA_Z_KZ", "CVA_Z_KXZ", "CVA_a_b", "CVA_a_KZ",
                                                "CVA_a_KXZ", "CVA_b_KZ", "CVA_b_KXZ", 
                                                "CVA_KZ_KXZ", "h2_Z", "h2_a", "h2_b", "h2_KZ",
                                                "h2_KXZ"))

d_h2_noZ_mkr$method <- "mkr"
d_h2_noZ_mrr$method <- "mrr"

d_h2 <- rbind(d_h2_noZ_mkr, d_h2_noZ_mrr)

# Clean data, remove additive
d_h2 <- d_h2 %>%
  distinct(gen, seed, modelindex, method, .keep_all = T) %>%
  dplyr::mutate(modelindex = as.factor(modelindex),
                seed = as.factor(seed)) %>%
  filter(is.na(VA_Z)) %>% distinct()

# inner join optPerc
d_h2 <- left_join(d_h2, d_qg_optPerc, by = c("gen", "seed", "modelindex"))

d_h2 <- AddCombosToDF(d_h2)

d_h2 <- d_h2 %>% filter(model != "Add")

# Counts for each model type: K+ harder to estimate than the other two
table(d_h2$model, d_h2$isAdapted)


# We have many recombination rates: choose a few
r_subsample <- c(1e-10, 1e-5, 1e-1)

# Distribution, how different are the estimates
ggplot(d_h2 %>% filter(isAdapted) %>%
         dplyr::select(gen, seed, modelindex, optPerc, h2_a, method, model) %>%
         distinct() %>%
         pivot_wider(names_from = method, values_from = h2_a), 
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
       aes(x = method, y = h2_a)) +
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

boxplot(d_h2$VA_a)
d_h2_all <- d_h2

# Detect outliers: Hampel filter
# variance depends on the tau group mainly - check outliers within groups
library(DMwR2)

lofscores <- lofactor(scale(d_h2$VA_a), 10)
threshold <- 1.5
outliers <- (lofscores > threshold)
outliers[is.na(outliers)] <- F

plot(density(lofscores[lofscores < 4 & !is.nan(lofscores)]))

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score")
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)
boxplot(d_h2[!outliers,]$VA_a)

# filter out outliers
d_h2 <- d_h2[!outliers,]
write.csv(d_h2, "d_h2_outliersremoved_noZ.csv")
d_h2 <- read.csv("d_h2_outliersremoved_noZ.csv")
d_h2 <- d_h2[,-1]


d_h2_sum <- d_h2 %>% 
  filter(r %in% r_subsample) %>%
  group_by(optPerc, model, tau, r, method, isAdapted) %>%
  dplyr::summarise(meanH2a = mean(h2_a, na.rm = T),
                   seH2a = se(h2_a, na.rm = T),
                   meanVAa = mean(VA_a, na.rm = T),
                   seVAa = se(VA_a, na.rm = T))
d_h2_sum$model <- as.factor(d_h2_sum$model)

# Number of loci doesn't seem to affect it too much, average across
ggplot(d_h2 %>% filter(isAdapted) %>%
         filter(method == "mkr", r %in% r_subsample) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance"),
       aes(x = optPerc, y = h2_a, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_h2_sum %>% ungroup() %>%
               filter(method == "mkr", r %in% r_subsample, isAdapted) %>% 
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance"),
             aes(x = optPerc, y = meanH2a, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Narrow-sense heritability $(h^2)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-")) +
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
       aes(x = optPerc, y = VA_a, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.0125),
             aes(x = optPerc, y = meanVAa, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_sml

ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 0.125),
       aes(x = optPerc, y = VA_a, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 0.125),
             aes(x = optPerc, y = meanVAa, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 0.7)) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes=list(shape = 15, size = 5))) +
  theme(text = element_text(size = 14),
        legend.position = "bottom") -> plt_add_va_med

ggplot(d_h2 %>% filter(isAdapted) %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci",
                tau_title = "Mutational effect size variance") %>%
         filter(method == "mkr", r %in% r_subsample, tau == 1.25),
       aes(x = optPerc, y = VA_a, colour = model)) +
  facet_nested(r_title + log10(r) ~ tau_title + tau) +
  geom_quasirandom(dodge.width = 0.9) +
  geom_point(data = d_h2_sum %>% filter(isAdapted) %>%
               mutate(r_title = "Recombination rate (log10)",
                      nloci_title = "Number of loci",
                      tau_title = "Mutational effect size variance") %>%
               filter(method == "mkr", r %in% r_subsample, tau == 1.25),
             aes(x = optPerc, y = meanVAa, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  #coord_cartesian(ylim = c(0, 2.5)) +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Additive variance $(V_A)$"),
       colour = "Model") +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-")) +
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
ggsave("plt_va_az_scaled.png", device = png, bg = "white",
       width = 560*4, height = 980*4, units = "px")

# G matrix analysis
d_h2$optPerc <- factor(d_h2$optPerc)

d_h2 %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 1) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc1

d_h2 %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 2) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc2

d_h2 %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 3) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc3

d_h2 %>% filter(isAdapted) %>%
  filter(model != "Add", tau == 0.0125, r %in% r_subsample) %>%
  filter(as.numeric(optPerc) == 4) %>%
  group_by(modelindex, optPerc, method, isAdapted) %>%
  group_split(.) -> split_h2_optPerc4


# Separate into model indices
# each sublist is replicates of a model index
Rcpp::sourceCpp("./getCovarianceMatrices.cpp")
lapply(split_h2_optPerc1, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op1
lapply(split_h2_optPerc2, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op2
lapply(split_h2_optPerc3, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op3
lapply(split_h2_optPerc4, function(x) {extractCovarianceMatrices(as.data.frame(x))}) -> cov_matrices_op4

lapply(split_h2_optPerc1, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op1
lapply(split_h2_optPerc2, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op2
lapply(split_h2_optPerc3, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op3
lapply(split_h2_optPerc4, function(x) {data.frame(optPerc = x$optPerc, seed = x$seed, modelindex = x$modelindex, isAdapted = x$isAdapted)}) -> cov_matrix_modelindex_op4

# We want to know if certain architectures are more/less important for describing
# variation between simulations and which components are most important for describing
# those differences

h2_mat_op1 <- unlist(cov_matrices_op1, recursive = F)
h2_mat_op2 <- unlist(cov_matrices_op2, recursive = F)
h2_mat_op3 <- unlist(cov_matrices_op3, recursive = F)
h2_mat_op4 <- unlist(cov_matrices_op4, recursive = F)

# get ids
GetMatrixIDs <- function(matList) {
  lapply(matList, function(x) {
    data.frame(optPerc = x$optPerc, 
               seed = x$seed, 
               modelindex = x$modelindex, 
               isAdapted = x$isAdapted)}) -> matList
  
  
  lapply(matList, function(x) {
    split(x, seq(nrow(x)))
  }) -> matList
  # unlist to full form
  matList <- unlist(matList, recursive = F)
  return(matList)
}

cov_matrix_modelindex_op1 <- GetMatrixIDs(split_h2_optPerc1)
cov_matrix_modelindex_op2 <- GetMatrixIDs(split_h2_optPerc2)
cov_matrix_modelindex_op3 <- GetMatrixIDs(split_h2_optPerc3)
cov_matrix_modelindex_op4 <- GetMatrixIDs(split_h2_optPerc4)


# Distance between G matrices
# dist matrix for optperc 1
Rcpp::sourceCpp("./distanceFunctions.cpp")

dist_matrix_op1 <- distanceMatrix(h2_mat_op1)
colnames(dist_matrix_op1) <- paste("Matrix", 1:nrow(dist_matrix_op1))
rownames(dist_matrix_op1) <- colnames(dist_matrix_op1)

#fviz_dist(as.dist(dist_matrix), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

hc <- hclust(as.dist(dist_matrix_op1), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 2 seems to be the best
library(factoextra)
# elbow plot
fviz_nbclust(dist_matrix_op1, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 2, border = 2)

clus <- cutree(hc, 2)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id_op1 <- rbindlist(cov_matrix_modelindex_op1, fill = T)
id_op1$label <- as.character(1:nrow(id_op1))
id_op1$modelindex <- as.factor(id_op1$modelindex)
id_op1 <- AddCombosToDF(id_op1)
id_op1$nloci_group <- "[4, 64)"
id_op1$nloci_group[id_op1$nloci >= 64 & id_op1$nloci < 1024] <- "[64, 256]"
id_op1$nloci_group[id_op1$nloci == 1024] <- "[1024]"
id_op1$nloci_group <- factor(id_op1$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id_op1$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id_op1[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id_op1$clus, id_op1$nloci_group, id_op1$r, id_op1$model, id_op1$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

# Model describes the clustering
model <- glm(Freq~model+nloci,family=poisson(),data=tab)
summary(model)

id_op1 %>% ungroup() %>%
  group_by(r, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_r

id_op1 %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_nloci

id_op1 %>% ungroup() %>%
  group_by(optPerc, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_optPerc

phylo <- full_join(as.phylo(phylo), id_op1, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op1_noZ.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id_op1$clus)) {
  if(length(id_op1$clus[id_op1$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id_op1$label[id_op1$clus == i]), 
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
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id_op1$clus)) {
  if(length(id_op1$clus[id_op1$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id_op1$label[id_op1$clus == i]), 
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
ggsave("plt_tree_gmatrix_optperc1_noZ.png", device = png, bg = "white",
       width = 12, height = 6)

# Progress to optimum <25%: very distinct clustering by model type
# early in the walk variance-covariance is different between the models
# outlier K+ clusters have lower recombination maybe




# Repeat for opt perc 4
dist_matrix_op4 <- distanceMatrix(h2_mat_op4)
colnames(dist_matrix_op4) <- paste("Matrix", 1:nrow(dist_matrix_op4))
rownames(dist_matrix_op4) <- colnames(dist_matrix_op4)

#fviz_dist(as.dist(dist_matrix), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

hc <- hclust(as.dist(dist_matrix_op4), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 2 seems to be the best again
# elbow plot
fviz_nbclust(dist_matrix_op4, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc)
rect.hclust(hc, 2, border = 2)

clus <- cutree(hc, 2)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id_op4 <- rbindlist(cov_matrix_modelindex_op4, fill = T)
id_op4$label <- as.character(1:nrow(id_op4))
id_op4$modelindex <- as.factor(id_op4$modelindex)
id_op4 <- AddCombosToDF(id_op4)
id_op4$nloci_group <- "[4, 64)"
id_op4$nloci_group[id_op4$nloci >= 64 & id_op4$nloci < 1024] <- "[64, 256]"
id_op4$nloci_group[id_op4$nloci == 1024] <- "[1024]"
id_op4$nloci_group <- factor(id_op4$nloci_group, levels = c("[4, 64)", "[64, 256]", "[1024]"))

id_op4$clus <- -1
# add cluster
for (i in 1:length(g)) {
  idx <- g[[i]]
  id_op4[idx,"clus"] <- i
}

# with id, check how frequent genetic architectures are with the clusters
tab <- table(id_op4$clus, id_op4$nloci_group, id_op4$r, id_op4$model, id_op4$isAdapted)
names(dimnames(tab)) <- c("cluster", "nloci", "r", "model", "isAdapted")
tab <- as.data.frame(tab)

model <- glm(Freq~model+nloci,family=poisson(),data=tab)
summary(model)

id_op4 %>% ungroup() %>%
  group_by(r, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_r

id_op4 %>% ungroup() %>%
  group_by(nloci_group, model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs_nloci

phylo <- full_join(as.phylo(phylo), id_op4, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op4_noZ.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(aes(shape = as.factor(log10(r))), size = 3) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", shape = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_r

# add clusters + proportions
for (i in unique(id_op4$clus)) {
  if(length(id_op4$clus[id_op4$clus == i]) < 2) next
  lab_dat <- cluster_percs_r[cluster_percs_r$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("r: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_r <- tree_r + geom_hilight(node = MRCA(phylo, id_op4$label[id_op4$clus == i]), 
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
  ggtitle("Progress to the optimum: >=75%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) -> tree_nloci

# add clusters + proportions
for (i in unique(id_op4$clus)) {
  if(length(id_op4$clus[id_op4$clus == i]) < 2) next
  lab_dat <- cluster_percs_nloci[cluster_percs_nloci$clus == i,]
  cluster_labels <- apply(lab_dat, 1, function(x) {
    sprintf("nloci: %s, model: %s = %.1f%%",
            x[1], x[2], as.numeric(x[5]) * 100)})
  tree_nloci <- tree_nloci + geom_hilight(node = MRCA(phylo, id_op4$label[id_op4$clus == i]), 
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
ggsave("plt_tree_gmatrix_optperc4_noZ.png", device = png, bg = "white",
       width = 12, height = 6)


# Analysis across all timepoints, time doesn't seem to matter too much
dist_matrix <- distanceMatrix(c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4))
colnames(dist_matrix) <- paste("Matrix", 1:nrow(dist_matrix))
rownames(dist_matrix) <- colnames(dist_matrix)

hc <- hclust(as.dist(dist_matrix), method="average")
#plot(as.phylo(hc), type="phylogram", main="Phylogenetic Tree of G Matrices")

# number of clusters: 2 seems to be the best
library(factoextra)
# elbow plot
fviz_nbclust(dist_matrix, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

# dendrogram
plot(hc, main = "Power Euclidean distances between molecular G matrices", labels = F)
rect.hclust(hc, 2, border = 2)

clus <- cutree(hc, 2)
g <- split(names(clus), clus)
g <- lapply(g, function(x) as.numeric(substring(x, 8)))

phylo <- as.phylo(hc)
phylo <- as_tibble(phylo)
phylo$label <- as.numeric(substring(phylo$label, 8))
phylo <- as.phylo(phylo)
id <- rbindlist(c(cov_matrix_modelindex_op1,
                  cov_matrix_modelindex_op2,
                  cov_matrix_modelindex_op3,
                  cov_matrix_modelindex_op4), 
                fill = T)
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
tab <- table(id$clus, id$model)
names(dimnames(tab)) <- c("cluster", "model")
tab <- as.data.frame(tab)

# Model describes the clustering
glm.clus <- glm(Freq~model,family=poisson(),data=tab)
summary(glm.clus)
report::report(glm.clus)

id %>% ungroup() %>%
  group_by(model, clus) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(clus) %>%
  dplyr::mutate(prop = n/sum(n)) -> cluster_percs

phylo <- full_join(as.phylo(phylo), id, by = "label")

clus_palette <- paletteer_d("ggsci::nrc_npg", 3)

ggtree(phylo, aes(colour = as.factor(clus)), layout="equal_angle") +
  #geom_text(aes(label=node)) +
  # scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
  #                     labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model", size = "Recombination rate (log10)") +
  ggtitle("Progress to the optimum: <25%") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) -> tree_clus
tree_clus

ggsave("tree_clus_op1_noZ.png", device = png, width = 4, height = 4)

ggtree(phylo, aes(colour = as.factor(model)), layout="equal_angle") +
  geom_tippoint(size = 2) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(colour = "Model") +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(shape=16, size = 5,
                                                   linetype = 0))) -> tree_full

tree_full
ggsave("plt_tree_gmatrix_full_noZ.png", device = png, bg = "white",
       width = 7/2, height = 9/2)


CalcECRA <- function(matList, id, noZ = F) {
  require(matrixcalc)
  require(Matrix)
  
  PCAdata <- data.frame(
    ev = numeric(length(matList)),
    cev = numeric(length(matList)),
    res = numeric(length(matList)),
    aut = numeric(length(matList)),
    cev_a = numeric(length(matList)),
    cev_b = numeric(length(matList)),
    cev_Z = numeric(length(matList)),
    cev_KXZ = numeric(length(matList)),
    cev_KZ = numeric(length(matList))
  )
  
  PCAdata <- PCAdata %>%
    mutate(optPerc = id$optPerc,
           seed = id$seed,
           modelindex = id$modelindex,
           clus = id$clus)
  
  Hx <- function(x) {
    1/mean(1/x)
  }
  
  Ix <- function(x) {
    var(x)/mean(x)^2
  }
  
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    
    # If its ODE, remove K parameters
    if (all(g[1:2,] < 1e-5)) {
      g <- g[3:5, 3:5]
    }
    
    # If we are a no Z run, remove the Z
    if (noZ) {
      # Z index changes depending on K+/K-
      if (nrow(g) == 5) {
        g <- g[-3, -3]
      } else {
        g <- g[-1, -1]
      }
    }
    
    
    # If the matrix isn't positive semi-definite, find the nearest PD
    if (!is.positive.semi.definite(g)) {
      g <- as.matrix(nearPD(g)$mat)
    }
    
    pca <- eigen(g)
    k <- length(pca$values)
    
    PCAdata$ev[i] <- mean(pca$values) #e
    PCAdata$cev[i] <- Hx(pca$values) * (1 + (2*Ix(1/pca$values)) / (k+2) ) #c
    PCAdata$res[i] <- sqrt(mean(pca$values^2)) * (1 - (Ix(pca$values^2) / (4*k+2) ) ) #r
    PCAdata$aut[i] <- (Hx(pca$values) / mean(pca$values)) * (1 + 2 * (Ix(pca$values) + Ix(1/pca$values) - 1 + Hx(pca$values)/mean(pca$values) + 2 * Ix(pca$values) * Ix(1/pca$values)/(k+2))/(k+2)) #a
    
    # What if we look at conditional evolvability of alpha and beta alone?
    # Maybe certain values of KZ and KXZ are more conducive to producing 
    # phenotype changes via alpha and beta mutations?
    # Rearrange
    if (nrow(g) > 3 && !noZ) {
      g <- g[c(3:5, 1:2), c(3:5, 1:2)]
    }
    
    if (nrow(g) > 3 && noZ) {
      g <- g[c(3:4, 1:2), c(3:4, 1:2)]
    }
    
    
    # conditional of trait y on the other traits x
    if (!noZ) {
      PCAdata$cev_z[i] <- g[1,1] - g[1,-1] %*% solve(g[-1, -1]) %*% g[-1,1]
      PCAdata$cev_a[i] <- g[2,2] - g[2,-2] %*% solve(g[-2, -2]) %*% g[-2,2]
      PCAdata$cev_b[i] <- g[3,3] - g[3,-3] %*% solve(g[-3, -3]) %*% g[-3,3]
    }
    
    if (noZ) {
      PCAdata$cev_a[i] <- g[1,1] - g[1,-1] %*% solve(g[-1, -1]) %*% g[-1,1]
      PCAdata$cev_b[i] <- g[2,2] - g[2,-2] %*% solve(g[-2, -2]) %*% g[-2,2]
    }
    
    if (nrow(g) > 3 && !noZ) {
      PCAdata$cev_KXZ[i] <- g[4,4] - g[4,-4] %*% solve(g[-4, -4]) %*% g[-4,4]
      PCAdata$cev_KZ[i] <- g[5,5] - g[5,-5] %*% solve(g[-5, -5]) %*% g[-5,5]
    }
    
    if (nrow(g) > 3 && noZ) {
      PCAdata$cev_KXZ[i] <- g[3,3] - g[3,-3] %*% solve(g[-3, -3]) %*% g[-3,3]
      PCAdata$cev_KZ[i] <- g[4,4] - g[4,-4] %*% solve(g[-4, -4]) %*% g[-4,4]
    }
  }
  
  return(PCAdata)
}

id <- rbindlist(c(cov_matrix_modelindex_op1,
                  cov_matrix_modelindex_op2,
                  cov_matrix_modelindex_op3,
                  cov_matrix_modelindex_op4), 
                fill = T)
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

d_ecr <- CalcECRA(h2_pd, 
                  #c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4), 
                  id, noZ = T)
d_ecr <- AddCombosToDF(d_ecr)

# Remove Z -> not included in this model
d_ecr <- d_ecr %>% select(-cev_Z)

# Outliers
library(DMwR2)
lofscores <- lofactor(d_ecr$cev, 20)
threshold <- 1.5
outliers <- lofscores > threshold

plot(density(lofscores[lofscores < 4]))

plot(lofscores, pch = 1, col = ifelse(outliers, "red", "blue"),
     main = "LOF Outlier Detection (k = 15)", xlab = "Data Point", 
     ylab = "LOF Score")
legend("topright", legend = c("Outlier", "Inlier"), col = c("red", "blue"), 
       pch = 1)
boxplot(d_ecr[!outliers,]$cev)

# filter out outliers
d_ecr <- d_ecr[!outliers,]

# Need to calculate cev means separately for the K- and K+ models
# K- shouldn't mean over cev_KZ and KXZ

d_ecr_sum <- d_ecr %>%
  group_by(optPerc, model, r) %>%
  dplyr::summarise_if(is.numeric, list(mean = mean, se = se))


ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean conditional evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_cev

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = res, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = res_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean respondability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_res

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = aut, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = aut_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean autonomy",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_aut

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = ev, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = ev_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Mean evolvability",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 14)) -> plt_ev

leg <- get_legend(plt_ev)

plt_evol <- plot_grid(plt_ev + theme(legend.position = "none"),
                      plt_cev + theme(legend.position = "none"),
                      plt_res + theme(legend.position = "none"),
                      plt_aut + theme(legend.position = "none"),
                      ncol = 2, labels = "AUTO")

plt_evol <- plot_grid(plt_evol,
                      leg, nrow = 2, rel_heights = c(1, 0.05))
plt_evol
ggsave("plt_evol_noZ.png", device = png, bg = "white",
       width = 10, height = 7)

# What if we look at conditional evolvability of alpha and beta alone?
# Maybe certain values of KZ and KXZ are more conducive to producing 
# phenotype changes via alpha and beta mutations?
ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev_a, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_a_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Mean conditional evolvability on $\\alpha_Z$"),
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev_a

ggplot(d_ecr %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev_b, colour = model)) +
  facet_nested(r_title + log10(r)~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_b_mean, group = model), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", 
       y = TeX("Mean conditional evolvability on $\\beta_Z$"),
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev_b

leg <- get_legend(plt_cev_a)

plt_cev_ab <- plot_grid(plt_cev_a + theme(legend.position = "none"),
                        plt_cev_b + theme(legend.position = "none"),
                        ncol = 2, labels = "AUTO")

plt_cev_ab <- plot_grid(plt_cev_ab,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_cev_ab
ggsave("plt_cev_ab_noZ.png", device = png, bg = "white",
       width = 9, height = 4)

# How about on KZ and KXZ, does that allow the populations to escape from
# low-fitness combinations?
ggplot(d_ecr %>% filter(model == "K") %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev_KZ, colour = as.factor(log10(r)))) +
  #facet_nested(r_title + r~.) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>% filter(model == "K") %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_KZ_mean, group = as.factor(log10(r))), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("rcartocolor::TealGrn", 6)[c(1,4,6)]) +
  labs(x = "Progress to the optimum", 
       y = TeX("Mean conditional evolvability on $K_Z$"),
       colour = "Recombination rate (log10)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev_KZ

ggplot(d_ecr %>% filter(model == "K") %>%
         mutate(r_title = "Recombination rate (log10)",
                nloci_title = "Number of loci"), 
       aes(x = optPerc, y = cev_KXZ, colour = as.factor(log10(r)))) +
  geom_quasirandom(shape = 1, dodge.width = 0.9, na.rm = F) +
  geom_point(data = d_ecr_sum %>% ungroup() %>% filter(model == "K") %>%
               mutate(r_title = "Recombination rate (log10)"),
             aes(x = optPerc, y = cev_KXZ_mean, group = as.factor(log10(r))), colour = "black",
             shape = 3, size = 2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("rcartocolor::TealGrn", 6)[c(1,4,6)]) +
  labs(x = "Progress to the optimum", 
       y = TeX("Mean conditional evolvability on $K_{XZ}$"),
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12)) -> plt_cev_KXZ

leg <- get_legend(plt_cev_KXZ)

plt_cev_K <- plot_grid(plt_cev_KZ + theme(legend.position = "none"),
                       plt_cev_KXZ + theme(legend.position = "none"),
                        ncol = 2, labels = "AUTO")

plt_cev_K <- plot_grid(plt_cev_K,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_cev_K
ggsave("plt_cev_K_noZ.png", device = png, bg = "white",
       width = 9, height = 4)

# So what we see is that there aren't many differences in the amount of available 
# variation for models. mean evolvability is similar between K+ and K-
# conditional evolvability is slightly lower in K+, suggesting stronger correlations
# which makes sense, there are more molecular traits - redcombination removes this
# There is marginally higher respondability under low recombination in many K+ 
# models, which might explain the rapid increase relative to K- models.
# autonomy suggests K+ are more constrained and recombination alleviates this
# In K+ all the traits suffer from less conditional evolvability, alleviated
# by increased recombination rate.

# So despite increased constraint on the amount of variance, the expected response
# to directional selection is still greater in K+ models
# Could be that although much of the variance is restricted, the additional
# axes of variation allow for more directions that can be explored, for less
# redundancy and more options to reach the optimal phenotype - i.e. more
# beneficial mutations

# Stats -> compare K+ and K- among recombination rates 
# (nloci doesn't affect, neither does optPerc)
hist(d_ecr$cev)
summary(ols <- lm(cev ~ model * as.factor(r), d_ecr))
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)

ggplot(d_ecr, 
       aes(x = model, y = cev, colour = as.factor(r))) +
  geom_quasirandom(dodge.width = 0.9)

# Variance differs between groups, use gls to account for unequal variance
library(nlme)
summary(gls.cev <- gls(cev ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
summary(ols.cev <- gls(cev ~ model * as.factor(r), d_ecr))
plot(gls.cev)
report(gls.cev)

# the variance between groups does improve the model, should use gls
anova(ols.cev, gls.cev)
anova(gls.cev)

# Marginal means
library(emmeans)
library(xtable)
em.cev <- emmeans(gls.cev, ~ model * r)
pairs(em.cev, simple = "model")
pairs(em.cev, simple = "r")
plot(em.cev, comparisons = T)

xtable(em.cev)


# Format to a nice table
library(stargazer)
stargazer(gls.cev)

# Repeat for autonomy, respondability, evolvability
summary(gls.aut <- gls(aut ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))

anova(gls.aut)
plot(gls.aut)
em.aut <- emmeans(gls.aut, ~ model * r)
pairs(em.aut, simple = "model")
pairs(em.aut, simple = "r")
plot(em.aut) # one of the comparisons has negative length because error is so small (ODE 0.1)

xtable(em.aut)

stargazer(gls.aut)

summary(gls.res <- gls(res ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.res)
anova(gls.res)

em.res <- emmeans(gls.res, ~ model * r)
pairs(em.res, simple = "model")
pairs(em.res, simple = "r")
plot(em.res, comparisons = T)
stargazer(gls.res)

xtable(em.res)


summary(gls.ev <- gls(ev ~ model * as.factor(r), d_ecr, 
                       weights = varIdent(form = ~ 1 | model * as.factor(r))))
plot(gls.ev)
anova(gls.ev)
em.ev <- emmeans(gls.ev, ~ model * r, type = "response")
pairs(em.ev, simple = "model")
pairs(em.ev, simple = "r")
plot(em.ev, comparisons = T)

stargazer(gls.ev)
xtable(em.ev)




# Will need to look at the shape of variance: Krzanowski correlation

# Bootstrap krzanowski correlation/subspace test:
# sample two matrices at random, get correlation
# compare to sampling two within r
# data input: dataframe with ids and a column with the matrix

# Nearest positive definite matrix
library(matrixcalc)
library(Matrix)
h2_pd <- lapply(c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4), function(x) {
  if (!is.positive.definite(x)) {as.matrix(nearPD(x)$mat)}
})


krz_in <- id %>%
  dplyr::mutate(g = h2_pd, #c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4),
         group = interaction(model, r))


# Remove null matrices (no nearest matrix found)
krz_in <- krz_in[!sapply(krz_in$g,is.null)];

bootKrzCorFn <- function(x, group = "", PCASim = F) {
  require(evolqg)
  require(dplyr)
  
  fn <- ifelse(PCASim, evolqg::PCAsimilarity, evolqg::KrzCor)
  
  if (group != "") {
    grps <- unique(x[,group])
    nGrps <- length(grps)
  
    
    # output data frame
    res <- data.frame(group1 = character(length(grps)^2),
                      group2 = character(length(grps)^2),
                      krzCor = numeric(length(grps)^2))
    
    # Temporary data frame for filling inner loop
    res_tmp <- data.frame(group1 = character(length(grps)),
                          group2 = character(length(grps)),
                          krzCor = numeric(length(grps)))
    
    for (i in seq_along(grps)) {
      for (j in seq_along(grps)) {
        # Sample matrices in different groups
        g_1 <- slice_sample(x[group == grps[i]], n = 1)
        g_2 <- slice_sample(x[group == grps[j]], n = 1)
        res_tmp$group1[j] <- as.character(g_1[1,group])
        res_tmp$group2[j] <- as.character(g_2[1,group])
        res_tmp$krzCor[j] <- fn(g_1$g[[1]], g_2$g[[1]])
      }
      indices <- (nGrps*(i-1) + 1):(nGrps*i)
      res[indices,] <- res_tmp
    }
    return(res)
  }
  
  # If group is "", sample two random matrices and return that
  g1 <- slice_sample(x, n = 1)
  g2 <- slice_sample(x, n = 1)
  return(fn(g1$g[[1]], g2$g[[1]]))
}

# Bootstrap in ten parts for RAM reasons
# exclude r -> doesn't contribute much variance at all
library(mcreplicate)
# newseed <- sample(1:.Machine$integer.max, 10)
# 1360932387 1900268993  991875895 1523108407  197897667  199526283 2070940443
# 128221287 1383031956  970870370
newseed <- c(1360932387, 1900268993,  991875895, 1523108407, 197897667, 199526283, 
             2070940443, 128221287, 1383031956, 970870370)
set.seed(1360932387)

bootKrzCor <- vector(mode = "list", length = 10)

for (i in seq_along(newseed)) {
  # Set seed
  set.seed(newseed[i])
  
  # Run replicate
  res <- mc_replicate(1000, bootKrzCorFn(krz_in, "model"))
  bootKrzCor[[i]] <- unnest(as.data.frame(t(res)), cols = everything())
}

# Output list into combined df
bootKrzCor2 <- bind_rows(bootKrzCor)

bootKrzCor <- bootKrzCor2 %>%
  separate(group1, c("model1", "r1"), "\\.",
           extra = "merge") %>%
  separate(group2, c("model2", "r2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = log10(as.numeric(r1)),
         r2 = log10(as.numeric(r2)),
         model1 = factor(model1, levels = c("ODE", "K")),
         model2 = factor(model2, levels = c("ODE", "K")))

bootKrzCor_sum <- bootKrzCor %>%
  group_by(model1, r1, model2, r2) %>%
  dplyr::summarise(meanKrzCor = mean(krzCor),
                   ciKrzCor = CI(krzCor))

ggplot(bootKrzCor_sum, aes(
  x = model1, y = model2
)) +
  facet_nested("Recombination rate 2 (log10))" + 
                 r2 ~ "Recombination rate 1 (log10))" + r1) +
  geom_tile(aes(fill = meanKrzCor)) +
  theme_bw() +
  geom_jitter(data = bootKrzCor, mapping = aes(fill = krzCor),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("K-", "K+")) +
  scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Model 1", y = "Model 2", fill = "Krzanowski Correlation") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))

ggsave("krzcor_r_model_noZ.png", device = png, width = 7, height = 5)

# Look at just recombination rate - is there any effect?
bootKrzCor_sum <- bootKrzCor %>%
  group_by(r1, r2) %>%
  dplyr::summarise(meanKrzCor = mean(krzCor),
                   ciKrzCor = CI(krzCor))

ggplot(bootKrzCor_sum, aes(
  x = as.factor(r1), y = as.factor(r2)
)) +
  geom_tile(aes(fill = meanKrzCor)) +
  theme_bw() +
  geom_jitter(data = bootKrzCor, mapping = aes(fill = krzCor),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #scale_x_discrete(labels = c("K-", "K+")) +
  #scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", fill = "Krzanowski Correlation") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("krzcor_r_noZ.png", device = png, width = 7, height = 5)

# Effect seems to be that decreasing recombination in both models
# makes the major axes of variation less similar, but this is a small effect
# Overall pretty close to 50% mean correlation (expected value)
# Increasing recombination in both models reduces constraint slightly?

# Now just model
bootKrzCor_sum <- bootKrzCor %>%
  group_by(model1, model2) %>%
  dplyr::summarise(meanKrzCor = mean(krzCor),
                   ciKrzCor = CI(krzCor))

ggplot(bootKrzCor_sum, aes(
  x = model1, y = model2
)) +
  geom_tile(aes(fill = meanKrzCor)) +
  theme_bw() +
  geom_jitter(data = bootKrzCor, mapping = aes(fill = krzCor),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("K-", "K+")) +
  scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Model 1", y = "Model 2", fill = "Krzanowski Correlation") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("krzcor_model_noZ.png", device = png, width = 7, height = 5)

# This is the strong effect: K- almost always have the same axes
# K+ often have the same axes
# comparing the two are most dissimilar

bootKrzCor <- bootKrzCor %>%
  mutate(modelCombo = ifelse(model1 != model2, "mix",
                             # paste(as.character(model1), 
                             #       as.character(model2), 
                             #      sep = "_"), 
                             as.character(model1)),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)))

# R doesn't explain much variance, removed
# beta regression for 0-1 data
library(report)
library(betareg)

# There is no variance in the ODE models: all have the same eigenstructure
# Makes sense, only two components to alter
ggplot(bootKrzCor, 
       aes(x = modelCombo, y = krzCor)) +
  geom_quasirandom(dodge.width = 0.9)

# 0 variance
var(bootKrzCor[bootKrzCor$modelCombo == "ODE",]$krzCor)
mean(bootKrzCor[bootKrzCor$modelCombo == "ODE",]$krzCor)

# t test is probably the better estimate here, compare K and mix to 1, to test
# if they are different from having at least one K+

# First test assumptions: normality, equal variance
# Plot looks like the variance is slightly different -> Welch's t test
# roughly normal? a fw bumps here and there but it is symmetric
# fat tailed

krz.KK <- bootKrzCor[bootKrzCor$modelCombo == "K", "krzCor"]
krz.KO <- bootKrzCor[bootKrzCor$modelCombo == "mix", "krzCor"]
krz.OO <- mean(bootKrzCor[bootKrzCor$modelCombo == "ODE",]$krzCor)

t.KK <- t.test(krz.KK, mu = krz.OO)
# One Sample t-test
# 
# data:  krz.KK
# t = -265.16, df = 8999, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 1
# 95 percent confidence interval:
#   0.5173946 0.5244776
# sample estimates:
#   mean of x 
# 0.5209361 

t.KO <- t.test(krz.KO, mu = krz.OO)
# One Sample t-test
# 
# data:  krz.KO
# t = -486.5, df = 17999, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 1
# 95 percent confidence interval:
#   0.4942841 0.4983428
# sample estimates:
#   mean of x 
# 0.4963134 

# Compare the two groups
t.KK.KO <- t.test(krz.KK, krz.KO)

# Welch Two Sample t-test
# 
# data:  krz.KK and krz.KO
# t = 11.825, df = 15067, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02054108 0.02870425
# sample estimates:
#   mean of x mean of y 
# 0.5209361 0.4963134 

krz.p <- c(t.KK$p.value, t.KO$p.value, t.KK.KO$p.value)

# Correct for multiple comparisons
p.adjust(krz.p, method = "bonferroni")

# Essentially, yes both groups are different from 1, so the K+ models use different
# combinations of molecular components to adapt. The K- are limited to aZ and bZ,
# so of course there is no opportunity to change that, but they can change the 
# amount of variance each component uses, which PCA similarity can investigate


# Repeat using PCA similarity to include amount of variation in morphospace
# newseed <- sample(1:.Machine$integer.max, 10)
# 1360932387 1900268993  991875895 1523108407  197897667  199526283 2070940443
# 128221287 1383031956  970870370
newseed <- c(1360932387, 1900268993,  991875895, 1523108407, 197897667, 199526283, 
             2070940443, 128221287, 1383031956, 970870370)

bootPCASim <- vector(mode = "list", length = 10)

for (i in seq_along(newseed)) {
  # Set seed
  set.seed(newseed[i])
  
  # Run replicate
  res <- mcreplicate::mc_replicate(1000, bootKrzCorFn(krz_in, "group", T))
  bootPCASim[[i]] <- unnest(as.data.frame(t(res)), cols = everything())
}

# Output list into combined df
bootPCASim2 <- bind_rows(bootPCASim)

# Null distribution
bootPCASim_null <- mcreplicate::mc_replicate(10000, bootKrzCor(krz_in, PCASim = T))

hist(bootPCASim_null)

bootPCASim <- bootPCASim2 %>%
  separate(group1, c("model1", "r1"), "\\.",
           extra = "merge") %>%
  separate(group2, c("model2", "r2"), "\\.",
           extra = "merge") %>%
  mutate(r1 = log10(as.numeric(r1)),
         r2 = log10(as.numeric(r2)),
         model1 = factor(model1, levels = c("ODE", "K")),
         model2 = factor(model2, levels = c("ODE", "K"))) %>%
  dplyr::rename(PCASim = krzCor)

bootPCASim_sum <- bootPCASim %>%
  group_by(model1, r1, model2, r2) %>%
  dplyr::summarise(meanPCASim = mean(PCASim),
                   ciPCASim = CI(PCASim))

ggplot(bootPCASim_sum, aes(
  x = model1, y = model2
)) +
  facet_nested("Recombination rate 2 (log10))" + 
                 r2 ~ "Recombination rate 1 (log10))" + r1) +
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("K-", "K+")) +
  scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Model 1", y = "Model 2", fill = "PCA Similarity") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))

ggsave("PCASim_r_model_noZ.png", device = png, width = 7, height = 5)

# Look at just recombination rate - is there any effect?
bootPCASim_sum <- bootPCASim %>%
  group_by(r1, r2) %>%
  dplyr::summarise(meanPCASim = mean(PCASim),
                   ciPCASim = CI(PCASim))

ggplot(bootPCASim_sum, aes(
  x = as.factor(r1), y = as.factor(r2)
)) +
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_r_noZ.png", device = png, width = 7, height = 5)

# Effect seems to be that decreasing recombination in both models
# makes the major axes of variation less similar, but this is a small effect
# Increasing recombination in both models reduces constraint slightly?
# Makes the path to adaptation more similar - finding beneficial combinations?
# If this is the case, then should be more pronounced in K- than K+

# Split into three by model
bootPCASim <- bootPCASim %>%
  mutate(modelCombo = ifelse(model1 != model2, "Mix",
                             as.character(model1)),
         rCombo = ifelse(r1 != r2, 
                         paste(as.character(r1), 
                               as.character(r2), sep = "_"), 
                         as.character(r1)))

# recomb by modelCombo
bootPCASim_sum <- bootPCASim %>%
  group_by(r1, r2, modelCombo) %>%
  dplyr::summarise(meanPCASim = mean(PCASim),
                   ciPCASim = CI(PCASim))


ggplot(bootPCASim_sum, aes(
  x = as.factor(r1), y = as.factor(r2)
)) +
  facet_nested(. ~ "Model comparison" + modelCombo,
               labeller = labeller(modelCombo = as_labeller(c("K" = "K+ vs K+",
                                                              "ODE" = "K- vs K-",
                                                              "Mix" = "K+ vs K-")))) + 
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Recombination rate 1 (log10)", y = "Recombination rate 2 (log10)", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 12), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_r_modelCombo_noZ.png", device = png, width = 7, height = 5)


ggplot(bootPCASim_sum %>% filter(r1 != -5,
                                 r2 != -5), 
       aes(x = as.factor(r1), y = as.factor(r2)
)) +
  facet_nested(. ~ "Model comparison" + modelCombo,
               labeller = labeller(modelCombo = as_labeller(c("K" = "K+ vs K+",
                                                              "ODE" = "K- vs K-",
                                                              "Mix" = "K+ vs K-")))) + 
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim %>% filter(r1 != -5,
                                           r2 != -5), 
              mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("Low", "High")) +
  scale_y_discrete(labels = c("Low", "High")) +
  labs(x = "Recombination rate 1", y = "Recombination rate 2", 
       fill = "PCA Similarity") +
  theme(text = element_text(size = 12), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_r_modelCombo_noZ_pres.png", device = png, width = 10, height = 4)

# Not the case: K+ most affected, but I guess it has the most to change
# K- pretty similar across all recombination rates, does get slightly
# more similar though, K+ is strongly affected - very dissimilar with low
# recombination, CORRELATES WITH MORE VA!!!

# Compare to null
bootPCASim_null_sum <- bootPCASim_null %>% as_tibble() %>%
  dplyr::summarise(meanPCASim = mean(bootPCASim_null),
                   ciPCASim = CI(bootPCASim_null))

ggplot(bootPCASim_null_sum, aes(x = 0, y = 0)) +
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim_null %>% as_tibble(), mapping = aes(fill = bootPCASim_null),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(guide = "none") +
  scale_y_continuous(guide = "none") +
  labs(x = "", y = "", fill = "PCA Similarity") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_nulldist_noZ.png", device = png, width = 5, height = 4)

# model only
bootPCASim_sum <- bootPCASim %>%
  group_by(model1, model2) %>%
  dplyr::summarise(meanPCASim = mean(PCASim),
                   ciPCASim = CI(PCASim))

ggplot(bootPCASim_sum, aes(
  x = model1, y = model2
)) +
  geom_tile(aes(fill = meanPCASim)) +
  theme_bw() +
  geom_jitter(data = bootPCASim, mapping = aes(fill = PCASim),
              shape = 21, size = 1) +
  scale_fill_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = c("K-", "K+")) +
  scale_y_discrete(labels = c("K-", "K+")) +
  labs(x = "Model 1", y = "Model 2", fill = "PCA Similarity") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = 10))
ggsave("PCASim_model_noZ.png", device = png, width = 7, height = 5)


# R doesn't explain much variance, removed
# beta regression for 0-1 data

# Distributions
ggplot(bootPCASim, 
       aes(x = modelCombo, y = PCASim)) +
  geom_quasirandom(dodge.width = 0.9)

# Definitely different between models, none are normally distributed
# beta regression
library(betareg)

# Floating point error: clamp to 1
bootPCASim <- bootPCASim %>%
  mutate(PCASim = raster::clamp(PCASim, 0, 1))

# Run regression: this is slow!
br.pcasim <- betareg(PCASim ~ modelCombo * as.factor(rCombo), bootPCASim)

# Save output
saveRDS(br.pcasim, "betareg_pcaSim_big.RDS")
br.pcasim <- readRDS("betareg_pcaSim_big.RDS")
summary(br.pcasim)
plot(br.pcasim)

# Significance test
library(lmtest)
lrtest(br.pcasim)

car::Anova(br.pcasim, type = 3, test.statistic = "F")

em.pcasim <- emmeans(br.pcasim, ~modelCombo * rCombo)
pairs(em.pcasim, simple = "modelCombo")
pairs(em.pcasim, simple = "rCombo")
plot(em.pcasim, comparisons = T)
pwpp(em.pcasim, by = "modelCombo", type = "response")
emmip(br.pcasim,  ~ modelCombo | rCombo)
joint_tests(em.pcasim)
xtable(em.pcasim)


# ODE models are always the ones with the highest PCA similarity, recombination
# doesn't really affect them mostly
# mix of K- vs K+ has the largest difference
# K+ also very different
# emmip shows the interaction is similar across all recombination pairs
# the largest difference between pairs is when at least one of the matrices has
# low recombination




# So what are the PC contributions in the models then?
# k = number of PCs
PCAContribs <- function(matList, id, k = 2) {
  singleTraitsMap <- c(
    "KXZ",
    "KZ",
    "Z",
    "a",
    "b"
  )
  
  nTraits <- length(singleTraitsMap)
  nG <- length(matList)
  
  PCAdata <- data.frame(
    trait = rep(singleTraitsMap, each = k),
    totalVariation = numeric(nG * nTraits * k),
    pc = rep(1:k, times = nG * nTraits),
    pc_prop = numeric(nG * nTraits * k),
    contrib = numeric(nG * nTraits * k)
  )
  

  # Iterate through matrices for PCA
  for (i in seq_along(matList)) {
    # Run PCA
    g <- matList[[i]]
    pca <- eigen(g)
    
    # Rows to fill for all traits
    i_range <- ( (i-1) * k * nTraits + 1 ):(i * k * nTraits)
    
    # Total variation
    totalVar <- sum(pca$values)
    PCAdata[i_range,]$totalVariation <- totalVar
    
    # Proportion contributed by each PC
    PCAdata[i_range,4] <- rep((pca$values/totalVar)[1:k], 
                                       length.out = length(i_range))
    
    pc_sqr <- pca$vectors^2
    
    pc_contrib <- sweep(pc_sqr, 2, colSums(pc_sqr), FUN="/") * 100
    
    # Contributions
    PCAdata[i_range,5] <- c(t(pc_contrib[,1:2]))

  }
  
  PCAdata <- PCAdata %>%
    mutate(optPerc = rep(id$optPerc, each = k * nTraits),
           seed = rep(id$seed, each = k * nTraits),
           modelindex = rep(id$modelindex, each = k * nTraits),
           clus = rep(id$clus, each = k * nTraits))
  
  return(PCAdata)
}



covpca <- PCAContribs(h2_pd,
  #c(h2_mat_op1, h2_mat_op2, h2_mat_op3, h2_mat_op4), 
  id)
covpca <- AddCombosToDF(covpca)

covpca_sum <- covpca %>%
  group_by(optPerc, r, model, trait, pc) %>%
  summarise_if(is.numeric, list(mean = mean, se = se))

traitLabels <- c(TeX("$K_{XZ}$", output = "character"), 
                 TeX("$K_Z$", output = "character"),
                 TeX("$Z$", output = "character"),
                 TeX("$\\alpha_Z$", output = "character"),
                 TeX("$\\beta_Z$", output = "character")
)

# Plot contributions of each to PC1 and 2

ggplot(covpca_sum %>% filter(trait != "Z", as.numeric(optPerc) == 1) %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = trait, y = contrib_mean, colour = model)) +
  facet_nested(r_title + log10(r) ~ "Principal component" + pc) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = contrib_mean - contrib_se,
                    ymax = contrib_mean + contrib_se),
                width = 0.2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Molecular Trait", y = "Contribution to PC (%)",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12))

# Total variation
ggplot(covpca_sum %>% 
         filter(trait == "Z", pc == 1) %>%
         mutate(r_title = "Recombination rate (log10)"), 
       aes(x = optPerc, y = totalVariation_mean, colour = model)) +
  facet_nested(r_title + log10(r) ~ .) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = totalVariation_mean - totalVariation_se,
                    ymax = totalVariation_mean + totalVariation_se),
                width = 0.2, position = position_dodge(0.9)) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, direction = -1)[2:3],
                      labels = c("K+", "K-"), breaks = c("K", "ODE")) +
  labs(x = "Progress to the optimum", y = "Total variation",
       colour = "Model") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin(-5, 0, 0, 0),
        text = element_text(size = 12))  

