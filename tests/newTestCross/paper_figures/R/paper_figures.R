library(tidyverse)
library(grid)
library(gridExtra)
library(latex2exp)
library(cowplot)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

# Filepath
path_add <- "/mnt/d/SLiMTests/tests/newTestCross/additive/getH2_newTestCross/data/"
path_net <- "/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/"

# Colour palette
cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")


# Pheno time add/net
# d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/addNetCombined/d_com_add+net_after.RDS")
# including moreReps2 data
d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/moreReps2/getH2_newTestCross/data/d_com_add+net_after_prefiltered.RDS")


# Some simulations ended up with populations accumulating massive alpha 
# values/really small beta values which are difficult to overcome
# in this case there is unlikely to be any phenotypic variance, meaning
# no h2 estimate. This was extremely common in the sigma = 0.1, nloci = 1000 case
# to the point where it masked the response of the cases where this didn't happen
d_com %>% filter(phenomean < 10, abs(estR) < 5 | is.na(estR),
                 !(is.na(AIC) & gen >= 50000),
                 aZ < 100 | is.na(aZ), bZ < 100 | is.na(bZ), 
                 KZ < 100 | is.na(KZ), KXZ < 100 | is.na(KXZ)) -> d_com

# Filter out crazy variances
d_com %>% filter((VarA >= 0 & VarA < 10) | is.na(VarA),
       (VarD >= 0 & VarD < 10) | is.na(VarD),
       (VarAA >= 0 & VarAA < 10) | is.na(VarAA),
       (VarR >= 0) | is.na(VarR)) -> d_com

d_qg_sum <- d_com %>% group_by(gen, model, nloci, sigma) %>%
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            ci95Pheno = qnorm(0.975)*sePheno,
            meanDist = mean(dist),
            seDist = se(dist))

pheno_time <- ggplot(d_qg_sum %>% mutate(gen = gen - 50000), 
                     aes(gen, meanPheno, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  scale_fill_manual(values = cc_ibm[c(1,3)], guide = "none") +
  scale_colour_manual(values = cc_ibm[c(1,3)]) +
  geom_ribbon(aes(ymin = (meanPheno - ci95Pheno), ymax = (meanPheno + ci95Pheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Model") +
  theme_bw() +
  theme(text = element_text(size=20), panel.spacing = unit(1, "lines"), legend.position = "bottom")
pheno_time

ggsave("pheno_time.png", pheno_time, width = 10, height = 7, bg = "white")


# GAM model
pred_gam <- readRDS("./pred_gam.RDS")
pred_gam %>% mutate(gen = gen - 50000) %>%
  ggplot(aes(gen, fit, colour = model)) +
  facet_grid(nloci~sigma) +
  scale_fill_manual(values = cc_ibm[c(1, 3)], guide = "none") +
  scale_color_manual(values = cc_ibm[c(1, 3)]) +
  geom_line() +
  #geom_point(data = d_qg_sbst %>% filter(phenomean < 3), aes(y = phenomean), size = 0.1) +
  geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit), 
                  fill = model), color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = TeX("$E\\[\\mu_{\\bar{z}}\\]$"), color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom")

ggsave("pheno_GAM.png", width = 8, height = 6, bg = "white")


# SFS: Calculate the mean allele frequency over time: means per frequency bin
library(gganimate)
library(transformr)
# Split the data set into a list of groups - each is a collection of allele frequencies
d_com %>%
  mutate(Freq_bin = cut(Freq,
                        breaks = seq(from = 0, to = 1, by = 0.05))) %>% # categorise freq
  group_by(gen, seed, model, nloci, sigma) %>%
  mutate(binCount = n()) %>%
  group_by(gen, seed, model, nloci, sigma, Freq_bin) %>%
  mutate(freqBinCount = n()/binCount) %>% # Get the number in each bin - normalise!
  ungroup() -> d_freqs

# Save data for stats
saveRDS(d_freqs, "d_freqs.RDS")

# Summarise for plots
d_freqs %>%
  group_by(gen, model, nloci, sigma, Freq_bin) %>%
  summarise(meanFreqBinCount = mean(freqBinCount),
            seFreqBinCount = se(freqBinCount),
            ci95FreqBinCount = qnorm(0.975)*seFreqBinCount) %>%
  ungroup() %>%
  complete(gen, model, nloci, sigma, Freq_bin, 
           fill = list(meanFreqBinCount = 0)) %>%
  mutate(Freq = as.numeric(Freq_bin) * 0.05, # convert freq back to continuous
         gen = as.integer(gen - 50000)) -> d_freqs 

plt_freqs <- ggplot(d_freqs %>% filter(gen < 2000), 
       aes(x = Freq, y = meanFreqBinCount, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_errorbar(aes(ymin = (meanFreqBinCount - ci95FreqBinCount), ymax = (meanFreqBinCount + ci95FreqBinCount), 
                    color = model)) +
  scale_fill_manual(values = cc_ibm[c(1, 3)], guide = "none") +
  scale_color_manual(values = cc_ibm[c(1, 3)]) +
  scale_y_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency", y = "Mean proportion of mutations", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom")

plt_freqs

# gganimate stuff
plt_freqs <- plt_freqs + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_freqs, nframes = 40, duration = 20, width = 720, height = 720, renderer = ffmpeg_renderer())
anim_save("freq_sfs.mp4", last_animation())

# Plot just the last generation
plt_freqs <- ggplot(d_freqs %>% filter(gen == 1950), 
                    aes(x = Freq, y = meanFreqBinCount, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm[c(1, 3)], guide = "none") +
  scale_color_manual(values = cc_ibm[c(1, 3)]) +
  geom_errorbar(aes(ymin = (meanFreqBinCount - ci95FreqBinCount), ymax = (meanFreqBinCount + ci95FreqBinCount), 
                    color = model)) +
  scale_y_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                           breaks = NULL, labels = NULL)) +
  scale_x_continuous(limits = c(0, 1), sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                           breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency", y = "Mean proportion of mutations", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_freqs
ggsave("sfs_end.png", plt_freqs, width = 7, height = 6, bg = "white")

# Allelic effects

# First add a molTrait value for the NAR models
d_com %>%
  mutate(molTrait = recode_factor(mutType, 
                                  `1`="Neutral", `2`="Del", `3`="$$\\alpha_Z$$", `4`="$$\\beta_Z$$", 
                                  `5`="$$K_Z$$", `6`="$$K_{XZ}$$")) -> d_com

# Lets see the total distribution to see how big these values get
d_com %>% filter(model == "NAR", gen == 51950, Freq == 1) %>%
  ggplot(aes(x = log(value), fill = molTrait)) + 
  geom_histogram(bins = 100)

d_com %>% filter(model == "Additive", gen == 51950) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50)


# summary
d_com %>% filter(gen == 51950) %>%
  mutate(value = case_when(model == "NAR" ~ log(value),
                           model == "Additive" ~ value)) %>%
  group_by(model, molTrait, nloci, sigma) %>%
  summarise(meanVal = mean(value),
            seVal = se(value),
            ci95Val = qnorm(0.975)*seVal) -> d_meanStats
  
# export table to latex
library(stargazer)
stargazer(as.data.frame(d_meanStats), summary = FALSE, rownames = F)

# Split the data set into a list of groups - each is a collection of allelic effects
MIN_VAL <- min(log(d_com[d_com$model == "NAR",]$value))
MAX_VAL <- max(log(d_com[d_com$model == "NAR",]$value))

d_com %>% filter(model == "NAR") %>%
  mutate(value = log(value)) %>%
  mutate(val_bin = cut(value,
                        breaks = seq(from = MIN_VAL, to = MAX_VAL, by = 0.05))) %>% # categorise freq
  group_by(gen, seed, nloci, sigma) %>%
  mutate(binCount = n()) %>%
  group_by(gen, seed, molTrait, nloci, sigma, val_bin) %>%
  mutate(valBinCount = n()/binCount) %>% # Get the number in each bin - normalise!
  ungroup(seed) -> d_values_NAR

# Save data for stats
saveRDS(d_values_NAR, "d_values_NAR.RDS")

d_values_NAR <- readRDS("d_values_NAR.RDS")

# Summarise for plots
MIN_VAL <- min(d_values_NAR$value)
MAX_VAL <- max(d_values_NAR$value)

VALUE_BIN_SEQ <- seq(from = MIN_VAL, to = MAX_VAL, by = 0.05)

d_values_NAR %>%
  group_by(gen, molTrait, nloci, sigma, val_bin) %>%
  summarise(meanValBinCount = mean(valBinCount),
            seValBinCount = se(valBinCount),
            ci95ValBinCount = qnorm(0.975)*seValBinCount,
            value = mean(VALUE_BIN_SEQ[as.numeric(val_bin)])) %>%
  ungroup() %>%
  mutate(gen = as.integer(gen - 50000)) -> d_values_NAR

mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$"),
  TeX("$K_Z$"),
  TeX("$K_{XZ}$")
)


plt_val_NAR <- ggplot(d_values_NAR %>% filter(gen == 1950), 
                    aes(x = value, y = meanValBinCount, color = molTrait)) +
  facet_grid(nloci~sigma, scales = "free") +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanValBinCount - ci95ValBinCount), ymax = (meanValBinCount + ci95ValBinCount), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                           breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                                           breaks = NULL, labels = NULL)) +
  labs(x = "Allelic effect", y = "Mean proportion of mutations", color = "Molecular trait") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))

plt_val_NAR

# Additive

MIN_VAL <- min(d_com[d_com$model == "Additive",]$value)
MAX_VAL <- max(d_com[d_com$model == "Additive",]$value)

d_com %>% filter(model == "Additive") %>%
  mutate(val_bin = cut(value,
                       breaks = seq(from = MIN_VAL, to = MAX_VAL, by = 0.05))) %>% # categorise freq
  group_by(gen, seed, nloci, sigma) %>%
  mutate(binCount = n()) %>%
  group_by(gen, seed, nloci, sigma, val_bin) %>%
  mutate(valBinCount = n()/binCount) %>% # Get the number in each bin - normalise!
  ungroup(seed) -> d_values_add

# Save data for stats
saveRDS(d_values_add, "d_values_add.RDS")

d_values_add <- readRDS("d_values_add.RDS")

# Summarise for plots
MIN_VAL <- min(d_values_add$value)
MAX_VAL <- max(d_values_add$value)

VALUE_BIN_SEQ <- seq(from = MIN_VAL, to = MAX_VAL, by = 0.05)

d_values_add %>%
  group_by(gen, nloci, sigma, val_bin) %>%
  summarise(meanValBinCount = mean(valBinCount),
            seValBinCount = se(valBinCount),
            ci95ValBinCount = qnorm(0.975)*se(valBinCount),
            value = mean(VALUE_BIN_SEQ[as.numeric(val_bin)])) %>%
  ungroup() %>%
  mutate(gen = as.integer(gen - 50000)) -> d_values_add

plt_val_add <- ggplot(d_values_add %>% filter(gen == 1950), 
                      aes(x = value, y = meanValBinCount)) +
  facet_grid(nloci~sigma, scales = "free") +
  geom_line() +
  scale_color_manual(values = cc_ibm) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanValBinCount - seValBinCount), ymax = (meanValBinCount + seValBinCount)
                  ), color = NA, alpha = 0.2) +
  scale_y_continuous(limits = c(0, 0.35), sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Allelic effect", y = "Mean proportion of mutations") +
  theme_bw() +
  theme(text=element_text(size=16), panel.spacing.x = unit(1, "lines"))

plt_val_add

plt_val <- plot_grid(plt_val_add, plt_val_NAR, 
                     ncol = 1, rel_heights = c(1, 1.2),
                     labels = "AUTO")

ggsave("val_end.png", plt_val, width = 7, height = 10, bg = "white")

# Mean effect and freq over time
d_com %>%
  filter(model == "NAR") %>%
  mutate(gen = gen - 50000,
         value = log(value)) %>%
  group_by(gen, molTrait, nloci, sigma) %>%
  summarise(meanValue = mean(value),
            seValue = se(value),
            ci95Value = qnorm(0.975)*seValue,
            meanFreq = mean(Freq),
            seFreq = se(Freq),
            ci95Freq = qnorm(0.975)*seFreq) -> d_meanfreq_NAR

# Summary
View(d_meanfreq_NAR %>% filter(gen == 0 | gen == 100, nloci == 10, sigma == 0.1))


plt_meanVal_NAR_s01 <- ggplot(d_meanfreq_NAR %>% filter(sigma == 0.1), 
                      aes(x = gen, y = meanValue, color = molTrait)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanValue - ci95Value), ymax = (meanValue + ci95Value), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  labs(x = NULL, y = "Mean mutational effect", color = "Molecular trait") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal_NAR_s01

leg_meanVal <- get_legend(plt_meanVal_NAR_s01)


plt_meanVal_NAR_s1 <- ggplot(d_meanfreq_NAR %>% filter(sigma == 1), 
                              aes(x = gen, y = meanValue, color = molTrait)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanValue - ci95Value), ymax = (meanValue + ci95Value), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  labs(x = NULL, y = "Mean mutational effect", color = "Molecular trait") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal_NAR_s1

plt_meanVal_NAR <- plot_grid(plt_meanVal_NAR_s01 +
            theme(strip.background.y = element_blank(),
                  strip.text.y = element_blank(),
                  plot.margin = margin(5.5, 12, 5.5, 5.5),
                  legend.position = "none"), 
          plt_meanVal_NAR_s1 + scale_y_continuous(
            sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                breaks = NULL, labels = NULL)) +
            labs(y = "") +
            theme(legend.position = "none"),
          ncol = 2,
          rel_widths = c(1.1, 1))

top.grob <- textGrob("Mutational effect size variance", 
                     gp = gpar(fontsize = 16))

x.grob <- textGrob("Generations after optimum shift", 
                   gp = gpar(fontsize = 16))


g <- arrangeGrob(plt_meanVal_NAR, bottom = x.grob, top = top.grob)
plt_meanVal_NAR <- plot_grid(g, leg_meanVal, ncol = 1, rel_heights = c(1, 0.1))

# Additive
d_com %>%
  filter(model == "Additive") %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, nloci, sigma) %>%
  summarise(meanValue = mean(value),
            seValue = se(value),
            ci95Value = qnorm(0.975)*seValue,
            meanFreq = mean(Freq),
            seFreq = se(Freq),
            ci95Freq = qnorm(0.975)*seFreq) -> d_meanfreq_add

plt_meanVal_add <- ggplot(d_meanfreq_add, 
                              aes(x = gen, y = meanValue)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_ribbon(aes(ymin = (meanValue - ci95Value), ymax = (meanValue + ci95Value)), 
              color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                                              breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean mutational effect") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal_add

plt_meanVal <- plot_grid(plt_meanVal_add,
                         plt_meanVal_NAR,
                         ncol = 1, 
                         rel_heights = c(0.9, 1),
                         labels = "AUTO")
plt_meanVal

ggsave("meanVal.png", plt_meanVal, width = 10, height = 10, bg = "white")

# Mean allele frequency
plt_meanFreq_NAR <- ggplot(d_meanfreq_NAR, 
                              aes(x = gen, y = meanFreq, color = molTrait)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanFreq - ci95Freq), ymax = (meanFreq + ci95Freq), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  scale_y_continuous(limits = c(-0.2398877, 1.0654877), 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean frequency", color = "Molecular trait") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanFreq_NAR

plt_meanFreq_add <- ggplot(d_meanfreq_add, 
                           aes(x = gen, y = meanFreq)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_ribbon(aes(ymin = (meanFreq - ci95Freq), ymax = (meanFreq + ci95Freq)), 
              color = NA, alpha = 0.2) +
  scale_y_continuous(limits = c(-0.2398877, 1.0654877), sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean frequency") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanFreq_add

plt_meanFreq <- plot_grid(plt_meanFreq_add,
                          plt_meanFreq_NAR,
                         ncol = 1, 
                         rel_heights = c(0.9, 1),
                         labels = "AUTO")
plt_meanFreq

ggsave("meanFreq.png", plt_meanFreq, width = 7, height = 10, bg = "white")

# Heritability and predictions
###################################################
d_qg_sum <- d_com %>% group_by(gen, model, nloci, sigma) %>%
  drop_na(AIC) %>%
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  mutate(gen = gen - 50000) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist),
            meanH2A = mean(H2.A.Estimate),
            seH2A = se(H2.A.Estimate),
            ci95H2A = qnorm(0.975)*seH2A,
            meanH2D = mean(H2.D.Estimate), 
            seH2D = se(H2.D.Estimate),
            ci95H2D = qnorm(0.975)*seH2D,
            meanH2AA = mean(H2.AA.Estimate),
            seH2AA = se(H2.AA.Estimate),
            ci95H2AA = qnorm(0.975)*seH2AA) %>%
  pivot_longer(cols = c(meanH2A, meanH2D, meanH2AA),
               names_to = "varComp", values_to = "prop") %>%
  mutate(varComp = factor(varComp, c("meanH2A", "meanH2D", "meanH2AA")),
         propCI95 = case_when(varComp == "meanH2A" ~ ci95H2A,
                              varComp == "meanH2D" ~ ci95H2D,
                              varComp == "meanH2AA" ~ ci95H2AA))

# stacked percent bar chart
plt_h2_add <- ggplot(d_qg_sum %>% filter(model == "Additive"),
                        aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Epistatic (AxA)")) +
  scale_color_viridis_d(guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~sigma) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom", 
        strip.background.y = element_blank(),
        strip.text.y = element_blank())

plt_h2_NAR <- ggplot(d_qg_sum %>% filter(model == "NAR"),
                     aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Epistatic (AxA)")) +
  scale_y_continuous(labels = scales::percent, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~sigma) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(plt_h2_add)

plt_h2 <- plot_grid(plt_h2_add + theme(legend.position = "none") + labs(x = NULL),
                     plt_h2_NAR + theme(legend.position = "none") + labs(y = NULL, x = NULL),
                     ncol = 2,
                     labels = "AUTO")

x.grob <- textGrob("Generations after optimum shift", 
                   gp = gpar(fontsize = 16))

g <- arrangeGrob(plt_h2, bottom = x.grob)

plt_h2 <- plot_grid(g, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_h2

# Predicted responses
d_resp <- d_com %>%
  distinct(gen, seed, modelindex, .keep_all = T) %>%
  group_by(gen, model, nloci, sigma) %>%
  mutate(expPheno = lag(phenomean, default = 0) + lag(estR, default = 0)) %>%
  ungroup() %>%
  group_by(gen, model, nloci, sigma) %>%
  summarise(meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            ci95Pheno = qnorm(0.975)*sePheno,
            seExpPheno = se(expPheno),
            ci95ExpPheno = qnorm(0.975)*seExpPheno,
            expPheno = mean(expPheno)
  ) %>%
  mutate(expPheno = case_when(gen == 49500 ~ meanPheno,
                              .default = expPheno)) %>%
  ungroup() %>%
  group_by(model, nloci, sigma) %>%
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000,
         sePhenoValue = case_when(phenoType == "Estimated" ~ seExpPheno,
                                  phenoType == "Observed" ~ sePheno),
         ci95PhenoValue = case_when(phenoType == "Estimated" ~ ci95ExpPheno,
                                  phenoType == "Observed" ~ ci95Pheno)) %>%
  ungroup()

plt_resp_add <- ggplot(d_resp %>% filter(model == "Additive"), 
                   aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~sigma) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  scale_fill_manual(values = c(cc_ibm[1], cc_ibm[4]), guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = phenoValue - ci95PhenoValue, ymax = phenoValue + ci95PhenoValue,
                  fill = phenoType), color = NA, alpha = 0.2) +
  coord_cartesian(ylim = c(-0.1, 3.5)) +
  scale_y_continuous() +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(5.5, 12, 5.5, 5.5))

plt_resp_NAR <- ggplot(d_resp %>% filter(model == "NAR"), 
                       aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~sigma) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  scale_fill_manual(values = c(cc_ibm[1], cc_ibm[4]), guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = phenoValue - ci95PhenoValue, ymax = phenoValue + ci95PhenoValue,
              fill = phenoType), colour = NA, alpha = 0.2) +
  coord_cartesian(ylim = c(-0.1, 3.5)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

h2_leg <- get_legend(plt_resp_NAR)

plt_resp <- plot_grid(plt_resp_add + theme(legend.position = "none") + labs(x = NULL),
                     plt_resp_NAR + theme(legend.position = "none") + labs(y = NULL, x = NULL),
                     ncol = 2,
                     labels = c("C", "D"))

x.grob <- textGrob("Generations after optimum shift", 
                   gp = gpar(fontsize = 16))

g <- arrangeGrob(plt_resp, bottom = x.grob)

plt_resp <- plot_grid(g, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_resp

plt_h2_resp <- plot_grid(plt_h2,
                         plt_resp,
                         nrow = 2)

plt_h2_resp
ggsave("h2_resp.png", plt_h2_resp, width = 16, height = 16, bg = "white")

# Allele age?

# https://stackoverflow.com/a/24241954
scientific_label <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00", "0", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  
  # return this as an expression
  parse(text=l)
}


d_com %>% 
  filter(model == "NAR") %>%
  group_by(gen, molTrait, nloci, sigma) %>%
  mutate(alleleAge = gen - originGen,
         gen = gen - 50000) %>%
  summarise(meanAge = mean(alleleAge),
            ci95Age = qnorm(0.975) * se(alleleAge)) -> d_age_NAR

plt_age_NAR <- ggplot(d_age_NAR, aes(x = gen, y = meanAge, colour = molTrait)) +
  facet_grid(nloci~sigma) +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = meanAge - ci95Age, ymax = meanAge + ci95Age,
                  fill = molTrait), colour = NA, alpha = 0.2) +
  scale_y_continuous(limits = c(-9609.96, 46877.33), labels = scientific_label,
    sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean allele age", colour = "Molecular trait") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))
  

d_com %>% 
  filter(model == "Additive") %>%
  group_by(gen, nloci, sigma) %>%
  mutate(alleleAge = gen - originGen,
         gen = gen - 50000) %>%
  summarise(meanAge = mean(alleleAge),
            ci95Age = qnorm(0.975) * se(alleleAge)) -> d_age_add

plt_age_add <- ggplot(d_age_add, aes(x = gen, y = meanAge)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_ribbon(aes(ymin = meanAge - ci95Age, ymax = meanAge + ci95Age), 
              alpha = 0.2) +
  scale_y_continuous(limits = c(-9609.96, 46877.33), labels = scientific_label,
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations after optimum shift", y = "Mean allele age") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))
plt_age_add


plt_age <- plot_grid(plt_age_add, plt_age_NAR, ncol = 1, labels = "AUTO")
plt_age

ggsave("meanAlleleAge.png", plt_age, width = 8, height = 10)


# Eigenvector of molecular traits
library(ggfortify)

d_com %>%
  filter(model == "NAR", gen == 51950, nloci == 100, sigma == 1) -> d_test

d_test %>%
  select(aZ, bZ, KZ, KXZ) %>%
  prcomp(., scale = T) -> pca_test

d_test %>%
  pivot_longer(c(aZ, bZ, KZ, KXZ), names_to = "moltrait_name", values_to = "moltrait_value") %>%
  autoplot(pca_test, data = ., colour = 'moltrait_name') + stat_ellipse()

# Eigenvector of pheno/moltrait/freq/value


# Data for groups that have/haven't adapted
path_both <- "/mnt/d/SLiMTests/tests/newTestCross/moreReps2/getH2_newTestCross/data/"

d_com_adapted_eg <- readRDS(paste0(path_both, "d_com_prefiltered_adapted_eg.RDS"))
d_com_maladapted_eg <- readRDS(paste0(path_both, "d_com_prefiltered_maladapted_eg.RDS"))
d_com_wasadapted_eg <- readRDS(paste0(path_both, "d_com_prefiltered_wasadapted_eg.RDS"))

d_com_adapted_eg %>% 
  filter(model == "NAR", nloci == 10, sigma == 0.1) %>%
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  View(.)
  mutate(outlier = phenomean > outlier_pheno_upper | phenomean < outlier_pheno_lower) %>% 
  select(outlier)

d_com_adapted_eg %>%
  filter(gen >= 49500) %>%
  #mutate(seed = as.character(seed)) %>%
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  ggplot(aes(x = gen, y = phenomean, color = nloci, group = seed)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 1.9, ymax = 2.1, alpha = 0.2) +
  facet_grid(model~sigma) +
  #coord_cartesian(ylim = c(0, 2.5)) +
  geom_line() +
  scale_color_manual(values = cc_ibm) +
  geom_point(size = 0.25) +
  labs(x = "Generations after optimum shift", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

d_com_maladapted_eg %>%
  filter(gen >= 49500) %>%
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  ggplot(aes(x = gen, y = phenomean, color = nloci, group = seed)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 1.9, ymax = 2.1, alpha = 0.2) +
  facet_grid(model~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm) +
  coord_cartesian(ylim = c(0, 2.5)) +
  geom_point(size = 0.25) +
  labs(x = "Generations after optimum shift", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

d_com_wasadapted_eg %>%
  filter(gen >= 49500) %>%
  mutate(seed = as.character(seed)) %>%
  distinct(gen, seed, model, nloci, sigma, .keep_all = T) %>%
  ggplot(aes(x = gen, y = phenomean, color = nloci, group = seed)) +
  facet_grid(model~sigma) +
  geom_line() +
  geom_point(size = 0.25) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 1.9, ymax = 2.1, alpha = 0.2) +
  scale_color_manual(values = cc_ibm) +
  labs(x = "Generations after optimum shift", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16), 
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))

View(d_com_adapted_eg %>% filter (gen >= 49500, model == "NAR", sigma == 0.1) %>% 
       distinct(gen, seed, model, nloci, sigma, .keep_all = T))
