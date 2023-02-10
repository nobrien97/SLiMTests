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
d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/addNetCombined/d_com_add+net_after.RDS")
d_com %>% filter(phenomean < 10) -> d_com

d_qg_sum <- d_com %>% group_by(gen, model, nloci, sigma) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
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
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Fixed effect") +
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
d_com <- readRDS("/mnt/d/SLiMTests/tests/newTestCross/addNetCombined/d_com_add+net_after.RDS")

d_com %>% filter(phenomean < 10) -> d_com
# Split the data set into a list of groups - each is a collection of allele frequencies
d_com %>%
  mutate(Freq_bin = cut(Freq,
                        breaks = seq(from = 0, to = 1, by = 0.05))) %>% # categorise freq
  group_by(gen, seed, model, nloci, sigma) %>%
  mutate(binCount = n()) %>%
  group_by(gen, seed, model, nloci, sigma, Freq_bin) %>%
  mutate(freqBinCount = n()/binCount) %>% # Get the number in each bin - normalise!
  ungroup(seed) -> d_freqs

# Save data for stats
saveRDS(d_freqs, "d_freqs.RDS")

# Summarise for plots
d_freqs %>%
  group_by(gen, model, nloci, sigma, Freq_bin) %>%
  summarise(meanFreqBinCount = mean(freqBinCount),
            seFreqBinCount = se(freqBinCount)) %>%
  ungroup() %>%
  mutate(Freq = as.numeric(Freq_bin) * 0.05, # convert freq back to continuous
         gen = as.integer(gen - 50000)) -> d_freqs 

plt_freqs <- ggplot(d_freqs %>% filter(gen < 2000), 
       aes(x = Freq, y = meanFreqBinCount, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
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

animate(plt_freqs, nframes = 40, duration = 20, width = 720, height = 720, renderer = av_renderer())
anim_save("freq_sfs.mp4", last_animation())

# Plot just the last generation
plt_freqs <- ggplot(d_freqs %>% filter(gen == 1950), 
                    aes(x = Freq, y = meanFreqBinCount, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm[c(1, 3)], guide = "none") +
  scale_color_manual(values = cc_ibm[c(1, 3)]) +
  geom_ribbon(aes(ymin = (meanFreqBinCount - seFreqBinCount), ymax = (meanFreqBinCount + seFreqBinCount), 
                  fill = model), color = NA, alpha = 0.2) +
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
  ggplot(aes(x = value)) + # Need to log twice because the value is exponentiated twice in this dataset... whoops
  geom_histogram(bins = 50)

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
  geom_ribbon(aes(ymin = (meanValBinCount - seValBinCount), ymax = (meanValBinCount + seValBinCount), 
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
library(ggh4x)
d_com %>%
  filter(model == "NAR") %>%
  mutate(gen = gen - 50000,
         value = log(value)) %>%
  group_by(gen, molTrait, nloci, sigma) %>%
  summarise(meanValue = mean(value),
            seValue = se(value),
            meanFreq = mean(Freq),
            seFreq = se(Freq)) -> d_meanfreq_NAR


plt_meanVal_NAR_s01 <- ggplot(d_meanfreq_NAR %>% filter(sigma == 0.1), 
                      aes(x = gen, y = meanValue, color = molTrait)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanValue - seValue), ymax = (meanValue + seValue), 
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
  geom_ribbon(aes(ymin = (meanValue - seValue), ymax = (meanValue + seValue), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  labs(x = NULL, y = "Mean mutational effect", color = "Molecular trait") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_meanVal_NAR_s1

plt_meanVal_NAR <- plot_grid(plt_meanVal_NAR_s01 +
            theme(strip.background.y = element_blank(),
                  strip.text.y = element_blank(),
                  legend.position = "none"), 
          plt_meanVal_NAR_s1 + scale_y_continuous(
            sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                breaks = NULL, labels = NULL)) +
            labs(y = "") +
            theme(legend.position = "none"),
          ncol = 2)

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
            meanFreq = mean(Freq),
            seFreq = se(Freq)) -> d_meanfreq_add

plt_meanVal_add <- ggplot(d_meanfreq_add, 
                              aes(x = gen, y = meanValue)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_ribbon(aes(ymin = (meanValue - seValue), ymax = (meanValue + seValue)), 
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

ggsave("meanVal.png", plt_meanVal, width = 10, height = 8, bg = "white")

# Mean allele frequency
plt_meanFreq_NAR <- ggplot(d_meanfreq_NAR, 
                              aes(x = gen, y = meanFreq, color = molTrait)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_color_manual(values = cc_ibm, labels = mutType_names) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  geom_ribbon(aes(ymin = (meanFreq - seFreq), ymax = (meanFreq + seFreq), 
                  fill = molTrait), color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
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
  geom_ribbon(aes(ymin = (meanFreq - seFreq), ymax = (meanFreq + seFreq)), 
              color = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
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
