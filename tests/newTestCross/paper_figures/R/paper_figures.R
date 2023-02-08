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
pheno_time_add <- readRDS(paste0(path_add, "pheno_time_add.RDS"))
pheno_time_net <- readRDS(paste0(path_net, "pheno_time_unfixed.RDS"))

pheno_time_leg <- get_legend(pheno_time_add)

pheno_time <- plot_grid(
  pheno_time_add + theme(legend.position = "none", panel.spacing = unit(1, "lines")),
  pheno_time_net + theme(legend.position = "none", panel.spacing = unit(1, "lines")),
  nrow = 2, labels = "AUTO"
)

pheno_time <- plot_grid(pheno_time, pheno_time_leg, nrow = 2, rel_heights = c(3, .4))
pheno_time

ggsave("pheno_time.png", pheno_time, width = 8, height = 10, bg = "white")


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
  labs(x = "Allele frequency", y = "Normalized count", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom")

plt_freqs

# gganimate stuff
plt_freqs <- plt_freqs + transition_states(gen) +
  labs(title = "Generation: {closest_state}")

animate(plt_freqs, nframes = 40, duration = 40, width = 720, height = 720, renderer = av_renderer())
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
  labs(x = "Allele frequency", y = "Normalized count", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom", panel.spacing.x = unit(1, "lines"))
plt_freqs
ggsave("sfs_end.png", plt_freqs, width = 7, height = 6, bg = "white")

