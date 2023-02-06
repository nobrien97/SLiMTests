# Plot heritability estimates
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
path <- "/mnt/d/SLiMTests/tests/newTestCross/additive/getH2_newTestCross/data/"

# Colour palette
cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")

d_h2 <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/additive/getH2_newTestCross/data/out_h2.csv", col_names = F)

names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
                 "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                 "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                 "H2.AA.SE", "H2.R.Estimate", "H2.R.SE", "AIC")

# Remove duplicates
d_h2 %>% distinct() -> d_h2

# Remove NAs
d_h2 %>% drop_na() -> d_h2

# Add our variables
combos <- read_delim("../../R/combos.csv", delim = " ", col_names = F)
names(combos) <- c("nloci", "locisigma")
d_h2 %>% mutate(nloci = combos$nloci[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_h2

# Remove weird values: negative variance, huge variance
d_h2 %>% filter(VarA >= 0, VarA < 10,
                VarD >= 0, VarAA >= 0,
                VarR > 0, sum(VarA, VarD, VarAA, VarR) > 0) -> d_h2

d_h2$nloci <- as.factor(d_h2$nloci)
d_h2$sigma <- as.factor(d_h2$sigma)


d_h2_sum <- d_h2 %>%
  group_by(gen, nloci, sigma) %>% 
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
ggplot(d_h2, aes(x = H2.A.Estimate)) +
  geom_density(alpha = 0.4)

ggplot(d_h2, aes(x = H2.AA.Estimate)) +
  geom_density(alpha = 0.4)

ggplot(d_h2, aes(x = H2.D.Estimate)) +
  geom_density(alpha = 0.4)


# Mean over time
h2A_time <- ggplot(d_h2_sum, aes(x = gen, y = meanH2A, color = nloci)) +
  facet_grid(.~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm, guide = guide_legend(override.aes = list(linewidth = 3))) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2A - seH2A), ymax = (meanH2A + seH2A), 
                  fill = nloci), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean heritability", color = "Number of loci") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "bottom")

# What about other variance components?
h2D_time <- ggplot(d_h2_sum, aes(x = gen, y = meanH2D, color = nloci)) +
  facet_grid(.~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2D - seH2D), ymax = (meanH2D + seH2D), 
                  fill = nloci), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average dominance variance proportion", color = "Number of\nloci") +
  theme_bw() +
  theme(text=element_text(size=16))


h2AA_time <- ggplot(d_h2_sum, aes(x = gen, y = meanH2AA, color = nloci)) +
  facet_grid(.~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2AA - seH2AA), ymax = (meanH2AA + seH2AA), 
                  fill = nloci), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average AxA variance proportion", color = "Number of\nloci") +
  theme_bw() +
  theme(text=element_text(size=16))

plt_row <- plot_grid(
  h2A_time + theme(legend.position = "none", panel.spacing = unit(1, "lines")) + 
    xlab(NULL) + draw_label(TeX("$h^2_{A}$"), x = 1800, y = 0.975),
  h2D_time + theme(legend.position = "none", panel.spacing = unit(1, "lines")) + 
    ylab(NULL) + xlab(NULL) +
    draw_label(TeX("$h^2_{D}$"), x = 1800, y = 0.975),
  h2AA_time + theme(legend.position = "none", panel.spacing = unit(1, "lines")) + 
    ylab(NULL) + xlab(NULL) +
    draw_label(TeX("$h^2_{AA}$"), x = 1800, y = 0.975),
  nrow = 1
)

x.grob <- textGrob("Generations after optimum shift", 
                   gp = gpar(fontsize = 16))

top.grob <- textGrob("Mutational effect size variance", 
                   gp = gpar(fontsize = 16))

g <- arrangeGrob(plt_row, bottom = x.grob, top = top.grob)

h2_leg <- get_legend(h2A_time)

plt_row_leg <- plot_grid(g, h2_leg, nrow = 2, rel_heights = c(3, .4))
plt_row_leg

ggsave("h2_time.png", plt_row_leg,  height = 4, width = 18, bg = "white")

# Now averaged over time
d_h2_sum_notime <- d_h2 %>%
  group_by(nloci, sigma) %>% 
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


h2_avg <- ggplot(d_h2_sum_notime %>% pivot_longer(
  cols = c(meanH2A, meanH2D, meanH2AA),
  names_to = "varComp", values_to = "prop") %>%
    mutate(propSE = case_when(varComp == "meanH2A" ~ seH2A,
                              varComp == "meanH2D" ~ seH2D,
                              varComp == "meanH2AA" ~ seH2AA),
           varComp = factor(varComp, c("meanH2A", "meanH2D", "meanH2AA"))), 
  aes(x = nloci, y = prop, fill = varComp)) +
  facet_grid(.~sigma) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = prop - propSE, ymax = prop + propSE), 
                position = position_dodge(0.9), width = 0.3) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Epistatic (AxA)")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Number of loci", 
         y = "Proportion of total\nphenotypic variance (%)",
         fill = "Variance component") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2_avg)

top.grob <- textGrob("Mutational effect size variance", 
                     gp = gpar(fontsize = 16))

g <- arrangeGrob(h2_avg + theme(legend.position = "none"), top = top.grob)


plt_row_leg <- plot_grid(g, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_row_leg

ggsave("h2_avg.png", plt_row_leg,  height = 6, width = 8, bg = "white")

###################################################
# stacked percent bar chart
h2_stk_var <- ggplot(d_h2_sum %>% 
                          pivot_longer(cols = c(meanH2A, meanH2D, meanH2AA),
                                       names_to = "varComp", values_to = "prop") %>%
                       mutate(varComp = factor(varComp, c("meanH2A", "meanH2D", "meanH2AA"))), 
  aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Epistatic (AxA)")) +
  scale_y_continuous(labels = scales::percent, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect size variance", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~sigma) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

h2_stk_var

ggsave("h2_stacked.png", h2_stk_var, width = 10, height = 8, bg = "white")


# Plot trait data over time
d_qg <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/additive/slim_qg.csv", col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")


d_qg %>% filter(phenomean < 10) -> d_qg

d_qg %>% mutate(nloci = combos$nloci[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_qg

d_qg$nloci <- as.factor(d_qg$nloci)
d_qg$sigma <- as.factor(d_qg$sigma)


d_qg_sum <- d_qg %>% group_by(gen, nloci, sigma) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

pheno_time <- ggplot(d_qg_sum %>% filter(gen > 49000) %>% mutate(gen = gen - 50000), 
       aes(gen, meanPheno, color = nloci)) +
  facet_grid(.~sigma) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm, guide = guide_legend(override.aes = list(linewidth = 3))) +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = nloci), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Number of loci") +
  theme_bw() +
  theme(text = element_text(size=20), panel.spacing = unit(2, "lines"))
pheno_time
ggsave("pheno_time.png", pheno_time, height = 5, width = 10)

# combined data: SFS
#TODO: There's too much information to unpack here, need a metric to describe the
# SFS, maybe just mean frequency? And similar deviation from model w/o fixed effect
d_com <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/additive/d_combined_after.csv")
d_com$fixTime <- d_com$fixGen - d_com$originGen

d_com_end <- d_com %>% filter(gen > 51000)
d_com_end %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                     h2D = cut(H2.D.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                     h2AA = cut(H2.AA.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                     fixTime = fixGen - originGen) -> d_com_end


d_com_end %>% mutate(nloci = combos$nloci[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_com_end

d_com_end$nloci <- as.factor(d_com_end$nloci)
d_com_end$sigma <- as.factor(d_com_end$sigma)
rm(d_com)

d_com_sum <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, h2A) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq)
            ) %>%
  complete(h2A, fill = list(meanFreq = 0,
                             seFreq = NA))

# Plot mean frequency for each group: averaged across sigma
h2A_freq <- ggplot(d_com_sum, 
       aes(x = h2A, y = meanFreq, fill = nloci)) +
  facet_grid(.~sigma) +
  scale_fill_manual(values = cc_ibm) +
  ylim(0, 1) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanFreq - seFreq, ymax = meanFreq + seFreq), 
                 position = position_dodge(0.9, preserve = "single"), width = 0.3) +
  labs(x = TeX("Heritability ($h^2_A$)"), y = "Mean frequency", fill = "Number of loci") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")
h2A_freq

h2_leg <- get_legend(h2A_freq)

top.grob <- ggdraw() + draw_label("Mutational effect size variance", size = 16) +
  theme(plot.margin = margin(0, 0, 0, 7))

g_h2A <- arrangeGrob(h2A_freq + theme(legend.position = "none"))


d_com_sum <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, h2AA) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq)
  ) %>%
  complete(h2AA, fill = list(meanFreq = 0,
                             seFreq = NA))

h2AA_freq <- ggplot(d_com_sum, 
                   aes(x = h2AA, y = meanFreq, fill = nloci)) +
  facet_grid(.~sigma) +
  scale_fill_manual(values = cc_ibm) +
  ylim(0, 1) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanFreq - seFreq, ymax = meanFreq + seFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3) +
  labs(x = TeX("Heritability ($h^2_{AA}$)"), y = "Mean frequency", fill = "Number of\nloci") +
  theme_bw() + theme(text = element_text(size = 16))
h2AA_freq

g_h2AA <- arrangeGrob(h2AA_freq + theme(legend.position = "none"))


d_com_sum <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, h2D) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq)
  ) %>%
  complete(h2D, fill = list(meanFreq = 0,
                             seFreq = NA))

h2D_freq <- ggplot(d_com_sum, 
                    aes(x = h2D, y = meanFreq, fill = nloci)) +
  facet_grid(.~sigma) +  
  ylim(0, 1) +
  scale_fill_manual(values = cc_ibm) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanFreq - seFreq, ymax = meanFreq + seFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3) +
  labs(x = TeX("Heritability ($h^2_D$)"), y = "Mean frequency", fill = "Number of\nloci") +
  theme_bw() + theme(text = element_text(size = 16))
h2D_freq

g_h2D <- arrangeGrob(h2D_freq + theme(legend.position = "none"))


plt_freq <- plot_grid(g_h2A, g_h2D, g_h2AA, ncol = 3, labels = "AUTO")
plt_freq <- plot_grid(top.grob, plt_freq, h2_leg, ncol = 1, rel_heights = c(.1, 1, .1))
plt_freq

ggsave("h2_freq.png", plt_freq,  height = 5, width = 24, bg = "white")

# Effect size
# Plot each heritability
y_lims <- c(-0.05, 0.6)

# summarise
d_val_sum_notime <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, h2A) %>%
  summarise(meanVal = mean(value),
            seVal = se(value)
  ) %>%
  complete(h2A, fill = list(meanVal = 0,
                            meanVal = NA))

h2A_val <- ggplot(d_val_sum_notime,
                   aes(x = h2A, y = meanVal, fill = nloci)) +
  facet_grid(.~sigma, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_fill_manual(values = cc_ibm) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanVal - seVal, ymax = meanVal + seVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_A$)"), y = "Mean molecular effect", fill = "Number of loci") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2A_val)

g_h2A <- arrangeGrob(h2A_val + theme(legend.position = "none"))


d_val_sum_notime <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, h2D) %>%
  summarise(meanVal = mean(value),
            seVal = se(value)
  ) %>%
  complete(h2D, fill = list(meanVal = 0,
                            meanVal = NA))


h2D_val <- ggplot(d_val_sum_notime, 
                  aes(x = h2D, y = meanVal, fill = nloci)) +
  facet_grid(.~sigma, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_fill_manual(values = cc_ibm) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanVal - seVal, ymax = meanVal + seVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_D$)"), y = "Mean molecular effect", fill = "Number of loci") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2D <- arrangeGrob(h2D_val + theme(legend.position = "none"))

d_val_sum_notime <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(nloci, sigma, h2AA) %>%
  summarise(meanVal = mean(value),
            seVal = se(value)
  ) %>%
  complete(h2AA, fill = list(meanVal = 0,
                            meanVal = NA))

h2AA_val <- ggplot(d_val_sum_notime, 
                  aes(x = h2AA, y = meanVal, fill = nloci)) +
  facet_grid(.~sigma, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_fill_manual(values = cc_ibm) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanVal - seVal, ymax = meanVal + seVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_{AA}$)"), y = "Mean molecular effect", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2AA <- arrangeGrob(h2AA_val + theme(legend.position = "none"))


plt_val <- plot_grid(g_h2A, g_h2D, g_h2AA, ncol = 3, labels = "AUTO")
plt_val <- plot_grid(top.grob, plt_val, h2_leg, ncol = 1, rel_heights = c(.1, 1, .1))
plt_val

ggsave("h2_val.png", plt_val,  height = 5, width = 24, bg = "white")

# Now plot the response to selection per generation vs the prediction
d_com_resp_h2A <- read_csv(paste0(path, "d_com_resp_h2A.csv"))
d_com_resp_h2A <- d_com_resp_h2A %>%
  mutate(expPheno = 1 + (expPheno/50)) %>% # This undoes the correction I did for there being 50 generations between samples
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000) %>%
  complete(h2A, fill = list(phenoValue = 0))

h2A_resp <- ggplot(d_com_resp_h2A %>% mutate(h2A = fct_relevel(h2A,
                                                               "[0,0.25]", "(0.25, 0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~h2A) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")
h2A_resp


d_com_resp_h2D <- read_csv(paste0(path, "d_com_resp_h2D.csv"))
d_com_resp_h2D <- d_com_resp_h2D %>%
  mutate(expPheno = 1 + (expPheno/50)) %>% # This undoes the correction I did for there being 50 generations between samples
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000) %>%
  complete(h2D, fill = list(phenoValue = 0))

h2D_resp <- ggplot(d_com_resp_h2D %>% mutate(h2D = fct_relevel(h2D,
                                                               "[0,0.25]", "(0.25, 0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~h2D) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")
h2D_resp

d_com_resp_h2AA <- read_csv(paste0(path, "d_com_resp_h2AA.csv"))
d_com_resp_h2AA <- d_com_resp_h2AA %>%
  mutate(expPheno = 1 + (expPheno/50)) %>% # This undoes the correction I did for there being 50 generations between samples
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000) %>%
  complete(h2AA, fill = list(phenoValue = 0))

h2AA_resp <- ggplot(d_com_resp_h2AA %>% mutate(h2AA = fct_relevel(h2AA,
                                                               "[0,0.25]", "(0.25, 0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~h2AA) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")
h2AA_resp

