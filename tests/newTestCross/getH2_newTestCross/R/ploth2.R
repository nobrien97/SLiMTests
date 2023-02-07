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
path <- "/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/"

# Colour palette
cc_ibm <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000")

d_h2 <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/getH2_newTestCross/data/out_h2_total.csv", col_names = F)

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
names(combos) <- c("model", "nloci", "locisigma")
d_h2 %>% mutate(fixedEffect = combos$model[.$modelindex],
                nloci = combos$nloci[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_h2

# Recode fixedEffect to mol trait name
d_h2 %>% mutate(fixedEffect = recode_factor(fixedEffect, `-1`="None", `0`="\u03B1", `1`="\u03B2", `2`="KZ", `3`="KXZ")) -> d_h2


# Remove weird values: negative variance, huge variance
d_h2 %>% filter(VarA >= 0, VarA < 10,
                VarD >= 0, VarAA >= 0,
                VarR > 0, sum(VarA, VarD, VarAA, VarR) > 0) -> d_h2

d_h2$nloci <- d_h2$nloci + 2


d_h2_sum <- d_h2 %>%
  group_by(gen, fixedEffect, nloci, sigma) %>% 
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
ggplot(d_h2, aes(x = H2.A.Estimate, fill = as.factor(fixedEffect))) +
  geom_density(alpha = 0.4)

ggplot(d_h2, aes(x = H2.AA.Estimate, fill = as.factor(fixedEffect))) +
  geom_density(alpha = 0.4)

ggplot(d_h2, aes(x = H2.D.Estimate, fill = as.factor(fixedEffect))) +
  geom_density(alpha = 0.4)


# Mean over time
h2A_time <- ggplot(d_h2_sum %>% filter(nloci != 100), aes(x = gen, y = meanH2A, color = fixedEffect)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2A - seH2A), ymax = (meanH2A + seH2A), 
                  fill = fixedEffect), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean heritability", color = "Fixed effect") +
  theme_bw() +
  theme(text=element_text(size=16))

# What about other variance components?
h2D_time <- ggplot(d_h2_sum %>% filter(nloci != 100), aes(x = gen, y = meanH2D, color = fixedEffect)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2D - seH2D), ymax = (meanH2D + seH2D), 
                  fill = fixedEffect), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average dominance variance proportion", color = "Fixed effect") +
  theme_bw() +
  theme(text=element_text(size=16))


h2AA_time <- ggplot(d_h2_sum %>% filter(nloci != 100), aes(x = gen, y = meanH2AA, color = fixedEffect)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2AA - seH2AA), ymax = (meanH2AA + seH2AA), 
                  fill = fixedEffect), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average AxA variance proportion", color = "Fixed effect") +
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

y.grob <- textGrob("Number of QTLs", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Mutational effect size variance", 
                   gp = gpar(fontsize = 16))

g <- arrangeGrob(plt_row, bottom = x.grob, right = y.grob, top = top.grob)

h2_leg <- get_legend(h2A_time)

plt_row_leg <- plot_grid(g, h2_leg, rel_widths = c(3, .4))
plt_row_leg

ggsave("h2_time.png", plt_row_leg,  height = 6, width = 18, bg = "white")

# Now averaged over time
d_h2_sum_notime <- d_h2 %>%
  group_by(fixedEffect, nloci, sigma) %>% 
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


h2_avg <- ggplot(d_h2_sum_notime %>% filter(nloci != 100) %>% pivot_longer(
  cols = c(meanH2A, meanH2D, meanH2AA),
  names_to = "varComp", values_to = "prop") %>%
    mutate(propSE = case_when(varComp == "meanH2A" ~ seH2A,
                              varComp == "meanH2D" ~ seH2D,
                              varComp == "meanH2AA" ~ seH2AA)), 
  aes(x = fixedEffect, y = prop, fill = varComp)) +
  facet_grid(nloci~sigma) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = prop - propSE, ymax = prop + propSE), 
                position = position_dodge(0.9), width = 0.3) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Fixed molecular trait", 
         y = "Proportion of total\nphenotypic variance (%)",
         fill = "Variance component") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2_avg)
y.grob <- textGrob("Number of QTLs", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Mutational effect size variance", 
                     gp = gpar(fontsize = 16))

g <- arrangeGrob(h2_avg + theme(legend.position = "none"), right = y.grob, top = top.grob)


plt_row_leg <- plot_grid(g, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_row_leg

ggsave("h2_avg.png", plt_row_leg,  height = 8, width = 8, bg = "white")

# As a deviation

# calc means
d_h2_ctrl <- d_h2 %>% filter(fixedEffect == "None") %>%
  group_by(nloci, sigma) %>% 
  summarise(meanH2A_ctrl = mean(H2.A.Estimate),
            seH2A_ctrl = se(H2.A.Estimate),
            meanH2D_ctrl = mean(H2.D.Estimate), 
            seH2D_ctrl = se(H2.D.Estimate),
            meanH2AA_ctrl = mean(H2.AA.Estimate),
            seH2AA_ctrl = se(H2.AA.Estimate),
            meanVarA_ctrl = mean(VarA),
            seVarA_ctrl = se(VarA),
            meanVarD_ctrl = mean(VarD),
            seVarD_ctrl = se(VarD),
            meanVarAA_ctrl = mean(VarAA),
            seVarAA_ctrl = se(VarAA)
  )

# Add to main dataframe
d_h2_total <- inner_join(d_h2, d_h2_ctrl, by = c("nloci", "sigma"))

# We do want directionality, so not squared deviation
d_h2_total <- d_h2_total %>% 
  group_by(nloci, sigma) %>% 
  mutate(devA = (H2.A.Estimate - meanH2A_ctrl),
         devD = (H2.D.Estimate - meanH2D_ctrl),
         devAA = (H2.AA.Estimate - meanH2AA_ctrl))

# summarise
d_h2_sum_notime <- d_h2_total %>%
  group_by(fixedEffect, nloci, sigma) %>% 
  summarise(meanDevH2A = mean(devA),
            seDevH2A = se(devA),
            meanDevH2D = mean(devD), 
            seDevH2D = se(devD),
            meanDevH2AA = mean(devAA),
            seDevH2AA = se(devAA),
  )

# plot
h2_avg <- ggplot(d_h2_sum_notime %>% filter(nloci != 100, fixedEffect != "None") %>% pivot_longer(
  cols = c(meanDevH2A, meanDevH2D, meanDevH2AA),
  names_to = "varComp", values_to = "prop") %>%
    mutate(propSE = case_when(varComp == "meanDevH2A" ~ seDevH2A,
                              varComp == "meanDevH2D" ~ seDevH2D,
                              varComp == "meanDevH2AA" ~ seDevH2AA)), 
  aes(x = fixedEffect, y = prop, fill = varComp)) +
  facet_grid(nloci~sigma) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = prop - propSE, ymax = prop + propSE), 
                position = position_dodge(0.9), width = 0.3, linewidth = 0.25) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  labs(x = "Fixed molecular trait", 
       y = "Deviation from no fixed effect",
       fill = "Variance component") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2_avg)
y.grob <- textGrob("Number of QTLs", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Mutational effect size variance", 
                     gp = gpar(fontsize = 16))

g <- arrangeGrob(h2_avg + theme(legend.position = "none"), right = y.grob, top = top.grob)


plt_row_leg <- plot_grid(g, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_row_leg

ggsave("h2_deviation.png", plt_row_leg,  height = 6, width = 8, bg = "white")



###################################################
# stacked percent bar chart
h2_stk_var0.1 <- ggplot(d_h2_sum %>% filter(sigma == 0.1) %>% 
                          pivot_longer(cols = c(meanH2A, meanH2D, meanH2AA),
                                       names_to = "varComp", values_to = "prop") %>%
                          mutate(varComp = factor(varComp, c("meanH2A", "meanH2D", "meanH2AA"))), 
  aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Epistatic (AxA)")) +
  scale_y_continuous(labels = scales::percent, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Fixed molecular trait", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~fixedEffect) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

h2_stk_var1 <- ggplot(d_h2_sum %>% filter(sigma == 1) %>% 
                          pivot_longer(cols = c(meanH2A, meanH2D, meanH2AA),
                                       names_to = "varComp", values_to = "prop") %>%
                        mutate(varComp = factor(varComp, c("meanH2A", "meanH2D", "meanH2AA"))),
                        aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Dominant", "Epistatic (AxA)")) +
  scale_y_continuous(labels = scales::percent, 
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Fixed molecular trait", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~fixedEffect) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2_stk_var0.1)

plt_stk <- plot_grid(h2_stk_var0.1 + theme(legend.position = "none"),
                     h2_stk_var1 + theme(legend.position = "none"),
                     nrow = 2,
                     labels = "AUTO")

plt_stk <- plot_grid(plt_stk, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_stk

ggsave("h2_stacked.png", plt_stk, width = 20, height = 10, bg = "white")

# Now just the unfixed data
h2_stk_var <- ggplot(d_h2_sum %>% filter(fixedEffect == "None") %>% 
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

ggsave("h2_stacked_unfixed.png", h2_stk_var, width = 10, height = 8, bg = "white")


# Plot trait data over time
d_qg <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/getH2_newTestCross/data/slim_qg_total.csv", col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw", "aZ", "bZ", "KZ", "KXZ")

d_qg %>% filter(phenomean < 10) -> d_qg

d_qg %>% mutate(fixedEffect = combos$model[.$modelindex],
                nloci = combos$nloci[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_qg

# Recode fixedEffect to mol trait name
d_qg %>% mutate(fixedEffect = recode_factor(fixedEffect, `-1`="None", `0`="\u03B1", `1`="\u03B2", `2`="KZ", `3`="KXZ")) -> d_qg

d_qg$nloci <- d_qg$nloci + 2

d_qg_sum <- d_qg %>% group_by(gen, fixedEffect, nloci, sigma) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

pheno_time <- ggplot(d_qg_sum %>% filter(gen > 49000) %>% mutate(gen = gen - 50000), 
       aes(gen, meanPheno, color = fixedEffect)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_colour_manual(values = cc_ibm) +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = fixedEffect), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Fixed effect") +
  theme_bw() +
  theme(text = element_text(size=20), panel.spacing = unit(1, "lines"))

ggsave("pheno_time.png", pheno_time, height = 9, width = 12)

# Now with only unfixed models
pheno_time <- ggplot(d_qg_sum %>% filter(gen > 49000, fixedEffect == "None") %>% mutate(gen = gen - 50000), 
                     aes(gen, meanPheno, color = as.factor(nloci))) +
  facet_grid(.~sigma) +
  geom_line() +
  ylim(0, 2.5) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  scale_fill_manual(values = cc_ibm, guide = "none") +
  scale_color_manual(values = cc_ibm, guide = guide_legend(override.aes = list(linewidth = 3))) +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = as.factor(nloci)), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = TeX("$\\mu_{\\bar{z}}$"), color = "Number of loci") +
  theme_bw() +
  theme(text = element_text(size=20), panel.spacing = unit(2, "lines"), legend.position = "bottom")
pheno_time

saveRDS(pheno_time, paste0(path, "pheno_time_unfixed.RDS"))
ggsave("pheno_time_unfixed.png", pheno_time, height = 5, width = 10)



d_com <- full_join(d_h2, d_qg, by = c("gen", "seed", "fixedEffect", "nloci", "sigma"))
d_com <- d_com %>% filter(phenomean < 10)

d_com_sum <- d_com %>% group_by(gen, fixedEffect, nloci, sigma) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

d_com_sum$gen <- d_com_sum$gen - 50000

ggplot(d_com_sum %>% filter(gen > -1000), aes(gen, meanPheno, color = fixedEffect)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_fill_viridis_d(guide = "none") +
  geom_hline(yintercept = 2, linetype = "dashed") + 
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = fixedEffect), color = NA, alpha = 0.2) +
  scale_color_viridis_d() +
  labs(x = "Generations after optimum shift", y = "Mean of population\nmean phenotypes", color = "Fixed\nmolecular\ntrait") +
  theme_bw() +
  theme(text = element_text(size=16))

ggsave("pheno_time.png", height = 4, width = 6)


# combined data: SFS
#TODO: There's too much information to unpack here, need a metric to describe the
# SFS, maybe just mean frequency? And similar deviation from model w/o fixed effect
d_com <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/d_combined_after.csv")
d_com$fixTime <- d_com$fixGen - d_com$originGen

d_com_end <- d_com %>% filter(gen > 51000)
d_com_end %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                     h2D = cut(H2.D.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                     h2AA = cut(H2.AA.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1), include.lowest = T),
                     fixTime = fixGen - originGen) -> d_com_end

d_com_end %>% mutate(molTrait = recode_factor(mutType, `1`="Neutral", `2`="Del", `3`="\u03B1", `4`="\u03B2", `5`="KZ", `6`="KXZ")) -> d_com_end

d_com_end %>% mutate(fixedEffect = combos$model[.$modelindex],
                nloci = combos$nloci[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_com_end

d_com_end %>% mutate(fixedEffect = recode_factor(fixedEffect, `-1`="None", `0`="\u03B1", `1`="\u03B2", `2`="KZ", `3`="KXZ")) -> d_com_end
d_com_end$nloci <- d_com_end$nloci + 2

d_com_begin <- d_com %>% filter(gen < 51000)
d_com_begin %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                       h2D = cut(H2.D.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                       h2AA = cut(H2.AA.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                       fixTime = fixGen - originGen) -> d_com_begin
d_com_begin %>% mutate(molTrait = recode_factor(mutType, `3`="aZ", `4`="bZ", `5`="KZ", `6`="KXZ")) -> d_com_begin

rm(d_com)

d_com_sum <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(fixedEffect, nloci, molTrait, h2A) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq)
            )

# Plot mean frequency for each group: averaged across sigma
h2A_freq <- ggplot(d_com_sum %>% filter(nloci != 100), 
       aes(x = h2A, y = meanFreq, fill = fixedEffect)) +
  facet_grid(molTrait~nloci) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanFreq - seFreq, ymax = meanFreq + seFreq), 
                 position = position_dodge(0.9, preserve = "single"), width = 0.3) +
  labs(x = "Heritability", y = "Mean frequency", fill = "Fixed\nmolecular\ntrait") +
  theme_bw() + theme(text = element_text(size = 16))
h2A_freq

d_com_sum <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(fixedEffect, nloci, molTrait, h2AA) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq)
  )

h2AA_freq <- ggplot(d_com_sum %>% filter(nloci != 100), 
                   aes(x = h2AA, y = meanFreq, fill = fixedEffect)) +
  facet_grid(molTrait~nloci) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanFreq - seFreq, ymax = meanFreq + seFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3) +
  labs(x = "Heritability (AxA)", y = "Mean frequency", fill = "Fixed\nmolecular\ntrait") +
  theme_bw() + theme(text = element_text(size = 16))
h2AA_freq

d_com_sum <- d_com_end %>%
  select(-c(fixGen, fixTime)) %>%
  drop_na() %>%
  group_by(fixedEffect, nloci, molTrait, h2D) %>%
  summarise(meanFreq = mean(Freq),
            seFreq = se(Freq)
  )

h2D_freq <- ggplot(d_com_sum %>% filter(nloci != 100), 
                    aes(x = h2D, y = meanFreq, fill = fixedEffect)) +
  facet_grid(molTrait~nloci) +
  geom_col(position = position_dodge(0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = meanFreq - seFreq, ymax = meanFreq + seFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3) +
  labs(x = "Heritability (D)", y = "Mean frequency", fill = "Fixed\nmolecular\ntrait") +
  theme_bw() + theme(text = element_text(size = 16))
h2D_freq



# Deviations from control
d_com_ctrl_sum <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/d_com_sum.csv")
d_com_ctrl_sum$nloci <- d_com_ctrl_sum$nloci + 2
# Add to main dataframe
d_com_total <- inner_join(d_com_end, d_com_ctrl_sum, by = c("nloci", "molTrait", "h2A", "h2D", "h2AA"))


# We do want directionality, so not squared deviation
d_com_total <- d_com_total %>% 
  group_by(nloci, sigma, molTrait, fixedEffect) %>% 
  mutate(devFreq = (Freq - meanFreq),
         devFX = (value - meanFX))

d_com_total$molTrait <- factor(d_com_total$molTrait, ordered = T,
                               levels = c("\u03B1", "\u03B2", "KZ", "KXZ"))


d_com_total %>% 
  mutate(fixedEffect = recode_factor(fixedEffect, 
                                     `-1`="None", 
                                     `0`="\u03B1", 
                                     `1`="\u03B2", 
                                     `2`="KZ", 
                                     `3`="KXZ", .ordered = T)) -> d_com_total

# Plot each heritability
y_lims <- c(-0.55, 1)

# summarise
d_freq_sum_notime <- d_com_total %>%
  group_by(nloci, molTrait, fixedEffect, h2A) %>%
  summarise(meanDevFreq = mean(devFreq),
            seDevFreq = se(devFreq)
  )


h2A_freq <- ggplot(d_freq_sum_notime %>% filter(fixedEffect != "None", nloci != 100) %>%
         mutate(h2A = fct_relevel(h2A, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = h2A, y = meanDevFreq, fill = molTrait)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevFreq - seDevFreq, ymax = meanDevFreq + seDevFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_A$)"), y = "Deviation in mean frequency", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2A_freq)
y.grob <- textGrob("Fixed molecular trait", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Number of QTLs", 
                     gp = gpar(fontsize = 16))

g_h2A <- arrangeGrob(h2A_freq + theme(legend.position = "none"), right = y.grob, top = top.grob)

d_freq_sum_notime <- d_com_total %>%
  group_by(nloci, molTrait, fixedEffect, h2D) %>%
  summarise(meanDevFreq = mean(devFreq),
            seDevFreq = se(devFreq)
  )


h2D_freq <- ggplot(d_freq_sum_notime %>% filter(fixedEffect != "None", nloci != 100) %>%
                     mutate(h2D = fct_relevel(h2D, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = h2D, y = meanDevFreq, fill = molTrait)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevFreq - seDevFreq, ymax = meanDevFreq + seDevFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_D$)"), y = "Deviation in mean frequency", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2D <- arrangeGrob(h2D_freq + theme(legend.position = "none"), right = y.grob, top = top.grob)

d_freq_sum_notime <- d_com_total %>%
  group_by(nloci, molTrait, fixedEffect, h2AA) %>%
  summarise(meanDevFreq = mean(devFreq),
            seDevFreq = se(devFreq)
  )


h2AA_freq <- ggplot(d_freq_sum_notime %>% filter(fixedEffect != "None", nloci != 100) %>%
                     mutate(h2AA = fct_relevel(h2AA, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = h2AA, y = meanDevFreq, fill = molTrait)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevFreq - seDevFreq, ymax = meanDevFreq + seDevFreq), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_{AA}$)"), y = "Deviation in mean frequency", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2AA <- arrangeGrob(h2AA_freq + theme(legend.position = "none"), right = y.grob, top = top.grob)


plt_freq <- plot_grid(g_h2A, g_h2D, g_h2AA, ncol = 3, labels = "AUTO")
plt_freq <- plot_grid(plt_freq, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_freq

ggsave("h2_freq_deviation.png", plt_freq,  height = 7, width = 20, bg = "white")

# Effect size
# Plot each heritability
y_lims <- c(-3, 3)

# summarise
d_val_sum_notime <- d_com_total %>%
  group_by(nloci, molTrait, fixedEffect, h2A) %>%
  summarise(meanDevVal = mean(devFX),
            seDevVal = se(devFX)
  )

h2A_val <- ggplot(d_val_sum_notime %>% filter(fixedEffect != "None", nloci != 100) %>%
                     mutate(h2A = fct_relevel(h2A, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), 
                   aes(x = h2A, y = meanDevVal, fill = molTrait)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevVal - seDevVal, ymax = meanDevVal + seDevVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_A$)"), y = "Deviation in mean molecular effect", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2A_val)
y.grob <- textGrob("Fixed molecular trait", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Number of QTLs", 
                     gp = gpar(fontsize = 16))

g_h2A <- arrangeGrob(h2A_val + theme(legend.position = "none"), right = y.grob, top = top.grob)


d_val_sum_notime <- d_com_total %>%
  group_by(nloci, molTrait, fixedEffect, h2D) %>%
  summarise(meanDevVal = mean(devFX),
            seDevVal = se(devFX)
  )

h2D_val <- ggplot(d_val_sum_notime %>% filter(fixedEffect != "None", nloci != 100) %>%
                    mutate(h2D = fct_relevel(h2D, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), 
                  aes(x = h2D, y = meanDevVal, fill = molTrait)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevVal - seDevVal, ymax = meanDevVal + seDevVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_D$)"), y = "Deviation in mean molecular effect", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2D <- arrangeGrob(h2D_val + theme(legend.position = "none"), right = y.grob, top = top.grob)

d_val_sum_notime <- d_com_total %>%
  group_by(nloci, molTrait, fixedEffect, h2AA) %>%
  summarise(meanDevVal = mean(devFX),
            seDevVal = se(devFX)
  )

h2AA_val <- ggplot(d_val_sum_notime %>% filter(fixedEffect != "None", nloci != 100) %>%
                    mutate(h2AA = fct_relevel(h2AA, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), 
                  aes(x = h2AA, y = meanDevVal, fill = molTrait)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevVal - seDevVal, ymax = meanDevVal + seDevVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_{AA}$)"), y = "Deviation in mean molecular effect", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2AA <- arrangeGrob(h2AA_val + theme(legend.position = "none"), right = y.grob, top = top.grob)

plt_val <- plot_grid(g_h2A, g_h2D, g_h2AA, ncol = 3, labels = "AUTO")
plt_val <- plot_grid(plt_val, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_val

ggsave("h2_val_deviation.png", plt_val,  height = 7, width = 20, bg = "white")



# Molecular trait value deviation
d_com_ctrl_molTraitVal_sum <- read_csv("/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/d_com_sum_molTraitVal.csv")

# Pivot longer
d_com_end %>% pivot_longer(c(aZ, bZ, KZ, KXZ), names_to = "molTraitVal_name", 
                                values_to = "molTraitVal_value") -> d_com_end


d_com_ctrl_molTraitVal_sum$nloci <- d_com_ctrl_molTraitVal_sum$nloci + 2

# Add to main dataframe
d_com_total <- inner_join(d_com_end, d_com_ctrl_molTraitVal_sum, by = c("nloci", "molTraitVal_name", "h2A", "h2D", "h2AA"))


# We do want directionality, so not squared deviation
d_com_total <- d_com_total %>% 
  group_by(nloci, sigma, molTraitVal_name, fixedEffect) %>% 
  mutate(devMolTraitVal = (molTraitVal_value - meanMolTraitVal))

d_com_total %>% 
  mutate(fixedEffect = recode_factor(fixedEffect, 
                                     `-1`="None", 
                                     `0`="\u03B1", 
                                     `1`="\u03B2", 
                                     `2`="KZ", 
                                     `3`="KXZ"),
         molTraitVal_name = recode_factor(molTraitVal_name, 
                                     "aZ"="\u03B1", 
                                     "bZ"="\u03B2", 
                                     "KZ"="KZ", 
                                     "KXZ"="KXZ")) -> d_com_total


d_val_sum_notime <- d_com_total %>% filter(molTraitVal_value < 10) %>%
  group_by(nloci, molTraitVal_name, fixedEffect, h2A) %>%
  summarise(meanDevVal = mean(devMolTraitVal),
            seDevVal = se(devMolTraitVal)
  )


# Now plot


h2A_val <- ggplot(d_val_sum_notime %>% filter(nloci != 100, fixedEffect != "None") %>%
         mutate(h2A = fct_relevel(h2A, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")),
       aes(x = h2A, y = meanDevVal, fill = molTraitVal_name)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevVal - seDevVal, ymax = meanDevVal + seDevVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_{A}$)"), y = "Deviation in mean molecular trait value", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2A_val)
y.grob <- textGrob("Fixed molecular trait", 
                   gp = gpar(fontsize = 16), rot=270)

top.grob <- textGrob("Number of QTLs", 
                     gp = gpar(fontsize = 16))

g_h2A <- arrangeGrob(h2A_val + theme(legend.position = "none"), right = y.grob, top = top.grob)



d_val_sum_notime <- d_com_total %>% filter(molTraitVal_value < 10) %>%
  group_by(nloci, molTraitVal_name, fixedEffect, h2D) %>%
  summarise(meanDevVal = mean(devMolTraitVal),
            seDevVal = se(devMolTraitVal)
  )


h2D_val <- ggplot(d_val_sum_notime %>% filter(nloci != 100, fixedEffect != "None") %>%
                    mutate(h2D = fct_relevel(h2D, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")),
                  aes(x = h2D, y = meanDevVal, fill = molTraitVal_name)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevVal - seDevVal, ymax = meanDevVal + seDevVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_{D}$)"), y = "Deviation in mean molecular trait value", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2D <- arrangeGrob(h2D_val + theme(legend.position = "none"), right = y.grob, top = top.grob)


d_val_sum_notime <- d_com_total %>% filter(molTraitVal_value < 10) %>%
  group_by(nloci, molTraitVal_name, fixedEffect, h2AA) %>%
  summarise(meanDevVal = mean(devMolTraitVal),
            seDevVal = se(devMolTraitVal)
  )


h2AA_val <- ggplot(d_val_sum_notime %>% filter(nloci != 100, fixedEffect != "None") %>%
                    mutate(h2AA = fct_relevel(h2AA, "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")),
                  aes(x = h2AA, y = meanDevVal, fill = molTraitVal_name)) +
  facet_grid(fixedEffect~nloci, scales = "free_x") +
  scale_y_continuous(limits = y_lims) +
  scale_x_discrete(labels = c(TeX("$<0.25$"), TeX("$<0.5$"), TeX("$<0.75$"), TeX(r"($\leq 1$)"))) +
  scale_fill_viridis_d(drop = FALSE, labels = c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), TeX("$K_Z$"), TeX("$K_{XZ}$"))) +
  geom_bar(position = position_dodge(0.9, preserve = "single"), stat="identity") +
  geom_errorbar(aes(ymin = meanDevVal - seDevVal, ymax = meanDevVal + seDevVal), 
                position = position_dodge(0.9, preserve = "single"), width = 0.3, linewidth = 0.25) +
  labs(x = TeX("Heritability ($h^2_{AA}$)"), y = "Deviation in mean molecular trait value", fill = "Molecular trait") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2AA <- arrangeGrob(h2AA_val + theme(legend.position = "none"), right = y.grob, top = top.grob)

plt_val <- plot_grid(g_h2A, g_h2D, g_h2AA, ncol = 3, labels = "AUTO")
plt_val <- plot_grid(plt_val, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_val

ggsave("h2_molTraitVal_deviation.png", plt_val,  height = 7, width = 20, bg = "white")

# Now plot the response to selection per generation vs the prediction
path <- "/mnt/d/SLiMTests/tests/newTestCross/moreReps/getH2_newTestCross/data/"

d_com_resp_h2A <- read_csv(paste0(path, "d_com_resp_h2A.csv"))
d_com_resp_h2A <- d_com_resp_h2A %>%
  #mutate(expPheno = expPheno/50) %>% # This undoes the correction I did for there being 50 generations between samples
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000) %>%
  complete(h2A, fill = list(phenoValue = 0)) %>%
  mutate(h2A = fct_relevel(h2A, "≤0.5", ">0.5"))

h2A_resp <- ggplot(d_com_resp_h2A, 
                   aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~h2A) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2A_resp)
top.grob <- textGrob("Heritability", 
                   gp = gpar(fontsize = 16))

y.grob <- textGrob("Number of QTLs", 
                     gp = gpar(fontsize = 16), rot=270)

g_h2A <- arrangeGrob(h2A_resp + theme(legend.position = "none"), right = y.grob, top = top.grob)



d_com_resp_h2D <- read_csv(paste0(path, "d_com_resp_h2D.csv"))
d_com_resp_h2D <- d_com_resp_h2D %>%
  #mutate(expPheno = expPheno/50) %>% # This undoes the correction I did for there being 50 generations between samples
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000) %>%
  complete(h2D, fill = list(phenoValue = 0)) %>%
  mutate(h2D = fct_relevel(h2D, "≤0.5", ">0.5"))


h2D_resp <- ggplot(d_com_resp_h2D, 
                   aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~h2D) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2D <- arrangeGrob(h2D_resp + theme(legend.position = "none"), right = y.grob, top = top.grob)


d_com_resp_h2AA <- read_csv(paste0(path, "d_com_resp_h2AA.csv"))
d_com_resp_h2AA <- d_com_resp_h2AA %>%
  #mutate(expPheno = expPheno/50) %>% # This undoes the correction I did for there being 50 generations between samples
  pivot_longer(c(meanPheno, expPheno), names_to = "phenoType", values_to = "phenoValue") %>%
  mutate(phenoType = recode_factor(phenoType,
                                   "expPheno" = "Estimated",
                                   "meanPheno" = "Observed"),
         gen = gen - 50000) %>%
  complete(h2AA, fill = list(phenoValue = 0)) %>%
  mutate(h2AA = fct_relevel(h2AA, "≤0.5", ">0.5"))


h2AA_resp <- ggplot(d_com_resp_h2AA,
                    aes(x = gen, y = phenoValue, colour = phenoType)) +
  facet_grid(nloci~h2AA) +
  scale_colour_manual(values = c(cc_ibm[1], cc_ibm[4])) +
  geom_line() +
  labs(x = "Generations after optimum shift", y = "Mean trait value", colour = "Trait value origin") +
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom")

g_h2AA <- arrangeGrob(h2AA_resp + theme(legend.position = "none"), right = y.grob, top = top.grob)


plt_resp <- plot_grid(g_h2A, g_h2D, g_h2AA, ncol = 3, labels = "AUTO")
plt_resp <- plot_grid(plt_resp, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_resp

ggsave("h2_R.png", plt_resp, height = 7, width = 20, bg = "white")

# Plot only the additive one since the others don't really give us any more information
top.grob <- textGrob(TeX("Heritability $(\\sigma^2_A)$"), 
                     gp = gpar(fontsize = 16))
g_h2A <- arrangeGrob(h2A_resp + theme(legend.position = "none"), right = y.grob, top = top.grob)
plt_resp_A <- plot_grid(g_h2A, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_resp_A
ggsave("h2_R_Add.png", plt_resp_A, height = 7, width = 8, bg = "white")



# Time to plot molecular trait landscapes



