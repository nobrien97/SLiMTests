# Plot heritability estimates
library(tidyverse)
library(grid)
library(gridExtra)
library(latex2exp)
library(cowplot)

se <- function(x) {
  return(sd(x)/sqrt(length(x)))
}

d_h2 <- read_csv("../data/out_h2.csv", col_names = F)

names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
                 "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                 "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                 "H2.AA.SE", "AIC")

# Remove duplicates
d_h2 %>% distinct() -> d_h2

# Remove NAs
d_h2 %>% drop_na() -> d_h2

d_h2 %>% filter(VarA >= 0, VarA < 10,
                VarD >= 0, VarAA >= 0,
                VarR > 0, sum(VarA, VarD, VarAA, VarR) > 0) -> d_h2


# Add our other variables
combos <- read_delim("../../R/combos.csv", delim = " ", col_names = F)
names(combos) <- c("model", "nloci", "width", "locisigma")
d_h2 %>% mutate(model = combos$model[.$modelindex],
                nloci = combos$nloci[.$modelindex],
                width = combos$width[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_h2

# Recode modelindex to additive/network
d_h2 %>% mutate(model = recode_factor(model, `0`="Additive", `1`="Network")) -> d_h2

d_h2_sum <- d_h2 %>%
  group_by(gen, model, nloci, sigma) %>% 
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
ggplot(d_h2, aes(x = H2.A.Estimate, fill = model)) +
  geom_density(alpha = 0.4)



# Mean over time
h2A_time <- ggplot(d_h2_sum %>% filter(nloci != 100, sigma != 0.01), aes(x = gen, y = meanH2A, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2A - seH2A), ymax = (meanH2A + seH2A), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean heritability", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16))

# What about other variance components?
h2D_time <- ggplot(d_h2_sum %>% filter(nloci != 100, sigma != 0.01), aes(x = gen, y = meanH2D, color = model)) +
  facet_grid(nloci~sigma) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2D - seH2D), ymax = (meanH2D + seH2D), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average dominance variance proportion", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16))


h2AA_time <- ggplot(d_h2_sum %>% filter(nloci != 100, sigma != 0.01), aes(x = gen, y = meanH2AA, color = model)) +
  facet_grid(nloci~sigma) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanH2AA - seH2AA), ymax = (meanH2AA + seH2AA), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average AxA variance proportion", color = "Model") +
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


###################################################
# stacked percent bar chart
h2_stk_var0.1 <- ggplot(d_h2_sum %>% filter(nloci != 100, sigma == 0.1) %>% 
                          pivot_longer(cols = c(meanH2A, meanH2D, meanH2AA),
                                       names_to = "varComp", values_to = "prop"), 
                        aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Model", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~model) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

h2_stk_var1 <- ggplot(d_h2_sum %>% filter(nloci != 100, sigma == 1) %>% 
                        pivot_longer(cols = c(meanH2A, meanH2D, meanH2AA),
                                     names_to = "varComp", values_to = "prop"), 
                      aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Model", 
                                         breaks = NULL, labels = NULL)) +
  facet_grid(nloci~model) +
  geom_bar(position="stack", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

h2_leg <- get_legend(h2_stk_var0.1)

plt_stk <- plot_grid(h2_stk_var0.1 + theme(legend.position = "none"),
                     h2_stk_var1 + theme(legend.position = "none"),
                     labels = "AUTO")

plt_stk <- plot_grid(plt_stk, h2_leg, ncol = 1, rel_heights = c(1, .1))
plt_stk

ggsave("h2_stacked.png", plt_stk, width = 20, height = 5, bg = "white")


# Plot trait data over time
d_qg <- read_csv("../data/slim_qg.csv", col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")
d_com <- inner_join(d_h2, d_qg, by = c("gen", "seed", "modelindex"))
write_csv(d_com, "d_h2_qg.csv")

d_qg %>% filter(phenomean < 100) -> d_qg

d_qg %>% mutate(model = combos$model[.$modelindex],
                nloci = combos$nloci[.$modelindex],
                width = combos$width[.$modelindex],
                sigma = combos$locisigma[.$modelindex]) -> d_qg

# Recode modelindex to additive/network
d_qg %>% mutate(model = recode_factor(model, `0`="Additive", `1`="Network")) -> d_qg



d_qg_sum <- d_qg %>% group_by(gen, model, nloci, sigma) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

ggplot(d_qg_sum %>% filter(nloci < 100, sigma != 0.01), 
       aes(gen, meanPheno, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Model") +
  theme_bw() +
  theme(text = element_text(size=20), panel.spacing = unit(1, "lines"))



d_com <- d_com %>% filter(phenomean < 10)

d_com_sum <- d_com %>% group_by(gen, model, nloci, sigma) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist))

d_com_sum$gen <- d_com_sum$gen - 50000

ggplot(d_com_sum %>% filter(nloci != 100, sigma != 0.01), aes(gen, meanPheno, color = model)) +
  facet_grid(nloci~sigma) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  scale_fill_discrete(guide = "none") +
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population mean phenotypes", color = "Model") +
  theme_bw() +
  theme(text = element_text(size=20), panel.spacing = unit(1, "lines"))

ggsave("pheno_time.png", width = 8, height = 6, bg = "white")

# Mutation data
d_com <- read_csv("../data/d_com_sum.csv")

# Plot mean freq and value

ggplot(d_com %>% nloci != 100)

