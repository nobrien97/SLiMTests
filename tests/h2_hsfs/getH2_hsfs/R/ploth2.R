# Plot heritability estimates
library(tidyverse)
library(gridExtra)
library(latex2exp)
library(cowplot)

se <- function(x) {
  return(sd(x)/sqrt(length(x)))
}

d_h2 <- read_csv("../data/new/out_h2.csv", col_names = F)

names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
                 "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
                 "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
                 "H2.AA.SE", "AIC")

# Remove duplicates
d_h2 %>% distinct() -> d_h2

# Remove NAs
d_h2 %>% drop_na() -> d_h2

# Recode modelindex to additive/network
d_h2 %>% mutate(model = recode_factor(modelindex, `0`="Additive", `1`="Network")) -> d_h2

# Remove weird values: negative variance, huge variance
d_h2 %>% filter(VarA >= 0, VarA < 10,
                VarD >= 0, VarAA >= 0,
                VarR > 0, sum(VarA, VarD, VarAA, VarR) > 0) -> d_h2

d_h2$VarA_scl <- scale(d_h2$VarA)
d_h2$VarD_scl <- scale(d_h2$VarD)
d_h2$VarAA_scl <- scale(d_h2$VarAA)

d_h2_sum <- d_h2 %>%
  group_by(gen, model) %>% 
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
h2A_time <- ggplot(d_h2_sum, aes(x = gen, y = meanH2A, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2A - seH2A), ymax = (meanH2A + seH2A), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean heritability", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16))

# What about other variance components?
h2D_time <- ggplot(d_h2_sum, aes(x = gen, y = meanH2D, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2D - seH2D), ymax = (meanH2D + seH2D), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average dominance variance proportion", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16))


h2AA_time <- ggplot(d_h2_sum, aes(x = gen, y = meanH2AA, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  geom_ribbon(aes(ymin = (meanH2AA - seH2AA), ymax = (meanH2AA + seH2AA), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Average AxA variance proportion", color = "Model") +
  theme_bw() +
  theme(text=element_text(size=16))

plt_row <- plot_grid(
  h2A_time + theme(legend.position = "none") + xlab(NULL) + 
    draw_label(TeX("$h^2_{A}$"), x = 1950, y = 0.95),
  h2D_time + theme(legend.position = "none") + ylab(NULL) + xlab(NULL) +
    draw_label(TeX("$h^2_{D}$"), x = 1950, y = 0.95),
  h2AA_time + theme(legend.position = "none") + ylab(NULL) + xlab(NULL) +
    draw_label(TeX("$h^2_{AA}$"), x = 1950, y = 0.95),
  nrow = 1
)

h2_leg <- get_legend(h2A_time)

plt_row_leg <- plot_grid(plt_row, h2_leg, rel_widths = c(3, .4))

x.grob <- textGrob("Generations after optimum shift", 
                   gp = gpar(fontsize = 16))

g <- arrangeGrob(plt_row_leg, bottom = x.grob)

ggsave("h2_time.png", g,  height = 3, width = 12, bg = "white")

# stacked percent bar chart
ggplot(d_h2_sum %>% pivot_longer(
  cols = c(meanH2A, meanH2D, meanH2AA, meanH2R),
  names_to = "varComp", values_to = "prop"), 
  aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(.~model) +
  geom_bar(position="fill", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total\nphenotypic variance", 
       fill = "Variance\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggsave("h2_stacked.png", width = 8, height = 4)


# Plot trait data over time
d_qg <- read_csv("../data/new/slim_qg.csv", col_names = F)
names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
                 "phenovar", "dist", "w", "deltaPheno", "deltaw")
d_com <- inner_join(d_h2, d_qg, by = c("gen", "seed", "modelindex"))
d_com <- d_com %>% filter(phenomean < 10)

d_com_sum <- d_com %>% group_by(gen, model) %>%
  summarise(He = mean(meanH),
            seHe = se(meanH),
            meanPheno = mean(phenomean),
            sePheno = se(phenomean),
            meanDist = mean(dist),
            seDist = se(dist),
            meanH2A = mean(H2.A.Estimate),
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
            seVarAA = se(VarAA))

d_com_sum$gen <- d_com_sum$gen - 50000

ggplot(d_com_sum, aes(gen, meanPheno, color = model)) +
  geom_line() +
  scale_fill_discrete(guide = "none") +
  geom_hline(yintercept = 2, linetype = "dashed") + 
  geom_ribbon(aes(ymin = (meanPheno - sePheno), ymax = (meanPheno + sePheno), 
                  fill = model), color = NA, alpha = 0.2) +
  labs(x = "Generations after optimum shift", y = "Mean of population\nmean phenotypes", color = "Model") +
  theme_bw() +
  theme(text = element_text(size=16))

ggsave("pheno_time.png", height = 4, width = 6)

# plots of variance components per phenotype range
ggplot(d_com_sum %>% pivot_longer(
  cols = c(meanH2A, meanH2D, meanH2AA),
  names_to = "varComp", values_to = "prop") %>%
  aes(x = gen, y = prop, fill = varComp)) +
  scale_fill_viridis_d(labels = c("Additive", "Epistatic (AxA)", "Dominant")) +
  facet_grid(meanPheno~model) +
  geom_bar(position="fill", stat="identity", width = 50) +
  labs(x = "Generations after optimum shift", y = "Proportion of total phenotypic variance", 
       fill = "Variance component") +
  theme_bw() +
  theme(text = element_text(size = 20))


# combined data: SFS
d_com <- read_csv("../data/new/d_combined_after.csv")
d_com %>% mutate(model = recode_factor(modelindex, `0`="Additive", `1`="Network"),
                 molTrait = recode_factor(mutType, `1`="Neutral", `2`="Del", `3`="\u03B1", `4`="\u03B2", `5`="KZ", `6`="KXZ")) -> d_com
d_com$fixTime <- d_com$fixGen - d_com$originGen

d_com_end <- d_com %>% filter(gen > 51000)
d_com_end %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                     h2D = cut(H2.D.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                     h2AA = cut(H2.AA.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                     fixTime = fixGen - originGen) -> d_com_end

d_com_begin <- d_com %>% filter(gen < 51000)
d_com_begin %>% mutate(h2A = cut(H2.A.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                     h2D = cut(H2.D.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                     h2AA = cut(H2.AA.Estimate, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                     fixTime = fixGen - originGen) -> d_com_begin

rm(d_com)

ggplot(d_com_end %>% filter(modelindex == 1, !is.na(h2AA), !is.na(molTrait)) %>% group_by(h2AA, molTrait) %>%
         mutate(weight = 1 / n()), 
       aes(x = Freq)) +
  facet_grid(h2AA~molTrait) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{AA}$"), breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Molecular trait", breaks = NULL, labels = NULL)) +
  labs(y = "Relative frequency across replicates", x = "Allele frequency (p)") +
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

ggsave("h2AA_net.png", height = 6, width = 10)

ggplot(d_com_end %>% filter(modelindex == 1) %>% drop_na(h2A, molTrait) %>% group_by(h2A, molTrait) %>%
         mutate(weight = 1 / n()), 
       aes(x = Freq)) +
  facet_grid(h2A~molTrait) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{A}$"), breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Molecular trait", breaks = NULL, labels = NULL)) +
  labs(y = "Relative frequency\nacross replicates", x = "Allele frequency (p)") +
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

ggsave("h2A_net.png", height = 6, width = 10)

ggplot(d_com_end %>% filter(modelindex == 1) %>% drop_na(h2D, molTrait) %>% group_by(h2D, molTrait) %>%
         mutate(weight = 1 / n()), 
       aes(x = Freq)) +
  facet_grid(h2D~molTrait) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{D}$"), breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Molecular trait", breaks = NULL, labels = NULL)) +
  labs(y = "Relative frequency\nacross replicates", x = "Allele frequency (p)") +
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

ggsave("h2D_net.png", height = 6, width = 10)


# Additive time
h2A_add <- ggplot(d_com_end %>% filter(modelindex == 0, !is.na(h2A)) %>% group_by(h2A) %>%
         mutate(weight = 1 / n()), 
       aes(x = Freq)) +
  facet_grid(h2A~.) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{A}$"), breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency (p)") + 
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

h2D_add <- ggplot(d_com_end %>% filter(modelindex == 0, !is.na(h2D)) %>% group_by(h2D) %>%
         mutate(weight = 1 / n()), 
       aes(x = Freq)) +
  facet_grid(h2D~.) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{D}$"), breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency (p)") + 
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

h2AA_add <- ggplot(d_com_end %>% filter(modelindex == 0, !is.na(h2AA)) %>% group_by(h2AA) %>%
         mutate(weight = 1 / n()), 
       aes(x = Freq)) +
  facet_grid(h2AA~.) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{AA}$"), breaks = NULL, labels = NULL)) +
  labs(x = "Allele frequency (p)") + 
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

plot_add <- plot_grid(h2A_add + ylab(NULL), 
                      h2D_add + ylab(NULL), 
                      nrow = 1) 


y.grob <- textGrob("Relative frequency across replicates", 
                   gp = gpar(fontsize = 16), rot = 90)

grid.arrange(arrangeGrob(plot_add, left = y.grob))
ggsave("h2A_D_add.png", height = 7, width = 8)

h2AA_add + labs(x = "Allele frequency (p)", y = NULL)
ggsave("h2AA_add.png", height = 3, width = 4)

# Effect size distributions
ggplot(d_com_end %>% filter(modelindex == 1, !is.na(h2A), !is.na(molTrait)) 
       %>% group_by(h2A, molTrait) %>% mutate(weight = 1 / n()), 
       aes(x = log(value))) +
  facet_grid(h2A~molTrait) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{A}$"), breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Molecular trait", breaks = NULL, labels = NULL)) +
  labs(y = "Relative frequency across replicates", x = "Molecular effect") +
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

ggsave("h2A_val_net.png", height = 6, width = 10)

ggplot(d_com_end %>% filter(modelindex == 1, !is.na(h2D), !is.na(molTrait)) 
       %>% group_by(h2D, molTrait) %>% mutate(weight = 1 / n()), 
       aes(x = log(value))) +
  facet_grid(h2D~molTrait) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{D}$"), breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Molecular trait", breaks = NULL, labels = NULL)) +
  labs(y = "Relative frequency across replicates", x = "Molecular effect") +
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

ggsave("h2D_val_net.png", height = 6, width = 10)

ggplot(d_com_end %>% filter(modelindex == 1, !is.na(h2AA), !is.na(molTrait)) 
       %>% group_by(h2AA, molTrait) %>% mutate(weight = 1 / n()), 
       aes(x = log(value))) +
  facet_grid(h2AA~molTrait) +
  geom_histogram(aes(weight = weight), bins = 40) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     sec.axis = sec_axis(~ . , name = TeX("$h^2_{AA}$"), breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Molecular trait", breaks = NULL, labels = NULL)) +
  labs(y = "Relative frequency across replicates", x = "Molecular effect") +
  theme_bw() + theme(text = element_text(size = 16), panel.spacing = unit(1, "lines"))

ggsave("h2AA_val_net.png", height = 6, width = 10)
