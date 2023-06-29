library(tidyverse)
library(paletteer)
library(latex2exp)
library(GGally)
library(cowplot)
library(ggnewscale)
library(ggalt)
library(ggarrow)
library(gghalves)
library(ggstance)
library(ggridges)
library(ggpmisc)

setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("wrangle_data.R")

# Fig 2 - phenotype mean
ggplot(d_adapted_sum %>% 
         filter(gen > 49000) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = meanPheno, colour = modelindex)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno, 
                  fill = modelindex), alpha = 0.2, colour = NA
              ) +
  scale_y_continuous(limits = c(0.9, 2.15), breaks = seq(1, 2, by = 0.25)) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  scale_fill_paletteer_d("ggsci::nrc_npg", guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean population phenotype",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenomean
plt_phenomean
ggsave("phenomean.png", plt_phenomean, png)

# phenotype median
ggplot(d_adapted_sum %>% 
         filter(gen > 49000) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = medPheno, colour = modelindex)) +
  geom_line(linewidth = 0.7) +
  scale_y_continuous(limits = c(0.9, 2.15), breaks = seq(1, 2, by = 0.25)) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  scale_fill_paletteer_d("ggsci::nrc_npg", guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Median population phenotype",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenomed
plt_phenomed
ggsave("phenomed.png", plt_phenomed, png)

# phenotype mode
ggplot(d_adapted_sum %>% 
         filter(gen > 49000) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = modePheno, colour = modelindex)) +
  geom_line(linewidth = 0.7) +
  scale_y_continuous(limits = c(0.9, 2.15), breaks = seq(1, 2, by = 0.25)) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  scale_fill_paletteer_d("ggsci::nrc_npg", guide = "none") +
  labs(x = "Generations post-optimum shift", y = "Mode of population\nphenotype distribution",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenomode
plt_phenomode
ggsave("phenomode.png", plt_phenomode, png)



d_qg %>% group_by(modelindex) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted))


# Fig 3 - effect sizes

## Additive
ggplot(d_fix_adapted %>% filter(modelindex == 1), 
       aes(x = value)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", linewidth = 0.8) +
  labs(x = "Phenotypic effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_phenotype_additive
dfe_phenotype_additive

ggsave("dfe_phenotype_additive.png", dfe_phenotype_additive, png)


ggplot(d_fix_add %>% filter(modelindex == 1), 
       aes(x = s)) +
  geom_density(show.legend = FALSE, linewidth = 0) +
  stat_density(geom="line", position="identity", size = 0.8) +
  labs(x = "Fitness effect", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> dfe_fitness_additive
dfe_fitness_additive

ggsave("dfe_fitness_additive.png", dfe_fitness_additive, png)



ggplot(d_fix_nar, 
       aes(x = avFX, colour = mutType)) +
  geom_density(show.legend = FALSE, size = 0) +
  stat_density(geom="line", position="identity", linewidth = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Phenotypic effect", y = "Density",
       colour = "Molecular component") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5), 
        legend.position = "bottom") -> dfe_phenotype_nar
dfe_phenotype_nar
ggsave("dfe_phenotype_nar.png", dfe_phenotype_nar, png)

ggplot(d_fix_nar, 
       aes(x = s, colour = mutType)) +
  geom_density(show.legend = FALSE, linewidth = 0) +
  stat_density(geom="line", position="identity", linewidth = 0.8) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Fitness effect", y = "Density",
       colour = "Molecular\ncomponent") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = margin(5.5, 10, 5.5, 5.5),
        legend.position = "none") -> dfe_fitness_nar
dfe_fitness_nar
ggsave("dfe_fitness_nar.png", dfe_fitness_nar, png)

# Combine the above figures
leg <- get_legend(dfe_phenotype_nar)

fig3 <- plot_grid(dfe_fitness_additive + ylim(c(0, 31.37901)),
          dfe_fitness_nar + theme(legend.position = "none") + ylim(c(0, 31.37901)),
          ncol = 2, labels = "AUTO")

# phenotype in supp figures
sfig1 <- plot_grid(dfe_phenotype_additive,
                   dfe_phenotype_nar + theme(legend.position = "none"),
                   ncol = 2, labels = "AUTO")

fig3 <- plot_grid(fig3, get_legend(dfe_phenotype_nar), ncol = 1, rel_heights = c(1, 0.1))
fig3
ggsave("dfe_combined_fitness.png", fig3, width = 8, height = 5, png, bg = "white")

# Fig 4: Plot adaptive walks - use d_fix_ranked 

plotPairwiseScatter <- function(dat, x, y, labels) {
  GRID_RES <- 200
  d_grid <- expand.grid(seq(from = min(dat[[x]]), 
                                   to = max(dat[[x]]), length.out = GRID_RES), 
                        seq(from = min(dat[[y]]), 
                            to = max(dat[[y]]), length.out = GRID_RES), 1, 1)
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  cc_arrow <- paletteer_d("RColorBrewer::Reds", 9)[c(3, 6, 9)]
  
  minFit <- 0.9 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  dat_arrow <- dat %>%
    filter(seed %in% sampled_seed) %>%
    arrange(rank) %>%
    group_by(seed) %>% mutate(xend = lead(aZ),
                              yend = lead(bZ))

  suppressWarnings(
  ggplot(dat %>% filter(modelindex == 2), aes(x = .data[[x]], y = .data[[y]])) +
    geom_raster(data = d_landscape, mapping = 
                     aes(x = .data[[x]], y = .data[[y]], fill = fitness)) +
    # geom_point(aes(colour = as.factor(rank)), size = 2) +
    #geom_encircle(s_shape = 0.2, expand = 0.2, colour = "black", mapping = aes(group = seed)) +
    scale_fill_gradientn(colors = c(cc[1], cc), 
                           limits = c(ifelse(minFit < 0.8, 0.8, minFit), 1),
                           values = wValues, na.value = cc[1]) +
    labs(fill = "Fitness (w)") +
    # geom_segment(data = dat_arrow,
    #              mapping = aes(x = .data[[x]], y = .data[[y]],
    #                            xend = xend, yend = yend,
    #                            group = seed,
    #                            colour = as.factor(rank)),
    #              arrow = arrow(length = unit(0.5, "cm")), linetype = "dashed",
    #              size = 1, show.legend = F) +
    scale_colour_manual(values = cc_arrow) +
    scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2], colour = "Adaptive step") +
    theme_bw() + 
    theme(legend.position = "bottom", text = element_text(size = 14)) +
    guides(
           fill = guide_colourbar(barwidth = 10, title.vjust = 0.87)) # i love magic numbers
  )
}

# d_qg_adapting <- d_adapted %>% filter(gen >= 50000)
# d_qg_adapting$gen_width <- scales::rescale(d_qg_adapting$gen, to = c(0.1, 1))

seed <- sample(0:.Machine$integer.max, 1)
set.seed(seed)
#set.seed(1098129538)
sampled_seed <- sample(d_fix_ranked[d_fix_ranked$rank == 2,]$seed, 3)
walk <- plotPairwiseScatter(d_fix_ranked %>% filter(seed %in% sampled_seed), 
                    "aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
suppressWarnings(walk)
ggsave("walk_landscape.png", walk, png)

# Ranked fixed effects
d_rank_av_nar <- d_fix_ranked %>%
  group_by(rank) %>%
  summarise(CIFit = CI(s),
            meanFit = mean(s),
            medFit = median(s),
            modeFit = estimate_mode(s, T),
            CIPheno = CI(phenomean),
            meanPheno = mean(phenomean),
            medPheno = median(phenomean),
            modePheno = estimate_mode(phenomean, T),
            meanDom = mean(h),
            meanGen = mean(gen),
            medGen = median(gen),
            modeGen = estimate_mode(gen, T),
            CIGen = CI(gen),
            CIDom = CI(h),
            nInRank = n())

ggplot(d_rank_av_nar %>% filter(rank > 0), aes(x = rank, y = meanFit)) +
  geom_point(position = position_dodge(1)) +
  geom_line(position = position_dodge(1)) +
  scale_y_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by = 0.005)) +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4, position = position_dodge(1)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Mean fitness effect", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize
plt_adaptivestepsize

d_rank_av_add <- d_fix_ranked_add %>%
  group_by(rank) %>%
  summarise(CIFit = CI(s),
            meanFit = mean(s),
            medFit = median(s),
            modeFit = estimate_mode(s, T),
            CIPheno = CI(phenomean),
            meanPheno = mean(phenomean),
            medPheno = median(phenomean),
            modePheno = estimate_mode(phenomean, T),
            meanGen = mean(gen),
            medGen = median(gen),
            modeGen = estimate_mode(gen, T),
            CIGen = CI(gen),
            CIDom = CI(h),
            meanDom = mean(h),
            nInRank = n())

ggplot(d_rank_av_add %>% filter(rank > 0), aes(x = rank, y = meanFit)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4) +
  labs(x = "Adaptive step", y = "Mean fitness effect") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize_add
plt_adaptivestepsize_add

# diminishing returns of each step in the walk
ggplot(d_rank_av_nar, aes(x = rank, y = meanPheno)) +
  geom_point(position = position_dodge(1)) +
  geom_line(position = position_dodge(1)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.9, 2.1), breaks = seq(1, 2, by = 0.25)) +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4, position = position_dodge(1)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Mean phenotype", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno

ggplot(d_rank_av_add, aes(x = rank, y = meanPheno)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.9, 2.1), breaks = seq(1, 2, by = 0.25)) +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4) +
  labs(x = "Adaptive step", y = "Mean phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno_add

# combined diminishing returns plot
d_rank_av <- rbind(d_rank_av_nar %>% mutate(model = "NAR"), 
                   d_rank_av_add %>% mutate(model = "Additive"))

ggplot(d_rank_av, aes(x = rank, y = meanPheno, colour = model)) +
  geom_point(position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_text(data = d_rank_av %>% filter(rank > 0), mapping = aes(x = rank, y = 2.1, colour = model,
                          label = paste0("n = ", nInRank)), size = 5,
            position = position_dodgev(height = 0.1), show.legend = F) +
  scale_y_continuous(limits = c(0.9, 2.15), breaks = seq(1, 2, by = 0.25)) +
  geom_errorbar(aes(ymin = meanPheno - CIPheno, ymax = meanPheno + CIPheno), 
                width = 0.4, position = position_dodge(0.3)) +
  labs(x = "Adaptive step", y = "Mean population phenotype", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostep
plt_phenostep
ggsave("phenostep.png", plt_phenostep, device = png, bg = "white")

ggplot(d_rank_av %>% filter(rank > 0), aes(x = rank, y = meanFit, colour = model)) +
  geom_point(position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = meanFit - CIFit, ymax = meanFit + CIFit), 
                width = 0.4, position = position_dodge(0.3)) +
  labs(x = "Adaptive step", y = "Mean fitness effect (s)", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_fitstep
plt_fitstep
ggsave("fitstep.png", plt_fitstep, device = png, bg = "white")

leg <- get_legend(plt_phenomean + theme(legend.position = "bottom"))

# grid with phenomean
plot_grid(plt_phenostep + theme(legend.position = "none",
                                plot.margin = margin(10, 10, 10, 10)),
          plt_phenomean + theme(legend.position = "none",
                                plot.margin = margin(10, 10, 10, 10)),
          ncol = 2, labels = "AUTO") -> plt_meanstep 

plot_grid(plt_meanstep,
          leg, ncol = 1, rel_heights = c(1, 0.1))
ggsave("phenomean_step.png", width = 10, height = 5.625, device = png, bg = "white")


plot_grid(plt_adaptivestepsize_add,
          plt_adaptivestepsize,
          plt_phenostepPheno_add,
          plt_phenostepPheno, 
          nrow = 2, labels = "AUTO") -> plt_steps
plt_steps
ggsave("stepsize_combined.png", plt_steps, png, bg = "white")

# Stick mean fixation times onto plt_phenomean
plt_phenomean +
  geom_point(data = d_rank_av %>% mutate(meanGen = meanGen - 50000,
                                         modelindex = ifelse(model == "NAR", 2, 1)), 
             mapping = aes(x = meanGen, y = meanPheno, colour = as.factor(modelindex)),
             size = 2) +
  geom_text(data = d_rank_av %>% filter(model == "Additive") %>% 
              mutate(meanGen = meanGen - 50000,
                     modelindex = ifelse(model == "NAR", 2, 1)), 
            mapping = aes(x = meanGen, y = meanPheno, colour = NULL,
                          label = paste("Step", as.factor(rank))),
            vjust = 0, nudge_y = -0.09, show.legend = F) +
  geom_errorbarh(data = d_rank_av %>% mutate(meanGen = meanGen - 50000,
                                            modelindex = ifelse(model == "NAR", 2, 1)),
                mapping = aes(x = NULL, y = meanPheno, xmin = meanGen - CIGen,
                              xmax = meanGen + CIGen,
                              colour = as.factor(modelindex)), width = 0.05) +
  geom_errorbar(data = d_rank_av %>% mutate(meanGen = meanGen - 50000,
                                            modelindex = ifelse(model == "NAR", 2, 1)),
                mapping = aes(x = meanGen, y = meanPheno, ymin = meanPheno - CIPheno,
                              ymax = meanPheno + CIPheno,
                              colour = as.factor(modelindex))) +
  theme(text = element_text(size = 12)) -> plt_phenomeanwalk

plt_phenomed +
  geom_point(data = d_rank_av %>% mutate(medGen = medGen - 50000,
                                         modelindex = ifelse(model == "NAR", 2, 1)), 
             mapping = aes(x = medGen, y = medPheno, colour = as.factor(modelindex)),
             size = 2) +
  geom_text(data = d_rank_av %>% filter(model == "Additive") %>% 
              mutate(medGen = medGen - 50000,
                     modelindex = ifelse(model == "NAR", 2, 1)), 
            mapping = aes(x = medGen, y = medPheno, colour = NULL,
                          label = paste("Step", as.factor(rank))),
            vjust = 0, nudge_y = -0.09, show.legend = F) +
  theme(text = element_text(size = 12)) -> plt_phenomedwalk


plt_phenomode +
  geom_point(data = d_rank_av %>% mutate(modeGen = modeGen - 50000,
                                         modelindex = ifelse(model == "NAR", 2, 1)), 
             mapping = aes(x = modeGen, y = modePheno, colour = as.factor(modelindex)),
             size = 2) +
  geom_text(data = d_rank_av %>% filter(model == "Additive") %>% 
              mutate(modeGen = modeGen - 50000,
                     modelindex = ifelse(model == "NAR", 2, 1)), 
             mapping = aes(x = modeGen, y = modePheno, colour = NULL,
                           label = paste("Step", as.factor(rank))),
            vjust = 0, nudge_y = -0.09, show.legend = F) +
  theme(text = element_text(size = 12)) -> plt_phenomodewalk
plt_phenomeanwalk
plt_phenomedwalk
plt_phenomodewalk
ggsave("phenomeanwalk.png", plt_phenomeanwalk, device = png)
ggsave("phenomedwalk.png", plt_phenomedwalk, device = png)
ggsave("phenomodewalk.png", plt_phenomodewalk, device = png)

# combine them
plot_grid(plt_phenomeanwalk, plt_phenomedwalk, plt_phenomodewalk,
          nrow = 3, labels = "AUTO")
ggsave("combinedphenowalks.png", height = 10, device = png)

ggplot(d_rank_av_stat %>% mutate(medGen = medGen - 50000,
                            modelindex = ifelse(model == "NAR", 2, 1)), 
       mapping = aes(x = gen, y = medPheno, colour = as.factor(modelindex)),) +
  geom_point(size = 2) +
  geom_text(data = d_rank_av %>% filter(model == "Additive") %>% 
              mutate(medGen = medGen - 50000,
                     modelindex = ifelse(model == "NAR", 2, 1)), 
            mapping = aes(x = medGen, y = medPheno, colour = NULL,
                          label = paste("Step", as.factor(rank))),
            vjust = 0, nudge_y = -0.09, show.legend = F) -> plt_phenomedwalk


# Now plot all the walks
plt_phenomean +
  geom_point(data = d_fix_ranked_combined %>% filter(rank > 0) %>% 
               mutate(gen = gen - 50000), 
             mapping = aes(x = gen, y = phenomean, colour = modelindex, 
                           shape = as.factor(rank)),
             size = 1.7) +
  geom_line(data = d_fix_ranked_combined %>% mutate(gen = gen - 50000), 
            mapping = aes(x = gen, y = phenomean, 
                          group = interaction(seed, modelindex),
                          colour = modelindex),
            linewidth = 0.1) +
  labs(shape = "Adaptive step") -> plt_phenowalks
ggsave("phenowalks.png", plt_phenowalks, device = png)

# as box plots
d_fix_ranked_combined$gen_group <- cut(d_fix_ranked_combined$gen - 50000, breaks = 10)

ggplot(d_fix_ranked_combined %>% filter(rank > 0),
       aes(x = as.factor(gen_group), y = phenomean, colour = model)) +
  geom_boxplot(position = position_dodge(1)) +
  labs(x = "Generations post-optimum shift", y = "Phenotype mean", 
       shape = "Adaptive step", colour = "Model") +
  theme_bw()

# Ridgeline works better for showing change in the distribution I think
d_adapted_walk <- d_adapted %>% filter(gen > 49000)
breaks <- seq(min(d_adapted_walk$gen - 50000), max(d_adapted_walk$gen - 50000), by = 1000)

d_adapted_walk$gen_group <- breaks[findInterval(d_adapted_walk$gen - 50000, breaks, rightmost.closed = TRUE)]

ggplot(d_adapted_walk,
       aes(y = as.factor(gen_group), x = phenomean, fill = modelindex)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(y = "Generations post-optimum shift", x = "Phenotype mean", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenomean_dist
plt_phenomean_dist
ggsave("phenomean_dist.png", plt_phenomean_dist, device = png)


ggplot(d_fix_ranked_combined %>% filter(rank > 0),
       aes(y = as.factor(rank), x = s, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivewalk_dist
plt_adaptivewalk_dist

ggplot(d_fix_ranked_combined %>% filter(rank > 0),
       aes(y = as.factor(rank), x = gen - 50000, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Generations post-optimum shift", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepgen_dist
plt_adaptivestepgen_dist


# distributions of phenotypes and fitness effects - departures from normality?
ggplot(d_fix_ranked_combined,
       aes(x = phenomean, fill = model)) +
  facet_grid(.~as.factor(rank)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Adaptive step", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Mean phenotype", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_pheno_dist
         
ggplot(d_fix_ranked_combined,
       aes(x = s, fill = model)) +
  facet_grid(.~as.factor(rank)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Adaptive step", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_s_dist

plot_grid(plt_pheno_dist, plt_s_dist, nrow = 2, labels = "AUTO") -> plt_pheno_s_dist
ggsave("pheno_s_dist.png", plt_pheno_s_dist, device = png)

cc3 <- paletteer_d("ggprism::viridis", 4)[c(1, 4)]
cc3

# Box plots with segregating sites as well
d_segFixRat_sum$model <- "NAR"
d_segFixRat_add_sum$model <- "Additive"
d_segFixRat_sum_combined <- rbind(d_segFixRat_sum, d_segFixRat_add_sum)
d_seg_ranked$model <- "NAR"
d_seg_ranked_add$model <- "Additive"
d_seg_ranked_combined <- rbind(d_seg_ranked, d_seg_ranked_add)

ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = as.factor(rank), y = s)) +
  facet_grid(.~model) +
  geom_half_boxplot(side = "l", center = T, width = 0.5, colour = cc3[1]) +
  geom_half_boxplot(side = "r", center = T, width = 0.5, 
                    data = d_seg_ranked_combined %>% distinct() %>% filter(rank > 0), 
                   mapping = aes(x = as.factor(rank), y = s), colour = cc3[2]) +
  geom_text(d_segFixRat_sum_combined, 
            mapping = aes(x = as.factor(rank), y = -1,
                          label = paste0(signif(meanPercFix * 100, 3), 
                                         " ± ", signif(CIPercFix * 100, 3), "%")),
            size = 3) +
  labs(x = "Adaptive step", y = "Fitness effect (s)") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_adaptivestepsize_bp
plt_adaptivestepsize_bp

# Distribution of all effects - should be normal! (Orr 2006)
d_ranked_combined <- rbind(d_fix_ranked_combined, d_seg_ranked_combined) %>% filter(!is.na(s))
d_ranked_combined$model <- if_else(d_ranked_combined$modelindex == 1, "Additive", "NAR")
ggplot(d_ranked_combined, aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_distallfx
plt_distallfx
ggsave("dist_allfx.png", device = png)

# Distribution of all beneficial effects - should be exponential (if Gumbel)
ggplot(d_ranked_combined %>% filter(s > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_distbenfx_alltimes
plt_distbenfx_alltimes
ggsave("dist_benfx_alltimes.png", device = png)

# Dist of fixed muts - should become less exponential over time
d_fix_ranked_combined$model <- if_else(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")
ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = s, y = as.factor(rank), fill = model)) +
  geom_density_ridges(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_distfixed_time
plt_distfixed_time
ggsave("dist_fixed_time.png", device = png)

ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_distfixed
plt_distfixed
ggsave("dist_fixed.png", device = png)



# distribution of generations the step happens at

d_fix_ranked_combined$model <- ifelse(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")

ggplot(d_fix_ranked_combined %>% filter(rank > 0), 
       aes(x = as.factor(rank), y = (gen - 50000), 
           colour = model)) +
  geom_boxplot(position = position_dodge2(0.8, preserve = "single"), alpha = 0.1) +
  geom_point(data = d_rank_av %>% mutate(meanGen = meanGen - 50000) %>% filter(rank > 0), 
             mapping = aes(y = meanGen), size = 2,
             position = position_dodge(0.8)) +
  geom_errorbar(data = d_rank_av %>% mutate(meanGen = meanGen - 50000) %>% filter(rank > 0),
                mapping = aes(y = meanGen, ymin = meanGen - CIGen, ymax = meanGen + CIGen), 
                position = position_dodge(0.8), width = 0.3) +
  geom_line(data = d_rank_av %>% mutate(meanGen = meanGen - 50000) %>% filter(rank > 0), 
            mapping = aes(y = meanGen, group = model), 
            position = position_dodge(0.8), linewidth = 1) +
  labs(x = "Adaptive step", y = "Average fixation generation", colour = "Model") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  scale_y_continuous(label = scales::comma) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.position = "bottom") -> plt_gensteptime
plt_gensteptime

# Hack together a legend for the grid
plt_legend <- ggplot(data.frame(x = 1,2, 
                                y = 1,2,
                                mutation = c("Fixed","Segregating")), 
                     aes(x = x, y = y, colour = mutation)) +
  geom_boxplot() +
  scale_colour_manual(values = cc3) +
  labs(colour = "Mutation type") + 
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
plt_legend

leg <- get_legend(plt_legend)

plot_grid(plt_adaptivestepsize_bp + ylim(c(-1, 0.45)), 
          leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_steps_fit

plot_grid(plt_steps_fit, 
          plt_gensteptime,
          nrow = 1, labels = "AUTO", rel_widths = c(1.5, 1)) -> plt_steps_fit
plt_steps_fit

ggsave("plt_steps_fit.png", plt_steps_fit, 
       width = 14, height = 5, png, bg = "white")

# diminishing returns of each step in the walk
ggplot(d_fix_ranked, aes(x = as.factor(rank), y = phenomean)) +
  geom_boxplot() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.64, 2.6), breaks = seq(1, 2.5, by = 0.5)) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Mean phenotype", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno_bp
plt_phenostepPheno_bp

ggplot(d_fix_ranked_add, aes(x = as.factor(rank), y = phenomean)) +
  geom_boxplot() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_y_continuous(limits = c(0.64, 2.6), breaks = seq(1, 2.5, by = 0.5)) +
  labs(x = "Adaptive step", y = "Mean phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_phenostepPheno_add_bp
plt_phenostepPheno_add_bp

plot_grid(plt_adaptivestepsize_add_bp,
          plt_adaptivestepsize_bp,
          plt_phenostepPheno_add_bp,
          plt_phenostepPheno_bp, 
          nrow = 2, labels = "AUTO") -> plt_steps_bp
plt_steps_bp
ggsave("stepsize_combined_boxplot.png", plt_steps_bp, png, bg = "white")

# Dominance
rbind(d_fix_ranked, d_fix_ranked_add) %>%
  group_by(rank, modelindex) %>%
  summarise(meanh = mean(h),
            CIh = CI(h))
ggplot(d_fix_ranked %>% filter(rank > 0, h < 1e+9), aes(x = as.factor(rank), y = h)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Dominance coefficient (h)", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_dom_fixed_bp
plt_dom_fixed_bp

ggplot(d_seg_ranked %>% filter(rank > 0, h < 1e+9, h > -1000), aes(x = as.factor(rank), y = h)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = mutType_names) +
  labs(x = "Adaptive step", y = "Dominance coefficient (h)", colour = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_dom_seg_bp
plt_dom_seg_bp


# Cumulative fitness effects for a few replicates
d_com_nar_sample %>% filter(gen > 49000) %>%
  group_by(gen, seed) %>%
  mutate(sumFit = sum(avFit)) -> d_com_nar_sample

d_com_add_sample %>% filter(gen > 49000) %>%
  group_by(gen, seed) %>%
  mutate(sumFit = sum(avFit)) -> d_com_add_sample

d_com_add_sample %>%
  ungroup(gen) %>%
  distinct(gen, .keep_all = T) %>%
  mutate(cumSumFit = cumsum(sumFit)) -> d_add_sample_sum

ggplot(d_com_sample_sum %>% distinct(), aes(x = gen, y = cumSumFit, colour = seed)) +
  geom_line() +
  theme_bw()

ggplot(d_add_sample_sum %>% distinct(), aes(x = gen, y = cumSumFit, colour = seed)) +
  geom_line() +
  theme_bw()

# combined ranks table
ggplot(d_rank_combined_tbl, 
       aes(x = rank, y = n, fill = mutType)) +
  facet_grid(.~isAdapted) +
  geom_col(position = position_dodge(0.9)) +
  geom_line(position = position_dodge(0.9)) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Adapted", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Adaptive step", y = "Count", fill = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16))


# effect size per step
ggplot(d_ranked_combined %>% drop_na(), 
       aes(x = as.factor(rank), y = value, fill = mutType)) +
  facet_grid(.~isAdapted) +
  geom_boxplot() +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Adaptive step", y = "Effect", fill = "Mutation type") +
  theme_bw() +
  theme(text = element_text(size = 16))

# total SFS at gen 50000
ggplot(d_muts %>% filter(gen == 50000), aes(x = Freq)) +
  facet_grid(.~modelindex) +
  geom_histogram() +
  theme_bw()

# Adaptive walk - ratio of alpha/beta
ggplot(d_ratio, aes(x = fixEffectSum_aZ, y = fixEffectSum_bZ, colour = aZbZ)) +
  geom_point() +
  scale_colour_paletteer_c("viridis::viridis") +
  theme_bw()

ggplot(d_ratio, aes(x = aZbZ, y = phenomean, colour = wAA)) +
  geom_point() +
  scale_colour_paletteer_c("viridis::viridis") +
  theme_bw()

# Fitness landscape aZbZ
plotRatioLandscape <- function(minRatio, maxRatio) {
  GRID_RES <- 50000
  rational <- MASS:::.rat(seq(minRatio, maxRatio, length.out = GRID_RES),
                          max.denominator = 20)$rat
  d_grid <- data.frame(aZ = rational[,1], 
                       bZ = rational[,2],
                       KZ = 1, 
                       KXZ = 1) %>% distinct()
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8) %>%
    filter(pheno > 0)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape  %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZbZ, y = pheno, colour = fitness)) +
      geom_point() +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                             limits = c(ifelse(minFit < 0.8, 0.8, minFit), 1),
                           values = wValues) +
      labs(x = TeX("$\\frac{\\alpha_Z}{\\beta_Z} ratio$"), y = "Phenotype", 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(colour = guide_colourbar(barwidth = 20))
  )
}

plt <- plotRatioLandscape(0.5, 3)
plt

plotaZbZLandscape <- function(minVal, maxVal) {
  GRID_RES <- 400
  d_grid <- expand.grid(seq(from = minVal, 
                            to = maxVal, length.out = GRID_RES), 
                        seq(from = minVal, 
                            to = maxVal, length.out = GRID_RES), 1, 1)
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZ, y = bZ, fill = fitness)) +
      geom_tile() +
      #geom_abline(slope = 1/1.27) +
      # geom_point(data = d_landscape %>% mutate(aZbZ = aZ/bZ) 
      #            %>% filter(aZbZ > 1.25, aZbZ < 1.35), size = 0.1, shape = 4) +
      scale_fill_gradientn(colors = c(cc[1], cc), 
                           limits = c(minFit, 1),
                             values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\alpha_Z$"), y = TeX("$\\beta_Z$"), 
           fill = "Fitness (w)", size = TeX("$\\frac{\\alpha_Z}{\\beta_Z} ratio$")) +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(fill = guide_colourbar(barwidth = 20))
  )
}

plotaZbZLandscape(0, 3)

# Plot difference in movement between alpha and beta
ggplot(d_molCompDiff,
       aes(x = molCompDiff)) +
  geom_density() +
  labs(x = TeX("Difference in molecular component evolution"), y = "Density") +
  theme_bw() + 
  theme(text = element_text(size = 16)) -> plt_molCompDiff
plt_molCompDiff
ggsave("molCompDiff.png", plt_molCompDiff, device = png)

# correlation between phenotype of fixations only with mean phenotype
ggplot(d_fix_ranked_combined %>% filter(rank > 0), 
       aes(x = AA_pheno, y = phenomean, colour = modelindex)) +
  geom_point(size = 0.5, shape = 1) +
  geom_abline(linetype = "dashed", colour = "#AAAAAA") +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2"))) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Phenotype with only fixed effects", y = "Mean population phenotype") +
  theme_bw() + 
  theme(text = element_text(size = 16))

