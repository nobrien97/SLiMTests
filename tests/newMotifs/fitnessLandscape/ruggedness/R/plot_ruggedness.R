library(tidyverse)
library(paletteer)
library(ggbeeswarm)
library(xtable)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

model_levels <- c("NAR", "PAR", "FFLC1", 
                 "FFLI1", "FFBH")
model_labels <- c("NAR", "PAR", "FFL-C1", "FFL-I1", "FFBH")

DATA_PATH <- "/mnt/c/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/ruggedness/R/"

RUG_DATA_PATH <- "/mnt/j/SLiMTests/tests/newMotifs/ruggedness/"


# Load data
RugRes <- data.table::fread(paste0(RUG_DATA_PATH, "d_ruggedness_permolcomp.csv"), header = F,
                               col.names = c("model", "startW", "endW", "netChangeW",
                                             "sumChangeW", "numFitnessHoles", "molComp",
                                             "bkg"))
RugRes$ruggedness <- RugRes$netChangeW - RugRes$sumChangeW



# Average over genetic backgrounds
d_ruggedness <- RugRes %>%
  mutate(model = factor(model, levels = model_levels)) %>%
  group_by(model, molComp, bkg) %>%
  mutate(replicate = row_number()) %>% # create replicate row
  ungroup() %>%
# Filter out irrelevant rows
  filter( ( model == "NAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) |
          ( model == "PAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) | 
          ( model == "FFLC1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                              "KXZ", "KY", "Hilln", "base", "XMult") ) |
          ( model == "FFLI1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                              "KXZ", "KY", "Hilln", "base", "XMult") ) |
            model == "FFBH" )
            

d_ruggedness_sum <- d_ruggedness %>%
  group_by(model, molComp, replicate) %>% # replicate: walk direction is the same, averaged across backgrounds
  summarise(meanRuggedness = mean(ruggedness),
            CIRuggedness = CI(ruggedness),
            meanFitnessHoles = mean(numFitnessHoles),
            CIFitnessHoles = CI(numFitnessHoles))

# To table
print(xtable(d_ruggedness_sum), include.rownames = F)



# Plot ruggedness
ggplot(d_ruggedness_sum %>% 
         mutate(model = as_factor(model),
                molComp = as_factor(molComp),
                meanRuggedness = abs(meanRuggedness)), 
       aes(x = model, y = meanRuggedness, colour = model)) +
  facet_wrap(.~molComp) +
  geom_pointrange(aes(ymin = meanRuggedness - CIRuggedness, ymax = meanRuggedness + CIRuggedness),
                linetype = "dotted", position = position_jitter(width = 0.5),
                shape = 1, size = 0.2) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3.5)) +
  scale_colour_paletteer_d("nationalparkcolors::Everglades") +
  labs(x = "Model", y = "Landscape ruggedness",
       colour = "Model") +
  theme(text = element_text(size = 10), legend.position = "none",
        panel.spacing.x = unit(1.0, "lines"))
ggsave("plt_landscaperuggedness_permolcomp.png", device = png, width = 9, height = 4)

# Ruggedness driven by mutations in base

# Plot number of fitness holes in the walk
ggplot(RugRes %>% filter(numFitnessHoles > 0) %>%
         filter( ( model == "NAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) |
                   ( model == "PAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) | 
                   ( model == "FFLC1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                                       "KXZ", "KY", "Hilln", "base", "XMult") ) |
                   ( model == "FFLI1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                                       "KXZ", "KY", "Hilln", "base", "XMult") ) |
                   model == "FFBH" ) %>%
         mutate(model = as_factor(model)), 
       aes(x = numFitnessHoles)) +
  facet_grid(molComp~model) +
  geom_histogram(bins = 11) +
  theme_bw() +
  scale_colour_paletteer_c("grDevices::Viridis") +
  labs(x = "Number of fitness holes during the random walk", y = "Count",
       colour = "Landscape Ruggedness") +
  theme(text = element_text(size = 10), legend.position = "bottom",
        panel.spacing.x = unit(1.5, "lines"))
ggsave("plt_landscapeholeyness_permolcomp.png", device = png, width = 6, height = 6)

# Number of fitness holes driven by mutations in aZ for both FFL-I1 and FFBH
# But comparably few fitness holes for the other models

# What about proportion? How many mutations of each type are generating fitness holes?
d_holes <- RugRes %>%
  mutate(model = factor(model, levels = model_levels)) %>%
  group_by(model, molComp, bkg) %>%
  mutate(replicate = row_number()) %>% # create replicate row
  ungroup() %>%
  # Filter out irrelevant rows
  filter( ( model == "NAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) |
            ( model == "PAR" & molComp %in% c("aZ", "bZ", "KZ", "KXZ", "Hilln", "base", "XMult") ) | 
            ( model == "FFLC1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                                "KXZ", "KY", "Hilln", "base", "XMult") ) |
            ( model == "FFLI1" & molComp %in% c("aZ", "aY", "bY", "bZ", "KZ", "KY", 
                                                "KXZ", "KY", "Hilln", "base", "XMult") ) |
            model == "FFBH" )

numSteps <- 10
d_holes <- d_holes %>%
  group_by(model, bkg, replicate) %>%
  mutate(nHolesAcrossMolComps = sum(numFitnessHoles)) %>%
  ungroup() %>%
  mutate(propHoles = numFitnessHoles / nHolesAcrossMolComps)


#TODO: calculate number of fitness holes across all possible steps in the walk
#propHolesAmongAllSteps = propHoles * nHolesAcrossMolComps / (numSteps * n()

# How many walks had fitness holes?
holes_per_model_molcomp <- d_holes %>%
  group_by(model, molComp) %>%
  summarise(propHolesTotal = sum(numFitnessHoles) / n())

d_holes_sum <- d_holes %>%
# filter to only the models with holes  
  filter(nHolesAcrossMolComps > 0) %>%
  group_by(model, molComp, replicate) %>%
  summarise(meanPropHoles = mean(propHoles),
            CIPropHoles = CI(propHoles))


ggplot(d_holes_sum %>% 
         mutate(model = as_factor(model),
                molComp = as_factor(molComp)), 
       aes(x = model, y = meanPropHoles, colour = model)) +
  facet_wrap(.~molComp) +
  geom_quasirandom(shape = 1, size = 0.2) +
  # geom_pointrange(aes(ymin = meanPropHoles - CIPropHoles, 
  #                     ymax = meanPropHoles + CIPropHoles),
  #                 linetype = "dotted", position = position_jitter(width = 0.5),
  #                 shape = 1, size = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_colour_paletteer_d("nationalparkcolors::Everglades") +
  labs(x = "Model", y = "Proportion of total fitness holes",
       colour = "Model") +
  theme(text = element_text(size = 10), legend.position = "none",
        panel.spacing.x = unit(1.0, "lines"))
ggsave("plt_landscape_holeyness_permolcomp_prop.png", device = png, width = 9, height = 4)

# Per model (average over molcomps and over genetic backgrounds)
d_holes_sum <- d_holes %>%
  group_by(model, replicate) %>%
  summarise(meanHoles = mean(numFitnessHoles),
            CIHoles = CI(numFitnessHoles))

ggplot(d_holes_sum %>% 
         mutate(model = as_factor(model)), 
       aes(x = model, y = meanHoles, colour = model)) +
  geom_quasirandom(shape = 1, size = 0.2) +
  theme_bw() +
  scale_colour_paletteer_d("nationalparkcolors::Everglades") +
  labs(x = "Model", y = "Mean number of fitness holes",
       colour = "Model") +
  theme(text = element_text(size = 12), legend.position = "none",
        panel.spacing.x = unit(1.0, "lines")) -> plt_holeyness
plt_holeyness
ggsave("plt_landscape_holeyness.png", device = png, width = 9, height = 4)

d_ruggedness_sum <- d_ruggedness %>%
  group_by(model, replicate) %>% # replicate: walk direction is the same, averaged across backgrounds
  summarise(meanRuggedness = mean(ruggedness),
            CIRuggedness = CI(ruggedness),
            meanFitnessHoles = mean(numFitnessHoles),
            CIFitnessHoles = CI(numFitnessHoles))

# Plot ruggedness
ggplot(d_ruggedness_sum %>% 
         mutate(model = as_factor(model),
                meanRuggedness = abs(meanRuggedness)), 
       aes(x = model, y = meanRuggedness, colour = model)) +
  geom_quasirandom(shape = 1, size = 0.2) +
  theme_bw() +
  scale_colour_paletteer_d("nationalparkcolors::Everglades") +
  labs(x = "Model", y = "Mean landscape ruggedness",
       colour = "Model") +
  theme(text = element_text(size = 12), legend.position = "none",
        panel.spacing.x = unit(1.0, "lines")) -> plt_ruggedness
plt_ruggedness
ggsave("plt_landscaperuggedness.png", device = png, width = 9, height = 4)

leg <- get_legend(plt_ruggedness)

plt_landscape <- plot_grid(plt_ruggedness + theme(legend.position = "none"), 
                        plt_holeyness + theme(legend.position = "none"),
                        nrow = 2, labels = "AUTO", label_size = 12)

plt_landscape <- plot_grid(plt_landscape,
                        leg, nrow = 2, rel_heights = c(1, 0.05))
plt_landscape
ggsave("plt_landscape.png", device = png, bg = "white",
       width = 5, height = 9)

# Tables

# Average across replicates
print(xtable(d_ruggedness %>%
               group_by(model, molComp) %>%
               summarise(meanRuggedness = mean(abs(ruggedness)),
                         CIRuggedness = CI(ruggedness),
                         meanFitnessHoles = mean(numFitnessHoles),
                         CIFitnessHoles = CI(numFitnessHoles),
                         count = n()), digits = 5
), include.rownames = F)

# Per model
print(xtable(d_ruggedness %>%
               group_by(model) %>%
               summarise(meanRuggedness = mean(abs(ruggedness)),
                         CIRuggedness = CI(ruggedness),
                         meanFitnessHoles = mean(numFitnessHoles),
                         CIFitnessHoles = CI(numFitnessHoles),
                         n = n()), digits = 5
), include.rownames = F)

