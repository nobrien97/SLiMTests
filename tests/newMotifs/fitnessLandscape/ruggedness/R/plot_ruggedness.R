library(tidyverse)
library(paletteer)
library(ggbeeswarm)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}


DATA_PATH <- "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/fitnessLandscape/ruggedness/R/"
# Load data
RugRes <- data.table::fread(paste0(DATA_PATH, "d_ruggedness_permolcomp.csv"), header = F,
                               col.names = c("model", "startW", "endW", "netChangeW",
                                             "sumChangeW", "numFitnessHoles", "molComp",
                                             "bkg"))
RugRes$ruggedness <- RugRes$netChangeW - RugRes$sumChangeW



# Average over genetic backgrounds
d_ruggedness <- RugRes %>%
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
            

d_ruggedness <- d_ruggedness %>%
  group_by(model, molComp, replicate) %>% # replicate: walk direction is the same, averaged across backgrounds
  summarise(meanRuggedness = mean(ruggedness),
            seRuggedness = se(ruggedness),
            meanFitnessHoles = mean(numFitnessHoles),
            seFitnessHoles = se(numFitnessHoles))

# Plot ruggedness
ggplot(d_ruggedness %>% 
         mutate(model = as_factor(model),
                molComp = as_factor(molComp),
                meanRuggedness = abs(meanRuggedness)), 
       aes(x = model, y = meanRuggedness, colour = model)) +
  facet_wrap(.~molComp) +
  geom_pointrange(aes(ymin = meanRuggedness - seRuggedness, ymax = meanRuggedness + seRuggedness),
                linetype = "dotted", position = position_jitter(width = 0.5),
                shape = 1, size = 0.2) +
  theme_bw() +
  scale_colour_paletteer_d("nationalparkcolors::Everglades") +
  labs(x = "Model", y = "Landscape ruggedness",
       colour = "Model") +
  theme(text = element_text(size = 10), legend.position = "none",
        panel.spacing.x = unit(1.0, "lines"))
ggsave("plt_landscaperuggedness_permolcomp.png", device = png, width = 9, height = 4)

# Ruggedness driven by mutations in base

# Plot number of fitness holes in the walk
ggplot(RugRes %>% filter(numFitnessHoles > 0) %>%
         mutate(model = as_factor(model)), 
       aes(x = numFitnessHoles)) +
  facet_grid(molComp~model) +
  geom_histogram(bins = 11) +
  theme_bw() +
  scale_colour_paletteer_c("grDevices::Viridis") +
  labs(x = "Number of fitness holes during the random walk", y = "Count",
       colour = "Landscape Ruggedness") +
  theme(text = element_text(size = 14), legend.position = "bottom",
        panel.spacing.x = unit(1.5, "lines"))
ggsave("plt_landscapeholeyness_permolcomp.png", device = png)

# Number of fitness holes driven by mutations in aZ for both FFL-I1 and FFBH

