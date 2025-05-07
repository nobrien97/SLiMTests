library(tidyverse)
library(ggh4x)

setwd("/mnt/c/GitHub/SLiMTests/tests/standingVar/mutVar/R")
DATA_PATH <- "/mnt/d/SLiMTests/tests/standingVar/mutVar/"
COMBOS_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/R/"
R_PATH <- "/mnt/c/GitHub/SLiMTests/tests/standingVar/calcMutationStats/R/"
source(paste0(R_PATH, "helperFunctionsAndSetup.R"))


# Read data
d_mutvar <- data.table::fread(paste0(DATA_PATH, "slim_mutvar.csv"), header = F,
                            col.names = c("replicate", "seed", "modelindex", "variance"))


d_combos <- read.table(paste0(COMBOS_PATH, "combos.csv"), header = F,
                       col.names = c("nloci", "tau", "r", "model"))


d_mutvar <- d_mutvar %>%
  mutate(replicate = replicate %/% 2 + 1, # convert from 1 3 5 to 1 2 3 
         seed = as.factor(seed),
         modelindex = as.factor(modelindex))

d_mutvar <- AddCombosToDF(d_mutvar)

d_mutvar_sum <- d_mutvar %>%
  group_by(model, r) %>%
  summarise(meanVar = mean(log10(variance)),
            CIVar = CI(log10(variance)))

# Plot mutational variance
ggplot(d_mutvar, aes(x = model, y = log10(variance), colour = model)) +
  facet_nested("Recombination rate (log10)" + log10(r) ~ .) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_point(data = d_mutvar_sum, 
             aes(x = model, y = meanVar),
             shape = 3, size = 2, colour = "black", inherit.aes = F) +
  scale_colour_manual(values = paletteer_d("nationalparkcolors::Everglades", 3, 
                                           direction = -1),
                    labels = c("Additive", "K+", "K-")) +
  scale_x_discrete(labels = c("Additive", "K+", "K-")) +
  labs(x = "Model", 
       y = "Mutational variance (log10)") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "none")
ggsave("plt_vm.png", device = png, width = 5, height = 9)  
