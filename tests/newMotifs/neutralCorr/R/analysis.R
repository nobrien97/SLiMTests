library(tidyverse)
library(ggh4x)

# Helper functions
ModelFromIndex <- function(id) {
    switch(id,
    1 = return("NAR"),
    2 = return("PAR"), 
    3 = return("FFLC1"),
    4 = return("FFLI1"),
    5 = return("FFBH"))
}

# Load in data
DATA_PATH <- "/mnt/d/SLiMTests/tests/newMotifs/neutralCorr/slim_qg.csv"
d_qg <- read_csv(DATA_PATH, col_Names = c("gen", "seed", "modelindex", "meanH", paste0("phenomean_", 1:4),
                                  paste0("phenovar_", 1:4), paste0("phenocor_", 1:4),
                                  paste0("phenocor_", c(12, 13, 14, 23, 24, 34))))


d_qg_sum <- d_qg %>%
    mutate(model = ModelFromIndex(modelindex)) %>%
    group_by(modelindex) %>%
    summarise_at(vars(starts_with("phenocor")), list(mean = mean, var = var)) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with("phenocor"),
                 names_to = "traitCombo",
                 values_to = "cor", names_prefix = "phenocor_")

ggplot(d_qg_sum, 
    aes(x = traitCombo, y = cor, colour = model)) +
    geom_point(position = position_dodge(0.9)) +
    geom_errorbar(aes(ymin = cor - cor_var, ymax = cor + cor_var), width = 0.2) +
    labs(x = "Trait combination", y = "Mean correlation", colour = "Model") +
    theme_bw() +
    theme(text = element_text(size = 12))
