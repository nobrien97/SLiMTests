library(tidyverse)
library(paletteer)

setwd("/mnt/c/GitHub/SLiMTests/tests/newMotifs/R")

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

d_times <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/slim_time.csv", col_names = F)
d_combos <- read_csv("/mnt/c/GitHub/SLiMTests/tests/newMotifs/R/combos.csv", col_names = F)

colnames(d_times) <- c("gen", "seed", "modelindex", "time")
colnames(d_combos) <- "model"

d_times$modelindex <- factor(d_times$modelindex)

AddCombosToDF <- function(df) {
  df %>% ungroup() %>%
    mutate(model = d_combos$model[as.numeric(levels(modelindex))[modelindex]])
}

d_times <- AddCombosToDF(d_times)


d_times %>%
  group_by(seed, modelindex, model) %>%
  mutate(delta = time - lag(time)) -> d_times

d_times %>%
  group_by(gen, model) %>%
  summarise(meanTime = mean(time),
            seTime = se(time),
            meanDelta = mean(delta),
            seDelta = se(delta)) -> time_summary

ggplot(time_summary, aes(x = gen, y = meanTime, colour = as.factor(model))) +
  geom_line() +
  geom_ribbon(aes(ymin = meanTime - seTime, ymax = meanTime + seTime, 
                  fill = as.factor(model)),
              colour = NA, alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Model", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  guides(fill = F) +
  labs(x = "Generation", y = "Time (s)", colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
ggsave("times.png")

# extrapolate to 60000 generations
time_summary %>%
  group_by(model) %>%
  filter(gen == 1000) %>%
  select(!c(meanDelta, seDelta, seTime, seDelta)) %>%
  mutate(wholeSimEstimate = paste(((meanTime * (60000/gen))/3600), "hrs"))
