library(tidyverse)
library(paletteer)

setwd("/mnt/c/GitHub/SLiMTests/tests/parallelTest/R")

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

d_times <- read_csv("/mnt/d/SLiMTests/tests/parallelTest/slim_time.csv", col_names = F)



colnames(d_times) <- c("gen", "seed", "modelindex", "time")
cores_treatments <- c(1, 2, 4, 6)
d_times %>% mutate(cores = cores_treatments[modelindex]) -> d_times

d_times %>%
  group_by(gen, modelindex, cores) %>%
  summarise(meanTime = mean(time),
            seTime = se(time)) -> time_summary

ggplot(time_summary, aes(x = gen, y = meanTime, colour = as.factor(cores))) +
  geom_line() +
  geom_ribbon(aes(ymin = meanTime - seTime, ymax = meanTime + seTime, 
                  fill = as.factor(cores)),
              colour = NA, alpha = 0.2) +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  guides(fill = F) +
  labs(x = "Generation", y = "Time (s)", colour = "Cores") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

ggplot(d_times, aes(x = gen, y = time, colour = as.factor(cores)), group = interaction(as.factor(cores), seed)) +
  geom_line() +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Generation", y = "Time (s)", colour = "Cores") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
