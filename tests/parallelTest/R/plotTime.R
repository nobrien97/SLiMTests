library(tidyverse)
library(paletteer)

setwd("/mnt/c/GitHub/SLiMTests/tests/parallelTest/R")

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

d_times <- read_csv("/mnt/d/SLiMTests/tests/parallelTest/slim_time.csv", col_names = F)
d_times_add <- read_csv("/mnt/d/SLiMTests/tests/parallelTest/slim_time_add.csv", col_names = F)


colnames(d_times) <- c("gen", "seed", "modelindex", "time")
colnames(d_times_add) <- c("gen", "seed", "modelindex", "time")

d_times$model <- "Network"
d_times_add$model <- "Additive"

d_times <- rbind(d_times, d_times_add)
rm(d_times_add)

cores_treatments <- c(1, 2, 4, 6)
d_times %>% mutate(cores = cores_treatments[modelindex]) -> d_times

d_times %>%
  group_by(seed, modelindex, model) %>%
  mutate(delta = time - lag(time)) -> d_times

d_times %>%
  group_by(gen, model, cores) %>%
  summarise(meanTime = mean(time),
            seTime = se(time),
            meanDelta = mean(delta),
            seDelta = se(delta)) -> time_summary

ggplot(time_summary, aes(x = gen, y = meanTime, colour = as.factor(cores))) +
  facet_grid(model~.) +
  geom_line() +
  geom_ribbon(aes(ymin = meanTime - seTime, ymax = meanTime + seTime, 
                  fill = as.factor(cores)),
              colour = NA, alpha = 0.2) +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Model", 
                                         breaks = NULL, labels = NULL)) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  guides(fill = F) +
  labs(x = "Generation", y = "Time (s)", colour = "Cores") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
ggsave("times.png")

# extrapolate to 60000 generations
time_summary %>%
  group_by(model, cores) %>%
  slice_tail(n = 1) %>%
  select(!c(meanTime, seTime, gen, seDelta)) %>%
  mutate(wholeSimEstimate = (meanDelta * (60000/50))/3600)
