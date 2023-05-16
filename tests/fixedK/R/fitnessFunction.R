library(tidyverse)
library(latex2exp)
times <- seq(-8, 8, length.out = 1000)

fitnessFn <- function(x, width, optimum) {
  dists <- (x - optimum)^2
  dists[dists < 0] <- 0.0
  
  exp(-(dists * width))
}

dat <- data.frame(
  timepoint = times,
  pow_05 = fitnessFn(times, 0.05, 0)#,
  #pow_7 = fitnessFn(times, 0.7, 0)
)


dat %>% 
  pivot_longer(cols = -timepoint, names_to = "power", values_to = "val") -> dat

cc <- c("#648FFF", "#DC267F")

ggplot(dat, aes(x = timepoint, y = val)) +
  geom_line(size = 1.1) +
  scale_colour_manual(values = cc, name = TeX("Width = $\\frac{1}{2{\\sigma_w}^2}$"), labels = c("Width = 0.05", 
                                                                  "Width = 0.7")) +
  labs(x = "Deviation from phenotypic optimum", y = "Fitness") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = c(-5, 0, 1, 5)) +
  theme_classic() +
  theme(text = element_text(size = 16), legend.text.align = 0)

ggsave("fitnessFunction.png", width = 8, height = 6)
