library(tidyverse)
library(paletteer)
library(latex2exp)
library(ggarrow)
library(GGally)
library(ggalt)
library(ggnewscale)

# Functions
se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}


# Load data
d_qg_add <- read_csv("/mnt/d/SLiMTests/tests/fixedK/additive/slim_qg.csv", col_names = F)
d_qg_net <- read_csv("/mnt/d/SLiMTests/tests/fixedK/slim_qg.csv", col_names = F)

colnames(d_qg_add) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "phenovar",
                        "dist", "w", "deltaPheno", "deltaw")

colnames(d_qg_net) <- c(colnames(d_qg_add), "aZ", "bZ", "KZ", "KXZ")

d_qg_add$model <- "Additive"
d_qg_net$model <- "NAR"
d_qg <- bind_rows(d_qg_add, d_qg_net)

rm(d_qg_net, d_qg_add)

d_combos <- read_delim("./combos.csv", delim = " ", col_names = F)
names(d_combos) <- c("nloci")

d_qg <- d_qg %>% mutate(model = as_factor(model),
                        nloci = as_factor(d_combos$nloci[.$modelindex]),
                        id = paste(seed, modelindex, model, sep = "_"))

# Mutation data
d_mut_add <- read_csv("/mnt/d/SLiMTests/tests/fixedK/additive/slim_muts.csv", col_names = F)
d_mut_net <- read_csv("/mnt/d/SLiMTests/tests/fixedK/slim_muts.csv", col_names = F)

colnames(d_mut_add) <- c("gen", "seed", "modelindex", "mutType", "id", "pos", "constraint",
                        "originGen", "value", "chi", "freq", "count", "fixGen")
colnames(d_mut_net) <- colnames(d_mut_add) 

d_mut_add$model <- "Additive"
d_mut_net$model <- "NAR"
d_mut <- bind_rows(d_mut_add, d_mut_net)

rm(d_mut_net, d_mut_add)

d_mut <- d_mut %>% mutate(model = as_factor(model),
                        nloci = as_factor(d_combos$nloci[.$modelindex]),
                        id = paste(seed, modelindex, model, sep = "_"))


# Figure 1: adaptive walk
d_qg_mean <- d_qg %>%
  group_by(gen, model, nloci) %>%
  summarise(phenoCI = CI(phenomean), phenomean = mean(phenomean))

ggplot(d_qg_mean %>% filter(gen > 49000) %>% 
         mutate(gen = gen - 50000), 
       aes(x = gen, y = phenomean, colour = model)) +
  facet_grid(nloci~.) +
  geom_line() +
  geom_ribbon(aes(ymin = phenomean - phenoCI, ymax = phenomean + phenoCI,
                  fill = model),
              colour = NA, alpha = 0.2) +
  scale_colour_paletteer_d("ggsci::default_aaas") +
  scale_fill_paletteer_d("ggsci::default_aaas", guide = "none") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Population phenotype mean", 
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16))

# Figure 2: distribution of fixed effects
d_fixed <- d_mut %>% filter(count >= 10000, originGen > 49500)

d_qg %>%
  group_by(seed, modelindex) %>%
  filter(any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_adapted

d_fixed <- d_fixed %>% 
  filter(interaction(seed, modelindex, model) %in% 
           interaction(d_adapted$seed, d_adapted$modelindex, d_adapted$model))

d_fixed <- d_fixed %>% mutate(molTrait = recode_factor(mutType, 
                                            `1`="Neutral", `2`="Del", `3`="aZ", `4`="bZ", 
                                            `5`="KZ", `6`="KXZ"))

ggplot(d_fixed,
       aes(x = value)) +
  facet_grid(nloci~model) +
  geom_histogram(bins = 100) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Model", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Effect size", y = "Count",
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16))

# Fig 3 - distribution among molecular traits

mutType_names <- as_labeller(c(
  "aZ" = TeX("$\\alpha_Z$"),
  "bZ" = TeX("$\\beta_Z$"),
  "KZ" = TeX("$K_Z$"),
  "KXZ" = TeX("$K_{XZ}$")
))

levels(d_fixed$molTrait) <- c("aZ" = TeX("$\\alpha_Z$"),
                              "bZ" = TeX("$\\beta_Z$"),
                              "KZ" = TeX("$K_Z$"),
                              "KXZ" = TeX("$K_{XZ}$"))


ggplot(d_fixed %>% filter(model == "NAR"),
       aes(x = value)) +
  facet_grid(nloci~molTrait, labeller = label_parsed) +
  geom_histogram(bins = 100) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Molecular component", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Effect size", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 16))

# Fig 4
cc <- paletteer_c("grDevices::Burg", 3)
cc2 <- paletteer_c("grDevices::Blues", 50, -1)

d_qg$gen_width <- scales::rescale(d_qg$gen, to = c(0.001, 1))

d_qg_sampled <- d_qg %>% filter(model == "NAR", nloci == 4) %>% filter(id %in% sample(id, 1))
sampled_id <- unique(d_qg_sampled$id)


plotPairwiseScatter <- function(x, y, xylabs) {
  ggplot(d_qg_sampled, aes(x = .data[[x]], y = .data[[y]], colour = phenomean)) +
    geom_point() +
    geom_encircle(colour = "black", mapping = aes(group = id)) +
    scale_colour_gradientn(colors = cc) +
    labs(colour = "Phenotype (Z)") +
    
    new_scale_colour() +
    geom_arrow_segment(data = d_qg_sampled, 
                       mapping = aes(x = lag(.data[[x]]), y = lag(.data[[y]]), 
                                     xend = .data[[x]], yend = .data[[y]], 
                                     group = id, linewidth_head = gen_width, 
                                     linewidth_fins = gen_width * 0.8,
                                     colour = gen_width), 
                       arrow_head = arrow_head_line()) +
    scale_colour_gradientn(colors = cc2, labels = c(0, 0.25*10500, 0.5*10500, 
                                                    0.75*10500, 10500)) +
    scale_linewidth(guide = "none") +
    labs(x = xylabs[1], y = xylabs[2], colour = "Generation") +
    theme_bw() + 
    theme(legend.position = "bottom") +
    guides(colour=guide_colourbar(barwidth=15))
}

plotPairwiseScatter("aZ", "KZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))

ggplot(d_qg_sampled, aes(x = aZ, y = bZ, colour = gen_width)) +
  geom_arrow_segment(data = d_qg_sampled, 
                     mapping = aes(x = lag(aZ), y = lag(bZ), 
                                   xend = aZ, yend = bZ, 
                                   group = id, linewidth_head = gen_width, 
                                   linewidth_fins = gen_width * 0.8,
                                   colour = gen_width), 
                     arrow_head = arrow_head_line()) +
  scale_colour_gradientn(colors = cc2, labels = c(0, 0.25*10500, 0.5*10500, 
                                                  0.75*10500, 10500)) +
  labs(colour = "Phenotype (Z)")

