library(tidyverse)
library(factoextra)

d_qg <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/checkpoint/d_qg.RDS")

d_qg$id <- paste(d_qg$seed, d_qg$modelindex, sep = "_")

# scale variables
d_qg_scaled <- d_qg %>% mutate_at(c("aZ", "bZ", "KZ", "KXZ"), ~(scale(.) %>% as.vector()))

res <- numeric(length(unique(d_qg$id)))

for (i in seq_along(res)) {
  # Get appropriate subset to construct a 140x4 matrix
  d_new <- d_qg[((i-1)*140+1):(i*140),] 
  id <- d_new[1, 18]
  # calculate l2 norm
  res[i] <- norm(as.matrix(d_new %>% select(aZ, bZ, KZ, KXZ)))
  names(res)[i] <- as.character(id)
}

res2 <- as.matrix(res)

# calculate distance matrix - log10 scale because we have those 
# giant distances from KZ
# problem: KZ is contributing a lot to the distance, much more so than the other components
# even though it doesn't contribute very much to the trait value when it's large
distmat <- dist(log10(res2))

# calculate hierarchical cluster
clust <- hclust(distmat, method = "ward.D2")
plot(clust)
rect.hclust(clust, k = 4, border = "blue")

# optimal clusters
fviz_nbclust(log10(res2), FUN = hcut, method = "wss")
ggsave("molTraitGen_elbow.png")

fviz_nbclust(log10(res2), FUN = hcut, method = "silhouette")
ggsave("molTraitGen_silhouette.png")

gap_stat <- cluster::clusGap(log10(res2), FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
ggsave("molTraitGen_gapstatistic.png")

## Looks like its either 2 or 3... we'll go with 3
sub_grp <- cutree(clust, k = 3)

# Add cluster number to original dataset
d_qg %>% distinct(seed, modelindex, .keep_all = T) %>% 
  mutate(cluster = as_factor(sub_grp)) -> d_qg_mdl

# Plot
ggplot(d_qg_mdl, aes(x = nloci, y = sigma, colour = cluster)) +
  geom_point() +
  stat_ellipse(geom = "polygon", aes(fill = cluster), alpha = 0.3) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d(guide = "none") +
  labs(x = "Number of loci", y = "Mutational effect variance", colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggsave("molTraitGen_cluster.png", width = 7, height = 5)

se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

# Time to reach the optimum
d_qg %>% filter(gen >= 49500) %>%
  mutate(isAdapted = between(phenomean, 1.9, 2.1),
         isAdapting = between(phenomean, 1.2, 1.9)) %>%
  group_by(seed, modelindex) %>%
  mutate(adaptTime = ifelse(any(isAdapted), min(gen[isAdapted]) - 50000, -1),
         initRespTime = ifelse(any(isAdapting), min(gen[isAdapting]) - 50000, -1)) %>%
  distinct(seed, modelindex, .keep_all = T) %>%
  ungroup() -> d_qg_adapttime #%>%
  filter(adaptTime != -1, initRespTime != -1) %>%
  group_by(nloci, sigma) %>%
  summarise(meanAdaptTime = mean(adaptTime), 
            CIAdaptTime = CI(adaptTime),
            meanInitRespTime = mean(initRespTime),
            CIInitRespTime = CI(initRespTime)) -> d_qg_meanadapttime

  
d_qg_adapttime$nloci_cat <- cut(d_qg_adapttime$nloci,
                            breaks=c(-Inf, 10, 50, 250, Inf))
  
d_qg_adapttime$sigma_cat <- cut(d_qg_adapttime$sigma,
                            breaks=c(-Inf, 0.05, 0.5, 1, Inf))
  
ggplot(d_qg_adapttime, aes(x = adaptTime)) + 
  facet_grid(nloci_cat~sigma_cat) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of QTLs", 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Mutational effect variance", 
                                         breaks = NULL, labels = NULL)) +
  geom_density() +
  labs(x = "Time to adaptation (generations)", y = "Density") +
  theme_bw()

# group by adaptation time - slow( >750 gens), fast (<=750 gens)
d_qg_adapttime$adaptSpeed <- ifelse(d_qg_adapttime$adaptTime > 750 | 
                                      d_qg_adapttime$adaptTime == -1, "slow", "fast")

d_qg_adapttime %>% mutate(adaptSpeed = as_factor(adaptSpeed)) -> d_qg_adapttime
d_qg_adapttime %>% group_by(adaptSpeed) %>% summarise(nRows = n())

# map clusters to adaptSpeed
d_qg_adapttime$cluster <- d_qg_mdl$cluster

plotAdaptSpeedCluster <- function(x, y, labels) {
  ggplot(d_qg_adapttime, aes(x = .data[[x]], y = .data[[y]], colour = adaptSpeed)) +
    geom_point() +
    scale_colour_viridis_d() +
    labs(colour = "Adaptation speed") +
    
    new_scale_colour() +
    stat_ellipse(mapping = aes(colour = cluster)) +
    scale_colour_paletteer_d("ggsci::lanonc_lancet", 1) +
    scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2], colour = "Cluster") +
    theme_bw() + 
    theme(legend.position = "bottom")
}

plotAdaptSpeedCluster("nloci", "sigma", c("Number of loci", "Mutational effect variance"))

d_qg_mdl %>% 
  group_by(cluster) %>% 
  summarise(meanaZ = mean(aZ), 
            CIaZ = CI(aZ),
            meanbZ = mean(bZ), 
            CIbZ = CI(bZ),
            meanKZ = mean(KZ), 
            CIKZ = CI(KZ),
            meanKXZ = mean(KXZ),
            CIKXZ = CI(KXZ))

CI(d_qg_mdl[d_qg_mdl$cluster == 1,]$aZ)
# Sample three examples from each cluster to look at more closely
seed <- sample(1:.Machine$integer.max, 1)
# 751129649
set.seed(seed)
set.seed(751129649)

d_sampled <- lapply(1:3, function(i) {
  d_qg %>% filter(id %in% sample(d_qg_mdl$id[d_qg_mdl$cluster == i], 3))
})

d_sampled <- bind_rows(d_sampled, .id = "cluster")

# Plot those phenomeans
ggplot(d_sampled %>% filter(gen > 49000), aes(x = gen, y = phenomean, colour = cluster, group = id)) +
  geom_line() +
  scale_fill_viridis_d(guide = "none") +
  scale_colour_viridis_d() +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype", colour = "Simulation parameters") +
  theme_bw()

# molecular trait trajectories
d_sampled$gen_width <- rescale(d_sampled$gen, to = c(0.001, 1))
cc2 <- paletteer_c("grDevices::Emrld", 50, -1)


plotPairwiseScatter <- function(x, y, labels) {
  ggplot(d_sampled %>% mutate(KZ = log10(KZ)), aes(x = .data[[x]], y = .data[[y]], colour = cluster)) +
    geom_point() +
    scale_colour_viridis_d() +
    labs(colour = "Cluster") +
    
    # new_scale_colour() +
    # geom_arrow_segment(data = d_sampled %>% mutate(KZ = log10(KZ)), 
    #                    mapping = aes(x = lag(.data[[x]]), y = lag(.data[[y]]), 
    #                                  xend = .data[[x]], yend = .data[[y]], 
    #                                  group = id, linewidth_head = gen_width, 
    #                                  linewidth_fins = gen_width * 0.8,
    #                                  colour = gen_width), 
    #                    arrow_head = arrow_head_line()) +
    # scale_colour_gradientn(colors = cc2, labels = c(0, 0.25*2500, 0.5*2500, 0.75*2500, 2500)) +
    # scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2])+#, colour = "Generation") +
    theme_bw() + 
    theme(legend.position = "bottom") #+
    #guides(colour=guide_colourbar(barwidth=20))
}

aZbZScatter <- plotPairwiseScatter("aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
aZKZScatter <- plotPairwiseScatter("aZ", "KZ", c(TeX("$\\alpha_Z$"), TeX("$K_Z$")))
aZKXZScatter <- plotPairwiseScatter("aZ", "KXZ", c(TeX("$\\alpha_Z$"), TeX("$K_{XZ}$")))
bZKZScatter <- plotPairwiseScatter("bZ", "KZ", c(TeX("$\\beta_Z$"), TeX("$K_Z$")))
bZKXZScatter <- plotPairwiseScatter("bZ", "KXZ", c(TeX("$\\beta_Z$"), TeX("$K_{XZ}$")))
KZKXZScatter <- plotPairwiseScatter("KZ", "KXZ", c(TeX("$K_Z$"), TeX("$K_{XZ}$")))

plot_list <- list(ggally_densityDiag(d_sampled, mapping = aes(x = aZ)), 
                  ggally_cor(d_sampled, mapping = aes(x = aZ, y = bZ)),
                  ggally_cor(d_sampled, mapping = aes(x = aZ, y = KZ)),
                  ggally_cor(d_sampled, mapping = aes(x = aZ, y = KXZ)),
                  
                  aZbZScatter, 
                  ggally_densityDiag(d_sampled, mapping = aes(x = bZ)), 
                  ggally_cor(d_sampled, mapping = aes(x = bZ, y = KZ)),
                  ggally_cor(d_sampled, mapping = aes(x = bZ, y = KXZ)),
                  
                  aZKZScatter,
                  bZKZScatter,
                  ggally_densityDiag(d_sampled, mapping = aes(x = KZ)), 
                  ggally_cor(d_sampled, mapping = aes(x = KZ, y = KXZ)),
                  
                  aZKXZScatter,
                  bZKXZScatter,
                  KZKXZScatter,
                  ggally_densityDiag(d_sampled, mapping = aes(x = KXZ)))

xlabs <- c(TeX("$\\alpha_Z$", output = "character"), TeX("$\\beta_Z$", output = "character"), 
           TeX("$log_{10}(K_Z)$", output = "character"), TeX("$K_{XZ}$", output = "character"))

ggmatrix(plot_list, nrow = 4, ncol = 4, xAxisLabels = xlabs, yAxisLabels = xlabs,
         progress = T, byrow = T, labeller = "label_parsed", 
         legend = grab_legend(plot_list[[5]])) + theme_bw() + 
  theme(legend.position = "bottom") -> pair_mat
pair_mat

ggsave("molTrait_landscape_singleWalk_clusters.png", pair_mat, width = 11, height = 8)

