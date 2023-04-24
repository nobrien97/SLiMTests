library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggnewscale)
library(paletteer)
library(cowplot)
library(latex2exp)
library(lattice)
library(latticeExtra)
library(ggarrow)
library(GGally)
library(ggalt)

d_qg <- readRDS("/mnt/d/SLiMTests/tests/equalK_sweep/data/checkpoint/d_qg.RDS")

d_qg$id <- paste(d_qg$seed, d_qg$modelindex, sep = "_")

# Plot all the phenotype curves - we need to cluster these responses
d_pheno <- d_qg %>% mutate(KZ = log10(KZ)) %>%
  pivot_longer(cols = c(phenomean, aZ, bZ, KZ, KXZ), names_to = "trait", values_to = "value")

molTrait_names <- c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$"), 
                    TeX("$K_{XZ}$"), TeX("$log_{10}(K_Z)$"), TeX("Z"))

ggplot(d_pheno %>% filter(gen > 49000) %>% mutate(gen = gen - 50000),
       aes(x = gen, y = value, colour = trait, group = interaction(trait,id))) +
  geom_line() +
  scale_colour_paletteer_d("ggsci::nrc_npg", 1, labels = molTrait_names) +
  labs(x = "Generations post-optimum shift", y = "Mean trait/\ncomponent value",
       colour = "Trait/component") +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 3)) +
  theme(text = element_text(size = 16), legend.position = "bottom") -> singlewalk2_pheno
singlewalk2_pheno


# First we need to pivot_wider so each simulation has a series of phenotype values for time
d_qg_wide <- d_qg %>% filter(gen > 49000) %>% 
  pivot_wider(names_from = gen, values_from = phenomean, values_fill = 0) %>%
  group_by(id) %>%
  mutate(across(`49500`:`52000`, sum)) %>%
  distinct(id, .keep_all = T) %>%
  select(`49500`:`52000`) %>% ungroup()

res <- numeric(length(unique(d_qg_wide$id)))

for (i in seq_along(res)) {
  # Get appropriate subset to construct a 140x4 matrix
  d_new <- d_qg_wide[i,]
  id <- d_new$id
  # calculate l2 norm
  res[i] <- norm(as.matrix(d_new %>% select(-id)))
  names(res)[i] <- as.character(id)
}

res2 <- as.matrix(res)

# Do PCA instead
res.pca <- PCA(d_qg_wide %>% select(-id), ncp = 2)
fviz <- fviz_eig(res.pca, addlabels = TRUE)
ggsave("fviz_phenotime.png", fviz, bg = "white")
var <- get_pca_var(res.pca)
head(var$contrib)

res2 <- data.frame(id = d_qg_wide$id,
                   pc1 = res.pca$ind$coord[, 1],
                   pc2 = res.pca$ind$coord[, 2])

pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)


# calculate distance matrix - log10 scale because we have those 
# giant distances from KZ
# problem: KZ is contributing a lot to the distance, much more so than the other components
# even though it doesn't contribute very much to the trait value when it's large
distmat <- dist(res2 %>% select(-id))

# calculate hierarchical cluster
clust <- hclust(distmat, method = "ward.D2")
plot(clust)
rect.hclust(clust, k = 2, border = "blue")

# optimal clusters
fviz_nbclust(res2 %>% select(-id), FUN = hcut, method = "wss")
ggsave("phenoGen_elbow.png")

fviz_nbclust(res2 %>% select(-id), FUN = hcut, method = "silhouette")
ggsave("phenoGen_silhouette.png")

gap_stat <- cluster::clusGap(res2 %>% select(-id), FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
ggsave("phenoGen_gapstatistic.png")

# https://ethen8181.github.io/machine-learning/clustering_old/clustering/clustering.html
Distance <- function(cluster)
{
  # the center of the cluster, mean of all the points
  center <- colMeans(cluster)
  
  # calculate the summed squared error between every point and 
  # the center of that cluster 
  distance <- apply( cluster, 1, function(row)
  {
    sum( ( row - center )^2 )
  }) %>% sum()
  
  return(distance)
}
CHCriterion <- function( data, kmax, clustermethod, ...  )
{
  if( !clustermethod %in% c( "kmeanspp", "hclust" ) )
    stop( "method must be one of 'kmeanspp' or 'hclust'" )
  
  # total sum squared error (independent with the number of cluster k)
  tss <- Distance( cluster = data )
  
  # initialize a numeric vector storing the score
  wss <- numeric(kmax)
  
  # k starts from 2, cluster 1 is meaningless
  if( clustermethod == "kmeanspp" )
  {
    for( k in 2:kmax )
    {
      results <- Kmeanspp( data, k, ... )
      wss[k]  <- results$tot.withinss 
    }		
  }else # "hclust"
  {
    d <- dist( data, method = "euclidean" )
    clustering <- hclust( d, ... )
    for( k in 2:kmax )
    {
      groups <- cutree( clustering, k )
      wss[k] <- WSS( data = data, groups =  groups )
    }
  }		
  
  # between sum of square
  bss <- tss - wss[-1]
  
  # cluster count start from 2! 
  numerator <- bss / ( 1:(kmax-1) )
  denominator <- wss[-1] / ( nrow(data) - 2:kmax )
  
  criteria <- data.frame( k = 2:kmax,
                          CHIndex = numerator / denominator,
                          wss = wss[-1] )
  
  # convert to long format for plotting 
  criteria_long <- gather( criteria, "index", "value", -1 )
  
  plot <- ggplot( criteria_long, aes( k, value, color = index ) ) + 
    geom_line() + geom_point( aes( shape = index ), size = 3 ) +
    facet_wrap( ~ index, scale = "free_y" ) + 
    guides( color = FALSE, shape = FALSE )
  
  return( list( data = criteria, 
                plot = plot ) )
}
WSS <- function( data, groups )
{
  k <- max(groups)
  
  # loop through each groups (clusters) and obtain its 
  # within sum squared error 
  total <- lapply( 1:k, function(k)
  {
    # extract the data point within the cluster
    cluster <- subset( data, groups == k )
    
    distance <- Distance(cluster)
    return(distance)
  }) %>% unlist()
  
  return( sum(total) )
}

criteria <- CHCriterion(res2 %>% select(-id), kmax = 10, clustermethod = "hclust",
                        method = "ward.D2")

criteria$plot + theme_bw()
ggsave("CHCriterion_phenotime.png")

## Looks like 2 clusters is enough
sub_grp <- cutree(clust, k = 2)
names(sub_grp) <- res2$id

# Add cluster number to original dataset
d_qg_adapting <- d_qg %>% filter(gen > 49000)
d_qg_adapting$cluster <- rep(sub_grp, each = 42)
d_qg_adapting$cluster <- as.factor(d_qg_adapting$cluster)

clustermeans <- d_qg_adapting %>%
  mutate(gen = gen - 50000) %>%
  group_by(gen, cluster) %>%
  summarise(phenomean = mean(phenomean))

# Plot walks
ggplot(d_qg_adapting %>% mutate(gen = gen - 50000),
       aes(x = gen, y = phenomean, group = id)) +
  facet_grid(.~cluster) +
  geom_line(colour = "grey") +
  geom_line(data = clustermeans, mapping = aes(group = NA), colour = "red") +
  scale_colour_paletteer_d("ggsci::nrc_npg", 1) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Cluster", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype",
       colour = "Cluster") +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 3)) +
  theme(text = element_text(size = 16), legend.position = "bottom") -> singlewalk2_pheno
singlewalk2_pheno
ggsave("phenotime_clusterwalks.png", singlewalk2_pheno)


ggplot(d_qg_adapting %>% mutate(gen = gen - 50000) %>%
         filter(id %in% sample(unique(id), min(length(unique(id)), 20))),
       aes(x = gen, y = phenomean, group = id)) +
  facet_grid(.~cluster) +
  geom_line(colour = "grey") +
  geom_line(data = clustermeans, mapping = aes(group = NA), colour = "red") +
  scale_colour_paletteer_d("ggsci::nrc_npg", 1) +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Cluster", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean phenotype",
     colour = "Cluster") +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 3)) +
  theme(text = element_text(size = 16), legend.position = "bottom") -> pheno_sampled
pheno_sampled
ggsave("phenotime_clusterwalks_sampled.png", pheno_sampled)

# How many in each cluster?
d_qg_adapting %>% group_by(cluster) %>% summarise(n = n()/42)


# Plot as a levelplot
d_qg_level <- d_qg_adapting %>% 
  distinct(seed, modelindex, .keep_all = T) %>%
  mutate(cluster = as.numeric(cluster)) %>%
  select(cluster, nloci, sigma)

# Draw each replicate in a circle around the point
calcPosFromAngle <- function(r, angle) {
  r <- runif(1, max = r)
  c(x = r * sin(angle * pi/180),
    y = r * cos(angle * pi/180))
}

angles <- d_qg_level %>%
  group_by(nloci, sigma) %>%
  mutate(row_num = row_number(), n = n()) %>%
  group_by(nloci, sigma, row_num) %>%
  mutate(posX = calcPosFromAngle(0.01, (360/n) * row_num)[1],
         posY = calcPosFromAngle(0.01, (360/n) * row_num)[2])


#d_qg_level[,2] <- jitter(as.matrix(d_qg_level[,2]), factor = 0.05 * max(d_qg_adapting$nloci))
#d_qg_level[,3] <- jitter(as.matrix(d_qg_level[,3]), factor = 0.05 * max(d_qg_adapting$nloci))

d_qg_level[,2] <- d_qg_level[,2] + angles$posX * max(d_qg_adapting$nloci)
d_qg_level[,3] <- d_qg_level[,3] + angles$posY * max(d_qg_adapting$sigma)


plt_level_phenocluster <- levelplot(cluster ~ nloci*sigma, d_qg_level, 
                          col.regions = paletteer_c("viridis::viridis", 100, 1),
                          panel = panel.levelplot.points,
                          xlab = list(label = "Number of loci", cex = 1.2), 
                          ylab = list(label = "Mutational effect variance", cex = 1.2), 
                          colorkey = list(space = "bottom",
                                          labels = list(cex = 1.2)),
                          par.settings = list(layout.heights = list(xlab.key.padding = 4),
                                              par.main.text = list(just = "left",
                                                                   x = grid::unit(5, "mm"))),
                          scales = list(tck = c(1, 0),
                                        cex = 1.2), cex = 0.4,
                          pretty = T) + 
  layer_(panel.2dsmoother(..., n = 200))
plt_level_phenocluster
png("phenocluster_nlocisigma.png",
    type = "cairo",
    units = "in", 
    pointsize = 12,
    res = 300, 
    bg = "white",
    width = 8.5, 
    height = 5.56)
plt_level_phenocluster
dev.off()

# What molecular trait behaviour describes each cluster?
# Are there particular combinations of molecular traits that drive
# these two different phenotypic responses?

# Sample some examples from each cluster - from the most extreme cases
d_qg_adapting %>% 
  distinct(id, .keep_all = T) %>%
  group_by(modelindex, cluster) %>%
  summarise(clusterPerc = n()/48) %>%
  ungroup() -> clusterPercs

inner_join(clusterPercs, d_qg_adapting %>% 
             distinct(modelindex, .keep_all = T) %>% 
             select(modelindex, nloci, sigma), by = "modelindex") -> clusterPercs

ggplot(clusterPercs, aes(x = nloci, y = clusterPerc)) +
  facet_grid(.~cluster) +
  geom_line() +
  theme_bw()

print(clusterPercs %>% filter(clusterPerc == min(clusterPerc, na.rm = T) |
                          clusterPerc == max(clusterPerc, na.rm = T) |
                          clusterPerc == median(clusterPerc, na.rm = T)), n = 25)

ggplot(clusterPercs, aes(x = sigma, y = clusterPerc)) +
  facet_grid(.~cluster) +
  geom_line() +
  theme_bw()


# Mostly cluster 2
sampled_ids <- clusterPercs %>% 
  pivot_wider(names_prefix = "clusterPerc_",
              names_from = "cluster", values_from = "clusterPerc") %>%
  filter(clusterPerc_1 == max(clusterPerc_1, rm.na = T) | 
           clusterPerc_1 == min(clusterPerc_1, rm.na = T))

sampled_ids <- pivot_longer(sampled_ids, cols = c("clusterPerc_1", "clusterPerc_2"),
                            names_to = "cluster", values_to = "perc") %>%
  group_by(cluster) %>%
  drop_na() %>%
  filter(modelindex %in% sample(modelindex, 1))

d_qg_sampled <- d_qg_adapting %>%
                 filter(modelindex %in% sample(sampled_ids$modelindex)) %>%
  group_by(modelindex) %>%
  filter(seed %in% sample(seed, 1))


d_qg_sampled <- d_qg_sampled %>% mutate(KZ = log10(KZ)) %>%
  pivot_longer(cols = c(phenomean, aZ, bZ, KZ, KXZ), names_to = "trait", values_to = "value")


# Plot molecular traits
ggplot(d_qg_sampled %>% mutate(gen = gen - 50000),
       aes(x = gen, y = value, colour = trait, group = interaction(trait, id))) +
  facet_grid(.~cluster) +
  geom_line() +
  scale_colour_viridis_d() +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Cluster", 
                                         breaks = NULL, labels = NULL)) +
  labs(x = "Generations post-optimum shift", y = "Mean trait/component value",
       colour = "Trait/component") +
  theme_bw()

# plot example walks from both clusters (pairwise)
cc <- paletteer_c("grDevices::Burg", 3)
cc2 <- paletteer_c("grDevices::Blues", 50, -1)

d_qg_sampled <- d_qg_sampled %>%
  pivot_wider(names_from = "trait", values_from = "value")


sampled_id <- unique(d_qg_sampled$id)

d_qg_sampled$gen_width <- scales::rescale(d_qg_sampled$gen, to = c(0.001, 1))

plotPairwiseScatter <- function(x, y, labels) {
  ggplot(d_qg_sampled %>% mutate(KZ = log10(KZ)), aes(x = .data[[x]], y = .data[[y]], colour = phenomean)) +
    geom_point() +
    geom_encircle(colour = "black", mapping = aes(group = id, linetype = cluster)) +
    scale_colour_gradientn(colors = cc) +
    labs(colour = "Phenotype (Z)", linetype = "Cluster") +

    new_scale_colour() +
    geom_arrow_segment(data = d_qg_sampled %>% mutate(KZ = log10(KZ)) %>% filter(id %in% sampled_id), 
                       mapping = aes(x = lag(.data[[x]]), y = lag(.data[[y]]), 
                                     xend = .data[[x]], yend = .data[[y]], 
                                     group = id, linewidth_head = gen_width, 
                                     linewidth_fins = gen_width * 0.8,
                                     colour = gen_width), 
                       arrow_head = arrow_head_line()) +
    scale_colour_gradientn(colors = cc2, labels = c(0, 0.25*2500, 0.5*2500, 0.75*2500, 2500)) +
    scale_linewidth(guide = "none") +
    labs(x = labels[1], y = labels[2], colour = "Generation") +
    theme_bw() + 
    theme(legend.position = "bottom") +
    guides(colour=guide_colourbar(barwidth=15))
}

aZbZScatter <- plotPairwiseScatter("aZ", "bZ", c(TeX("$\\alpha_Z$"), TeX("$\\beta_Z$")))
aZKZScatter <- plotPairwiseScatter("aZ", "KZ", c(TeX("$\\alpha_Z$"), TeX("$K_Z$")))
aZKXZScatter <- plotPairwiseScatter("aZ", "KXZ", c(TeX("$\\alpha_Z$"), TeX("$K_{XZ}$")))
bZKZScatter <- plotPairwiseScatter("bZ", "KZ", c(TeX("$\\beta_Z$"), TeX("$K_Z$")))
bZKXZScatter <- plotPairwiseScatter("bZ", "KXZ", c(TeX("$\\beta_Z$"), TeX("$K_{XZ}$")))
KZKXZScatter <- plotPairwiseScatter("KZ", "KXZ", c(TeX("$K_Z$"), TeX("$K_{XZ}$")))

plot_list <- list(ggally_densityDiag(d_qg_adapting, mapping = aes(x = aZ)), 
                  ggally_cor(d_qg_adapting, mapping = aes(x = aZ, y = bZ)),
                  ggally_cor(d_qg_adapting, mapping = aes(x = aZ, y = KZ)),
                  ggally_cor(d_qg_adapting, mapping = aes(x = aZ, y = KXZ)),
                  
                  aZbZScatter, 
                  ggally_densityDiag(d_qg_adapting, mapping = aes(x = bZ)), 
                  ggally_cor(d_qg_adapting, mapping = aes(x = bZ, y = KZ)),
                  ggally_cor(d_qg_adapting, mapping = aes(x = bZ, y = KXZ)),
                  
                  aZKZScatter,
                  bZKZScatter,
                  ggally_densityDiag(d_qg_adapting, mapping = aes(x = KZ)), 
                  ggally_cor(d_qg_adapting, mapping = aes(x = KZ, y = KXZ)),
                  
                  aZKXZScatter,
                  bZKXZScatter,
                  KZKXZScatter,
                  ggally_densityDiag(d_qg_adapting, mapping = aes(x = KXZ)))

xlabs <- c(TeX("$\\alpha_Z$", output = "character"), TeX("$\\beta_Z$", output = "character"), 
           TeX("$log_{10}(K_Z)$", output = "character"), TeX("$K_{XZ}$", output = "character"))

ggmatrix(plot_list, nrow = 4, ncol = 4, xAxisLabels = xlabs, yAxisLabels = xlabs,
         progress = T, byrow = T, labeller = "label_parsed", 
         legend = grab_legend(plot_list[[5]])) + theme_bw() + 
  theme(legend.position = "bottom") -> pair_mat
pair_mat
ggsave("molTrait_landscape_clusterExtremes.png", pair_mat, width = 11, height = 8)

# Recluster by molecular traits? How many combinations of ways do they respond?



# Timee to reach optimum distribution
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
    stat_ellipse(geom = "polygon", mapping = aes(colour = cluster), fill = NA) +
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

