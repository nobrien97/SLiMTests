# Generates combinations of molecular component values from a variety of distributions
# to see how the bug in twice-exponentiating molecular trait values affected things
library(SimDesign)
library(tidyverse)
library(factoextra)
library(FactoMineR)

sigma <- diag(4) * 0.1
mu <- c(0, 0, 0, 0)

# Set up a dataframe
d_normal <- data.frame(aZ = rep(-1, 1024), 
                       bZ = rep(-1, 1024), 
                       KZ = rep(-1, 1024), 
                       KXZ = rep(-1, 1024))

# each sample is the sum of 5 sampled mutations from a normal distribution, exponentiated
for (i in seq_len(nrow(d_normal))) {
  d_normal[i,] <- colSums(exp(rmvnorm(1, mu, sigma)))
}

write.table(d_normal, "exp.csv", sep = ",", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# Twice exponentiated data (bug)

write.table(exp(d_normal), "expexp.csv", sep = ",", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

system("ODELandscaper -i exp.csv -o ./exp_out.csv -t 10 -p 2 -w 0.05")
system("ODELandscaper -i expexp.csv -o ./expexp_out.csv -t 10 -p 2 -w 0.05")

d_exp <- read_csv("exp_out.csv", col_names = F)
d_doubleexp <- read_csv("expexp_out.csv", col_names = F)

colnames(d_exp) <- c("w", "z", "aZ", "bZ", "KZ", "KXZ")
colnames(d_doubleexp) <- c("w", "z", "aZ", "bZ", "KZ", "KXZ")

dist1D <- function(x, y) {
  return(abs(x - y))
} 

d_dist <- data.frame(dist_z = dist1D(d_exp$z, d_doubleexp$z),
                     dist_w = dist1D(d_exp$w, d_doubleexp$w))

ggplot(d_dist, aes(x = dist_w, y = dist_z)) + 
  geom_point() + 
  labs(x = "Distance between fitnesses",
       y = "Distance between phenotypes") +
  theme_bw()

d_exp$model <- "Allelic effect exponentiation"
d_doubleexp$model <- "Molecular component + Allelic effect exponentiation"

d_com <- rbind(d_exp, d_doubleexp)

res.pca <- PCA(d_com %>% 
                 select(aZ, bZ, KZ, KXZ), 
               scale.unit = T, graph = F)
fviz_eig(res.pca, addlabels = TRUE)
ggsave("scree_phenofreq.png", bg = "white")
var <- get_pca_var(res.pca)
var$contrib


# https://tem11010.github.io/Plotting-PCAs/
d_com$pc1 <- res.pca$ind$coord[, 1]
d_com$pc2 <- res.pca$ind$coord[, 2]
pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars 

ggplot(d_com, aes(x = pc1, y = pc2, colour = z)) +
  facet_grid(.~model) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(shape = 1, size = 2) +
  scale_colour_viridis_c() +
  labs(x = sprintf("PC 1 (%.2f%%)", res.pca$eig[1,2]), y = sprintf("PC 2 (%.2f%%)", res.pca$eig[2,2]), colour = "Phenotypic value") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_com, aes(x = aZ, y = bZ)) +
  geom_contour(aes(z = z, colour = after_stat(level))) +
  scale_colour_viridis_c() +
  labs(x = "aZ", y = "bZ", fill = "Phenotype") +
  theme_bw() +
  theme(text = element_text(size = 16))
