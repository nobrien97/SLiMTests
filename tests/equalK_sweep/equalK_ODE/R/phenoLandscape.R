# plot molecular components vs phenotype w/ samples overlaid 
# as points connected by arrows (for time)
library(tidyverse)
library(gganimate)
library(latex2exp)
library(lattice)
library(latticeExtra)
library(paletteer)

# Load data
d_ode_phasemeasures <- readRDS("d_phaseDiagramMeasures.RDS")

# Get range of each mol trait
molTraitRange <- list(aZ = quantile(d_ode_phasemeasures$aZ, c(0.01, 0.99)),
                      bZ = quantile(d_ode_phasemeasures$bZ, c(0.01, 0.99)),
                      KZ = quantile(d_ode_phasemeasures$KZ, c(0.01, 0.99)),
                      KXZ = quantile(d_ode_phasemeasures$KXZ, c(0.01, 0.99)))

# Get sequences between that range
molTraitSeq <- lapply(molTraitRange, function(x) {
  seq(from = x[1], to = x[2], length.out = 100)
})

# Generate data frame of mol trait combos
d_molTrait <- expand.grid(molTraitSeq)
write.table(samples, "./samples.csv", sep = ",", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# Use ODE Landscaper
runLandscaper <- function(df_path, width, optimum, threads, multi=FALSE) {
  if (multi) {
    system(sprintf("ODELandscaper -i %s* -o ./landscape.csv -w %f -p %f -t %i",
                   df_path, width, optimum, threads))
    result <- read_csv("./landscape.csv", col_names = F, col_types = "d")
    fig <- plotFigure(result, input)
    
    return(fig)
  }
  
  system(sprintf("ODELandscaper -i %s -o ./landscape.csv -w %f -p %f -t %i",
                 df_path, width, optimum, threads))
  result <- read_csv("./landscape.csv", col_names = F, col_types = "d")
  names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  return(result)
}
