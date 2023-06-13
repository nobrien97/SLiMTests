# From each adaptive step, generate some mutations and get the distribution of 
# fitness effects
library(tidyverse)
library(future)
library(doParallel)
library(foreach)
source("wrangle_data.R")

MutationScreenExp <- function(fixed, n) {
  # at each step in the walk, sample n mutations for each molecular component
  # and add them to the ancestor (the previous step) - then measure fitness
  foreach (curSeed = unique(fixed$seed)) %dopar%  {
    sampled_aZ = rnorm(n)
    sampled_bZ = rnorm(n)
    fixed[fixed$seed == curSeed,] 

  }
  fixed %>%
    group_by(rank, seed, modelindex) %>%
    {
      newRow = tibble_row(rank = rank, seed = seed, modelindex = modelindex, 
                          aZ = fixEffectSum_aZ + sampled_aZ, 
                          bZ = fixEffectSum_bZ + sampled_bZ)
    }
}