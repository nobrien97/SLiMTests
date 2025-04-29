# Calculates density and histogram counts for the massive pairwise epistasis dataset

library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(data.table)

# Get command line arguments
## 1: modelindex
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

# Paths
GDATA_PATH <- "/g/data/ht96/nb9894/newMotifs/randomisedStarts/calcMutationStats/"
WRITE_PATH <- "/scratch/ht96/nb9894/newMotifs/randomisedStarts/epistasisDensity/"
EPISTASIS_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_density", model, ".csv")
EPISTASIS_WEIGHTED_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_density", model, ".csv")
EPISTASIS_NOMOLCOMP_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_nomolcomp_density", model, ".csv")
EPISTASIS_WEIGHTED_NOMOLCOMP_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_nomolcomp_density", model, ".csv")
EPISTASIS_MEAN_FILE <- paste0(WRITE_PATH, "d_epi_mean", model, ".csv")
EPISTASIS_WEIGHTED_MEAN_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_mean", model, ".csv")
EPISTASIS_NOMOLCOMP_MEAN_FILE <- paste0(WRITE_PATH, "d_epi_nomolcomp_mean", model, ".csv")
EPISTASIS_WEIGHTED_NOMOLCOMP_MEAN_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_nomolcomp_mean", model, ".csv")

# extract the model from the database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "epistasis.db"))

d_epistasis <- tbl(con, "tab_epistasis") %>% 
  filter(modelindex == model) %>%
  collect()

# Make sure it has all worked properly
if (nrow(d_epistasis) < 1) {
  print(paste("Model ", run, " couldn't get epistasis data (likely no adapted pops)."))
  q(save = "no")
}

# Add optPerc labels to use instead of gen
dpdt <- read.csv(paste0(GDATA_PATH, "d_dpdt.csv"), header = F)
dpdt <- dpdt %>% filter(V2 == model)
d_epistasis$timePoint <- rep(dpdt$V1, each = nrow(d_epistasis)/nrow(dpdt))

# Calculate density for each optPerc across Pa - Pwt, Pb - Pwt, Pab - Pwt, ew, ep
# and across log fitness
d_epistasis %>% 
  nest_by(timePoint, modelindex, mutType_ab) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         wwt = list(data.frame(density(log(data$wwt))[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density

# No molcomp
d_epistasis %>% 
  nest_by(timePoint, modelindex) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         wwt = list(data.frame(density(log(data$wwt))[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density_nomolcomp


# Mean: compare the mean epistasis between models to pick the models to compare
# in density curves

d_epistasis %>%
  group_by(timePoint, modelindex, mutType_ab) %>%
  summarise(meanEP = mean(ep),
            sdEP = sd(ep),
            meanEW = mean(ew),
            sdEW = sd(ew),
            meanEW_s = mean(ew_s),
            sdEW_s = sd(ew_s),
            n = n()) -> d_epistasis_mean

d_epistasis %>%
  group_by(timePoint, modelindex) %>%
  summarise(meanEP = mean(ep),
            sdEP = sd(ep),
            meanEW = mean(ew),
            sdEW = sd(ew),
            meanEW_s = mean(ew_s),
            sdEW_s = sd(ew_s),
            n = n()) -> d_epistasis_mean_nomolcomp

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_DENSITY_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_mean,
                   EPISTASIS_MEAN_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_density_nomolcomp,
                   EPISTASIS_NOMOLCOMP_DENSITY_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_mean_nomolcomp,
                   EPISTASIS_NOMOLCOMP_MEAN_FILE, sep = ",", col.names = F, row.names = F)


rm(d_epistasis)

# Repeat for frequency-weighted data
d_epistasis_freq <- tbl(con, "tab_epistasis_freq") %>% 
  filter(modelindex == model) %>%
  collect()

d_epistasis_freq$timePoint <- rep(dpdt$V1, each = nrow(d_epistasis_freq)/nrow(dpdt))

d_epistasis_freq %>% 
  nest_by(timePoint, modelindex, mutType_ab) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         wwt = list(data.frame(density(log(data$wwt))[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density

# No molcomp
d_epistasis_freq %>% 
  nest_by(timePoint, modelindex) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         wwt = list(data.frame(density(log(data$wwt))[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density_nomolcomp


d_epistasis_freq %>%
  group_by(timePoint, modelindex- mutType_ab) %>%
  summarise(meanEP = mean(ep),
            sdEP = sd(ep),
            meanEW = mean(ew),
            sdEW = sd(ew),
            meanEW_s = mean(ew_s),
            sdEW_s = sd(ew_s),
            n = n()) -> d_epistasis_mean_nomolcomp

d_epistasis_freq %>%
  group_by(timePoint, modelindex) %>%
  summarise(meanEP = mean(ep),
            sdEP = sd(ep),
            meanEW = mean(ew),
            sdEW = sd(ew),
            meanEW_s = mean(ew_s),
            sdEW_s = sd(ew_s),
            n = n()) -> d_epistasis_mean

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_WEIGHTED_DENSITY_FILE, sep = ",", 
                   col.names = F, row.names = F)

data.table::fwrite(d_epistasis_mean,
                   EPISTASIS_WEIGHTED_MEAN_FILE, sep = ",", 
                   col.names = F, row.names = F)

data.table::fwrite(d_epistasis_density_nomolcomp,
                   EPISTASIS_WEIGHTED_NOMOLCOMP_DENSITY_FILE, sep = ",", col.names = F, row.names = F)

data.table::fwrite(d_epistasis_mean_nomolcomp,
                   EPISTASIS_WEIGHTED_NOMOLCOMP_MEAN_FILE, sep = ",", col.names = F, row.names = F)
