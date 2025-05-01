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
GDATA_PATH <- "/g/data/ht96/nb9894/newMotifs/randomisedStarts/"
CALCMUTATIONSTATS_PATH <- paste0(GDATA_PATH, "calcMutationStats/")
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
                      dbname = paste0(CALCMUTATIONSTATS_PATH, "epistasis.db"))

d_epistasis <- tbl(con, "tab_epistasis") %>% 
  filter(modelindex == model) %>%
  collect()

# Read in qg
# Load mutation data from database
con_qg <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "newMotifsMuts.db"))

# Quantitative data
d_qg <- tbl(con_qg, "slim_qg") %>%
  filter(gen >= 49500, modelindex == model) %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & mean_w > 0.95)) %>%
  ungroup()
d_qg <- d_qg %>% collect()


# d_epistasis <- data.table::fread(paste0(CALCMUTATIONSTATS_PATH, "d_epistasis.csv"), header = F, 
#                           sep = ",", colClasses = c("integer", "factor", "factor", 
#                                                     rep("numeric", times = 7)), 
#                           col.names = c("gen", "seed", "modelindex", "mutType_ab",
#                                         "wa", "wb", "wab",
#                                         "wwt", "ew", "ew_s"), 
#                           fill = T)



# Make sure it has all worked properly
if (nrow(d_epistasis) < 1) {
  print(paste("Model ", model, " couldn't get epistasis data (likely no adapted pops)."))
  q(save = "no")
}

# Filter out rows with invalid values (e.g. fitness holes)
d_epistasis <- d_epistasis %>%
  mutate(ew = as.numeric(ew),
         ew_s = as.numeric(ew_s)) %>%
  filter(!is.na(ew) & !is.na(ew_s) & !is.infinite(ew) & !is.infinite(ew_s))

# Calculate density for each optPerc across Pa - Pwt, Pb - Pwt, Pab - Pwt, ew, ep
# and across log fitness
d_epistasis %>% 
  nest_by(gen, modelindex, mutType_ab) %>%
  mutate(ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(ew, ew_s), names_sep = "_") -> d_epistasis_density

# No molcomp
d_epistasis %>% 
  nest_by(gen, modelindex) %>%
  mutate(ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(ew, ew_s), names_sep = "_") -> d_epistasis_density_nomolcomp


# Mean: compare the mean epistasis between models to pick the models to compare
# in density curves

# Combine
d_epistasis <- inner_join(d_qg, d_epistasis, 
                            by = c("gen", "seed", "modelindex"))

d_epistasis %>%
  group_by(gen, modelindex, isAdapted, mutType_ab) %>%
  summarise(meanEW = mean(ew),
            sdEW = sd(ew),
            meanEW_s = mean(ew_s),
            sdEW_s = sd(ew_s),
            n = n()) -> d_epistasis_mean

d_epistasis %>%
  group_by(gen, modelindex, isAdapted) %>%
  summarise(meanEW = mean(ew),
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

if (nrow(d_epistasis) < 1) {
  print(paste("Model ", model, " couldn't get epistasis data (likely no adapted pops)."))
  q(save = "no")
}

# Filter out rows with invalid values (e.g. fitness holes)
d_epistasis_freq <- d_epistasis_freq %>%
  mutate(ew = as.numeric(ew),
         ew_s = as.numeric(ew_s)) %>%
  filter(!is.na(ew) & !is.na(ew_s) & !is.infinite(ew) & !is.infinite(ew_s))

d_epistasis_freq %>% 
  nest_by(gen, modelindex, mutType_ab) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         wwt = list(data.frame(density(log(data$wwt))[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, wwt,
                  ew, ew_s), names_sep = "_") -> d_epistasis_density

# No molcomp
d_epistasis_freq %>% 
  nest_by(gen, modelindex) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         wwt = list(data.frame(density(log(data$wwt))[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ew_s = list(data.frame(density(data$ew_s)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, wwt,
                  ew, ew_s), names_sep = "_") -> d_epistasis_density_nomolcomp


d_epistasis_freq %>%
  group_by(gen, modelindex, mutType_ab) %>%
  summarise(meanEW = mean(ew),
            sdEW = sd(ew),
            meanEW_s = mean(ew_s),
            sdEW_s = sd(ew_s),
            n = n()) -> d_epistasis_mean_nomolcomp

d_epistasis_freq %>%
  group_by(gen, modelindex) %>%
  summarise(meanEW = mean(ew),
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
