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

EPISTASIS_SIGN_CHANGES_NOMOLCOMP_FILE <- paste0(WRITE_PATH, "d_epi_sign_nomolcomp", model, ".csv")

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
  summarise(meanEW = mean(ew_s),
            sdEW = sd(ew_s),
            minEW = min(ew_s),
            maxEW = max(ew_s),
            q025EW = quantile(ew_s, probs = 0.025, na.rm = T, names = F),
            q25EW = quantile(ew_s, probs = 0.25, na.rm = T, names = F),
            q50EW = quantile(ew_s, probs = 0.5, na.rm = T, names = F),
            q75EW = quantile(ew_s, probs = 0.75, na.rm = T, names = F),
            q975EW = quantile(ew_s, probs = 0.975, na.rm = T, names = F),
            n = n(),  
            fracOverDriftBarrier = sum(abs(ew_s) > 1e-4) / n) -> d_epistasis_mean

d_epistasis %>%
  group_by(gen, modelindex, isAdapted) %>%
  summarise(meanEW = mean(ew_s),
            sdEW = sd(ew_s),
            minEW = min(ew_s),
            maxEW = max(ew_s),
            q025EW = quantile(ew_s, probs = 0.025, na.rm = T, names = F),
            q25EW = quantile(ew_s, probs = 0.25, na.rm = T, names = F),
            q50EW = quantile(ew_s, probs = 0.5, na.rm = T, names = F),
            q75EW = quantile(ew_s, probs = 0.75, na.rm = T, names = F),
            q975EW = quantile(ew_s, probs = 0.975, na.rm = T, names = F),
            n = n(),  
            fracOverDriftBarrier = sum(abs(ew_s) > 1e-4) / n) -> d_epistasis_mean_nomolcomp

# Get sign changes: change in average epistasis across mutations in a run
eps <- 1e-6

d_epistasis %>%
  group_by(gen, seed, modelindex) %>%
  summarise(ew = mean(ew), ew_s = mean(ew_s)) %>%
  ungroup() %>%
  arrange(gen) %>%
  group_by(seed, modelindex) %>%
  mutate(
    sign_EW = case_when(
      ew > eps ~ 1,
      ew < -eps ~ -1,
      TRUE ~ 0
    ),
    sign_EW_s = case_when(
      ew_s > eps ~ 1,
      ew_s < -eps ~ -1,
      TRUE ~ 0
    ),
    # Compare with previous non-zero sign
    prev_sign_EW = lag(sign_EW),
    prev_sign_EW_s = lag(sign_EW_s),
    sign_EW_change = sign_EW != prev_sign_EW & !is.na(prev_sign_EW),
    sign_EW_s_change = sign_EW_s != prev_sign_EW_s & !is.na(prev_sign_EW_s)) %>%
  summarise(n= n(), n_sign_EW_changes = sum(sign_EW_change, na.rm = TRUE),
            n_sign_EW_s_changes = sum(sign_EW_s_change, na.rm = TRUE)) -> d_signchanges

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_DENSITY_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_mean,
                   EPISTASIS_MEAN_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_density_nomolcomp,
                   EPISTASIS_NOMOLCOMP_DENSITY_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_mean_nomolcomp,
                   EPISTASIS_NOMOLCOMP_MEAN_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_signchanges,
                   EPISTASIS_SIGN_CHANGES_NOMOLCOMP_FILE, sep = ",", col.names = F, row.names = F)


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
  summarise(meanEW = mean(ew_s),
            sdEW = sd(ew_s),
            minEW = min(ew_s),
            maxEW = max(ew_s),
            q025EW = quantile(ew_s, probs = 0.025, na.rm = T, names = F),
            q25EW = quantile(ew_s, probs = 0.25, na.rm = T, names = F),
            q50EW = quantile(ew_s, probs = 0.5, na.rm = T, names = F),
            q75EW = quantile(ew_s, probs = 0.75, na.rm = T, names = F),
            q975EW = quantile(ew_s, probs = 0.975, na.rm = T, names = F),
            n = n(),  
            fracOverDriftBarrier = sum(abs(ew_s) > 1e-4) / n) -> d_epistasis_mean_nomolcomp

d_epistasis_freq %>%
  group_by(gen, modelindex) %>%
  summarise(meanEW = mean(ew_s),
            sdEW = sd(ew_s),
            minEW = min(ew_s),
            maxEW = max(ew_s),
            q025EW = quantile(ew_s, probs = 0.025, na.rm = T, names = F),
            q25EW = quantile(ew_s, probs = 0.25, na.rm = T, names = F),
            q50EW = quantile(ew_s, probs = 0.5, na.rm = T, names = F),
            q75EW = quantile(ew_s, probs = 0.75, na.rm = T, names = F),
            q975EW = quantile(ew_s, probs = 0.975, na.rm = T, names = F),
            n = n(),  
            fracOverDriftBarrier = sum(abs(ew_s) > 1e-4) / n) -> d_epistasis_mean

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


# extract the model from the database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "epistasis_newEpi.db"))

d_epistasis <- tbl(con, "tab_epistasis_new_epi") %>% 
  filter(modelindex == model) %>%
  collect()

# Make sure it has all worked properly
if (nrow(d_epistasis) < 1) {
  print(paste("Model ", run, " couldn't get epistasis data (likely no adapted pops)."))
  q(save = "no")
}

# if additive, we need to rearrange some of the data since it wasn't stored correctly:
# mutType_ab is not right and everything needs to be shifted
isAdditive <- ((model - 1) %% 3 == 0)

if (isAdditive) {
  d_epistasis <- d_epistasis %>%
    mutate(ep = ew,
           ew = Pab,
           Pab = Pb,
           Pb = Pa,
           Pa = Pwt,
           Pwt = wab,
           wab = wb,
           wb = wa,
           wa = as.numeric(mutType_ab),
           mutType_ab = "3_3")
}

# Add optPerc labels to use instead of gen
dpdt <- read.csv(paste0(GDATA_PATH, "d_dpdt.csv"), header = F)
dpdt <- dpdt %>% filter(V2 == model)
d_epistasis$optPerc <- rep(dpdt$V1, each = nrow(d_epistasis)/nrow(dpdt))

# Calculate density for each optPerc across Pa - Pwt, Pb - Pwt, Pab - Pwt, ew, ep
# and across log fitness
d_epistasis %>% 
  nest_by(optPerc, modelindex, mutType_ab) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         Pa = list(data.frame(density(data$Pa - data$Pwt)[c("x", "y")])),
         Pb = list(data.frame(density(data$Pb - data$Pwt)[c("x", "y")])),
         Pab = list(data.frame(density(data$Pab - data$Pwt)[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ep = list(data.frame(density(data$ep)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density

# Mean: compare the mean epistasis between models to pick the models to compare
# in density curves
d_epistasis %>%
  group_by(optPerc, modelindex) %>%
  summarise(meanEP = mean(ep),
            sdEP = sd(ep),
            meanEW = mean(ew),
            sdEW = sd(ew),
            minEW = min(ew),
            maxEW = max(ew),
            q025EW = quantile(ew, probs = 0.025, na.rm = T, names = F),
            q25EW = quantile(ew, probs = 0.25, na.rm = T, names = F),
            q50EW = quantile(ew, probs = 0.5, na.rm = T, names = F),
            q75EW = quantile(ew, probs = 0.75, na.rm = T, names = F),
            q975EW = quantile(ew, probs = 0.975, na.rm = T, names = F),
            n = n(),  
            fracOverDriftBarrier = sum(abs(ew) > 1e-4) / n) -> d_epistasis_mean

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_DENSITY_FILE, sep = ",", col.names = F, row.names = F)
data.table::fwrite(d_epistasis_mean,
                   EPISTASIS_MEAN_FILE, sep = ",", col.names = F, row.names = F)


rm(d_epistasis)

# Repeat for frequency-weighted data
d_epistasis_freq <- tbl(con, "tab_epistasis_freq_new_epi") %>% 
  filter(modelindex == model) %>%
  collect()

# Same problem with mutType_ab being wrong
if (isAdditive) {
  d_epistasis_freq <- d_epistasis_freq %>%
    mutate(ep = ew,
           ew = Pab,
           Pab = Pb,
           Pb = Pa,
           Pa = Pwt,
           Pwt = wab,
           wab = wb,
           wb = wa,
           wa = as.numeric(mutType_ab),
           mutType_ab = "3_3")
}

d_epistasis_freq$optPerc <- rep(dpdt$V1, each = nrow(d_epistasis_freq)/nrow(dpdt))

d_epistasis_freq %>% 
  nest_by(optPerc, modelindex, mutType_ab) %>%
  mutate(wa = list(data.frame(density(log(data$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(data$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(data$wab))[c("x", "y")])),
         Pa = list(data.frame(density(data$Pa - data$Pwt)[c("x", "y")])),
         Pb = list(data.frame(density(data$Pb - data$Pwt)[c("x", "y")])),
         Pab = list(data.frame(density(data$Pab - data$Pwt)[c("x", "y")])),
         ew = list(data.frame(density(data$ew)[c("x", "y")])),
         ep = list(data.frame(density(data$ep)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density

d_epistasis_freq %>%
  group_by(optPerc, modelindex) %>%
  summarise(meanEP = mean(ep),
            sdEP = sd(ep),
            meanEW = mean(ew),
            sdEW = sd(ew),
            minEW = min(ew),
            maxEW = max(ew),
            q025EW = quantile(ew, probs = 0.025, na.rm = T, names = F),
            q25EW = quantile(ew, probs = 0.25, na.rm = T, names = F),
            q50EW = quantile(ew, probs = 0.5, na.rm = T, names = F),
            q75EW = quantile(ew, probs = 0.75, na.rm = T, names = F),
            q975EW = quantile(ew, probs = 0.975, na.rm = T, names = F),
            n = n(),  
            fracOverDriftBarrier = sum(abs(ew) > 1e-4) / n) -> d_epistasis_mean

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_WEIGHTED_DENSITY_FILE, sep = ",", 
                   col.names = F, row.names = F)

data.table::fwrite(d_epistasis_mean,
                   EPISTASIS_WEIGHTED_MEAN_FILE, sep = ",", 
                   col.names = F, row.names = F)
