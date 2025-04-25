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
GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/calcMutationStats/"
WRITE_PATH <- "/scratch/ht96/nb9894/standingVar/epistasisDensity/"
EPISTASIS_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_density_nomut", model, ".csv")
EPISTASIS_WEIGHTED_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_density_nomut", model, ".csv")

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
  nest_by(timePoint, modelindex) %>%
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

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_DENSITY_FILE, sep = ",", col.names = F, row.names = F)


rm(d_epistasis)

# Repeat for frequency-weighted data
d_epistasis_freq <- tbl(con, "tab_epistasis_freq") %>% 
  filter(modelindex == model) %>%
  collect()

d_epistasis_freq$timePoint <- rep(dpdt$V1, each = nrow(d_epistasis_freq)/nrow(dpdt))

d_epistasis_freq %>% 
  nest_by(timePoint, modelindex) %>%
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

# write
data.table::fwrite(d_epistasis_density,
                   EPISTASIS_WEIGHTED_DENSITY_FILE, sep = ",", 
                   col.names = F, row.names = F)
