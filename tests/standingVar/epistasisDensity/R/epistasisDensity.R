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
GDATA_PATH <- "/g/data/ht96/nb9894/standingVar/"
WRITE_PATH <- "/scratch/ht96/nb9894/standingVar/epistasisDensity/"
EPISTASIS_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_density", model, ".csv")
EPISTASIS_WEIGHTED_DENSITY_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_density", model, ".csv")
EPISTASIS_HIST_FILE <- paste0(WRITE_PATH, "d_epi_hist", model, ".csv")
EPISTASIS_WEIGHTED_HIST_FILE <- paste0(WRITE_PATH, "d_epi_freqweight_hist", model, ".csv")

# extract the model from the database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(GDATA_PATH, "epistasis.db"))

d_epistasis <- tbl(con, "tab_epistasis") %>% 
  filter(modelindex == model) %>%
  collect()

# Add optPerc labels to use instead of gen
dpdt <- read.csv(paste0(GDATA_PATH, "d_dpdt.csv"), header = F)
dpdt <- dpdt %>% filter(V2 == model)
d_epistasis$optPerc <- rep(dpdt$V1, each = nrow(d_epistasis)/nrow(dpdt))

# Calculate density for each optPerc across Pa - Pwt, Pb - Pwt, Pab - Pwt, ew, ep
# and across log fitness

# mock data
# d_epistasis <- data.frame(
#   optPerc = rep(dpdt$V1, each = 4800/5),
#   wa = runif(4800),
#   wb = runif(4800),
#   wab = runif(4800),
#   Pwt = rnorm(4800),
#   Pa = rnorm(4800),
#   Pb = rnorm(4800),
#   Pab = rnorm(4800)
# )
# d_epistasis$ew <- log(d_epistasis$wab) - (log(d_epistasis$wa) + log(d_epistasis$wb))
# d_epistasis$ep <- (d_epistasis$Pab - d_epistasis$Pwt) - 
#   ((d_epistasis$Pa - d_epistasis$Pwt) + (d_epistasis$Pb - d_epistasis$Pwt))


d_epistasis %>% 
  nest_by(optPerc) %>%
  mutate(wa = list(data.frame(density(log(d_epistasis$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(d_epistasis$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(d_epistasis$wab))[c("x", "y")])),
         Pa = list(data.frame(density(d_epistasis$Pa - d_epistasis$Pwt)[c("x", "y")])),
         Pb = list(data.frame(density(d_epistasis$Pb - d_epistasis$Pwt)[c("x", "y")])),
         Pab = list(data.frame(density(d_epistasis$Pab - d_epistasis$Pwt)[c("x", "y")])),
         ew = list(data.frame(density(d_epistasis$ew)[c("x", "y")])),
         ep = list(data.frame(density(d_epistasis$ep)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_density

# ggplot(d_epistasis, aes(x = Pab_x, y = Pab_y)) +
#   geom_line()

# write
data.table::fwrite(d_epistasis_density, 
                   EPISTASIS_DENSITY_FILE, sep = ",", col.names = F, row.names = F)

# Repeat for histograms: count 
d_epistasis %>% 
  group_by(optPerc) %>%
  { breaks <- hist() }
  mutate()
d_epistasis %>% 
  nest_by(optPerc) %>%
  mutate(wa = list(data.frame(hist(log(d_epistasis$wa))[c("x", "y")])),
         wb = list(data.frame(density(log(d_epistasis$wb))[c("x", "y")])),
         wab = list(data.frame(density(log(d_epistasis$wab))[c("x", "y")])),
         Pa = list(data.frame(density(d_epistasis$Pa - d_epistasis$Pwt)[c("x", "y")])),
         Pb = list(data.frame(density(d_epistasis$Pb - d_epistasis$Pwt)[c("x", "y")])),
         Pab = list(data.frame(density(d_epistasis$Pab - d_epistasis$Pwt)[c("x", "y")])),
         ew = list(data.frame(density(d_epistasis$ew)[c("x", "y")])),
         ep = list(data.frame(density(d_epistasis$ep)[c("x", "y")]))) %>%
  select(-data) %>% 
  unnest(cols = c(wa, wb, wab, 
                  Pa, Pb, Pab,
                  ew, ep), names_sep = "_") -> d_epistasis_hist


rm(d_epistasis)

# Repeat for frequency-weighted data
d_epistasis_freq <- tbl(con, "tab_epistasis_freq") %>% 
  filter(modelindex == model)
