# Calculates statistics involving mutation data
# Runs per model across all seeds of that model
library(dplyr)
library(tibble)

# Get command line arguments
## 1: model (ODE, K, or Additive)
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])

# Path to write output
WRITE_PATH <- "/scratch/ht96/nb9894/standingVar/calcMutationStats/"
EPISTASIS_FILE <- paste0(WRITE_PATH, "out_e_", model, ".csv")
EFFECTS_FILE <- paste0(WRITE_PATH, "out_fx_", model, ".csv")

# Load mutation data from database
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      dbname = paste0(dataPath, "standingVarMuts.db"))

# Quantitative data
d_adapted <- tbl(con, "slim_qg") %>%
  filter(gen >= 49500, modelindex == model) %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  filter(any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup()
d_adapted <- d_adapted %>% collect()

# If the model has 0 adapted runs, end early
if (nrow(d_adapted) == 0) {
  print(paste("Model ", run, " never adapted across all seeds, closing R."))
  q(save = "no")
  
}

# Mutation data
adapted_seeds <- unique(d_adapted$seed)
d_muts_adapted <- tbl(con, "slim_muts") %>% 
  filter(modelindex == model, gen >= 49500, seed %in% adapted_seeds) %>%
  select(!c(pos, chi))
d_muts_adapted <- d_muts_adapted %>% collect()


# Set data types for incorrect columns
d_muts_adapted$seed <- as.factor(d_muts_adapted$seed)
d_muts_adapted$modelindex <- as.factor(d_muts_adapted$modelindex)
d_muts_adapted$mutType <- as.factor(d_muts_adapted$mutType)
d_muts_adapted$fixGen <- as.numeric(d_muts_adapted$fixGen)
d_muts_adapted <- d_muts_adapted %>% 
  rename(value = effect)


d_adapted$seed <- as.factor(d_adapted$seed)
d_adapted$modelindex <- as.factor(d_adapted$modelindex)

# inner join the mutation w/ quantitative data
# Filter to include only adapted populations
d_com_adapted <- inner_join(d_adapted, d_muts_adapted, 
                            by = c("gen", "seed", "modelindex"))
d_com_adapted <- AddCombosToDF(d_com_adapted) 

d_fixed_adapted <- d_com_adapted %>% filter(!is.na(fixGen))

# Calculate phenotype effects
d_phenofx <- CalcPhenotypeEffects(d_com_adapted %>% filter(is.na(fixGen)),
                                  d_fixed_adapted)

d_epistasis <- PairwiseEpistasis(d_fixed_adapted,
                                 d_com_adapted %>% 
                                   filter(is.na(fixGen)) %>%
                                   select(gen, seed, modelindex, mutType, value),
                                 m = 10, n = 10, T)

d_epistasis <- AddCombosToDF(d_epistasis) 

# Write to table
RSQLite::dbAppendTable()

# dP/dt
sampleRate <- 50 # sample every 50 generations, so divide deltaP by 50
d_adapted %>%
  mutate(dPdT = deltaPheno / 50) -> d_adapted

d_dpdt <- d_adapted %>%
  filter(gen >= 50000) %>%
  mutate(gen = gen - 50000,
         optPerc = (phenomean - 1))    # percent to optimum

d_dpdt$optPerc <- cut(d_dpdt$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))
d_dpdt <- AddCombosToDF(d_dpdt)

d_dpdt %>%
  group_by(optPerc, modelindex, model, nloci, r, tau) %>%
  mutate(nloci = as.factor(nloci),
         r = as.factor(r),
         tau = as.factor(tau)) %>%
  summarise(meandPdT = mean(dPdT),
            CIdPdT = CI(dPdT)) -> d_dpdt_sum

ggplot(d_dpdt,
       aes(x = dPdT, y = optPerc, fill = model)) +
  facet_wrap(vars(nloci)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Change in phenotype per generation ($\\frac{\\delta P}{\\delta t}$)"), 
       y = "Closeness to optimum (%)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")


d_SFS <- CalcSFS(d_com_adapted)

ggplot(d_SFS,
       aes(x = freq, y = optPerc, height = after_stat(density), fill = model)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Allele frequency"), 
       y = "Closeness to optimum (%)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

