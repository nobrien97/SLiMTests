# Transform/summarise data - 
# calculates fitness effects, ranks fixations in the adaptive walk,
# tidies data for plotting

# load combo info
d_combos <- read.table("../../R/combos.csv", header = F,
                       col.names = c("nloci", "tau", "r", "model"))

# load trait evolution data
d_qg <- read.table(paste0(dataPath, "slim_qg.csv"), header = F, 
                   sep = ",", colClasses = c("integer", "factor", "factor", 
                                             rep("numeric", times = 12)), 
                   col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                 "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                 "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                   fill = T)

 # Add predictors
d_qg <- AddCombosToDF(d_qg) 

# Filter for models which adapted (within 10% of optimum by the end of the simulation)
d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg
d_adapted <- d_qg %>% filter(isAdapted)

# Check how many there are
View(d_qg %>% group_by(model, nloci, tau, r) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted)))


# load in observed heterozygosity data
d_het <- read.table(paste0(dataPath, "slim_locusHo.csv"), 
                     header = F, sep = ",", 
                     colClasses = c("integer", "factor", "factor", rep("numeric", 256)), 
                                    col.names = c("gen", "seed", "modelindex", paste0("Ho_l", 1:256)), 
                                    fill = T)

# Pivot to calculate mean heterozygosity across the loci
d_het %>% pivot_longer(cols = -c(gen, seed, modelindex), names_to = "locus", values_to = "Ho") %>%
  group_by(modelindex) %>%
  summarise(meanHo = mean(Ho, na.rm = T),
            CIHo = CI(Ho, na.rm = T))

d_Ho <- d_het %>% 
  pivot_longer(cols = -c(gen, seed, modelindex), names_to = "locus", values_to = "Ho")

d_Ho <- AddCombosToDF(d_Ho) 

# Summarise, mean H_O over time
d_Ho %>%
  group_by(gen, model, nloci, tau, r) %>%
  summarise(meanHo = mean(Ho, na.rm = T),
            CIHo = CI(Ho, na.rm = T)) -> d_Ho_sum

# Load mutation data
d_muts <- read.table(paste0(dataPath, "slim_muts.csv"), header = F, 
                   sep = ",", colClasses = c("integer", "factor", "factor", 
                                             "factor", rep("integer", times = 4),
                                             rep("numeric", times = 3),
                                             rep("integer", times = 2)), 
                   col.names = c("gen", "seed", "modelindex", "mutType", "mutID",
                                 "pos", "constraint", "originGen", "value", "chi",
                                 "Freq", "Count", "fixGen"), 
                   fill = T)

# Filter to include only adapted populations
d_muts_adapted <- d_muts %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_adapted$seed, d_adapted$modelindex))

# Combine with d_qg
d_com_adapted <- inner_join(d_adapted, d_muts_adapted, by = c("gen", "seed", "modelindex"))

# Testing
d_com_adapted <- d_com_adapted %>%
  filter(modelindex == 25 | modelindex == 26)

d_fixed_adapted <- d_com_adapted %>% filter(!is.na(fixGen), gen >= 50000)
# Calculate the fitness effects in additive populations
d_add_fx <- CalcAddEffects(d_com_adapted %>% filter(model == "Add", 
                                                    is.na(fixGen),
                                                    gen >= 50000),
                           d_fixed_adapted %>% filter(model == "Add"))

# Calculate fitness effects for NAR populations
d_nar_fx <- CalcNARPhenotypeEffects(d_com_adapted %>% filter(model != "Add",
                                                              is.na(fixGen),
                                                              gen >= 50000),
                                     d_fixed_adapted %>% filter(model != "Add"))

# Calculate pairwise epistasis
d_add_e <- PairwiseEpistasisAdditive(d_fixed_adapted %>% filter(model == "Add"),
                          d_com_adapted %>% filter(model == "Add",
                                                   is.na(fixGen)) %>%
                            select(gen, seed, modelindex, value) %>%
                            rename(a = value),
                          d_com_adapted %>% filter(model == "Add",
                                                   is.na(fixGen)) %>%
                            select(gen, seed, modelindex, value) %>%
                            rename(b = value))

d_net_e <- PairwiseEpistasisNAR(d_fixed_adapted %>% filter(model != "Add"),
                                     d_com_adapted %>% filter(model != "Add",
                                                              is.na(fixGen)) %>%
                                       select(gen, seed, modelindex, mutType, value) %>%
                                       rename(a = value),
                                     d_com_adapted %>% filter(model != "Add",
                                                              is.na(fixGen)) %>%
                                       select(gen, seed, modelindex, mutType, value) %>%
                                       rename(b = value))

d_add_e_av <- d_add_e %>%
  group_by(gen, modelindex) %>%
  summarise(meanEW = mean(ew),
            CIEW = CI(ew),
            meanEP = mean(ep),
            CIEP = CI(ep))


d_net_e_av <- d_net_e %>%
  mutate(mutTypeAB = interaction(mutType_a, mutType_b)) %>%
  group_by(gen, modelindex, mutTypeAB) %>%
  summarise(meanEW = mean(ew),
            CIEW = CI(ew),
            meanEP = mean(ep),
            CIEP = CI(ep))

ggplot(d_net_e_av %>% mutate(gen = gen - 50000) %>%
         pivot_longer(c(meanEP, meanEW, CIEP, CIEW), 
                      names_to = c(".value", "EType"),
                      names_pattern = "(..*)(..$)"),
       aes(x = gen, y = mean, colour = EType)) +
  facet_wrap(vars(mutTypeAB)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - CI, ymax = mean + CI, colour = NULL,
                  fill = EType), alpha = 0.4) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Phenotype", "Fitness")) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  guides(fill = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean epistasis", colour = "Epistasis type") +
  theme_bw() +
  theme(text = element_text(size = 16))

ggplot(d_add_e_av %>% mutate(gen = gen - 50000) %>%
         pivot_longer(c(meanEP, meanEW, CIEP, CIEW), 
                      names_to = c(".value", "EType"),
                      names_pattern = "(..*)(..$)"),
       aes(x = gen, y = mean, colour = EType)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - CI, ymax = mean + CI, colour = NULL,
                  fill = EType), alpha = 0.4) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Phenotype", "Fitness")) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  guides(fill = "none") +
  labs(x = "Generations post-optimum shift", y = "Mean epistasis", 
       colour = "Epistasis type") +
  theme_bw() +
  theme(text = element_text(size = 16))

# dP/dt
sampleRate <- 50 # sample every 50 generations, so divide deltaP by 50
d_adapted %>%
  mutate(dPdT = deltaPheno / 50) -> d_adapted

d_dpdt <- d_adapted %>%
  filter(gen >= 50000) %>%
  mutate(gen = gen - 50000,
         optPerc = (phenomean - 1))    # percent to optimum

d_dpdt$optPerc <- cut(d_dpdt$optPerc, c(-Inf, 0.25, 0.5, 0.75, Inf))

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

# SFS
d_SFS <- CalcSFS(d_com_adapted)

ggplot(d_SFS,
       aes(x = Freq, y = optPerc, height = after_stat(density), fill = model)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Allele frequency"), 
       y = "Closeness to optimum (%)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

# Deviation plots
# Calculate control means
# deviation from grand mean across all groups (nloci, r, tau)
d_dpdt_ctrl <- d_dpdt %>%
  group_by(optPerc, model) %>% 
  summarise(mean_dpdt_ctrl = mean(dPdT),
            CI_dpdt_ctrl = se(dPdT),
  )

# Add to main dataframe
d_dpdt_total <- inner_join(d_dpdt, d_dpdt_ctrl, 
                           by = c("optPerc", "model"))

# We do want directionality, so not squared deviation
d_dpdt_total <- d_dpdt_total %>% 
  group_by(optPerc, model) %>% 
  mutate(devdPdT = (dPdT - mean_dpdt_ctrl))

# Plot deviations
ggplot(d_dpdt_total,
       aes(x = devdPdT, y = optPerc, colour = model)) +
  facet_grid(nloci~r) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Recombination rate",
                                         breaks = NULL, labels = NULL)) +
  guides(y.sec = guide_axis_manual(
    breaks = NULL, labels = NULL, title = "Number of loci"
  )) +
  geom_jitter(size = 0.5) +
  # geom_density_ridges(alpha = 0.4) +
  scale_colour_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Deviation from overall mean change in phenotype per generation ($\\frac{\\delta P}{\\delta t}$)"), 
       y = "Closeness to optimum (%)", 
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

# SFS Deviations
# Calculate control means: mean proportion of mutations within
# each bin
# deviation from grand mean across all groups (nloci, r, tau)
d_SFS <- d_SFS %>%
  group_by(optPerc, seed, modelindex) %>%
  mutate(totalMuts = n()) %>%
  group_by(optPerc, seed, modelindex, FreqBin) %>%
  mutate(binProp = n()/totalMuts) %>%
  ungroup()

d_SFS_ctrl <- d_SFS %>%
  group_by(optPerc, model, FreqBin) %>%
  summarise(mean_SFS_ctrl = mean(binProp),
            CI_SFS_ctrl = CI(binProp)
  )

# Add to main dataframe
d_SFS_total <- inner_join(d_SFS, d_SFS_ctrl, 
                           by = c("optPerc", "model", "FreqBin"))

# We do want directionality, so not squared deviation
d_SFS_total <- d_SFS_total %>% 
  group_by(optPerc, modelindex, FreqBin) %>% 
  mutate(devSFS = (binProp - mean_SFS_ctrl)) %>%
  ungroup() %>%
  select(optPerc, model, r, nloci, tau, devSFS, FreqBin) %>%
  distinct()

# Plot
ggplot(d_SFS_total,
       aes(x = FreqBin, y = optPerc, stat = devSFS, fill = model)) +
  facet_grid(nloci~tau) +
  geom_bar() +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "ODE", "K")) +
  labs(x = TeX("Allele frequency"), 
       y = "Change in Closeness to optimum (%)", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")



# Fitness calculations aren't correct for those before optimum shift: recalculate
d_seg_del[d_seg_del$gen < 50000,]$wAA <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$AA_pheno, 1, 0.05)
d_seg_del[d_seg_del$gen < 50000,]$wAa <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$Aa_pheno, 1, 0.05)
d_seg_del[d_seg_del$gen < 50000,]$waa <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$aa_pheno, 1, 0.05)
d_seg_del$avFit <- d_seg_del$wAa - d_seg_del$waa
d_seg_del$avFit_AA <- d_seg_del$wAA - d_seg_del$waa
d_seg_del$s <- d_seg_del$wAA - d_seg_del$waa

d_seg_del_add[d_seg_del_add$gen < 50000,]$wAA <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$AA_pheno, 1, 0.05)
d_seg_del_add[d_seg_del_add$gen < 50000,]$wAa <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$Aa_pheno, 1, 0.05)
d_seg_del_add[d_seg_del_add$gen < 50000,]$waa <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$aa_pheno, 1, 0.05)
d_seg_del_add$avFit <- d_seg_del_add$wAa - d_seg_del_add$waa
d_seg_del_add$avFit_AA <- d_seg_del_add$wAA - d_seg_del_add$waa
d_seg_del_add$s <- d_seg_del_add$wAA - d_seg_del_add$waa

d_seg_del <- rbind(d_seg_del, d_seg_del_add)
d_seg_del$model <- ifelse(d_seg_del$modelindex == 2, "NAR", "Additive")
rm(d_seg_del_add)

# Get the difference in selection coefficient after optimum shift
# as well as mean frequency just prior to optimum shift
d_seg_del %>%
  filter(gen == 49500 | gen == 50000) %>%
  arrange(gen) %>%
  group_by(seed, model, mutID) %>%
  filter(n() > 1) %>%
  summarise(diff_s = s[2] - s[1],
            Freq = Freq[1]) -> d_del_diffs

# Calculate the contributions of adaptive steps vs segregating mutations on phenotypic variation 
d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(molCompDiff = sum(abs(2 * value_aZ), na.rm = T) - sum(abs(2 * value_bZ), na.rm = T),
         molCompCorr = cor(abs(2 * value_aZ), abs(2 * value_bZ))) %>%
  ungroup() %>%
  dplyr::select(seed, molCompDiff, molCompCorr) %>%
  distinct(seed, .keep_all = T) -> d_molCompDiff
