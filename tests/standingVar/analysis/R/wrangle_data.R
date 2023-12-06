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


# Get adaptive step order and attach step 0 (phenotype from before the first step in the walk)
d_fix_ranked <- RankFixations(d_fix_nar, T, d_fix %>% filter(modelindex == 2))
d_fix_ranked_add <- RankFixations(d_fix_add, F, d_fix %>% filter(modelindex == 1))

d_seg_ranked <- CalcNARPhenotypeEffects(d_seg_ranked, F, d_fix_adapted)

d_seg_ranked_add <- CalcAddEffects(d_seg_ranked_add, F, d_fix_adapted)


d_fix_ranked_combined <- rbind(d_fix_ranked, d_fix_ranked_add)
d_fix_ranked_combined$model <- if_else(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")

# Calculate fitness effects when they arose, 
# and when they first reached 50% or greater freq - not fixed
d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_muts_del$gen, d_muts_del$seed, d_muts_del$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_del <- inner_join(d_muts_del, d_qg_matched_seg, 
                        by = c("gen", "seed", "modelindex"))

# Calculate phenotypes
d_seg_del_add <- CalcAddEffects(d_seg_del %>% filter(modelindex == 1), isFixed = F,
                                dat_fixed = d_fix_adapted %>% filter(modelindex == 1))
d_seg_del <- CalcNARPhenotypeEffects(d_seg_del %>% filter(modelindex == 2), isFixed = F, 
                        dat_fixed = d_fix_adapted %>% filter(modelindex == 2))

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
