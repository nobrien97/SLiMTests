library(AlphaSimR)

haplos <- read.csv("~/Desktop/haplos.csv", header = F)
names(haplos) <- paste0("q_", 1:ncol(haplos))
ped <- read.csv("~/Desktop/ped.csv")

genMap <- read.csv("~/Desktop/test.csv")

pheno_means <- read.csv("~/Desktop/mean.csv", header = F)
rel <- read.csv("~/Desktop/rel.csv", header = F)

founderPop <- importHaplo(haplos, genMap, ploidy = 2, ped = ped)
founderPop

SP <- SimParam$new(founderPop)

SP$addTraitADE(10, mean = pheno_means[[1]], var = pheno_means[[2]], useVarA = T)

pop = newPop(founderPop)
genMean = meanG(pop)

for(generation in 1:20){
  pop = setPheno(pop, varE = 0, simParam = SP)
  pop = selectCross(pop=pop, nInd = 500, use="pheno", nCrosses=1000)
  genMean = c(genMean, meanG(pop))
}

plot(0:20, genMean, xlab="Generation", ylab="Mean Genetic Value", type="l")

pop = setPheno(pop, varE = 0, simParam = SP)
ans = genParam(pop, SP)
blup = fastRRBLUP(pop, simParam = SP)
pop = setEBV(pop, blup)
ebv(pop)


# Estimate heritability using K
source("helpFns.R")
library(heritability)

genos <- hapToGen(haplos)
# kinship <- popkin(as.matrix(genos), loci_on_cols = T)
rel_k <- read.csv("../k_mat.csv", header = F)
rel_k <- as.matrix(rel_k)
colnames(rel_k) <- paste0("i_", 1:nrow(genos))
rownames(rel_k) <- paste0("i_", 1:nrow(genos))
# rel <- as.matrix(rel)
# colnames(rel) <- paste0("i_", 1:nrow(genos))
# rownames(rel) <- paste0("i_", 1:nrow(genos))
phenos <- unlist(read.csv("~/Desktop/phenos.csv", header=F))
names(phenos) <- NULL
geno_names <- paste0("i_", 1:nrow(genos))

h2 <- marker_h2(scale(phenos), geno_names, K = rel_k, max.iter = 1000)

# Estimate heritability using more complicated linear models
library(rrBLUP)
library(lme4GS)
library(pedigreemm)
library(AGHmatrix)
library(GENESIS)
library(GWASTools)

colnames(rel) <- 1:nrow(genos)
rownames(rel) <- 1:nrow(genos)


geno <- MatrixGenotypeReader(genotype = t(genos),
                             snpID = 1:100,
                             chromosome = rep(1L, times = 100),
                             position = 1:100,
                             scanID = 1:100)
pheno <- data.frame(scanID = geno@scanID, pheno = phenos)
# Fit mixed model
nullmod <- fitNullModel(pheno, cov.mat = rel, outcome = "pheno", family = "gaussian")
varCompCI(nullmod, prop= TRUE)


# Add parents to pedigree
parents <- unique(unlist(ped[,c(2,3)]))
ped <- rbind(data.frame(id = parents, mother = 0, father = 0), ped)
Amat <- Amatrix(ped)
pedEdit <- editPed(sire = ped$father, dam = ped$mother, label = ped$id, verbose = T)
pedFinal <- with(pedEdit, pedigree(label=label, sire=sire, dam=dam))
Amat <- A.mat(pedFinal)

phenoDat <- data.frame(pheno = phenos, 
                       id = as.factor(ped$id)
                       )

ans1 <- mmer(pheno~1, random = ~ id, nIters=3, data=phenoDat, verbose=F)
summary(ans1)$varcomp
vpredict(ans1, h2 ~ V1 / (V1 + V2))

Dmat
