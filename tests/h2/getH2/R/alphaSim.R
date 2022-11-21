library(AlphaSimR)

haplos <- read.csv("~/Desktop/slim_haplo2587819639178264576_0.csv", header = F)[, -(1:3)]
names(haplos) <- paste0("q_", seq_len(ncol(haplos)))

genMap <- data.frame(markerName=paste0("q_", seq_len(ncol(haplos))),
                     chromosome=seq(1, ncol(haplos)),
                     position=rep(1, times = ncol(haplos)))

ped <- read.csv("~/Desktop/slim_pedigree2587819639178264576_0.csv", header = F)[,-(1:3)]
names(ped) <- c("id", "mother", "father")
ped$id <- as.character(ped$id); ped$mother <- as.character(ped$mother); ped$father <- as.character(ped$father)

ped <- 

founderPop <- importHaplo(haplo=haplos, genMap=genMap, ploidy=2L, ped=ped)
