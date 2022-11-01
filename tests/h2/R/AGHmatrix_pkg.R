library(AGHmatrix)
source("helpFns.R")

haplos <- read.csv("~/Desktop/haplos.csv", header = F)
names(haplos) <- paste0("q_", 1:ncol(haplos))
genos <- hapToGen(haplos)

Amat <- Amatrix(ped)

A <- Gmatrix(as.matrix(genos), method = "VanRaden")
D <- Gmatrix(as.matrix(genos), method = "Vitezica")
E <- A*A

