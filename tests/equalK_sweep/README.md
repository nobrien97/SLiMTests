# Equal mutational target, individual tracking, hypercube sweep

This experiment runs a hypercube sweep across $n_{loci}$ and $\sigma_{loci}$ to create an adaptedness/fitness
surface based on molecular architecture.
The experiment also has equal mutational targets for each molecular trait - i.e. $\alpha_Z$, $\beta_Z$, $K_Z$, and $K_{XZ}$. It also tracks the mutations and phenotypes of 10 individuals per generation.
In a sense, it combines the newTestCross/equalK and indTrack experiments I previously ran, with hypercube sampling instead of factorial levels of $n_{loci}$ and $\sigma_{loci}$.

There are 128 total combinations in the hypercube and 48 replicates for a total 6144 simulations. This is too large to store all the haplotype information, so instead I'll store the population data from each simulation at burn-in so that it can be replicated when needed. If I need the heritability information from a particular simulation, then I can load in a given model, run only the 2000 generations of test time, get the haplotypes, and go from there. 