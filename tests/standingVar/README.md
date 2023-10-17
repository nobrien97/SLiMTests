# Chapter 2: Adaptation by standing variation

We're running some parameter sweeps here across three variables:
- Recombination rate [0 - 0.5]
- Number of loci [1 - 1000]
- Mutational effect size variance [0.1 - 1]

To sample this space well with <= 0.05 correlations between parameters, we will use 1024 parameters

The goal is to see how increasing the number of loci might influence the production of non-additive variance
The expectation is that it should decrease, as the mutational target is so large and mutational effects
should be favoured to be very small to compensate for the higher effective mutation rate. When mutational effect
size variance is small, this should be particularly apparent - we should approach an infinitesimal model.
Recombination rate will further this - with free recombination, non-additivity should be less frequent.

To test this we need to measure the standard outputs:
- Heterozygosity
- Phenotype means (and a sample of individuals as well to get an estimated distribution)
- Allele frequencies and effect sizes (Note - Fix the fixation generation)

And also
- Heritability
- Additive variance