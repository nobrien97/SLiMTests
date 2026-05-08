# Paper 1: Network topology constrains evolvability: developmental bias in the mutational landscape across gene regulatory motifs

These simulations extend chapter 3
5% drop in fitness now instead of 10%; this corrects the tau = 0.004 aligning with 100 additive beneficial mutations to get to the optimum again


- adjTau.slim runs an experiment for how much the results depend on mutational scaling. Mutation scale adjustment: 
from starting conditions measure M, then scale each tau by the covariance between each component and the trait on the total variance in the trait. Then continue burn-in as normal.
- randomisedStartsM.slim runs evolutionary simulations with randomised starting conditions, measuring the M matrix over the simulation
- parallelSel.slim runs simulations with selection directions chosen to be parallel to the trait correlations identified during neutral simulations
- orthSel.slim runs simulations with selection directions orthogonal to the trait correlations

parallelSel and orthSel were run for r = 0.1 as recombination was found to have similar effects across model treatments.
