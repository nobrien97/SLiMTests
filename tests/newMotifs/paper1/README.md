# Paper 1: Network topology constrains evolvability: developmental bias in the mutational landscape across gene regulatory motifs

These simulations extend chapter 3
Mutational variance is rescaled - not tau = 0.004 anymore, it correctly accounts for the 10% drop rather than 5% drop in fitness at the optimum shift
tau = 0.006


- adjTau.slim runs an experiment for how much the results depend on mutational scaling. Mutation scale adjustment: 
from starting conditions measure M, then scale each tau by the covariance between each component and the trait on the total variance in the trait. Then continue burn-in as normal.
- randomisedStartsM.slim runs evolutionary simulations with randomised starting conditions, measuring the M matrix over the simulation
- parallelSel.slim runs simulations with selection directions chosen to be parallel to the trait correlations identified during neutral simulations
- orthoSel.slim runs simulations with selection directions orthogonal to the trait correlations