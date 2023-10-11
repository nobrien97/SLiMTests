# ODE vs Additive vs Multiplicative

A few short test simulations to compare the effect of the ODE on the additive
We additionally run a multiplicative model which describes how much the allele
transformation affects adaptation. By comparing this to the ODE result we can
see how much the delayed response to selection is due to the ODE vs multiplicative
allelic effects

The simulations also test a combined script that combines the additive and network
script running, which should be far easier to maintain

## Analysis
From the simulations, we're measuring
- Trait evolution over time (phenomean)
- Adaptive walks (steps)
- Allelic effects on fitness
- Mutant screen results