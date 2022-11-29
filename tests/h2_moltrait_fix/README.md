# Molecular trait effects on heritability

In this experiment, I'll be exploring the effects of each molecular trait on estimates of variance components.
To do this, I'll be looking at three levels of nloci (10, 100, 1000) under two levels of selection (5%, 20%).
Mutational effect variance has two levels (0.1, 1).

The SLiM model considers three of the molecular traits as usual, but disables mutation in one of them, leading to the default value (1) being retained over time.

This gives:
- nloci = [ 10, 100, 1000 ]
- s = [ 0.05, 0.22 ]
- molecular trait fixed = [ aZ, bZ, KZ, KXZ ]
- effect variance = [ 0.1, 1 ]

Which is 24 combinations, we'll use 50 seeds to get 2400 simulations