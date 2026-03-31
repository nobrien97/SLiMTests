# Neutral trait correlations of motifs

The idea of this model is to measure the correlations between traits in populations under drift.
Selection directions will then be chosen to be either parallel to the trait direction/correlation or perpendicular
e.g. if response time and steady state have correlation 0.5, then a parallel selection regime would have an optimum
with response time and steady state following that correlation (e.g. a randomly sampled multivariate normal optimum, 
with correlation structure informed by the estimated correlations).
A perpendicular regime would multiply the correlations by -1 and use that to sample the multivariate optimum.
Somewhere in-between would multiply by a number between 0 and -1.

