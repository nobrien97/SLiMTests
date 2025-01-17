# Multi-Motif pilot simulations

I've run some initial simulations to test each motif's ability to adapt following an optimum shift.
Each population burnt in to a state 10% away from the initial state of the ODE system (when all molecular components = 1).
The optimum shift was then a further 10% away from the burn-in optimum.

![](./plt_adapt_mdist.png)

Fig. 1 - The mean of mean Mahalanobis distances from the optimum across individuals and replicates. Mahalanobis distance accounts for differences in selection strengths across traits (attributable to differences in trait scale).

![](./plt_adapt_w.png)

Fig. 2 - The mean of mean fitnesses across individuals and replicates. Fitness was calculated via a multivariate normal fitness function.

Figs 1 and 2 show how populations adapted over time. Despite starting at the same fitness deficit, only the NAR populations were able to perfectly match the burn-in optimum, with all other motifs reaching stable points away from the optimum. These might represent local optima on the fitness landscape. The ramifications when the optimum shifts are severe drops in fitness for maladapted motifs because the optimum shift is relative to the burn-in optimum, of which these motifs are not adapted to. Recombination created larger drops in fitness after the optimum shift for the FFBH motif.

![](./plt_adapt_wvar.png)

Fig. 3 - The mean variance in fitness across across replicates.

Fig 3 shows motifs maintained different variance in fitness: FFBH motifs were most variable, whereas the NAR was the least. Recombination tended to slightly increase variance in fitness, but not by a great deal.

![](./plt_adapt_dist_pertrait.png)

Fig. 4 - Per-trait mean distances from the optimum adjusted for selection strength. Top label is the log10 recombination rate, bottom label is the trait. Note that the traits differ between motifs:
NAR and PAR traits are (1) response time and (2) steady state concentration; FFL-C1 traits are (1) response time, (2) response delay, and (3) steady state concentration; FFL-I1 traits are
(1) time to half-maximum expression, (2) maximum expression concentration, and (3) time above half-maximum expression; and FFBH traits are (1) time to half-maximum expression, (2) maximum expression concentration, (3) response time to the second steady state, and (4) the second steady state concentration.

Figure 4 shows how adaptation progressed across all traits for the models. FFBH motifs were limited in reaching the optimum among all traits. The response time to the second steady state especially was limited. Fluctuations in maximum expression concentration show volatility in this trait even with small mutational effects. FFBH populations showed very slow movement to the optimum under low and intermediate recombination in their first trait, the time to half-maximum expression.

PAR populations adapted quickly in steady state concentration, but response time remained unoptimizable. FFL-C1 were the opposite, being able to adjust their response time rapidly but being limited in their ability to adjust steady state concentration beyond a certain limit: the initial change was rapid, but steady state change capped at around 0.1 distance from the optimum. A similar effect was observed for the response delay.

FFL-I1 populations quickly reached the burn-in optimum for time to half-maximum expression, but when the optimum shifted, they remained stable at a maladapted trait value. They were unable to optimise their time above half-maximum expression, despite quickly optimising maximum expression concentration.


## Trait correlations

The adaptation difficulty might be due to trait correlations. Here are figures of the pairwise correlations between trait means. Each point in the plot is a trait mean in a replicate (288 total points in each figure; 3 recombination rates and 96 replicates). Note that because these are trait measurements taking during adaptation, these are a combination of selection and network derived correlations:

### NAR
![](./plt_trait_corr_nar.png)

There is a strong correlation between response time and steady state concentration. The two densities might correspond to time differences (e.g. at the optimum during burn-in vs the shifted optimum), but I need to verify this.

### PAR
![](./plt_trait_corr_par.png)

The same traits, but they appear much less correlated. A weaker correlation would slow adaptation in this case because the optimum shifts in an increasing direction for both traits simultaneously: hence, a mutation that changes both traits positively will be strongly favoured. This type of mutation might be less likely under the PAR.

### FFL-C1
![](./plt_trait_corr_fflc1.png)

The FFL-C1 has some strange patterns going on. Especially the correlation between response time and steady state concentration, which is almost logistic, recapitulating the shape of the Hill activation function in its circuit. There are two plateaus where steady state concentration is maintained at a given response time, and a sharp exponential switch between the two states at response time ~ 1.42s. The other two traits don't really seem correlated either, with response delay being confined to very small values, and response time being largely confined to one of two specific values (although there is some variation around those points).

### FFL-I1
![](./plt_trait_corr_ffli1.png)

Expected correlations: large maximum concentrations are correlated with longer times to reach half-max concentration. Similarly, longer times to reach half-max concentration leaves less time to remain above the half-max concentration threshold. There is some strange sampling going on in the distributions though, especially in time to half max vs time above half-max. Not sure what to make of this one.

### FFBH
![](./plt_trait_corr_ffbh.png)

For the most complex network, we see relatively weak correlations among traits except for a -0.273 correlation between the response time to second steady state and the time to half maximum concentration. This means that network configurations that take longer to reach the maximum tend to decline rapidly to the second steady state or vice versa. Maximum concentration is heavily restricted, with some outliers producing quite large concentrations, but generally being quite modest.

## Trait distance correlations
In addition to the trait correlations, we can look at correlations among distances to the optimum on each axis. Strong correlations here mean that when one trait is far from the optimum, the other trait tends to be as well. Negative correlations imply that when one trait is maladapted, the other is usually adapted.

### NAR
![](./plt_dist_corr_nar.png)

As expected, distance to the optimum response time and steady state concentration are tightly coupled.

### PAR
![](./plt_dist_corr_par.png)

The correlation from the NAR model is still there, but weaker: being close to the optimum in one trait does not always guarantee being close to the optimum in the other.

### FFL-C1
![](./plt_dist_corr_fflc1.png)

There is no correlation between response delay and steady state concentration, meaning that mutations affecting one trait don't always affect the other. Note that these two traits are the ones that hovered slightly away from the optimum.


### FFL-I1
![](./plt_dist_corr_ffli1.png)

Weaker correlations again (for two of the three combinations) and strange distributions. 

### FFBH
![](./plt_dist_corr_ffbh.png)

A variety of weaker correlations among half of the pairs, and almost discrete distributions, with distances outside of certain ranges among all traits being inaccessible.

## Conclusions

Together these results show that adaptation in these networks is strongly dependent on initial network configurations, the burn-in optima, and the direction of selection (and how that corresponds to trait correlations). Some future ideas: 

- Maybe we should randomise the direction of selection amongst traits? This will allow us to separate trait correlations from selection.
- It might be a good idea to randomise the initial states of the networks to reduce the effect of this initial state on adaptation. The optima can still be relative to this state (i.e. 10% potential increase in fitness when reaching the optimum).
- It would be good to map out the ruggedness of the fitness landscape in higher dimensions. Perhaps a simple metric like the number of beneficial mutations randomly sampled from a given point in the multidimensional space could give us a simple visual for this: e.g. if there are many beneficial mutations relative to the total sampled, it is likely on the "side" of a fitness incline, especially if the beneficial effect is large. Problem would be effectively sampling this space, but perhaps in a small area around where populations reach/equilibriate would suffice.
