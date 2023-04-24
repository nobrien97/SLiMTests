# Fixed KZ/KXZ

Fixing K values for simplicity - allows us to just plot alpha/beta maps
Note: for the nloci argument, it is multiplied by 4 because each molecular trait is still active, its
just that after being sampled, KZ/KXZ mutations are disabled - they still exist in the QTL map, they just
don't do anything and never cause a m5/m6 mutation. Hence, to get 10 loci for aZ and bZ, you need nloci = 40

Also need to get haplotype data, sampling every 50 generations for 10k generations
144 models * 10000/50 = 28800 samples

96000 samples was 397.4 GB for haplotype files in the h2_moltrait_fix experiment
So, 28000/96000 * 397.4 = 119.22GB