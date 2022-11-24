# h2 experiment: altering nloci

When doing the half-sib/full-sib experiment with 1000 QTLs, the network model had quite some difficulty reaching the optimum. Is this because of mutational load? Selection strength? In this set of experiments, I'm altering some parameters:

- nloci: [10, 100, 1000]
- width: [ 0.05, 0.105, 0.22 ] (5% 10% 20% selection strength)
- locisigma: [0.01, 0.1, 1]

3^3 models = 27 models

27 * 24 replicates = 648 models

Additive + network = 648*2 = 1296

Gadi queue: 1440 cores, ~12 hours per run

In future we can add more replicates