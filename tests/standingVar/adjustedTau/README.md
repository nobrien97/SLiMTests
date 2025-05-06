# Equal mutational variance
Need to measure if mutational variance is equivalent between additive, K+ and K-
To fix this, adjust the tau parameter along each molecular component based on their sensitivities:

| Component | Sensitivity (Si = ∂Z∗/∂pi) | Squared S2 | Contribution to variance |
|-----------|----------------------------|------------|--------------------------|
| KXZ       | −11.60                     | 134.60     |   87.10%                 |
| αZ        | −3.40                      | 11.60      |   7.50%                  |
| βZ        | +2.90                      | 8.40       |  5.40%                   |
| KZ        | +0.09                      | 0.0081     | 0.005%                   |

Each tau value needs to be calibrated so they have equal phenotypic effects, based on tau = 0.0125
The raw weights are given by 
$$w_i = 1/S_i^2$$

where $S_i^2$ is the squared sensitivity of molecular component i.
From the squared sensitivities, need to normalise around tau = 0.0125 with

$$\lambda = \frac{0.0125 * 4}{\sum w_i} = \frac{0.05}{123.67} \approx 0.000404 $$

where 4 is the number of components and w_i is the weight of molecular component i

The final tau for each component is given by
$$\tau_i = \lambda w_i$$

| Component  | Weight   | Calibrated τi |
|------------|----------|---------------|
| KXZ        | 0.00743  | 3.005e-6      |
| αZ         | 0.0865   | 3.497e-5      |
| βZ         | 0.1189   | 4.807e-5      |
| KZ         | 123.457  | 4.991e-2      |
