# Equal mutational variance
Need to measure if mutational variance is equivalent between additive, K+ and K-
To fix this, adjust the tau parameter along each molecular component based on their sensitivities:

| Component | Sensitivity (Si = ∂Z∗/∂pi) | Squared S2 | Contribution to variance |
|-----------|----------------------------|------------|--------------------------|
| KXZ       | −11.60                     | 134.60     |   87.10%                 |
| αZ        | −3.40                      | 11.60      |   7.50%                  |
| βZ        | +2.90                      | 8.40       |  5.40%                   |
| KZ        | +0.09                      | 0.0081     | 0.005%                   |

Each tau value needs to be calibrated so they have equal phenotypic effects, based n tau = 0.0125

| Component  | Calibrated τi |
|------------|---------------|
| KXZ        | 2.97 × 10−7   |
| αZ         | 3.45 × 10−6   |
| βZ         | 4.76 × 10−6   |
| KZ         |  4.94 × 10−3  |
