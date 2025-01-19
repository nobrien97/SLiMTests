# Find functions to calculate new fitness optima
library(mvtnorm)

# FFBH initial value when all mol comps = 1
startTraits = c(1.54156, 0.179451, 1.03305, 0.00635902)

# First: selection strength found so that a 10% shift on any given trait axis reduces fitness by 10%
optSigma = c(0.335821, 0.0390923, 0.225045, 0.00138527)

# test that this is the case
dnorm(x = startTraits * 1.1, mean = startTraits, sd = optSigma) / 
  dnorm(x = startTraits, mean = startTraits, sd = optSigma) 


# Now we need to find a shift (x) in all traits so that fitness is reduced by 10%, given our selection strengths
# multivariate normal with diagonal Sigma is just a product of univariate normals
n = length(startTraits)
d = 0.1

# Distance formula

perTraitValue <- function(mu, sigma, p, n) {
  return(sqrt(-2 * sigma^2 * log(p)/n) + mu)
} 

newVals = perTraitValue(startTraits, optSigma, 1 - d, n)

dmvnorm(x = newVals, mean = startTraits, sigma = diag(optSigma^2)) / 
  dmvnorm(x = startTraits, mean = startTraits, sigma = diag(optSigma^2)) 

