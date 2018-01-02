# This script examplfies the way simulations were carried for the scenario of PH calibration model with covariates. 
# Check the 'Scripts' folder for the actual scripts used in the cluster. 
rm(list= ls())
source('IterSimCox.R', echo=F) 
mu = 0.2 # Exponential censoring paramter (mean=10)
lambda = 0.1 # Scale for the time-to-event gompartz distribution
alpha = 0.25 # Shape for the time-to-event gompartz distribution
gamma.q = c(log(0.75), log(2.5)) # Coefficients of Q in the Calibration model 
gamma.z = log(1.5) # Coefficient of Z in the Calibration model 
beta0 <- log(0.2) # For example
n.points <- 2 # Mstar
pts.for.ints <- seq(0,4.5,0.25) # values for grouping risk-sets for risk-set calibration
n.int <- 10 # Defaults for the I-splines
order <- 3 # Defaults for the I-splines
n.sample <- 200 # Sample size - 200 just for having a fast example

set.seed(114)

IterSimCox(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, beta0 = beta0, n.points = n.points, 
           gamma.z = gamma.z, gamma.q = gamma.q, pts.for.ints = pts.for.ints,  n.int = n.int, order = order)

# See IterSimCox for interptation of each number
# 200.00000000   0.20000000   0.10000000   0.25000000  -1.60943791   2.00000000  47.00000000  -1.41929056  -0.62790081   0.77259819   0.27009314  -2.08393665 
# Z1           Z2           Z3                                                                                                                      
# -0.55903948   0.69588310   0.31512152  -1.69014978  -0.52439421   0.67145697   0.27684345  -1.68523691  -0.53267907   0.68258125   0.28133657   0.32632343 
# 
# 0.09231453   0.08513218   0.02655686   0.31173906   0.09247648   0.08707242   0.02703852   0.40519246   0.09165485   0.08551875   0.02736641   0.41148317 
# 
# 0.09187881   0.08500212   0.02712825 
