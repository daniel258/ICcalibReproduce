# This script examplfies the way simulations were carried for the scenario of calibration model with no covariates. 
# Check the 'Scripts' folder for the actual scripts used in the cluster. 
# This scripts is for the case the true dsitribution of V was Weibull. Check the 'Scripts' folder for examples under piecewise exponential
rm(list= ls())
source('IterSimSingle.R', echo=F) 
# True Paramter Values
mu = 0.2 # Exponential censoring paramter (mean=10)
lambda = 0.1 # Scale for the time-to-event gompartz distribution
alpha = 0.25 # Shape for the time-to-event gompartz distribution
weib.shape = 1.5 # Shape for the weibull distribution for time to change
weib.scale = 7.5  # Scale for the weibull distribution for time to change
BS = 200 # Bootstrap replications for variance estimation
beta0 <- log(0.2) # For example
n.points <- 2 # Mstar
n.sample <- 200 # Sample size - 200 just for having a fast example

set.seed(411)

IterSimSingle(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, 
              weib.shape = weib.shape, weib.scale = weib.scale, 
              beta = beta0, n.points = n.points, BS = BS)

# See IterSimSingle.R for interptation of each number
# 200.0000000   0.2000000   0.1000000   0.2500000   1.5000000   7.5000000  -1.6094379   2.0000000  73.0000000 200.0000000   1.9202186   5.4726712  -0.6472622 
# X                                                                                                                                2.5%       97.5% 
#   -1.3774132  -0.7434752  -0.7891523  -1.9264580  -1.8233972   0.2921976   0.2810620   0.2796934   0.2702503   0.4916103   0.3975773  -2.7098727  -0.3838089 
# 2.5%       97.5% 
#   -2.7462596  -0.3027576 
# There were 50 or more warnings (use warnings() to see the first 50)
# > warnings()
# Warning messages:
#   1: In stats::optimize(f = CoxLogLikX, tm = tm.bs, event = event.bs,  ... :
#                           NA/Inf replaced by maximum positive value

# The warnings occuer during the bootstrap, when the optimzation of the likelihood is more delicate for some of the bootstrap repetitions 
# Nevertheless, this seems to no affect the obtained results too much and these are reasonable.
# This scripts surves as an exmaple for the simulation scripts, where the larger sample size  results in signficantly less issues.