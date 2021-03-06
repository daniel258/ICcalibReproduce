# Sept ,2017
# Daniel Nevo
rm(list=ls())
library(foreach)
library(doParallel)
library(doRNG)

setwd("/home/dn84/ChangePointBin/Sims/Raux")
source('IterSimSingle.R')


setwd("/home/dn84/ChangePointBin/Sims")

# True Paramter Values
mu = 0.2 # Exponential censoring paramter (mean=10)
lambda = 0.1 # Scale for the time-to-event gompartz distribution
alpha = 0.25 # Shape for the time-to-event gompartz distribution
weib.shape = 1.5 # Shape for the weibull distribution for time to change
weib.scale = 7.5  # Scale for the weibull distribution for time to change
BS = 200 # Bootstrap replications for variance estimation
beta0 <- log(5)
n.points <- 10
n.sample <- 1000
n.sim <- 1000

cl <- makeCluster(12,outfile="")
registerDoParallel(cl)

# Load prior seeds used for set.seed, and pick the one corresponding to beta0 and n.points
my.seeds <- read.table("Raux/seeds.simple.weib.txt",header = T, check.names = F)
seed <- my.seeds[which(as.numeric(rownames(my.seeds))==n.points), which(round(as.numeric(colnames(my.seeds)),1)==round(exp(beta0),1))] 
set.seed(seed)
Results <- foreach(i = 1:n.sim,.packages = c("ICcalib"),
                   .combine = "rbind")   %dorng% {
  cat("i = ", i, "\n")                                              
  IterSimSingle(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, 
                                               weib.shape = weib.shape, weib.scale = weib.scale, 
                                               beta = beta0, n.points = n.points, BS = BS)
}
stopCluster(cl)
save.image(paste0("CPBSingleWeibBeta",exp(beta0),"W",n.points,"N",n.sample,".RData"))
write.table(Results, paste0("CPBSingleWeibBeta",exp(beta0),"W",n.points,"N",n.sample,".txt"))
sessionInfo()

