# Sep 21 ,2017
# Daniel Nevo
rm(list=ls())
library(foreach)
library(doParallel)
library(doRNG)

setwd("/home/dn84/ChangePointBin/Sims/Raux")
source('IterSimCox.R')


setwd("/home/dn84/ChangePointBin/Sims")

# True Paramter Values
mu = 0.2 # Exponential censoring paramter (mean=10)
lambda = 0.1 # Scale for the time-to-event gompartz distribution
alpha = 0.25 # Shape for the time-to-event gompartz distribution
gamma.q = c(log(0.75), log(2.5))
gamma.z = log(1.5)
beta0 <- log(0.5)
n.points <- 5
pts.for.ints <- seq(0,4.5,0.25)
n.int <- 10
order <- 3
n.sample <- 1000
n.sim <- 1000



cl <- makeCluster(12,outfile="")
registerDoParallel(cl)

my.seeds <- read.table("Raux/seeds.cox.txt",header = T, check.names = F)
seed <- my.seeds[which(as.numeric(rownames(my.seeds))==n.points), which(round(as.numeric(colnames(my.seeds)),1)==round(exp(beta0),1))] 
set.seed(seed)

Results <- foreach(i = 1:n.sim,.packages = c("ICcalib"), .combine = "rbind")   %dorng% {
  cat("i = ", i, "\n")                                              
       IterSimCox(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, beta0 = beta0, n.points = n.points, 
                  gamma.z = gamma.z, gamma.q = gamma.q, pts.for.ints = pts.for.ints,  n.int = n.int, order = order)
}
stopCluster(cl)
save.image(paste0("CPBcoxBeta",exp(beta0),"W",n.points,"N",n.sample,".RData"))
write.table(Results, paste0("CPBcoxBeta",exp(beta0),"W",n.points,"N",n.sample,".txt"))
sessionInfo()


