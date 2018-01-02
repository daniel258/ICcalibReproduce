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
pexp.rates <- c(0.025, 0.3, 0.1, 0.025)
pexp.ts <- c(0, 2, 3, 4)
BS = 200 # Bootstrap replications for variance estimation
beta0 <- log(0.5)
n.points <- 10
n.sample <- 1000
n.sim <- 1000

cl <- makeCluster(12,outfile="")
registerDoParallel(cl)

# Load prior seeds used for set.seed, and pick the one corresponding to beta0 and n.points
my.seeds <- read.table("Raux/seeds.simple.pexp.txt",header = T, check.names = F)
seed <- my.seeds[which(as.numeric(rownames(my.seeds))==n.points), which(as.numeric(colnames(my.seeds))==exp(beta0))] 
set.seed(seed)
Results <- foreach(i = 1:n.sim,.packages = c("ICcalib"),
                   .combine = "rbind")   %dorng% {
  cat("i = ", i, "\n")                                              
  IterSimSinglePexp(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, 
                    beta = beta0, n.points = n.points, BS = BS, pexp.rates = pexp.rates, pexp.ts = pexp.ts)
                   }
stopCluster(cl)
save.image(paste0("CPBSinglePexpBeta",exp(beta0),"W",n.points,"N",n.sample,".RData"))
write.table(Results,paste0("CPBSinglePexpBeta",exp(beta0),"W",n.points,"N",n.sample,".txt"))
sessionInfo()


