
## This script reproduce the results in Table 2 of the paper

#### Preperations ####

## Clean workspace and load pacakges 
rm(list = ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(xtable)

## Set working directory
setwd("Simulations/PH/Results")


## load all txt files amd save in one table/data.frame


all.beta0 <- log(c(1/7, 0.2, 0.5, 1, 2, 5, 7))
all.n.points <- c(2,5,10)
my.n.sample <- 1000
keep.list <- c("all.files", "all.results.cox", "free", "all.beta0", "all.n.points", "my.n.sample","my.beta0",
               "my.n.points","j","keep.list")
all.files <- list.files()
all.results.cox <- matrix(nr = 21*1000, nc = 39) # 1000 iterations per simulation scenario, 21 simulations scenarios
free <- 1  
j <- 0

for (my.beta0 in all.beta0)
  {
  for (my.n.points in all.n.points)
{
  j <- j + 1
  cat("j = ", j, "\n")
  yes <- paste0("CPBcoxBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".txt") %in% all.files
  if (yes)
    {
    temp <- as.matrix(read.table(paste0("CPBcoxBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".txt")) )
    all.results.cox[free:(free+999), ] <- temp
    free <- free + 1000
    rm(list=setdiff(ls(), keep.list))
  }
  }}

all.results.df <- as.data.frame(all.results.cox) 
# return.vec <- c(n.sample, mu, lambda, alpha,  beta0, n.points, n.cases, 
#                 est.lvcf, est.midI, est.cox, est.cox.rs,
#                 var.lvcf, var.midI, var.cox, var.cox.rs)

# Check IterSimCox.R to see what is each column
colnames(all.results.df) <-c("n.sample", "mu", "lambda", "alpha", "beta0", "n.points", "n.cases", 
                              "lvcf.beta", paste0("lvcf.gamma",1:3),"midI.beta", paste0("midI.gamma",1:3),
                              "cox.beta", paste0("cox.gamma",1:3), "cox.rs.beta", paste0("cox.rs.gamma",1:3),
                              "var.beta.lvcf", paste0("var.gamma.lvcf",1:3), "var.beta.midI", paste0("var.gamma.midI",1:3),  
                              "var.beta.cox", paste0("var.gamma.vox",1:3), "var.beta.cox.rs", paste0("var.gamma.cox.rs",1:3))
#View(all.results.df) 

##################################################
# calculate coverage rates for CI in all methods
all.results.df$lvcf.in.ci <- (all.results.df$beta0 > all.results.df$lvcf.beta - sqrt(all.results.df$var.beta.lvcf)*1.96) & 
  (all.results.df$beta0 < all.results.df$lvcf.beta + sqrt(all.results.df$var.beta.lvcf)*1.96) 

all.results.df$cox.in.ci <- (all.results.df$beta0 > all.results.df$cox.beta - sqrt(all.results.df$var.beta.cox)*1.96) & 
                           (all.results.df$beta0 < all.results.df$cox.beta + sqrt(all.results.df$var.beta.cox)*1.96) 
all.results.df$cox.rs.in.ci <- (all.results.df$beta0 > all.results.df$cox.rs.beta - sqrt(all.results.df$var.beta.cox.rs)*1.96) & 
  (all.results.df$beta0 < all.results.df$cox.rs.beta + sqrt(all.results.df$var.beta.cox.rs)*1.96) 
mean(all.results.df$lvcf.in.ci, na.rm = T)
mean(all.results.df$cox.in.ci, na.rm = T)
mean(all.results.df$cox.rs.in.ci, na.rm = T)

##################################################

#### Summarize Results #####
# Calculate estimated standard erros from estimated variances
all.results.df$se.lvcf <- sqrt(all.results.df$var.beta.lvcf)
all.results.df$se.cox <- sqrt(all.results.df$var.beta.cox)
all.results.df$se.cox.rs <- sqrt(all.results.df$var.beta.cox.rs)

# Summarize seprately for mean, se, se.hat and ci
#  Use "melt()" to change the data format from wide to long
melted.est <- all.results.df %>% select(beta0, n.points, lvcf.beta, cox.beta, cox.rs.beta) %>%  melt(id.vars = c("beta0","n.points"))
colnames(melted.est)[3] <- "Method"
est <- ddply(melted.est, c("beta0", "n.points",  "Method"), summarise,
               MEAN = round(mean(value, na.rm = T),3), EMP.SD = round(sd(value, na.rm = T),3))

melted.se <- all.results.df %>% select(beta0, n.points, se.lvcf, se.cox, se.cox.rs)%>%  melt(id.vars = c("beta0","n.points"))
colnames(melted.se)[3] <- "Method"
se <- ddply(melted.se, c("beta0", "n.points",  "Method"), summarise, se = round(mean(value, na.rm = T), 3))

melted.ci <- all.results.df %>% select(beta0, n.points, lvcf.in.ci, cox.in.ci, cox.rs.in.ci) %>%  melt(id.vars = c("beta0","n.points"))
colnames(melted.ci)[3] <- "Method"
ci <- ddply(melted.ci, c("beta0", "n.points",  "Method"), summarise, cp95 = round(mean(value, na.rm = T), 3))

# Put all togather 
final <- cbind(est, SE = se$se, CP95 = ci$cp95)
levels(final$Method) <- c("LVCF", "PH-OC", "PH-RSC")
final$beta0 <- final$beta0 %>%  round(3)

# For Table 2, keep only certain values of beta0 and n.points (Mstar)
final.pos <- final %>% filter(beta0>=0)
final.pos.no.10 <- final.pos %>% filter(n.points<10)

# Table 2
print(xtable(final.pos.no.10, digits = 3), include.rownames =  F)






