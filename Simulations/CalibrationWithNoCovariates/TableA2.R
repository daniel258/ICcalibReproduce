## This script reproduce the results in Table A.2 in the Supplementary Materials
setwd("/Users/nhdne/Dropbox/ChangePointBin/Code/ICCalibReproduce/")
#### Preparations ####

## Clean workspace and load packages  
rm(list = ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(xtable)

## Set working directory
setwd("Simulations/CalibrationWithNoCovariates/Results")

## Load all txt files and save in one table/data.frame
all.beta0 <- log(c(0.2,0.5,1,2,5))
all.n.points <- c(2,5,10) # all.n.points are all possible values of M^\star
my.n.sample <- 1000
keep.list <- c("all.files", "all.results.weib","all.results.pexp", "free", "all.beta0", "all.n.points", "my.n.sample","my.beta0",
               "my.n.points","j","keep.list", "all.results.weib.df", "all.results.pexp.df")
# First Weibull distribution for V
setwd("Weib")
all.files <- list.files()
all.results.weib <- matrix(nr = 15*1000, nc = 29) 
all.results.weib[,1] <- 0
free <- 1  
j <- 0
for (my.beta0 in all.beta0)
    {
  for (my.n.points in all.n.points)
{
  j <- j + 1
  cat("j = ", j, "\n")
  yes <- paste0("CPBSingleWeibBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".txt") %in% all.files
  if (yes)
    {
    temp <- as.matrix(read.table(paste0("CPBSingleWeibBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".txt")) )
    all.results.weib[free:(free+999), 2:29] <- temp
    free <- free + 1000
    rm(list=setdiff(ls(), keep.list))
  }
  }}

all.results.weib.df <- as.data.frame(all.results.weib) 

colnames(all.results.weib.df) <- c("V.dist", "n.sample", "mu", "lambda", "alpha", "weib.shape", "weib.scale", 
                              "beta0", "n.points", "n.cases", "BS", 
                              "shape.est", "scale.est", 
                              "lvcf.beta", "midI.beta", "weib.calib.beta", "weib.rs.beta", "np.calib.beta", "np.rs.beta",
                              "var.beta.lvcf", "var.beta.midI", "var.beta.wb", "var.beta.wb.rs", "var.beta.np", "var.beta.np.rs", 
                              "ci.beta.np.l","ci.beta.np.r", "ci.beta.np.rs.l", "ci.beta.np.rs.r")


## Now Pexp (Piecewise Exponential)
setwd("../Pexp")

all.results.pexp <- matrix(nr = 15*1000, nc = 35) 
all.results.pexp[,1] <- 1

all.files <- list.files()
free <- 1  
j <- 0
for (my.beta0 in all.beta0)
{
  for (my.n.points in all.n.points)
  {
    j <- j + 1
    cat("j = ", j, "\n")
    yes <- paste0("CPBSinglePexpBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".txt") %in% all.files
    if (yes)
    {
#      load(paste0("CPBSinglePexpBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".RData"))
      temp <- as.matrix(read.table(paste0("CPBSinglePexpBeta",exp(my.beta0),"W",my.n.points,"N",my.n.sample,".txt")) )
      all.results.pexp[free:(free+999), 2:35] <- temp
      free <- free + 1000
      rm(list=setdiff(ls(), keep.list))
    }
  }}


all.results.pexp.df <- as.data.frame(all.results.pexp) 

colnames(all.results.pexp.df) <- c("V.dist", "n.sample", "mu", "lambda", "alpha", 
                                    paste0("pexp.rates",1:4), paste0("pexp.ts",1:4),
                                    "beta0", "n.points", "n.cases", "BS", 
                                    "shape.est", "scale.est", 
                                   "lvcf.beta", "midI.beta", "weib.calib.beta", "weib.rs.beta", "np.calib.beta", "np.rs.beta",
                                   "var.beta.lvcf", "var.beta.midI", "var.beta.wb", "var.beta.wb.rs", "var.beta.np", "var.beta.np.rs", 
                                   "ci.beta.np.l","ci.beta.np.r", "ci.beta.np.rs.l", "ci.beta.np.rs.r")

both.cols <- intersect(colnames(all.results.weib.df) , colnames(all.results.pexp.df))

all.results.df <- rbind(select(all.results.weib.df, match(both.cols,colnames(all.results.weib.df))), 
                        select(all.results.pexp.df, match(both.cols,colnames(all.results.pexp.df))))

##################################################
### calculate CI coverage rates for all methods ##
all.results.df$lvcf.in.ci <- (all.results.df$beta0 > all.results.df$lvcf.beta - sqrt(all.results.df$var.beta.lvcf)*1.96) & 
  (all.results.df$beta0 < all.results.df$lvcf.beta + sqrt(all.results.df$var.beta.lvcf)*1.96) 

all.results.df$midI.in.ci <- (all.results.df$beta0 > all.results.df$midI.beta - sqrt(all.results.df$var.beta.midI)*1.96) & 
  (all.results.df$beta0 < all.results.df$midI.beta + sqrt(all.results.df$var.beta.midI)*1.96) 

all.results.df$wb.in.ci <- (all.results.df$beta0 > all.results.df$weib.calib.beta - sqrt(all.results.df$var.beta.wb)*1.96) & 
                           (all.results.df$beta0 < all.results.df$weib.calib.beta + sqrt(all.results.df$var.beta.wb)*1.96) 

all.results.df$wb.rs.in.ci <- (all.results.df$beta0 > all.results.df$weib.rs.beta - sqrt(all.results.df$var.beta.wb.rs)*1.96) & 
  (all.results.df$beta0 < all.results.df$weib.rs.beta + sqrt(all.results.df$var.beta.wb.rs)*1.96) 

all.results.df$np.in.ci <- (all.results.df$beta0 > all.results.df$np.calib.beta - sqrt(all.results.df$var.beta.np)*1.96) & 
  (all.results.df$beta0 < all.results.df$np.calib.beta + sqrt(all.results.df$var.beta.np)*1.96) 

all.results.df$np.rs.in.ci <- (all.results.df$beta0 > all.results.df$np.rs.beta - sqrt(all.results.df$var.beta.np.rs)*1.96) & 
  (all.results.df$beta0 < all.results.df$np.rs.beta + sqrt(all.results.df$var.beta.np.rs)*1.96) 


all.results.df$boot.in.ci <- (all.results.df$beta0 > all.results.df$ci.beta.np.l) & 
  (all.results.df$beta0 < all.results.df$ci.beta.np.r) 

all.results.df$boot.rs.in.ci <- (all.results.df$beta0 > all.results.df$ci.beta.np.rs.l) & 
  (all.results.df$beta0 < all.results.df$ci.beta.np.rs.r) 

##################################################

#### Summarize Results

all.results.df %>% filter(beta0==log(2) & n.points==2)
all.results.df$sd.beta.lvcf <- sqrt(all.results.df$var.beta.lvcf)
all.results.df$sd.beta.midI <- sqrt(all.results.df$var.beta.midI)
all.results.df$sd.beta.wb <- sqrt(all.results.df$var.beta.wb)
all.results.df$sd.beta.np <- sqrt(all.results.df$var.beta.np)
all.results.df$sd.beta.wb.rs <- sqrt(all.results.df$var.beta.wb.rs)
all.results.df$sd.beta.np.rs <- sqrt(all.results.df$var.beta.np.rs)

melted.est <- all.results.df %>% select(V.dist, beta0, n.points, lvcf.beta, midI.beta, weib.calib.beta, weib.rs.beta, np.calib.beta, np.rs.beta) %>% 
  melt(id.vars = c("V.dist","beta0","n.points"))
melted.est.right <- melted.est %>% filter(abs(value)<5)
est.wb <- melted.est.right %>% filter(V.dist==0) %>% ddply(c("beta0", "n.points",  "variable"), summarise,
             mean.wb = round(mean(value, na.rm = T),3), 
#             median = round(median(value, na.rm = T),3), 
             emp.sd.pexp = round(sd(value, na.rm = T),3))
est.pexp <- melted.est.right %>% filter(V.dist==1) %>% ddply(c("beta0", "n.points",  "variable"), summarise,
                                                       mean.pexp = round(mean(value, na.rm = T),3), 
                                                       median = round(median(value, na.rm = T),3), 
                                                       emp.sd.pexp = round(sd(value, na.rm = T),3))

est <- cbind(est.wb, mean.pexp = est.pexp$mean.pexp, sd.pexp = est.pexp$emp.sd.pexp)

melted.se <- all.results.df %>% select(V.dist, beta0, n.points, sd.beta.lvcf, sd.beta.midI, sd.beta.wb, sd.beta.wb.rs ,sd.beta.np, sd.beta.np.rs) %>%  
  melt(id.vars = c("V.dist","beta0","n.points"))

se.wb <- melted.se %>% filter(V.dist==0) %>% ddply(c("beta0", "n.points",  "variable"), summarise, se = round(mean(value, na.rm = T), 3))
se.pexp <- melted.se %>% filter(V.dist==1) %>% ddply(c("beta0", "n.points",  "variable"), summarise, se = round(mean(value, na.rm = T), 3))
melted.ci <- all.results.df %>% select(V.dist, beta0, n.points, lvcf.in.ci, midI.in.ci, wb.in.ci, wb.rs.in.ci, np.in.ci, np.rs.in.ci) %>%  
  melt(id.vars = c("V.dist","beta0","n.points"))
ci.wb <- melted.ci %>% filter(V.dist==0) %>% ddply(c("beta0", "n.points",  "variable"), summarise, cp95.wb = round(mean(value, na.rm = T), 3))
ci.pexp <- melted.ci %>% filter(V.dist==1) %>% ddply(c("beta0", "n.points",  "variable"), summarise, cp95.pexp = round(mean(value, na.rm = T), 3))
final <- cbind(est.wb, emp.se.wb = se.wb$se, cp95.wb = ci.wb$cp95.wb,  
               mean.pexp = est.pexp$mean.pexp, emp.se.pexp = est.pexp$emp.sd.pexp, se.pexp = se.pexp$se, CP95 = ci.pexp$cp95.pexp)
levels(final$variable) <- c("LVCF","MidI", "WBcalib", "WBcalibRS", "NPcalib", "NPcalibRS")
final$beta0 <- final$beta0 %>%  round(3)
View(final)

print(xtable(final, digits = 3), include.rownames =  F)
# This is Table A.2



