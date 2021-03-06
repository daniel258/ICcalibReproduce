### Analyze cmv data used by Goggins et al. (1999)
# The data was downloaded from "http://hedwig.mgh.harvard.edu/biostatistics/node/32"
# interval_censr_data.zip contained cmvshed.dat and cmvshed.sas 
# We ran cmvshed.sas while adding a proc export to get the data at a .csv format
# The revised SAs script is called cmvshedDN.sas 
# Variable names are given in cmvshedDN.sas and VarNames.txt

rm(list = ls())

library(ICcalib)
packageVersion("ICcalib")
#[1] ‘1.0.5’
library(xtable)
library(dplyr)

a <- proc.time()

# For reproduciblity. We are going to use the bootstrap.
set.seed(314) 


# Change folder path to where cmvshed.csv is stored on your computer
setwd("/Users/nhdne/Dropbox/ChangePointBin/Data/cmvshed")
setwd("C:/Users/User/Dropbox/ChangePointBin/Data/cmvshed")

cmv.data.raw <- read.csv("cmvshed.csv", header = T)
colnames(cmv.data.raw) <- tolower(colnames(cmv.data.raw))
cmv.data.raw <- cmv.data.raw %>% dplyr::select(pid, sex, race, cmvind, sheddind, deathcen, deathdat, cmvdiag , digdate, startrx, offstdat,
                                sbas181d, blposind, bloodedt, bldnegdt, bldposdt,  urposind, urineedt, urnnegdt, urnposdt)

#### For all date variables, create new variables with R-friendly date formats
date.vars <- c("deathdat", "digdate", "startrx", "offstdat", "sbas181d",
               "bloodedt", "bldnegdt", "bldposdt",  "urineedt", "urnnegdt", "urnposdt")
cmv.data <- cmv.data.raw
for (j in 1:length(date.vars))
{
  cmv.data[date.vars[j]] <- cmv.data.raw[date.vars[j]] %>% unlist %>%  as.character %>% as.Date("%d-%b-%y")
}

### Keep only participants CMV-free at the beginning of the study  (As in Goggins et AL., 1999)
# We can use the observations with no shedding, unlike Goggins et AL., 1999. 
cmv.data <- cmv.data %>% filter(cmvind==0)

# One observation had study start time = study end time, so we took it out.
cmv.data  <- cmv.data  %>% filter(offstdat-startrx>0)

# Correct NA to 0 in the event indicator variable (cmvdiag)
cmv.data$cmvdiag[is.na(cmv.data$cmvdiag)] <- 0 

# We get 37 events over 221 observations.  
nrow(cmv.data) # [1] 221
cmv.data %>% select(cmvdiag) %>% sum  # [1] 37

#### Create actual times (in approximate months)

# For blood and urn shedding, create left and right time points.
# Replace NA with Inf for right points of the interval, and replace NA with zero for left side of the interval
# For time-to-cmv-diagnosis, use diagnoisis date if cmvdiag==1 (cases) and end of study if cmvdiag==0 (event-free participants)
cmv.data <- cmv.data %>% dplyr::mutate(bld.lft =  ifelse(is.na(bldnegdt), 0, ((bldnegdt - startrx)/30)),
                              bld.rgt =  ifelse(is.na(bldposdt), Inf ,((bldposdt - startrx)/30)),
                              urn.lft =  ifelse(is.na(urnnegdt), 0, ((urnnegdt - startrx)/30)),
                              urn.rgt =  ifelse(is.na(urnposdt), Inf, ((urnposdt- startrx)/30)),
                              obs.tm = ifelse(cmvdiag==1, ((digdate - startrx)/30)+1, ((offstdat - startrx)/30)))
# For some people with positive shedding status at time 0,  there is one day mismatch in the left and/or right point of the interval. 
# We change their intervals to [0,1/30]  
#That is, for people with positive shedding at time 0 [not positive diagnosis! those where excluded], change right time to epsilon.
cmv.data$bld.rgt[cmv.data$bld.rgt<=0] <- 1/30
cmv.data$urn.rgt[cmv.data$urn.rgt<=0] <- 1/30
cmv.data$bld.lft[cmv.data$bld.lft<=0] <- 0
cmv.data$urn.lft[cmv.data$urn.lft<=0] <- 0

# Format the data to fit to our R package
w.bld <- cmv.data %>% dplyr::select(bld.lft, bld.rgt) %>% as.matrix
w.urn <- cmv.data %>% dplyr::select(urn.lft, urn.rgt) %>% as.matrix
w.res.bld <- matrix(nr = nrow(w.bld), nc = 2, 0)
w.res.bld[, 2] <- ifelse(w.bld[, 2]==Inf, Inf, 1)
w.res.urn <- matrix(nr = nrow(w.urn), nc = 2, 0)
w.res.urn[, 2] <- ifelse(w.urn[, 2]==Inf, Inf, 1)

# For one observation (pid = 250237), right point of w.urn equal to obs.tm. That is, participant left the study on the day it urine shedding was discovered.
# This would cause problem later: I am adding 1/30 to the observed time to prevent problems for the naive methods later
sum(cmv.data$obs.tm==w.bld[, 2] | cmv.data$obs.tm==w.urn[, 2])

cmv.data$obs.tm[cmv.data$obs.tm==w.bld[, 2] | cmv.data$obs.tm==w.urn[, 2]] <- cmv.data$obs.tm[cmv.data$obs.tm==w.bld[, 2] | cmv.data$obs.tm==w.urn[, 2]] + 1
obs.tm <- cmv.data$obs.tm
delta <- cmv.data$cmvdiag

######### Analysis of  CMV **blood** shedding and CMV diagnosis ##########

# Create datasets for the naive methods
df.lvcf.bld <-   ICcalib:::LVCFdata(w = w.bld, w.res = w.res.bld, obs.tm = obs.tm, delta = delta)
df.midI.bld <-   ICcalib:::MidIdata(w = w.bld, w.res = w.res.bld, obs.tm = obs.tm, delta = delta)

# LVCF +  MidI model fitting
fit.lvcf.bld <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.lvcf.bld)
fit.midI.bld <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.midI.bld)

## Nonparametric calibration model
fit.npmle.bld <- ICcalib::FitCalibNpmle(w = w.bld, w.res = w.res.bld)

## Calculate probability matrices. These matrices represent P(X(t)=1|history) using different estimators.
# Number of rows: number of cases , number of columns: sample size. So these probabilities are calculated for all the sample at each case time
case.times <- obs.tm[delta]
px.np.rs.bld <- t(sapply(case.times, ICcalib:::CalcNpmleRSP, w = w.bld, w.res =  w.res.bld, obs.tm = obs.tm))
px.np.bld <- t(sapply(case.times, ICcalib:::CalcNpmleCalibP, w = w.bld, w.res =  w.res.bld, fit.npmle = fit.npmle.bld))

###  Calculate estimates
est.np.bld <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px.np.bld, interval = c(-50, 50), maximum = T)$maximum
est.np.calib.rs.bld <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px.np.rs.bld, interval = c(-50, 50), maximum = T)$maximum
est.lvcf.bld <- fit.lvcf.bld$coefficients
est.midI.bld <- fit.midI.bld$coefficients
all.est.bld <- c(est.lvcf.bld, est.midI.bld, est.np.bld, est.np.calib.rs.bld)
all.est.bld %>% round(2)

### SE, p-values, CIs, and final results' table 
boot.np.bld <- ICcalib::CalcVarNpmle(tm = obs.tm, event = delta, w = w.bld , w.res = w.res.bld, BS = 1000)
boot.np.rs.bld <- ICcalib::CalcVarNpmleRS(tm = obs.tm, event = delta, w = w.bld , w.res = w.res.bld, BS = 1000)
var.beta.np.bld <- boot.np.bld$v
var.beta.np.rs.bld <- boot.np.rs.bld$v
ci.beta.boot.np.bld <- boot.np.bld$ci
ci.beta.boot.np.rs.bld <- boot.np.rs.bld$ci
out.bld <- matrix(nr = 4, nc = 4)
rownames(out.bld) <- c("LVCF", "MidI", "NP.calib", "NP.RS")
colnames(out.bld) <- c("Est", "HR", "SE(Est)", "CI for HR")
out.bld[ , 1] <- all.est.bld  %>% round(2)
out.bld[ , 2] <- exp(all.est.bld) %>% round(2)
out.bld[ , 3] <- c(sqrt(vcov(fit.lvcf.bld)), sqrt(vcov(fit.midI.bld)), sqrt(var.beta.np.bld), sqrt(var.beta.np.rs.bld)) %>% round(3)
#out.bld[ , 4] <- 2 * pnorm(-abs(out.bld[,1]/out.bld[,3])) # decided not to show p-values, they are ridiculously small 
out.bld[c(1,2) , 4] <- apply(out.bld[c(1, 2), c(1, 3)], 1, 
                                   function(x) paste0("(", round(exp(x[1]-1.96*x[2]), 2), ", ", round(exp(x[1]+1.96*x[2]), 2), ")"))
out.bld[3 , 4] <-  paste0("(", round(exp(ci.beta.boot.np.bld[1]), 2), ", ", round(exp(ci.beta.boot.np.bld[2]), 2), ")")
out.bld[4 , 4] <-  paste0("(", round(exp(ci.beta.boot.np.rs.bld[1]), 2), ", ", round(exp(ci.beta.boot.np.rs.bld[2]), 2), ")")

######### Analysis of  CMV **urine** shedding and CMV diagnosis ##########

# Create datasets for the naive methods
df.lvcf.urn <-   ICcalib:::LVCFdata(w = w.urn, w.res = w.res.urn, obs.tm = obs.tm, delta = delta)
df.midI.urn <-   ICcalib:::MidIdata(w = w.urn, w.res = w.res.urn, obs.tm = obs.tm, delta = delta)

# LVCF +  MidI model fitting
fit.lvcf.urn <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.lvcf.urn)
fit.midI.urn <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.midI.urn)

## Calibration model
fit.npmle.urn <- ICcalib::FitCalibNpmle(w = w.urn, w.res = w.res.urn)

## Calculate probability matrices. These matrices represent P(X(t)=1|history) using different estimators.
# Number of rows: number of cases , number of columns: sample size. So these probabilities are calculated for all the sample at each case time
px.np.rs.urn <- t(sapply(case.times, ICcalib:::CalcNpmleRSP, w = w.urn, w.res =  w.res.urn, obs.tm = obs.tm))
px.np.urn <- t(sapply(case.times, ICcalib:::CalcNpmleCalibP, w = w.urn, w.res =  w.res.urn, fit.npmle = fit.npmle.urn))

###  Calculate estimates
est.np.urn <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px.np.urn, interval = c(-100,100), maximum = T)$maximum
est.np.calib.rs.urn <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px.np.rs.urn, interval = c(-50,50), maximum = T)$maximum
est.lvcf.urn <- fit.lvcf.urn$coefficients
est.midI.urn <- fit.midI.urn$coefficients

all.est.urn <- c(est.lvcf.urn, est.midI.urn, est.np.urn, est.np.calib.rs.urn)
all.est.urn %>% round(2)

### SE, p-values, CIs, and final results' table 
boot.np.urn <- ICcalib::CalcVarNpmle(tm = obs.tm, event = delta, w = w.urn , w.res = w.res.urn, BS = 1000)
boot.np.rs.urn <- ICcalib::CalcVarNpmleRS(tm = obs.tm, event = delta, w = w.urn , w.res = w.res.urn, BS = 1000)
var.beta.np.urn <- boot.np.urn$v
var.beta.np.rs.urn <- boot.np.rs.urn$v
ci.beta.boot.np.urn <- boot.np.urn$ci
ci.beta.boot.np.rs.urn <- boot.np.rs.urn$ci
out.urn <- matrix(nr = 4, nc = 4)
rownames(out.urn) <- c("LVCF", "MidI", "NP.calib", "NP.RS")
colnames(out.urn) <- c("Est", "HR", "SE(Est)", "CI for HR")
out.urn[ , 1] <- all.est.urn  %>% round(2)
out.urn[ , 2] <- exp(all.est.urn) %>% round(2)
out.urn[ , 3] <- c(sqrt(vcov(fit.lvcf.urn)), sqrt(vcov(fit.midI.urn)), sqrt(var.beta.np.urn), sqrt(var.beta.np.rs.urn)) %>% round(3)
#out.urn[ , 4] <- 2 * pnorm(-abs(out.urn[,1]/out.urn[,3])) # I decided not to show p-values, they are ridiculously small 
out.urn[c(1,2), 4] <- apply(out.urn[c(1, 2), c(1, 3)], 1, 
                                 function(x) paste0("(", round(exp(x[1]-1.96*x[2]), 2), ", ", round(exp(x[1]+1.96*x[2]), 2), ")"))
out.urn[3 , 4] <-  paste0("(", round(exp(ci.beta.boot.np.urn[1]), 2), ", ", round(exp(ci.beta.boot.np.urn[2]), 2), ")")
out.urn[4 , 4] <-  paste0("(", round(exp(ci.beta.boot.np.rs.urn[1]), 2), ", ", round(exp(ci.beta.boot.np.rs.urn[2]), 2), ")")

# First four rows in Table A.3 in the supplementary materials 
# The fifth row was taken from from Goggins et al 1999
xtable(cbind(out.urn, out.bld))
b <- proc.time()
b - a 

# On my iMac: 
# user     system elapsed 
# 583.294   2.492 608.385 
# ~10 minutes. Most of the effort is due to using the bootstrap
# MY iMac specs:
# Processor: 2.9GHz Intel Core i5
# Memory: 8 GB 1600 MHz DDR3
# If you have a "weaker" computer, you may change  BS=1000 to another number (BS=100, 200 or 500) 
# to make sure the code run on your computer with no errors. The SE(est) and CI columns may be different though


packageVersion("ICcalib")
