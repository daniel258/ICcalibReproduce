IterSimCox <- function(n.sample, mu, lambda, alpha, beta0, n.points, gamma.z, gamma.q, pts.for.ints, n.int = 5, order = 2)
{
  my.data <- ICcalib:::SimCoxIntervalCensCox(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, gamma.z = gamma.z, gamma.q = gamma.q,
                                   beta0 = beta0, n.points = n.points)
  Z <- my.data$Z
  Q <- my.data$Q
  obs.tm <- my.data$obs.tm
  delta <- my.data$delta
  case.times <- obs.tm[delta==1]
  n.cases <- length(case.times)
  w <- my.data$w
  w.res <- my.data$w.res
  ##### Naive methods ########
  #### Create data frame with time-dependent covariate suitable for coxph (for naive methods, see LVCFdata and MIdIdata) #######
  df.lvcf <- ICcalib:::LVCFdata(w = w, w.res = w.res, obs.tm = obs.tm, delta = delta, Z = Z)
  df.midI <- ICcalib:::MidIdata(w = w, w.res = w.res, obs.tm = obs.tm, delta = delta, Z = Z)
  fit.lvcf <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X + Z1 + Z2 + Z3, data = df.lvcf)
  fit.midI <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X + Z1 + Z2 + Z3, data = df.midI)
  est.lvcf <- coef(fit.lvcf)
  est.midI <- coef(fit.midI)
  var.lvcf <- diag(fit.lvcf$var)
  var.midI <- diag(fit.midI$var)
  ##### Fit calibration and risk set calibration models and then calculate P(X(t)=1|history)       ############
  ###### for each person in the risk set, for all risk sets                                        ############
  cox.hz.times <- sort(unique(obs.tm))
  fit.cox <- tryCatch(ICcalib:::FitCalibCox(w = w, w.res = w.res, Q = Q, hz.times = cox.hz.times, n.int = n.int, order = order), error = function(e){e})
  fit.cox.fail <- 0
  if(inherits(fit.cox, "error")){ fit.cox.fail <- 1 }
  n.fits <- length(pts.for.ints)
  fit.cox.rs.ints <- tryCatch(ICcalib:::FitCalibCoxRSInts(w = w, w.res = w.res, Q = Q, hz.times = cox.hz.times, tm = obs.tm, event = delta, 
                                                pts.for.ints = pts.for.ints, n.int = n.int, order = order),           error = function(e){e})
  fit.cox.rs.ints.fail <- 0 
  if(inherits(fit.cox.rs.ints, "error")){ fit.cox.rs.ints.fail <- 1 } else {
    n.etas.per.fit <- vector(length = n.fits)
    for (j in 1:n.fits)
    {
      n.etas.per.fit[j] <- length(fit.cox.rs.ints[[j]]$b) + length(fit.cox.rs.ints[[j]]$g)
    }
  }
  est.cox.fail <- 0
  if (fit.cox.fail==0)
  {
  px.cox <- t(sapply(case.times, ICcalib:::CalcCoxCalibP, w = w, w.res =  w.res, fit.cox, hz.times = cox.hz.times,  Q = Q))
  px.cox.deriv <- sapply(case.times, ICcalib:::CalcCoxCalibPderiv, w = w, w.res =  w.res, fit.cox, hz.times = cox.hz.times,  Q = Q, simplify = "array")
  est.cox <- tryCatch(optim(par = rep(0,ncol(Z)+1),fn = ICcalib:::CoxLogLik,  tm = obs.tm, event = delta, ps = px.cox, Z = Z,
                            control = list(fnscale=-1))$par, error = function(e){e})
  if(inherits(est.cox, "error")){
    est.cox.fail <- 1
  } else {
    var.cox <- diag(ICcalib:::CalcVarParam(theta = est.cox, tm = obs.tm, event = delta, Z = Z, Q = Q,
                            ps = px.cox, ps.deriv = px.cox.deriv, w = w, w.res = w.res, fit.cox = fit.cox))
    }}
  est.cox.rs.fail <- 0
  if (fit.cox.rs.ints.fail==0) {
  px.cox.rs.ints <- t(sapply(case.times, ICcalib:::CalcCoxCalibRSIntsP, w = w, w.res =  w.res, fit.cox.rs.ints = fit.cox.rs.ints, 
                             hz.times = cox.hz.times,  Q = Q, pts.for.ints = pts.for.ints))
  px.cox.rs.ints.deriv <- sapply(case.times, ICcalib:::CalcCoxCalibPderivRSInsts, w = w, w.res =  w.res, fit.cox.rs.ints = fit.cox.rs.ints, 
                                 hz.times = cox.hz.times,  Q = Q,  tm = obs.tm,  pts.for.ints = pts.for.ints, n.etas.per.fit = n.etas.per.fit,
                                 simplify = "array")
  est.cox.rs <- tryCatch(optim(par = rep(0,ncol(Z)+1),fn = ICcalib:::CoxLogLik,  tm = obs.tm, event = delta, ps = px.cox.rs.ints, Z = Z,
                            control = list(fnscale=-1))$par, error = function(e){e})
  if(inherits(est.cox.rs, "error")){
    est.cox.rs.fail <- 1
  } else {
    var.cox.rs.ints <- ICcalib:::CalcVarParamRSInts(theta = est.cox.rs, tm = obs.tm, event = delta, Z = Z, Q = Q,
                            ps = px.cox.rs.ints, ps.deriv = px.cox.rs.ints.deriv, w = w, w.res = w.res, fit.cox.rs.ints = fit.cox.rs.ints,
                            pts.for.ints = pts.for.ints, n.etas.per.fit = n.etas.per.fit)
    var.cox.rs <- diag(var.cox.rs.ints)
    }
  }
  # Returns a vector
  return.vec <- c(n.sample, mu, lambda, alpha,  beta0, n.points, n.cases, # design parameters
                  est.lvcf, est.midI, est.cox, est.cox.rs, # parameter estimates 
                  var.lvcf, var.midI, var.cox, var.cox.rs) #variance estimates
  #CI using estimated variance can be calculated when summarizing data
  return(return.vec)
}


