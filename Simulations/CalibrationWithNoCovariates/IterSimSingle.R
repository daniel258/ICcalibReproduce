IterSimSingle <- function(n.sample, mu, lambda, alpha, weib.shape, weib.scale, 
                          beta0, n.points, BS)
{
  my.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, weib.shape = weib.shape, weib.scale = weib.scale,
                                      beta0 = beta0, n.points = n.points)
  obs.tm <- my.data$obs.tm
  delta <- my.data$delta
  case.times <- obs.tm[delta==1]
  n.cases <- length(case.times)
  w <- my.data$w
  w.res <- my.data$w.res
  #### calculate P(X(t)=1|history) for each person in the risk set, for all risk sets for carry forward and midpoint #######
  df.lvcf <- ICcalib:::LVCFdata(w = w, w.res = w.res, obs.tm = obs.tm, delta = delta)
  df.midI <- ICcalib:::MidIdata(w = w, w.res = w.res, obs.tm = obs.tm, delta = delta)
  fit.lvcf <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.lvcf)
  fit.midI <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.midI)
  est.lvcf <- coef(fit.lvcf)
  est.midI <- coef(fit.midI)
  var.lvcf <- diag(fit.lvcf$var)
  var.midI <- diag(fit.midI$var)
  ##### Fit calibration and risk set calibration models and then calculate P(X(t)=1|history)       ############
  ###### for each person in the risk set, for all risk sets                                        ############
  ### Fit Weibull calibration model 
  calib.weib.params <- ICcalib:::FitCalibWeibull(w = w, w.res = w.res)
  px <- t(sapply(case.times,ICcalib:::CalcWeibullCalibP, w = w, w.res =  w.res, weib.params = calib.weib.params))
  # Calculate derivative matrices
  px.deriv.shape <- t(sapply(case.times, ICcalib:::CalcWeibullCalibPderivShape, w = w, w.res =  w.res, weib.params = calib.weib.params))
  px.deriv.scale <- t(sapply(case.times, ICcalib:::CalcWeibullCalibPderivScale, w = w, w.res =  w.res, weib.params = calib.weib.params))
  ### Fit Weibull Risk Set model 
  calib.rs.params <- ICcalib:::FitCalibWeibullRS(w = w, w.res = w.res, tm = obs.tm, event = delta)
  px.rs <- t(sapply(1:n.cases,function(m) {ICcalib:::CalcWeibullRSP(w = w, w.res = w.res, point = case.times[m], weib.params = calib.rs.params[m,])}))
  # Calculate derivative matrices
  px.rs.deriv.shape <- ICcalib:::CalcWeibullCalibPderivShapeRS(w = w, w.res = w.res, obs.tm = obs.tm, event = delta, weib.rs.params = calib.rs.params)
  px.rs.deriv.scale <- ICcalib:::CalcWeibullCalibPderivScaleRS(w = w, w.res = w.res, obs.tm = obs.tm, event = delta, weib.rs.params = calib.rs.params)
  ### Fit Nonparametric calibration model
  fit.npmle <- ICcalib:::FitCalibNpmle(w = w, w.res = w.res)
  px.np <- t(sapply(case.times, ICcalib:::CalcNpmleCalibP, w = w, w.res =  w.res, fit.npmle = fit.npmle))
  
  ### Fit Nonparametric risk set calibration model
  px.np.rs <- t(sapply(case.times, ICcalib:::CalcNpmleRSP, w = w, w.res =  w.res, obs.tm = obs.tm))
  
  
  ###### Calculate beta estimates #######
  est.weib.calib <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px, 
                                interval = c(-50,50), maximum = T)$maximum
  est.weib.calib.rs <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = my.data$delta, ps = px.rs, 
                                   interval = c(-50,50), maximum = T)$maximum
  est.np.calib <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = my.data$delta, ps = px.np, 
                              interval = c(-50,50), maximum = T)$maximum
  est.np.calib.rs <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px.np.rs, 
                                 interval = c(-50,50), maximum = T)$maximum
  # Estimate Var(\hat{beta}) of different estimators 
  var.beta.wb <- ICcalib:::CalcVarThetaWeib(beta = est.weib.calib, etas = calib.weib.params, tm = obs.tm, event = delta, ps = px, 
                                 ps.deriv.shape = px.deriv.shape, ps.deriv.scale =  px.deriv.scale, w = w, w.res = w.res)
  var.beta.wb.rs <- ICcalib:::CalcVarThetaWeibRS(beta = est.weib.calib.rs, etas.matrix =  calib.rs.params, tm = obs.tm, event = delta, ps = px.rs, 
                                      ps.deriv.shape.rs = px.rs.deriv.shape, ps.deriv.scale.rs = px.rs.deriv.scale, w = w, w.res = w.res)
  boot.np <- ICcalib:::CalcVarNpmle(tm = obs.tm, event = delta, w = w , w.res = w.res, BS = BS)
  boot.np.rs <- ICcalib:::CalcVarNpmleRS(tm = obs.tm, event = delta, w = w , w.res = w.res, BS = BS)
  var.beta.np <- boot.np$v
  var.beta.np.rs <- boot.np.rs$v
  boot.ci.np <- boot.np$ci
  boot.ci.np.rs <- boot.np.rs$ci
  # Returns a vector
  return.vec <- c(n.sample, mu, lambda, alpha, weib.shape, weib.scale, beta0, n.points, n.cases, BS, # design parameters
                  calib.weib.params, est.lvcf, est.midI, est.weib.calib, est.weib.calib.rs, est.np.calib, est.np.calib.rs, # parameter estimates 
                  var.lvcf, var.midI, var.beta.wb, var.beta.wb.rs, var.beta.np, var.beta.np.rs, # variance estimates 
                  boot.ci.np, boot.ci.np.rs) #  non parametric BS confidence intervals
  return(return.vec)
}
#CI using estimated variance can be calculated when summarizing data

IterSimSinglePexp <- function(n.sample, mu, lambda, alpha, beta, n.points, pexp.rates, pexp.ts, BS)
{
  my.data <- ICcalib:::SimCoxIntervalCensSinglePexp(n.sample = n.sample, mu = mu, lambda = lambda, alpha = alpha, 
                                                       beta0 = beta0, n.points = n.points, rates = pexp.rates, ts = pexp.ts)
  obs.tm <- my.data$obs.tm
  delta <- my.data$delta
  case.times <- obs.tm[delta==1]
  n.cases <- length(case.times)
  w <- my.data$w
  w.res <- my.data$w.res
  #### calculate P(X(t)=1|history) for each person in the risk set, for all risk sets for carry forward and midpoint #######
  df.lvcf <- ICcalib:::LVCFdata(w = w, w.res = w.res, obs.tm = obs.tm, delta = delta)
  df.midI <- ICcalib:::MidIdata(w = w, w.res = w.res, obs.tm = obs.tm, delta = delta)
  fit.lvcf <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.lvcf)
  fit.midI <- survival::coxph(survival::Surv(start.time, stop.time, delta) ~ X, data = df.midI)
  est.lvcf <- coef(fit.lvcf)
  est.midI <- coef(fit.midI)
  var.lvcf <- diag(fit.lvcf$var)
  var.midI <- diag(fit.midI$var)
  ##### Fit calibration and risk set calibration models and then calculate P(X(t)=1|history)       ############
  ###### for each person in the risk set, for all risk sets                                        ############
  ### Fit Weibull calibration model 
  calib.weib.params <- ICcalib:::FitCalibWeibull(w = w, w.res = w.res)
  px <- t(sapply(case.times,ICcalib:::CalcWeibullCalibP, w = w, w.res =  w.res, weib.params = calib.weib.params))
  # Calculate derivative matrices
  px.deriv.shape <- t(sapply(case.times,ICcalib:::CalcWeibullCalibPderivShape, w = w, w.res =  w.res, weib.params = calib.weib.params))
  px.deriv.scale <- t(sapply(case.times,ICcalib:::CalcWeibullCalibPderivScale, w = w, w.res =  w.res, weib.params = calib.weib.params))
  ### Fit Weibull Risk Set model 
  calib.rs.params <- ICcalib:::FitCalibWeibullRS(w = w, w.res = w.res, tm = obs.tm, event = delta)
  px.rs <- t(sapply(1:n.cases,function(m) {ICcalib:::CalcWeibullRSP(w = w, w.res = w.res, point = case.times[m], weib.params = calib.rs.params[m,])}))
  # Calculate derivative matrices
  px.rs.deriv.shape <- ICcalib:::CalcWeibullCalibPderivShapeRS(w = w, w.res = w.res, obs.tm = obs.tm, event = delta, weib.rs.params = calib.rs.params)
  px.rs.deriv.scale <- ICcalib:::CalcWeibullCalibPderivScaleRS(w = w, w.res = w.res, obs.tm = obs.tm, event = delta, weib.rs.params = calib.rs.params)
  ### Fit Nonparametric calibration model
  fit.npmle <- ICcalib:::FitCalibNpmle(w = w, w.res = w.res)
  px.np <- t(sapply(case.times, ICcalib:::CalcNpmleCalibP, w = w, w.res =  w.res, fit.npmle = fit.npmle))
  ### Fit Nonparametric risk set calibration model
  px.np.rs <- t(sapply(case.times, ICcalib:::CalcNpmleRSP, w = w, w.res =  w.res, obs.tm = obs.tm))
  ###### Calculate beta estimators #######
  est.weib.calib <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px, 
                             interval = c(-50,50), maximum = T)$maximum
  est.weib.calib.rs <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = my.data$delta, ps = px.rs, 
                                interval = c(-50,50), maximum = T)$maximum
  est.np.calib <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = my.data$delta, ps = px.np, 
                           interval = c(-50,50), maximum = T)$maximum
  est.np.calib.rs <- optimize(f = ICcalib:::CoxLogLikX,  tm = obs.tm, event = delta, ps = px.np.rs, 
                              interval = c(-50,50), maximum = T)$maximum
  # Estimate Var(\hat{beta}) of different estimators 
  var.beta.wb <- ICcalib:::CalcVarThetaWeib(beta = est.weib.calib, etas = calib.weib.params, tm = obs.tm, event = delta, ps = px, 
                              ps.deriv.shape = px.deriv.shape, ps.deriv.scale =  px.deriv.scale, w = w, w.res = w.res)
  var.beta.wb.rs <- ICcalib:::CalcVarThetaWeibRS(beta = est.weib.calib.rs, etas.matrix =  calib.rs.params, tm = obs.tm, event = delta, ps = px.rs, 
                                   ps.deriv.shape.rs = px.rs.deriv.shape, ps.deriv.scale.rs = px.rs.deriv.scale, w = w, w.res = w.res)
  boot.np <- ICcalib:::CalcVarNpmle(tm = obs.tm, event = delta, w = w , w.res = w.res, BS = BS)
  boot.np.rs <- ICcalib:::CalcVarNpmleRS(tm = obs.tm, event = delta, w = w , w.res = w.res, BS = BS)
  var.beta.np <- boot.np$v
  var.beta.np.rs <- boot.np.rs$v
  boot.ci.np <- boot.np$ci
  boot.ci.np.rs <- boot.np.rs$ci
  # Returns a vector
  return.vec <- c(n.sample, mu, lambda, alpha, pexp.rates, pexp.ts, beta0, n.points, n.cases, BS, # design parameters
                  calib.weib.params, est.lvcf, est.midI, est.weib.calib, est.weib.calib.rs, est.np.calib, est.np.calib.rs, # parameter estimates
                  var.lvcf, var.midI, var.beta.wb, var.beta.wb.rs, var.beta.np, var.beta.np.rs, # Variance estimates
                  boot.ci.np, boot.ci.np.rs)  #  non parametric BS confidence intervals
  return(return.vec)
}
#CI using estimated variance can be calculated when summarizing data
