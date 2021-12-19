misterigauss = function(Z=Z,A=A,Y=Y){
  ## Step 1: regress Y on 1, A, Z, AZ
  fit1 = lm(Y ~ A + Z + A:Z)
  fit1_coef = summary(fit1)$coef[,1]
  ## obtain E[Y|A=0,Z]
  theta_0 = fit1_coef[1]
  theta_z = fit1_coef[3]
  E.hat_Y_A0_Z =  theta_0 + theta_z*Z
  ## obtain squared residuals
  res.sq = (fit1$residuals)^2
  ### Step 2: estimate the residual variance functional form of Z
  fit2 = glm(res.sq ~  Z,family = Gamma(link="log"))
  eta = summary(fit2)$coef[,1]
  eta_0 = eta[1]
  eta_z = eta[2]
  sigma.hat.Z = exp( eta_0 + eta_z*Z)
  ### Step 3: refit the model
  y.new = Y - E.hat_Y_A0_Z
  A.sigma = A*sigma.hat.Z  ## interaction term
  fit3 = lm(y.new ~ A + A.sigma - 1 ) ## no intercept
  beta = summary(fit3)$coef[1,1]
  gamma = summary(fit3)$coef[2,1]
  par_est = c(beta, gamma, eta, theta_0,theta_z)
  
  neg.log.like = function(pars){
    beta = pars[1]
    gamma = pars[2]
    eta = pars[3:4]
    theta_0 = pars[5]
    theta_z = pars[6]
    ## the mean of the normal
    mu_AZ = beta*A + gamma*A*exp(eta[1] + eta[2]*Z) + theta_0 + theta_z*Z
    ## the SD of the normal
    sigma_AZ = exp(0.5*(eta[1]+eta[2]*Z))
    ## negative log normal likelihood
    nll = log(sigma_AZ) + (Y- mu_AZ)^2/(2*sigma_AZ^2)
    return(sum(nll))
  }
  
  MLE_opt = optim(par =par_est,fn = neg.log.like,hessian = TRUE)
  MLE_res = MLE_opt$par[c(1,2,4)] 
  MLE_se = sqrt(diag(solve(MLE_opt$hessian)))[c(1,2,4)]
  names(MLE_res) = c("beta","gamma","eta_z") ## eta_z is the IV strength
  names(MLE_se) = c("SE_beta","SE_gamma","SE_eta_z") ## eta_z is the IV strength
  tmp = list(MLE_res,MLE_se)
  names(tmp)= c("Estimates","SE")
  #tmp = c(MLE_res[1:2],MLE_se)
  return(tmp)
}

