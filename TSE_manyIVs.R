#################################
TSE_manyIVs = function(Z,A,Y){
  p<-dim(Z)[2] ## number of Markers
  ## Step 1: regress Y on 1, A, Z, AZ
  fit1 = lm(Y ~ A + Z + A:Z)
  fit1_coef = summary(fit1)$coef[,1]
  ## obtain E[Y|A=0,Z]
  theta_0 = fit1_coef[1]
  theta_z = fit1_coef[3:(2+p)]
  E.hat_Y_A0_Z =  theta_0 + Z%*%theta_z
  ## obtain squared residuals
  res.sq = (fit1$residuals)^2
  
  ### Step 2: estimate the residual variance functional form of Z
  fit2 = glm(res.sq ~ Z,family = Gamma(link="log")) # log transformation
  eta = summary(fit2)$coef[,1]
  eta_0 = eta[1]
  eta_z = eta[2:(1+p)]
  sigma.hat.Z = exp( eta_0 + Z%*% eta_z)
  ### Step 3: refit the model 
  y.new = Y - E.hat_Y_A0_Z
  A.sigma = A*sigma.hat.Z  ## interaction term
  fit3 = lm(y.new ~ A + A.sigma - 1 ) ## no intercept 
  beta = summary(fit3)$coef[1,1]  
  gamma = summary(fit3)$coef[2,1]
  
  par_est = c(beta, gamma, eta, theta_0,theta_z)
  beta_gamma_se = summary(fit3)$coef[1:2,2]
  tmp = list(par_est,beta_gamma_se)
  names(tmp) = c("par_est","beta_gamma_se")
  return(tmp)
}
