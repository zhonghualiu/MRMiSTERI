## simulation for MR-MiSTERI
#rm(list = ls())
#set.seed(123)
expit = function(x){
  exp(x)/(1+exp(x))
}

TSE = function(Z,A,Y){
  p<-dim(Z)[2]
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

CMLE<-function(Z,A,Y,par_start){
  p<-dim(Z)[2]
  neg.log.like = function(pars){
    beta = pars[1]
    gamma = pars[2]
    eta_0 = pars[3]
    eta_z = pars[4:(3+p)]
    theta_0 = pars[4+p]
    theta_z = pars[(5+p):(4+2*p)]
    
    ## the mean of the normal
    mu_AZ = beta*A + gamma*A*exp(eta_0 + Z%*% eta_z) + theta_0 + Z %*% theta_z
    ## the SD of the normal 
    sigma_AZ = exp(0.5*(eta_0+ Z %*% eta_z ))
    
    ## negative log normal likelihood 
    nll = log(sigma_AZ) + (Y- mu_AZ)^2/(2*sigma_AZ^2)
    return(sum(nll))
  }
  neg.log.like.deriv<-function(pars){
    beta = pars[1]
    gamma = pars[2]
    eta_0 = pars[3]
    eta_z = pars[4:(3+p)]
    theta_0 = pars[4+p]
    theta_z = pars[(5+p):(4+2*p)]
    
    res<-numeric(4+2*p)
    mu_AZ = beta*A + gamma*A*exp(eta_0 + Z%*% eta_z) + theta_0 + Z %*% theta_z
    sigma_AZ = exp(0.5*(eta_0+ Z %*% eta_z ))
    res[1]<-sum((Y-mu_AZ)*A/sigma_AZ^2)
    res[2]<-sum((Y-mu_AZ)*A)
    res[3]<-sum(-1/2+(Y-mu_AZ)*A*gamma+(Y-mu_AZ)^2/2/sigma_AZ^2)
    res[4:(3+p)]<-(-1/2)*colSums(Z)+ gamma*t((Y-mu_AZ)*A) %*% Z + t((Y-mu_AZ)^2/2/sigma_AZ^2) %*% Z
    res[4+p]<-sum((Y-mu_AZ)/sigma_AZ^2)
    res[(5+p):(4+2*p)]<-t((Y-mu_AZ)/sigma_AZ^2) %*% Z
    
    -res
  }
  tmp<-list()
  ntry<-10
  for (try in 1:ntry){
    print(try)
    if(try==1){
      par_start_random<-par_start
    }
    else{
      par_start_random<-rnorm(length(par_start),mean=par_start,sd=(par_start)^2)
    }
    try(tmp[[try]] <- optim(par_start_random,fn = neg.log.like,gr=neg.log.like.deriv,method = "BFGS"),
        silent = TRUE)
  }
  index <- which.min(sapply(1:ntry, function(i) {res <- tmp[[i]]$value; if(is.null(res)){Inf}
  else{ ifelse(tmp[[i]]$convergence!=0, Inf, res)}}))
  MLE_est = tmp[[index]]$par
  MLE_hessian =optimHess(MLE_est,fn = neg.log.like)
  criteria<-min(eigen(MLE_hessian)$values)/length(par_start)
  MLE_se = sqrt(diag(solve(MLE_hessian))[1:2])
  tmp = list(MLE_est,MLE_se,criteria)
  names(tmp)= c("MLE_est","MLE_se","criteria")
  return(tmp)
}

############################
runSim = function(B=100){

Z<-matrix(nrow=n,ncol=p)
for(j in 1:p){
  Z[,j]<-rbinom(n,size = 2,p=0.3)
}
####################
U = rnorm(n)
#mu = beta*A + gamma*A*exp(eta_0 + Z %*% eta_z) + theta0 + Z%*% thetaz
sd = sqrt(exp(eta_0 + Z %*% eta_z))
Y0 = theta0 + Z%*% thetaz + sd*U
### generate A, which depends on U through Y_0
### extended propensity score model
gamma0 = -0.2 ## to make Pr(A=1) reasonable 
gammaz = rep(0.01,p)
Pr.A = expit(gamma0 + Z%*%gammaz + gamma*Y0)
#summary(Pr.A)
A = rbinom(n,size = 1,p=Pr.A)
Y = beta*A + Y0
#################
## run three stage estimation
##################

three_stage_est<-tryCatch(TSE(Z=Z,A=A,Y=Y),error=function(e) NA)
if(!is.na(three_stage_est[1])){
  par_est = three_stage_est$par_est
  boot<-matrix(nrow=B,ncol=2)
  for(b in 1:B){
    boot_ind<-sample(n,size=n,replace = TRUE)
    boot[b,]<-TSE(Z=Z[boot_ind,],A=A[boot_ind],Y=Y[boot_ind])$par_est[1:2]
  }
  beta_gamma_se = apply(boot,2,sd)
  par_start<-par_est
}else{
  par_est<-NA
  beta_gamma_se<-NA
  par_start<-rep(0,4+2*dim(Z)[2])
}


#######################################################
## one-step estimator using Newton-Raphson
## numerical 1st and 2nd derivative for log-likelihood
#######################################################

## use negative log-like, so it is minimization problem.
## assume Z, A, Y are in the environment
CMLE_res<-CMLE(Z,A,Y,par_start)
tmp = list(par_est,beta_gamma_se,CMLE_res$MLE_est,CMLE_res$MLE_se,CMLE_res$criteria)
names(tmp)= c("par_est","beta_gamma_se","MLE_est","MLE_se","criteria")
return(tmp)
}
######
