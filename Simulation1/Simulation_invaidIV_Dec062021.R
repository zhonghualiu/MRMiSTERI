##############################
## Simulation for MR-MiSTERI
## Dec. 2021
## Methods: TSLS, MR GENIUS, MR MiSTERI
#############################
rm(list = ls()) ## remove everything
library(ivreg)
library(genius)
source("./misterigauss.R")
###########################
expit = function(x){
  exp(x)/(1+exp(x))
}
##############################################
##############################################
#beta = 0.8 #causal effect of interest
#gamma = 0.2 #selection bias
#thetaz = 0.3 #if thetaz=0, then valid. 
## the third number in eta measures the IV strength
#eta_0 = -0.1
##############################################
##############################################
runSim = function(n=1e5,eta_z=0.05,thetaz=0.3,gamma=0.2,gammaz=0.6,beta=0.8){
  theta0 = 1
  eta_0 = -0.1
  gamma0 = -0.2 ## to make Pr(A=1) reasonable 
  ##########################
  ### simulate U, unmeasured each individual's overall health index
  U = rnorm(n)
  ### simulate Z, MAF = 0.3
  Z = rbinom(n,size = 2, p=0.3)
  ### simulate Y0, the potential outcome, viewed as each individual's ultimate confounder 
  sd.z = sqrt(exp(eta_0 + + eta_z*Z))
  Y0 = theta0 + thetaz*Z + sd.z*U ## Y0 for everyone
  ### generate A, which depends on U through Y_0
  ### extended propensity score model
  Pr.A = expit(gamma0 + gammaz*Z + gamma*Y0)
  #summary(Pr.A)
  A = rbinom(n,size = 1,p=Pr.A)
  #cat(mean(A))
  Y = beta*A + Y0
  ##################
  ## 2SLS estimate
  ##################
  fit_2s = ivreg(Y~A|Z)
  beta_hat_2sls = summary(fit_2s)$coefficients[2,1]
  beta_se_2sls = summary(fit_2s)$coefficients[2,2]
  beta_2sls_lo = beta_hat_2sls - 1.96*beta_se_2sls
  beta_2sls_up = beta_hat_2sls + 1.96*beta_se_2sls
  cover_2sls = (beta_2sls_lo < beta & beta <  beta_2sls_up)
  res_2sls = c(beta_hat_2sls,beta_se_2sls, cover_2sls)
  ######################
  ## MR GENIUS method
  ######################
  fit_genius = genius_addY(Y=Y,A=A,G=Z)
  beta_genius = fit_genius$beta.est
  beta_genius_se = sqrt(fit_genius$beta.var)
  beta_genius_ci = fit_genius$ci 
  beta_genius_is_cover = (beta_genius_ci[1] < beta & beta <  beta_genius_ci[2])
  ##########################
  res_genius = c(beta_genius,beta_genius_se,beta_genius_is_cover)
  names(res_genius) = c("beta_genius","beta_genius_se","covered")
  ################################
  ## MR-MISTERI
  ################################
  est_misteri = misterigauss(Z=Z,A=A,Y=Y)
  beta_mi = est_misteri$Estimates[1]
  beta_mi_se = est_misteri$SE[1]
  gamma_mi = est_misteri$Estimates[2]
  gamma_mi_se = est_misteri$SE[2]
  
  ### 95% CI coverage check 
  beta_mi_lo = beta_mi - 1.96*beta_mi_se
  beta_mi_up = beta_mi + 1.96*beta_mi_se
  beta_mi_is_cover = (beta_mi_lo < beta & beta <  beta_mi_up) ## boolean
  
  gamma_mi_lo = gamma_mi - 1.96*gamma_mi_se
  gamma_mi_up = gamma_mi + 1.96*gamma_mi_se
  gamma_mi_is_cover = (gamma_mi_lo < gamma & gamma < gamma_mi_up) ## boolean
  ###############################
  res_misteri = c(beta_mi,beta_mi_se,beta_mi_is_cover,gamma_mi,gamma_mi_se,gamma_mi_is_cover)
  names(res_misteri) = c("beta","beta_se","beta_cover","gamma","gamma_se","gamma_cover")
  res = list(res_2sls,res_genius,res_misteri)
  names(res) = c("res_2sls","res_genius","res_misteri")
  return(res)
}

###########################
nsim = 2000
n.vec = c(1e4,1e4,1e4,3e4,3e4,1e5) 
eta_z.vec = c(0.2,0.15,0.1,0.1,0.05,0.05)
#### output results
mis_res = matrix(NA,ncol = 11,nrow = 6)
genius_res = matrix(NA,ncol = 5,nrow = 6)
tsls_res = matrix(NA,ncol = 5,nrow = 6)
###########################

for(k in 1:6){
  n = n.vec[k]
  eta_z = eta_z.vec[k]
  cat("n:",n,"eta_z:",eta_z,"\n")
  ###########################
  ########################
  res_2sls_mat = matrix(NA,ncol = 3,nrow = nsim)
  res_genius_mat = matrix(NA,ncol = 3,nrow = nsim)
  res_misteri_mat = matrix(NA,ncol = 6,nrow = nsim)
  ###########################
  
  for(j in 1:nsim){
    # if(j%%100==0){
    #   cat(j," ")
    # }
    res = runSim(n=n,eta_z=eta_z,thetaz=0.3,gamma=0.2,gammaz=0.6,beta=0.8)
    res_2sls_mat[j,] = res$res_2sls
    res_genius_mat[j,] = res$res_genius
    res_misteri_mat[j,] = res$res_misteri
  }
  
  ################################
  ## output simulation results 
  ################################
  ## for 2SLS
  ################################
  cat("2SLS results--","beta, SE, coverage:","\n")
  tsls_res_tmp = colMeans(res_2sls_mat)
  tsls_bias = (tsls_res_tmp[1] - beta)/beta*100
  tsls_sd = sd(res_2sls_mat[,1])
  tsls_res[k,] =c(tsls_res_tmp[1],tsls_bias,tsls_res_tmp[2],tsls_sd,tsls_res_tmp[3])
  print(tsls_res[k,])
  ###############################
  ## for GENIUS
  ## est, SE, and coverage
  ################################
  cat("GENIUS results--","beta, SE, coverage:","\n")
  
  genius_res_tmp = colMeans(res_genius_mat)
  genius_bias  = (genius_res_tmp[1] - beta)/beta*100
  genius_sd = sd(res_genius_mat[,1])
  
  genius_res[k,] = c(genius_res_tmp[1],genius_bias,genius_res_tmp[2],genius_sd,genius_res_tmp[3])
  print(genius_res[k,])
  
  #######################
  ## for MISTERI
  ### good to know the relationship between eta the IV strength and the SE relationship
  ## output simulation results 
  ################################
  MLE_est = res_misteri_mat
  beta_hat = mean(MLE_est[,1])
  gamma_hat = mean(MLE_est[,4])
  beta_bias = (beta_hat - beta)/beta*100 ## percentage
  gamma_bias = (gamma_hat - gamma)/gamma*100 ## percentage
  beta_se = mean(MLE_est[,2])
  gamma_se = mean(MLE_est[,5])
  beta_sd = sd(MLE_est[,1])
  gamma_sd = sd(MLE_est[,4])
  beta_coverage = mean(MLE_est[,3])*100 ## percentage
  gamma_coverage = mean(MLE_est[,6])*100 ## percentage
  ################
  cat("MISTERI results:","n:",n,"eta_z:",eta_z,"\n")
  mis_res[k,] = c(eta_z,beta_hat,beta_bias,beta_se,beta_sd,beta_coverage,gamma_hat,gamma_bias,gamma_se,gamma_sd,gamma_coverage)
  print(mis_res[k,])
}
#########################
### results for 2sls 
#######################
tsls_res
### results for genius 
genius_res
### results for misteri
mis_res
n.vec
#######################
#write.csv(cbind(n.vec,mis_res),file = "Table1_invalidIV2000_misteri.csv")
#write.csv(cbind(n.vec,tsls_res),file = "Table1_invalidIV2000_tsls.csv")
#write.csv(cbind(n.vec,genius_res),file = "Table1_invalidIV2000_genius.csv")

#### debug
boxplot(MLE_est[,c(2,5)])
colMeans(MLE_est[,c(2,5)])
apply(MLE_est[,c(2,5)],2,summary)
boxplot(MLE_est[,c(1,4)])
