##### Condor ######
#### many possibly invalid IV
source("MR-MiSTERI_function.R")

n_rep<-1000 ## number of simulations
n=1e5
p<-25 ## number of IVs
beta = 0.8
gamma = 0.2
theta0 = -0.5
thetaz = rep(0.1,p) ## pleiotropy  bias
## the third number in eta measures the IV strenth
eta_0 = -0.1
eta_z<-rep(0.1,p)
alpha<-0.05 ## sig level


res<-matrix(nrow=n_rep,ncol=5)
res_se<-matrix(nrow=n_rep,ncol=4)
res_CP<-matrix(0,nrow=n_rep,ncol=4)
colnames(res)<-c("3stage_beta","3stage_gamma","mle_beta","mle_gamma","criteria")
colnames(res_se)<-c("3stage_beta","3stage_gamma","mle_beta","mle_gamma")
colnames(res_CP)<-c("3stage_beta","3stage_gamma","mle_beta","mle_gamma")

for(k in 1:n_rep){
  cat("Iteration:",k,"\n")
  set.seed(k)
  tmp = runSim()
  res[k,1:2]<-tmp$par_est[1:2]
  res_se[k,1:2]<-tmp$beta_gamma_se
  res[k,3:4]<-tmp$MLE_est[1:2]
  res[k,5]<-tmp$criteria
  res_se[k,3:4]<-tmp$MLE_se
  for(i in 1:2){
    #print(c(beta, res[1,(i-1)*2+1],res_se[1,(i-1)*2+1]) )
    res_CP[k,(i-1)*2+1]<-1*(beta<res[k,(i-1)*2+1]+qnorm(1-alpha/2)*res_se[k,(i-1)*2+1] &
                            beta>res[k,(i-1)*2+1]-qnorm(1-alpha/2)*res_se[k,(i-1)*2+1] )
  }
  for(i in 1:2){
    #print(c(gamma, res[1,(i-1)*2+2],res_se[1,(i-1)*2+2]) )
    res_CP[k,(i-1)*2+2]<-1*(gamma<res[k,(i-1)*2+2]+qnorm(1-alpha/2)*res_se[k,(i-1)*2+2] &
                            gamma>res[k,(i-1)*2+2]-qnorm(1-alpha/2)*res_se[k,(i-1)*2+2] )
  }
  print("The results:")
  print(res[k,])
  print(res_se[k,])
  print(res_CP[k,])
}

#write.csv(res,file=paste("res_","p",p,"n",n,"eta_z",eta_z[1],".csv",sep=""),row.names = FALSE)
#write.csv(res_se,file=paste("res_se_","p",p,"n",n,"eta_z",eta_z[1],".csv",sep=""),row.names = FALSE)
#write.csv(res_CP,file=paste("res_CP_","p",p,"n",n,"eta_z",eta_z[1],".csv",sep=""),row.names = FALSE)
n
p
colMeans(res)
apply(res,2,sd)
colMeans(res_se)
colMeans(res_CP)
### point estimates
beta.tse = mean(res[,1])
beta.mle = mean(res[,3])
gamma.tse = mean(res[,2])
gamma.mle = mean(res[,4])
#### beta bias for tse and mle
(c(beta.tse,beta.mle) - beta)/beta*100
#### gamma bias for tse and mle
(c(gamma.tse,gamma.mle) - gamma)/gamma*100
####
### SE and SD
beta.se.tse = colMeans(res_se)[1]
beta.se.mle = colMeans(res_se)[3]
gamma.se.tse = colMeans(res_se)[2]
gamma.se.mle = colMeans(res_se)[4]
colMeans(res_se)
apply(res,2,sd)
### 95% coverage
## beta TSE
beta.tse.lo = res[,1] - 1.96*res_se[,1]
beta.tse.up =  res[,1] + 1.96*res_se[,1]
beta.tse.cover = (beta.tse.lo < beta & beta <  beta.tse.up) ## boolean
mean(beta.tse.cover)
## beta mle
beta.mle.lo = res[,3] - 1.96*res_se[,3]
beta.mle.up =  res[,3] + 1.96*res_se[,3]
beta.mle.cover = (beta.mle.lo < beta & beta <  beta.mle.up) ## boolean
mean(beta.mle.cover)

## gamma TSE
gamma.tse.lo = res[,2] - 1.96*res_se[,2]
gamma.tse.up =  res[,2] + 1.96*res_se[,2]
gamma.tse.cover = (gamma.tse.lo < gamma & gamma <  gamma.tse.up) ## boolean
mean(gamma.tse.cover)
## gamma mle
gamma.mle.lo = res[,4] - 1.96*res_se[,4]
gamma.mle.up =  res[,4] + 1.96*res_se[,4]
gamma.mle.cover = (gamma.mle.lo < gamma & gamma <  gamma.mle.up) ## boolean
mean(gamma.mle.cover)
