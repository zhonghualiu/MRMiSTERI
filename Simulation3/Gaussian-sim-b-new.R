source("Gaussian_mixture_fun-c.R")

#compute mixture weights
wgt <- function(A, s, gamma, pi, mu, se) {
  K=length(pi);
  w=matrix(NA,n,K);
  den=numeric(n);
  for (k in 1:K) {
    w[,k]= exp(gamma*A*s*mu[k]+se[k]^2*(gamma*A*s)^2/2);
  }
  w
}

expit = function(x){
  exp(x)/(1+exp(x))
}

# random generation from mixture normal
rmixnormal <- function(n, pi, mu, se) {
  k <- sample.int(length(pi), n, replace = TRUE, prob = pi)
  rnorm(n, mu[k], se[k])
}

#############################
### Simulation parameters ###
#############################

#Number of replicates
niter= 1000
#sample size per replicate
n= 1e5;
# normal mixture component proportion
pi<- c(0.4, 0.6);
# normal mixture component mean
mu<- c(0.6,-0.4)
# normal mixture component se
se<- c(0.5,sqrt((1-sum(pi*mu^2)-pi[1]*0.5^2)/pi[2]));

#sum(pi*se^2)+sum(pi*mu^2)
#sum(pi*mu)

beta = 0.8
gamma = 0.2
gamma0 = -0.2 ## to make Pr(A=1) reasonable
gammaz = 0.6
#theta = c(1,0.3); #invalid IV
theta = c(1,0); #valid IV

eta.z=c(0.1,0.25,0.5) ## is this part useful???--ZL
#the second number in eta measures the IV strength
#eta = c(-0.1,0.5)
#eta = c(-0.1,0.2)
#eta = c(-0.1,0.5)
eta = c(-0.1,0.4)
par_est = matrix(NA,nrow = niter,ncol = 6)
par_se = matrix(NA,nrow = niter,ncol = 6)

set.seed(1)

for (j in 1:niter) {

    #############################
    ###### Data generation ######
    #############################
    ###Gaussian mixture error####
    #############################
    U = rmixnormal(n, pi, mu, se);
    #var(eps)=1
    Z = rbinom(n,size=2, prob=0.3)
    #A = rnorm(n)
    s.z= exp(0.5*(eta[1] + eta[2]*Z ))

    #mu.az = beta*A  + theta[1] + theta[2]*Z + s.z*as.vector((w%*%(pi*mu) + gamma*A*s.z*w%*%(pi*se^2))/(w%*%pi));
    Y0 =  theta[1] + theta[2]*Z + s.z*U ## Y0 for everyone

    Pr.A = expit(gamma0 + gammaz*Z + gamma*Y0)
    #summary(Pr.A)
    A = rbinom(n,size = 1,p=Pr.A)
    Y = beta*A + Y0
    #w=wgt(A,s.z,gamma,pi,mu,se)
    out = Gaussian_mix(Z,A,Y, maxiter=100,tol=1e-3);

    par_est[j,] = out$est
    par_se[j,]  = out$se
    print(paste("******* ",j," *******",sep=""))
}


mean(par_est[,2]);
(mean(par_est[,2])-beta)/beta*100
(mean(par_se[,2]));
sd(par_est[,2]);
sum(par_est[,2]-1.96*par_se[,2]<beta & par_est[,2]+1.96*par_se[,2]>beta)/niter

mean(par_est[,4]);
(mean(par_est[,4])-gamma)/gamma*100
(mean(par_se[,4]));
sd(par_est[,4]);
sum(par_est[,4]-1.96*par_se[,4]<gamma & par_est[,4]+1.96*par_se[,4]>gamma)/niter


filename=paste("simulation_n",n,"_etaz",eta[2],"New.RData",sep="")
save.image(filename)


###################
