### MR-MiSTERI with Gaussian mixture error in 2 components ###

#Solver with contraints
library(alabama)

#set.seed(1)

#Gaussian_mix(Z,A,Y)

Gaussian_mix<- function(Z,A,Y, maxiter=100,tol=1e-3) {

    ########## Auxiliary functions ##########
    ####### Gaussian mixture density ########
    dmixnormal <- function(x, pi, mu, se) {
      k <- length(pi)
      n <- length(x)
      rowSums(vapply(1:k, function(i) pi[i] * dnorm(x, mu[i], se[i]), numeric(n)))
    }

    #Zhonghua's functions

    TSE = function(Z,A,Y){
  
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
      fit2 = glm(res.sq ~ Z,family = Gamma(link="log"))
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
      beta_gamma_se = summary(fit3)$coef[1:2,2]
      tmp = list(par_est,beta_gamma_se)
      names(tmp) = c("par_est","beta_gamma_se")
      return(tmp)
    }


    neg.log.like = function(pars){
      beta = pars[1]
      gamma = pars[2]
      eta = pars[3:4]
      theta_0 = pars[5]
      theta_z = pars[6]
  
      ## the mean of the normal
      mu_AZ = beta*A + gamma*A*exp(eta[1] + eta[2]*Z) + theta_0 + theta_z*Z
      ## the SD of the normal 
      sigma_AZ = exp(0.5*(eta[1]+ eta[2]*Z))
    
      ## negative log normal likelihood 
      nll = log(sigma_AZ) + (Y- mu_AZ)^2/(2*sigma_AZ^2)
      return(sum(nll))
    }


    #compute mixture weights
    wgt<- function(A, s, gamma, pi, mu, se) {
      K=length(pi);
      w=matrix(NA,n,K);
      den=numeric(n);
      for (k in 1:K) {
        w[,k]= exp(gamma*A*s*mu[k]+se[k]^2*(gamma*A*s)^2/2);
      }
      w
    }

    #constraints of mixture parameters
    eval_nll_eq<- function(g){
    h<-rep(NA,2);
    h[1]=g[1]*g[2]+(1-g[1])*g[3];
    h[2]=g[1]*g[4]^2+(1-g[1])*g[5]^2+g[1]*g[2]^2+(1-g[1])*g[3]^2-1
    h
    }

    eval_nll_ineq<- function(g){
    c(g[1],g[4],g[5]);
    }

    #negative log-likelihood of Gaussian mixture
    nll.mix<- function(g) {
      #mixture params
      pi=c(g[1],(1-g[1]));
      mu=g[2:3];
      se=g[4:5];
      return(-sum(log(dmixnormal(Y.tilde, pi,mu,se))))
    }

    #Stacked estimating functions

    est.fun <- function(g) {
      h<-rep(0,11);

      s2= exp(g[10]+g[11]*Z)
      w.tilde = wgt(A, sqrt(s2),g[9], c(g[1],(1-g[1])), g[2:3], g[4:5])
      d.tilde =as.vector(w.tilde%*%c(g[1],(1-g[1])))
      n1.tilde =as.vector(w.tilde%*%(c(g[1],(1-g[1]))*g[2:3]))
      n2.tilde =as.vector(w.tilde%*%(c(g[1],(1-g[1]))*g[4:5]^2))
      delta  = (sqrt(s2)*n1.tilde)/d.tilde
      V      = (A*n2.tilde*s2)/d.tilde

      res=Y-g[6]-g[7]*A-g[8]*Z-g[9]*V-delta;
      eps=res/sqrt(s2);
      wt =dmixnormal(eps, c(g[1],(1-g[1])), g[2:3], g[4:5])

      #scores for Gaussian mixture
      h[1]=sum((dnorm(eps,g[2],g[4])-dnorm(eps,g[3],g[5]))/wt);
      h[2]=sum((g[1]*dnorm(eps,g[2],g[4])/wt)*(eps-g[2])/g[4]);
      h[3]=sum(((1-g[1])*dnorm(eps,g[3],g[5])/wt)*(eps-g[3])/g[5]);
      h[4]=sum((g[1]*dnorm(eps,g[2],g[4])/wt)*((eps-g[2])^2-g[4]^2)/g[4]^3);
      h[5]=sum(((1-g[1])*dnorm(eps,g[3],g[5])/wt)*((eps-g[3])^2-g[5]^2)/g[5]^3);

      #regression parameters
      h[6]=sum(res);
      h[7]=sum(A*res);
      h[8]=sum(Z*res);
      h[9]=sum(V*res);

      #variance parameters
      h[10]=sum(res^2-s2);
      h[11]=sum(Z*(res^2-s2));
      h
    }

    est.mm <- function(g) {

      s2= exp(g[10]+g[11]*Z)
      w.tilde = wgt(A, sqrt(s2),g[9], c(g[1],(1-g[1])), g[2:3], g[4:5])
      d.tilde =as.vector(w.tilde%*%c(g[1],(1-g[1])))
      n1.tilde =as.vector(w.tilde%*%(c(g[1],(1-g[1]))*g[2:3]))
      n2.tilde =as.vector(w.tilde%*%(c(g[1],(1-g[1]))*g[4:5]^2))
      delta  = (sqrt(s2)*n1.tilde)/d.tilde
      V      = (A*n2.tilde*s2)/d.tilde

      res=Y-g[6]-g[7]*A-g[8]*Z-g[9]*V-delta;
      eps=res/sqrt(s2);
      wt =dmixnormal(eps, c(g[1],(1-g[1])), g[2:3], g[4:5])

      rbind(
      #scores for Gaussian mixture
      ((dnorm(eps,g[2],g[4])-dnorm(eps,g[3],g[5]))/wt),
      ((g[1]*dnorm(eps,g[2],g[4])/wt)*(eps-g[2])/g[4]),
      (((1-g[1])*dnorm(eps,g[3],g[5])/wt)*(eps-g[3])/g[5]),
      ((g[1]*dnorm(eps,g[2],g[4])/wt)*((eps-g[2])^2-g[4]^2)/g[4]^3),
      (((1-g[1])*dnorm(eps,g[3],g[5])/wt)*((eps-g[3])^2-g[5]^2)/g[5]^3),

      #regression parameters
      (res),
      (A*res),
      (Z*res),
      (V*res),

      #variance parameters
      (res^2-s2),
      (Z*(res^2-s2)))
    }



    ########## Alternating optimization procedure ###########
    ## Step 1: Initialize assuming standard Gaussian error ##

    three_stage_est = TSE(Z=Z,A=A,Y=Y)
    par.est = three_stage_est$par_est

    MLE_opt = optim(par =par.est,fn = neg.log.like,hessian = TRUE)
    MLE_est = MLE_opt$par

    s.z.tilde= exp(0.5*(MLE_est[3] + MLE_est[4]*Z));
    gamma.tilde= MLE_est[2];
    Y.tilde=(Y-(MLE_est[1]*A + gamma.tilde*A*s.z.tilde^2 + MLE_est[5] + MLE_est[6]*Z))/exp(0.5*(MLE_est[3]+MLE_est[4]*Z));

    ## Step 2: Alternate between maximum likelihood and regression ##
    loglike=numeric(maxiter);
    loglike[1]<-0;
    loglike[2]<-1;
    k=2;

    while (abs(loglike[k]-loglike[k-1])/abs(loglike[k-1]) >= tol) {
      k=k+1;
 
      #A set of initial values that fit the contraints
      out<-constrOptim.nl(par=c(0.5,-0.5,0.5,0.5,1.12),
      fn=nll.mix,hin=eval_nll_ineq,heq=eval_nll_eq,control.outer=list(trace=FALSE))

      loglike[k]=-nll.mix(out$par)
      pi.tilde=c(out$par[1],(1-out$par[1])); mu.tilde= out$par[2:3]; se.tilde= out$par[4:5];
      w.tilde = wgt(A, s.z.tilde,gamma.tilde,pi.tilde,mu.tilde,se.tilde)
      d.tilde =as.vector(w.tilde%*%pi.tilde)
      n1.tilde =as.vector(w.tilde%*%(pi.tilde*mu.tilde))
      n2.tilde =as.vector(w.tilde%*%(pi.tilde*se.tilde^2))

      delta  = (s.z.tilde*n1.tilde)/d.tilde
      V      = (A*n2.tilde*s.z.tilde^2)/d.tilde
      reg    = lm(Y~A+Z+V, offset=delta)
      res.sq = (reg$residuals)^2
      eta.fit= glm(res.sq ~ Z,family = Gamma(link="log"))
      eta.est= summary(eta.fit)$coef[,1]
      reg.est = summary(reg)$coef[,1]
      names(eta.est)<-c("eta0","eta1")
      names(reg.est)<-c("Intercept","beta","thetaz","gamma")
      s.z.tilde= exp(0.5*(eta.est[1] + eta.est[2]*Z));
      gamma.tilde= reg.est[4];
      Y.tilde=(Y-reg$fit)/s.z.tilde;
      print(paste("Iteration: ", (k-2),sep=""))
    }

    if ((k-2)<maxiter) {
      print("Convergence reached")
    } else {
      print("Maximum iteration reached")
    }
    par.est = c(reg.est,eta.est);
    ## compute se ##
    dM     <- jacobian(func=est.fun,x=c(out$par,reg.est,eta.est))/n
    m      <- est.mm(c(out$par,reg.est,eta.est))
    var.est<- diag(solve(dM)%*%(m%*%t(m)/n)%*%t(solve(dM))/n)
    par.se<- sqrt(var.est)[6:11]
    names(par.se)<- c("Intercept","beta","thetaz","gamma","eta0","eta1");
    return(list(est=par.est,se=par.se))
}










