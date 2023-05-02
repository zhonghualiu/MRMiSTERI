#' MR MiSTERI with many weak invalid IVs
#'
#' misterimawii combines many weak invalid IVs to reduce weak IV bias.
#' @param Z  an IV matrix with columns representing IVs
#' @param A  the exposure variable
#' @param Y  the continuous outcome variable
#' @return a list object that contains causal effect estimates and standard errors.
#' @references https://www.medrxiv.org/content/10.1101/2020.09.29.20204420v3
#' @import alabama
#' @export
misterimawii<-function(Z,A,Y){

  three.stage <-TSE_manyIVs(Z=Z,A=A,Y=Y)
  par_start = three.stage$par_est

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
    try(tmp[[try]] <- optim(par_start_random,fn = neg.log.like,method = "BFGS"),
        silent = TRUE)
  }
  index <- which.min(sapply(1:ntry, function(i) {res <- tmp[[i]]$value; if(is.null(res)){Inf}
  else{ ifelse(tmp[[i]]$convergence!=0, Inf, res)}}))
  MLE_est = tmp[[index]]$par
  MLE_hessian =optimHess(MLE_est,fn = neg.log.like)
  criteria<-min(eigen(MLE_hessian)$values)/length(par_start)
  MLE_se = sqrt(diag(solve(MLE_hessian))[1:2])
  Zvalue = MLE_est/MLE_se
  pvalue = 2*(1-pnorm(abs(Zvalue))) ## two-sided test pvalue
  CI.lower = MLE_est - 1.96*MLE_se
  CI.upper = MLE_est + 1.96*MLE_se
  tmp = list(MLE_est,MLE_se,pvalue,CI.lower, CI.upper, criteria)
  names(tmp)= c("effect_estimate","effect_SE","pvalue","95%CI.lower", "95%CI.upper","criteria")
  return(tmp)
}
