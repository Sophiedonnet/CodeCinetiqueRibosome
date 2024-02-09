#------------------------------- Fit log norm mixutre
fitMixtureCtr <- function(Y,Kmax){
  
  loglik <- c()
  allRes[[1]] <- list(mu = mean(log(Y)),sigma = sd(log(Y)),lambda = 1)
  loglik[1] = sum(dlnorm(Y,allRes[[1]]$mu,allRes[[1]]$sigma,log = TRUE))
  for (k in 2:Kmax){
    allRes[[k]] <- mixtools::normalmixEM(log(Y),maxit = 10000,k=k)
    loglik[k] <- allRes[[k]]$loglik
  }
  res <- allRes[[which.max(loglik)]]
  return(res)
}



#-------------------------  density log-Normal mixture on Ctrl 
dCtrData <- function(t,resFitMixture,log = FALSE){
  K <- length(resFitMixture[[k]])
  y <- apply(sapply(1:K,function(k){resFitMixture$lambda[k]*dlnorm(t,resFitMixture$mu[k],resFitMixture$sigma[k])}),1,sum)
  if(log){y <- log(y)}
  return(y)
}

#-------------------------  probability function log-Normal mixture on Ctrl 
pCtrData <- function(t,resFitMixture){
  K <- length(resFitMixture[[k]])
  p <- apply(sapply(1:K,function(k){resFitMixture$lambda[k]*plnorm(t,resFitMixture$mu[k],resFitMixture$sigma[k])}),1,sum)
  return(p)
}
#-------------------------- Probability function of min(Ctr,EMG)
pminCtrEMG <- function(t,resFitMixture,mu,sigma,lambda,Tmax){
  FU <- pemg(t,mu=mu,sigma,lambda) 
  FV <- pCtrData(t,resFitMixture)
  Fmin <- 1- (1-FU)*(1-FV)
  # normalisation sur T max
  FUTmax <- pemg(Tmax,mu=mu,sigma,lambda)
  FVTmax <- pCtrData(Tmax,resFitMixture)
  FminTmax <- 1- (1-FUTmax)*(1-FVTmax)
  return(Fmin/FminTmax)
}
#----------------------- 
dminCtrEMG <- function(t,resFitMixture,mu,sigma,lambda,Tmax,log = FALSE){
  FU <- pemg(t,mu=mu,sigma=sigma,lambda = lambda) 
  FV <- pCtrData(t,resFitMixture)
  fU <- demg(t,mu=mu,sigma,lambda)
  fV <- dCtrData(t,resFitMixture)
  fmin <- fV*(1-FU) + fU*(1-FV)
  
  #Normalisation 
  FUTmax <- pemg(Tmax,mu=mu,sigma=sigma,lambda = lambda)
  FVTmax <- pCtrData(Tmax,resFitMixture)
  FminTmax <- 1- (1-FUTmax)*(1-FVTmax)
  d <- fmin/FminTmax
  if(log){d <- log(fmin) - log(FminTmax) }
  return(d)
}

logLik_mixture <-function(paramUP,paramDN, resFitMixture,data,Tmax){
  
  LL <- 0
  #UP
  if(!is.null(data$Texp_UP)){
    logLUP <- sum(dminCtrEMG(data$Texp_UP,resFitMixture$UP,mu=paramUP[1],sigma = paramUP[2],lambda = paramUP[3],Tmax,log = TRUE))
    LL <- LL + logLUP
  }

  #DN
  if(!is.null(data$Texp_DN)){
    logLDN <- sum(dminCtrEMG(data$Texp_DN,resFitMixture$DN,mu=paramDN[1],sigma = paramDN[2],lambda = paramDN[3],Tmax,log = TRUE))
    LL <- LL + logLDN
  }
  return(LL)  
  
}

