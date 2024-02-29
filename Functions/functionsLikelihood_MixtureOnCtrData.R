#------------------------------- Fit log norm mixutre
fitMixtureCtr <- function(Y, Kmax, Kmin = 1,select = TRUE ){
  
  nbK <- Kmax - Kmin + 1 
  allRes <- vector("list", nbK)
  loglik <- rep(-Inf,nbK)
  pk <-  0
  for (k in Kmin:Kmax){
    pk <- pk  + 1

    skip_to_next <- FALSE
  
    if(k==1){
      allRes[[1]] <- list(mu = mean(log(Y)),sigma = sd(log(Y)),lambda = 1)
      loglik[1] = sum(dlnorm(Y,allRes[[1]]$mu,allRes[[1]]$sigma,log = TRUE))
    }
    else{
      tryCatch({
        init <- mixtools::normalmixEM(Y,mu = kmeans(log(Y),k)$centers, maxit = 10000,k=k,maxrestarts = 50)
        allRes[[pk]] <- mixtools::normalmixEM(log(Y),mu = log(init$mu), maxit = 20000,k=k,maxrestarts = 20,verb = FALSE)
        loglik[pk] <- allRes[[pk]]$loglik
        },
        error = function(e){skip_to_next <<-TRUE}
      )
      if(skip_to_next){ next } 
    }  
    
  }
  Crit <- loglik - 1/2*(3*(Kmin:Kmax)-1)*length(Y) 
  
  ##########################################
  if (select){
    res <- allRes[[which.max(loglik)]]
  }else{
    res <- list() 
    res$fit <- allRes
    res$Crit = Crit
    res$logik = loglik
  }
  return(res)
}

#------------------------------- Density log norm mixtre
dlnormMixture <- function(x, param,log = FALSE ){
  
  nbClass = length(param$mu)
  res <- rowSums(sapply(1:nbClass,function(k){dlnorm(x,param$mu[k],param$sigma[k])*param$lambda[k]}))
  if(log){res <- log(res)}
  return(res)
}

#------------------------------- Fit log norm mixutre
plnormMixture <- function(x, param){
  nbClass = length(param$mu)
  res <- rowSums(sapply(1:nbClass,function(k){plnorm(x,param$mu[k],param$sigma[k])*param$lambda[k]}))
  return(res)
}


#-------------------------  density log-Normal mixture on Ctrl 
dCtrData <- dlnormMixture

#-------------------------  probability function log-Normal mixture on Ctrl 
pCtrData <- plnormMixture
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

#----------------------- 
loglik_allData_mixture <-function(paramUP,paramDN, resFitMixture,data,Tmax){
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

