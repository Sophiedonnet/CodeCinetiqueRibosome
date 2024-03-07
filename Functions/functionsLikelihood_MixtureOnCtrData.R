dlnormTrunc <- function(x,mu,sigma,Tmax,log = FALSE){
  C <- plnorm(Tmax,mu,sigma)
  y <- (x<=Tmax)*dlnorm(x,mu,sigma)/C
  if(log){y <- log(y)}
  return(y)
}
#-------------------------------
dexpTrunc <- function(x,lambda,Tmax,log = FALSE){
  C <- pexp(Tmax,lambda)
  y <- (x<=Tmax)*dexp(x,lambda)/C
  if(log){y <- log(y)}
  return(y)
}
#-------------------------------
pexpTrunc <- function(x,lambda,Tmax){
  C <- pexp(Tmax,lambda)
  y <- (x<=Tmax)*pexp(x,lambda)/C
  return(y)
}

fitExpoTruncCtr  <- function(Tctr,Tmax){
  Y <- Tctr[(Tctr < Tmax) & (Tctr > 0)]
  myF <- function(lambda){
  U<- -sum(dexpTrunc(x=Y,lambda,Tmax=Tmax,log = TRUE))
  return(U)}
  res_optim <- optimize(f=myF, 1/mean(Y)*c(0.1 ,10))
  res_optim$objective <- - res_optim$objective
  names(res_optim) <- c('lambda_hat','loglik')
  return(res_optim)
} 
#------------------------------- Fit log norm mixutre
fitMixtureCtr <- function(Tctr,Tmax, Kmax, Kmin = 1,select = TRUE ){
  
  nbK <- Kmax - Kmin + 1 
  Tctr_trunc <- Tctr[(Tctr<Tmax)&(Tctr>0)]
  allRes <- vector("list", nbK)
  loglik <- rep(-Inf,nbK)
  pk <-  0
  for (k in Kmin:Kmax){
    pk <- pk  + 1

    skip_to_next <- FALSE
  
    if(k==1){
      
      myloglik <- function(param,Y,Tmax){
        mu <- param[1]
        sigma <- exp(param[2])
        return(-sum(dlnormTrunc(Y,mu,sigma,Tmax,log = TRUE)))
      } 
      res_optim <- optim(par = c(mean(log(Tctr_trunc)),log(sd(log(Tctr_trunc)))) ,fn=myloglik,Y = Tctr_trunc,Tmax = Tmax)
      allRes[[1]]<- list(mu=res_optim$par[1],sigma=exp(res_optim$par[2]),lambda=1)
      allRes[[1]]$loglik = sum(dlnormTrunc(Tctr_trunc,allRes[[1]]$mu,allRes[[1]]$sigma,Tmax = Tmax,log = TRUE))
      allRes[[1]]$loglik -> loglik[1]
    }
    else{
      tryCatch({
        init <- mixtools::normalmixEM(Tctr_trunc, maxit = 10000,k=k,maxrestarts = 50)
        allRes[[pk]] <- mixtools::normalmixEM(log(Tctr_trunc),mu = log(init$mu), maxit = 20000,k=k,maxrestarts = 20,verb = FALSE)
        loglik[pk] <- allRes[[pk]]$loglik
        },
        error = function(e){skip_to_next <<-TRUE}
      )
      if(skip_to_next){ next } 
    }  
    
  }
  Crit <- loglik - 1/2*(3*(Kmin:Kmax)-1)*length(Tctr) 
  
  ##########################################
  if (select){
    res <- allRes[[which.max(loglik)]]
  }else{
    res <- list() 
    res$fit <- allRes
    res$Crit = Crit
    res$loglik = loglik
  }
  return(res)
}

#------------------------------- Density log norm mixtre
dlnormMixture <- function(x, param,log = FALSE ){
  
  nbClass = length(param$mu)
  re = 0
  if(length(x)==1){re = 1;x = c(x,x)}
  res <- rowSums(sapply(1:nbClass,function(k){dlnorm(x,param$mu[k],param$sigma[k])*param$lambda[k]}))
  if(re == 1){res = res[1]}
  if(log){res <- log(res)}
  return(res)
}

#------------------------------- Fit log norm mixutre
plnormMixture <- function(x, param){
  nbClass = length(param$mu)

  re = 0
  if(length(x)==1){re = 1; x = c(x,x)}
  res <- rowSums(sapply(1:nbClass,function(k){plnorm(x,param$mu[k],param$sigma[k])*param$lambda[k]}))
  if(re == 1){res = res[1]}
  return(res)
}

#------------------------------- simulation log norm mixutre
rlnormMixture <- function(n, param){
  K = length(param$mu)
  Z <- sample(1:K,n,replace = TRUE,prob = param$lambda)
  muZ <- param$mu[Z]
  sigmaZ<- param$sigma[Z]
  return(rlnorm(n,muZ,sigmaZ))
}
#-------------------------  density log-Normal mixture on Ctrl 
pCtrData  <- function(t, paramCtr,Tmax){
  
  y <- plnormMixture(t,paramCtr)
  if(Tmax<Inf){C <- plnormMixture(Tmax,paramCtr)}else{C = 1}
  res <- y/C * (t<Tmax) + 1*(t>=Tmax)
  return(res)
}

#-------------------------  density log-Normal mixture on Ctrl 
rCtrData  <- function(t, paramCtr,Tmax){
  
  y <- rlnormMixture(n,paramCtr)
  w <- which(y > Tmax)
  for(i in w){
    cond = FALSE
    while(!cond){
      yi =  rlnormMixture(1,paramCtr)
      cond <- (yi<Tmax)
      }
    y[i] = yi
  }
  return(y)
}


#-------------------------  probability function log-Normal mixture on Ctrl 
dCtrData <- function(t, paramCtr,Tmax,log = FALSE){
  
  y <- dlnormMixture(t,paramCtr)
  if(Tmax<Inf){C <- plnormMixture(Tmax,paramCtr)}else{C = 1}
  res <- y/C * (t<Tmax) 
  if(log){res <- log(res)}
  return(res)
}

#-------------------------- Probability function of min(Ctr,EMG)
pminCtrEMG <- function(t,paramExp,paramCtr,Tmax){
  mu <- paramExp[1]
  sigma <- paramExp[2]
  lambda <-  paramExp[3]
  
  FU <- pemg(t,mu=mu,sigma,lambda) 
  FV <- pCtrData(t,paramCtr,Tmax = Inf)
  Fmin <- 1- (1-FU)*(1-FV)
  # normalisation sur [0,Tmax]
  FUTmax <- pemg(Tmax,mu=mu,sigma=sigma,lambda = lambda)
  FVTmax <- pCtrData(Tmax,paramCtr,Tmax = Inf)
  FminTmax <- 1- (1-FUTmax)*(1-FVTmax)
  FUT0 <- pemg(0,mu=mu,sigma=sigma,lambda = lambda)
  FVT0 <- pCtrData(0,paramCtr,Tmax = Inf)
  Fmin0 <- 1- (1-FUT0)*(1-FVT0)
  return(Fmin/(FminTmax-Fmin0))
}
#----------------------- 
dminCtrEMG <- function(t,paramExp,paramCtr,Tmax,log = FALSE){
  
  mu <- paramExp[1]
  sigma <- paramExp[2]
  lambda <- paramExp[3]
  FU <- pemg(t,mu=mu,sigma=sigma,lambda = lambda) 
  fU <- demg(t,mu=mu,sigma,lambda) 
  FV <- pCtrData(t,paramCtr,Tmax = Inf) 
  fV <- dCtrData(t,paramCtr,Tmax = Inf)
  fmin <- fV*(1-FU) + fU*(1-FV)
  
  #Normalisation 
  FUTmax <- pemg(Tmax,mu=mu,sigma=sigma,lambda = lambda)
  FVTmax <- pCtrData(Tmax,paramCtr,Tmax = Inf)
  FminTmax <- 1- (1-FUTmax)*(1-FVTmax)

  FUT0 <- pemg(0,mu=mu,sigma=sigma,lambda = lambda)
  FVT0 <- pCtrData(0,paramCtr,Tmax = Inf)
  Fmin0 <- 1- (1-FUT0)*(1-FVT0)
  d <- fmin/(FminTmax-Fmin0)
  if(log){d <- log(fmin) - log(FminTmax-Fmin0) }
  return(d)
}

#-------------------------- Probability function of min(Ctr,EMG)
rminCtrEMG <- function(n,paramExp,paramCtr,Tmax){
  yread <- remg(n,mu=mu,sigma,lambda) 
  y <- rlnormMixture(n,param = paramCtr) 
  w.min <- which(yread<y)
  y[w.min] <- yread[w.min]
  
  w <- which( (y<0) | (y>Tmax))
  for (i in w){
    cond = FALSE
    while (!cond){
      yi <- min(remg(1,mu=mu,sigma,lambda),rlnormMixture(1,param = paramCtr))
      cond <-  ((yi>=0) | (yi<=Tmax))
    }
    y[i] = yi
  }
  return(Z)
  
  
  
  return(Fmin/(FminTmax-Fmin0))
}
#-----------------

#----------------------- 
loglik_Texp_mixture <-function(logparamExp, paramCtr,Texp,Tmax){
  
  paramExp <- exp(logparamExp)
  LL <- -sum(dminCtrEMG(Texp[Texp<Tmax],paramExp,paramCtr,Tmax=Tmax,log = TRUE))
  return(LL)
}
#----------------------- 
fitRepEmp_Texp_mixture <-function(logparamExp, paramCtr,Texp,Tmax){
  
  paramExp <- exp(logparamExp)
  t <- sort(Texp[Texp<Tmax])
  FN_exp <- ecdf(t)
  yEmp <- FN_exp(t)
  yMod <- pminCtrEMG(t,paramExp,paramCtr,Tmax=Tmax)
  LL <- sum((yEmp-yMod)^2)
  return(LL)
}

#---------------------------- 
estimProc_Up_or_Dn <- function(Texp,Tmax, paramCtr){
  
  #---------- moment
  resMoment <- estim_param_moment(Texp, Tmax,paramCtr)
  paramExp_Moment <- resMoment$paramExpMoment
  
  #---------- Fit repartition function
  resFitRepEmp <- optim(par = log(paramExp_Moment), fn = fitRepEmp_Texp_mixture, paramCtr= paramCtr, Texp = Texp,Tmax = Tmax)
  paramExp_FitRepEmp <-exp(resFitRepEmp$par)

  #---------- Fit max likelihood 
  resMaxLoglik <- optim(par = log(paramExp_Moment), fn = loglik_Texp_mixture, paramCtr= paramCtr, Texp = Texp,Tmax = Tmax)
  estim = list(moment = paramExp_Moment, loglik =  exp(resMaxLoglik$par), fitFn  = exp(resFitRepEmp$par))
  return(estim)
}

