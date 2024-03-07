dlnormTrunc <- function(x,mu,sigma,Tmax,log = FALSE){
  C <- plnorm(Tmax,mu,sigma)
  y <- (x<=Tmax)*dlnorm(x,mu,sigma)/C
  if(log){y <- log(y)}
  return(y)
}
#-------------------------------

plnormTrunc <- function(x,mu,sigma,Tmax){
  C <- plnorm(Tmax,mu,sigma)
  y <- (x<=Tmax)*plnorm(x,mu,sigma)/C
  if(length(x>=Tmax)>0){
    y[x>=Tmax] = 1}
  
  return(y)
}

#-------------------------------

rlnormTrunc <- function(n,mu,sigma,Tmax){
  y <- rlnorm(n,mu,sigma)
  w <- which(y > Tmax)
  for(i in w){
    cond = FALSE
    while(!cond){
      yi =  rlnorm(1,mu,sigma)
      cond <- (yi<Tmax)
    }
    y[i] = yi
  }
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
#-------------------------------
fitExpoTruncCtr  <- function(Tctr,Tmax){
  Y <- Tctr[(Tctr < Tmax) & (Tctr > 0)]
  myF <- function(lambda,Y,Tmax){
    U <- -sum(dexpTrunc(Y,lambda,Tmax=Tmax,log = TRUE))
    return(U)
  }
  lambda_init <- 1/mean(Y)*c(0.1 ,10)
  res_optim <- optimize(f=myF, lambda_init,Y = Y,Tmax = Tmax)
  res_optim$objective <- - res_optim$objective
  names(res_optim) <- c('lambda_hat','loglik')
  return(res_optim)
} 
#------------------------------- Fit log norm mixutre
fitLogNormTruncCtr <- function(Tctr,Tmax){
  
  Y <- Tctr[(Tctr<Tmax)&(Tctr>0)]
  myF <- function(param,Y,Tmax){
    U  <- -sum(dlnormTrunc(Y,mu = param[1],sigma = exp(param[2]),Tmax,log = TRUE))
    return(U)
  } 
  param_init <- c(mean(log(Y)),log(sd(log(Y))))
  res_optim <- optim(par =  param_init ,fn=myF,Y = Y,Tmax = Tmax)
  paramCtr <- res_optim$par
  paramCtr[2] <- exp(paramCtr[2])
  paramCtr[3] <- mean(Tctr==Tmax)
  res <- list(paramCtr = paramCtr, loglik  = -res_optim$value)
  return(res)
}
  

#-------------------------  density log-Normal mixture on Ctrl 
pCtrData  <- function(t, paramCtr,Tmax){
  res <- (1-paramCtr[3])*plnormTrunc(t,paramCtr[1],paramCtr[2], Tmax = Tmax)
  res[t>=Tmax]=1
  return(res)
 }

#-------------------------  density log-Normal mixture on Ctrl 
dCtrData  <- function(t, paramCtr,Tmax,log = FALSE){
  res <- (1-paramCtr[3])*dlnormTrunc(t,paramCtr[1],paramCtr[2], Tmax = Tmax,log = log)
  res[t==Tmax] <- paramCtr[3]
  return(res)
}

#-------------------------  density log-Normal mixture on Ctrl 
rCtrData  <- function(n, paramCtr,Tmax){
  Z <- rbinom(n,paramCtr[3])
  res <- rlnormTrunc(n,paramCtr[1],paramCtr[2], Tmax = Tmax,log = log)
  res[Z==1] = Tmax
  return(res)
}


#-------------------------- Probability function of min(Ctr,EMG)
pExpData <- function(t,paramRead,paramCtr,Tmax){
  
  FU <- pemg(t,mu=paramRead[1],paramRead[2],paramRead[3]) 
  FV <- plnorm(t,paramCtr[1],paramCtr[2])
  Fmin <- 1- (1-FU)*(1-FV)
  # normalisation sur [0,Tmax]
  FUTmax <- pemg(Tmax,mu=paramRead[1],sigma=paramRead[2],lambda = paramRead[3])
  FVTmax <-  plnorm(Tmax,paramCtr[1],paramCtr[2])
  FminTmax <- 1- (1-FUTmax)*(1-FVTmax)
  FUT0 <- pemg(0,mu=paramRead[1],sigma=paramRead[2],lambda = paramRead[3])
  FVT0 <- plnorm(0,paramCtr[1],paramCtr[2])
  Fmin0 <- 1- (1-FUT0)*(1-FVT0)
  return((Fmin-Fmin0)/(FminTmax-Fmin0))
}
#----------------------- 
dExpData <- function(t,paramRead,paramCtr,Tmax,log = FALSE){
  
  FU <- pemg(t,mu=paramRead[1],paramRead[2],paramRead[3]) 
  fU <- demg(t,mu=paramRead[1],paramRead[2],paramRead[3])
  FV <- plnorm(t,paramCtr[1],paramCtr[2])
  fV <- dlnorm(t,paramCtr[1],paramCtr[2])
  fmin <- fV*(1-FU) + fU*(1-FV)
  
  # normalisation sur [0,Tmax]
  FUTmax <- pemg(Tmax,mu=paramRead[1],sigma=paramRead[2],lambda = paramRead[3])
  FVTmax <-  plnorm(Tmax,paramCtr[1],paramCtr[2])
  FminTmax <- 1- (1-FUTmax)*(1-FVTmax)
  FUT0 <- pemg(0,mu=paramRead[1],sigma=paramRead[2],lambda = paramRead[3])
  FVT0 <- plnorm(0,paramCtr[1],paramCtr[2])
  Fmin0 <- 1- (1-FUT0)*(1-FVT0)
  d <- fmin/(FminTmax-Fmin0)
  if(log){d <- log(fmin) - log(FminTmax-Fmin0) }
  return(d)
}

#-------------------------- Probability function of min(Ctr,EMG)
rExpData <- function(n,paramRead,paramCtr,Tmax){
  
  yread <- remg(n,mu = paramRead[1],sigma = paramRead[2], lambda = paramRead[3]) 
  y <- rCtrData(n,paramCtr,Tmax = Tmax ) 
  w.min <- which(yread<y)
  y[w.min] <- yread[w.min]
  
  w <- which( (y<0) | (y>Tmax))
  for (i in w){
    cond = FALSE
    while (!cond){
      yi <- min(remg(1,mu = paramRead[1],sigma = paramRead[2], lambda = paramRead[3]),rCtrData(1,paramCtr,Tmax= Inf))
      cond <-  ((yi>=0) | (yi<=Tmax))
    }
    y[i] = yi
  }
  return(y)
}

#----------------------- 
loglikExp <-function(logparamRead, paramCtr,Texp,Tmax){
  
  paramRead <- exp(logparamRead)
  U <- Texp[Texp<Tmax]
  LL <- -sum(dExpData(U,paramRead,paramCtr,Tmax=Tmax,log = TRUE))
  return(LL)
}
#----------------------- 
fitRepEmp_Texp_mixture <-function(logparamRead, paramCtr,Texp,Tmax){
  
  paramRead <- exp(logparamRead)
  U <- Texp[Texp<Tmax]
  t <- sort(U)
  FN_exp <- ecdf(t)
  yEmp <- FN_exp(t)
  yMod <- pExpData(t,paramRead,paramCtr,Tmax=Tmax)
  LL <- sum((yEmp-yMod)^2)
  return(LL)
}

#---------------------------- 
estimProc_Up_or_Dn <- function(Texp,Tctr,Tmax, paramCtr){
  
  #---------- moment
  resMoment <- estim_param_moment(Texp, Tctr,Tmax,paramCtr)
  paramRead_Moment <- resMoment$paramRead
  
  #---------- Fit repartition function
  #resFitRepEmp <- optim(par = log(paramRead_Moment), fn = fitRepEmp_Texp_mixture, paramCtr= paramCtr, Texp = Texp,Tmax = Tmax)
  #paramRead_FitRepEmp <-exp(resFitRepEmp$par)

  #---------- Fit max likelihood 
  #resMaxLoglik <- optim(par = resFitRepEmp$par, fn = loglikExp, paramCtr= paramCtr, Texp = Texp,Tmax = Tmax)
  estim = list(moment = paramExp_Moment) #, loglik =  exp(resMaxLoglik$par), fitFn  = exp(resFitRepEmp$par))
  return(estim)
}

