library(emg)
library(moments)


#----------------------------------------------------------------
#-----------  : simulation
#----------------------------------------------------------------
rExpCensored = function(n,lambda,Tmax = Inf, piTrunc=0){
  
  if(piTrunc==0){Tmax=Inf}
  if(Tmax==Inf){piTrunc=0}
  
  Z <- rbinom(n,1,piTrunc)
  Y <- -1/lambda*log(1-runif(n)*(1-exp(-lambda*Tmax)))
  e <- (1-Z)*Y
  if(Tmax<Inf){e <- e + Z*Tmax} 
  return(e)
}


remgCensored = function(n,mu, sigma,lambda,Tmax = Inf, piTrunc=0){
  
  if(piTrunc==0){Tmax=Inf}
  if(Tmax==Inf){piTrunc=0}
  
  
  Z <- rbinom(n,1,piTrunc)
  N <- floor(n*1/pemg(Tmax,mu, sigma,lambda))
  Y <- remg(N,mu,sigma,lambda)
  Y <- Y[Y<Tmax]; 
  if(length(Y)>n){Y = Y[1:n]}
  while(length(Y)<n){
    y <- remg(1,mu,sigma,lambda)
    if(y<Tmax){Y <- c(Y,y)}
  }
  e <- (1-Z)*Y
  if(Tmax<Inf){e <- e + Z*Tmax} 
  return(e)
}


# 
# #----------------------------------------------------------------
# #----------- min(Exponentielle, Tmax)  : simulation
# #----------------------------------------------------------------
# pExpCensored = function(x,lambda,Tmax = Inf){
#   e = pexp(x,lambda)
#   e[x>=Tmax] = 1
#   return(e)
# }

#----------------------------------------------------------------
#----------- min(Exponentielle, Tmax)  : simulation
#----------------------------------------------------------------
pExpCensored = function(x,lambda,Tmax = Inf,piTrunc = 0){
  
  e = (1-piTrunc)*pexp(x,lambda)/pexp(Tmax,lambda)
  e[x>=Tmax] = 1
  return(e)
}

pemgCensored = function(x,mu, sigma,lambda,Tmax = Inf, piTrunc=0){
  e = (1-piTrunc)/pemg(Tmax,mu,sigma,lambda)*pemg(x,mu,sigma,lambda)
  e[x>=Tmax] = 1
  return(e)
}

# #----------------------------------------------------------------
# #----------- min(Exponentielle, Tmax)  : densité
# #----------------------------------------------------------------
# dExpCensored = function(x,lambda,Tmax = Inf,log = FALSE){
#   
#   if(length(x)==0){
#     return(0)
#     }else{
#     if(Tmax==Inf){
#       d = dexp(x,lambda,log)
#       }else{
#       d = 0*x
#       d[x==Tmax] = 1-pexp(Tmax,lambda)
#       d[x<Tmax] = dexp(x[x<Tmax],lambda)
#       if(log){d <- log(d)}
#   }
#   return(d)
#   }
#}

#----------------------------------------------------------------
#----------- min(Exponentielle, Tmax)  : densité
#----------------------------------------------------------------
dExpCensored  = function(x,lambda,Tmax = Inf,piTrunc=0,log = FALSE){
  
  if(length(x)==0){
    return(0)
  }else{
    if((Tmax==Inf) | (piTrunc==0)){
      d = dexp(x,lambda,log)
    }else{
      d = 0*x
      d[x==Tmax] = piTrunc
      d[x<Tmax] = (1-piTrunc)/pexp(Tmax,lambda) * dexp(x[x<Tmax],lambda)
      if(log){d <- log(d)}
    }
    return(d)
  }
}
#----------------------------------------------------------------
demgCensored  = function(x,mu,sigma,lambda,Tmax = Inf,piTrunc=0,log = FALSE){
  
  if(is.null(Tmax)){Tmax==Inf}
  if(Tmax == -Inf){Tmax = Inf}
  if(length(x)==0){
    return(0)
  }else{
    if((Tmax==Inf) | (piTrunc==0)){
      d = demg(x,mu,sigma,lambda,log)
    }else{
      d = 0*x
      d[x==Tmax] = piTrunc
      d[x<Tmax] = (1-piTrunc)/pemg(Tmax,mu,sigma,lambda) * demg(x[x<Tmax],mu,sigma,lambda)
      if(log){d <- log(d)}
    }
    return(d)
  }
}


 
#----------------------------------------------------------------
#--------------------------  densité min(Gamma,Exponentially Modified Gaussian) 
#-----------------------------------------------------------------

# TV_LD ~ min( UV,VV) with UV ~ Exp(lamda_ND), VV ~ Exp(lambda) + N(mu,sigma)   
dminExpExpplusGaussian <- function(x,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax= Inf,log = FALSE){
 
  # x: vector of observations
  # lambda: positive real number. Param of the exponential distribution
  # mu:  real number. Mean of the gaussian
  # sigma:  positive real number. SD of the gaussian
  # theta:  vector of 2 positive real number. Param of the gamma  distribution
  # output : vector of length x
  
  
  if(length(x) == 0){
    return(0)
  }else{
    F1z <- pemgCensored(x,mu=mu,sigma=sigma,lambda=lambda,Tmax,piTrunc)
    f1z <- demgCensored(x,mu=mu,sigma=sigma,lambda=lambda,Tmax,piTrunc)
    F2z <- pExpCensored(x,lambda_ND,Tmax,piTrunc_ND)
    f2z <- dExpCensored(x,lambda_ND,Tmax,piTrunc_ND)
    y <- f1z*(1-F2z) + f2z*(1-F1z)
    y[x==Tmax] = piTrunc * piTrunc_ND
    if(log){y <- log(y)}
    return(y)
  }
}



#----------------------------------------------------------------
#--------------------------  Probability distribution min(Gamma,EMG) 
#----------------------------------------------------------------

# TV  ~ min( UV,VV) with UV ~ Gamma(theta[1], theta[2]), VV ~ Exp(lambda) + N(mu,sigma)   
pminExpExpplusGaussian <- function(x,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax= Inf){
  # x: vector of real number
  # lambda: positive real number. Param of the exponential distribution
  # mu:  real number. Mean of the gaussian
  # sigma:  positive real number. SD of the gaussian
  # theta:  vector of 2 positive real number. Param of the gamma  distribution
  # output : vector of length x
  F1z <- pemgCensored(x,mu=mu,sigma=sigma,lambda=lambda,Tmax,piTrunc)
  F2z <- pExpCensored(x,lambda_ND,Tmax,piTrunc_ND)
  1-(1-F1z)*(1-F2z)
  
  
  
}

#----------------------------------------------------------------
#--------------------------  Simulation of a sample from  min(Gamma,EMGamma + delta) 
#-----------------------------------------------------------------

# TV_LD ~ min( UV,VV) with UV ~ Exp(lambda_ND), VV ~ Exp(lambda_c) + Normal(k/lambda_e,sqrt(k)/lambda)   


rminExpExpplusGaussian <- function(n,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax){
  
  U <- rExpCensored(n,lambda_ND,Tmax,piTrunc_ND)    # mort naturelle lente
  V <- remgCensored(n,mu, sigma, lambda=lambda,Tmax,piTrunc) # extinction by reading
  Y <- apply(cbind(U,V),1,min)
  return(Y)
}



#---------------------------------------------------------------------------------------------------
#----------Simulation of a sample from  pi_QD Gamma() + (1-pi_QD)*min(Exp,Exp+Gaussian)  
#----------------------------------------------------------------------------------------------------

rminExpExpplusGaussianwithQD <- function(n,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,pi_QD,lambda_QD,Tmax){
  
  # n: number of simulation 
  # output : vector of length n
  Z_QD <- sample(c(1,0),n,replace=TRUE,prob = c(pi_QD, 1-pi_QD)) 
  T_QD <- rexp(n,lambda_QD)
  T_LD <- rminExpExpplusGaussian(n,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax)
  Y <- Z_QD * T_QD  + (1-Z_QD)*T_LD
  return(list(Y = Y,Z  = Z_QD))
}

#---------------------------------------------------------------------------------------------------
#----------density of a sample from  pi_QD Gamma() + (1-pi_QD)*min(Exp,Exp+Gaussian)  
#----------------------------------------------------------------------------------------------------

dminExpExpplusGaussianwithQD <- function(x,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,pi_QD,lambda_QD,Tmax,log = FALSE){
  
  # x: values 
  # output : vector of length n
  if(length(x) == 0){
    return(0)
  }else{
    if(pi_QD>0){
      d_QD <- dexp(x,lambda_QD)
      d_LD <- dminExpExpplusGaussian(x,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax)
      y <- pi_QD* d_QD  + (1-pi_QD)*d_LD
      if(log){y <- log(y)}
      }else{
      y <- dminExpExpplusGaussian(x,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax,log)
    }
    return(y)
  }
  
}

#---------------------------------------------------------------------------------------------------
#----------Density of  pi_QD Gamma() + (1-pi_QD)*min(Gamma,EMGamma)  
#----------------------------------------------------------------------------------------------------
rOurModelExp <- function(n,param,k,kprime,color='red',Tmax = Inf){
  
  lambda_ND <- ifelse(color=='red',param[3],param[1])
  piTrunc_ND <- ifelse(color=='red',param[4],param[2])
  nbcodons <- ifelse(color=='red',k+kprime,k)
  lambda_c <- param[5]
  lambda_e <- param[6]
  piTrunc_Read <- param[7]
  lambda_QD <- param[8]
  pi_QD <- param[9]
  
  mu= nbcodons/lambda_e
  sigma=sqrt(nbcodons)/lambda_e
  
  E <- rminExpExpplusGaussianwithQD(n,mu,sigma,lambda_c,piTrunc_Read,lambda_ND,piTrunc_ND,pi_QD,lambda_QD,Tmax) 
  
  return(E)
}




dOurModelExp <- function(x,param,k,kprime,color='red',Tmax = Inf,log = FALSE){
  
  # param = (lambda_ND_V,piTrunc_ND_V,lambda_ND_R,piTrunc_ND_R , lambda_c,lambda_e,piTrunc_Read, lambda_QD,pi_QD 
  
  
  #: Param of the exponential distribution for the time to go to the ADN
  # k,kprime:  Param of the Gamma. Number of codon
  # Tmax : truncature of data
  # output : vector of same length as x
  
  
  lambda_ND <- ifelse(color=='red',param[3],param[1])
  piTrunc_ND <- ifelse(color=='red',param[4],param[2])
  nbcodons <- ifelse(color=='red',k+kprime,k)
  
  
  lambda_c <- param[5]
  lambda_e <- param[6]
  piTrunc_Read <- param[7]
  lambda_QD <- param[8]
  pi_QD <- param[9]
  
  mu= nbcodons/lambda_e
  sigma=sqrt(nbcodons)/lambda_e
  
  f <- dminExpExpplusGaussianwithQD(x,mu,sigma,lambda_c,piTrunc_Read,lambda_ND,piTrunc_ND,pi_QD,lambda_QD,Tmax,log) 
  return(f)
} 

# #========================================================================
# #----------------- Likelihood (TV, TR,ZR,ZV)
# #========================================================================
# log_lik <- function(log_param,data,ZQD){
# 
#   # log_param: log(lambda_ND_V),log(lambda_ND_R), log(lambda_c),log(lambda_e), log(lambda_QD) pi_QD
# 
#   Tmax <- data$Tmax
#   k <- data$k
#   kprime <- data$kprime
#   
# 
#   d = 0
#   #--------------------------------------
#   lambda_ND_V <- exp(log_param[1])
#  
#   
#   lambda_ND_R <- exp(log_param[2])
#   lambda_c <- exp(log_param[3])
#   lambda_e <- exp(log_param[4])      
#   lambda_QD <- exp(log_param[5])
#   #--------------------------------------
#   d1 <- dExpCensored(data$T_Contr_V, lambda_ND_V,Tmax,log  = TRUE)
#   d  <- sum(d1)
# 
#   #--------------------------------------
#   d2 <- dExpCensored(data$T_Contr_R, lambda_ND_R,Tmax,log  = TRUE)
#   d  <- d + sum(d2)
#   #--------------------------------------
#   ZV <- Z$ZV
#   mu_V <- k/lambda_e
#   sigma_V <- sqrt(k)/lambda_e
#   U3_NQD <- dminExpExpplusGaussian(data$T_Exp_V[ZV==0], lambda = lambda_c,mu = mu_V,sigma = sigma_V,lambda_ND_V,Tmax,log = TRUE)
#   d <- d + sum(U3_NQD)
#   U3_QD <- dExpCensored(data$T_Exp_V[ZV==1],lambda_QD,Tmax,log=TRUE)
#   d <- d  + sum(U3_QD)
#   #--------------------------------------
#   
#   ZR <- Z$ZR
#   mu_R <- (k+kprime)/lambda_e
#   sigma_R <- sqrt(k+kprime)/lambda_e
#   U4_NQD <-  dminExpExpplusGaussian(data$T_Exp_R[ZR==0], lambda = lambda_c,mu = mu_R,sigma = sigma_R,lambda_ND_R,Tmax,log = TRUE)
#   d <- d + sum(U4_NQD) 
#   U4_QD <- dExpCensored(data$T_Exp_R[ZR==1],lambda_QD,Tmax,log=TRUE)
#   d <- d + sum(U4_QD)
#   
#   return(d)
# }
# 
# #========================================================================
# log_prior = function(log_param,hyperparams,withQD){
# #========================================================================
#   if('mean' %in% names(hyperparams)){
#     d <- sum(dnorm(log_param[1:(4+withQD)], hyperparams$mean[1:(4+withQD)], hyperparams$sd[1:(4+withQD)],log = TRUE))
#   }
#   if('upperbound' %in% names(hyperparams)){
#   d <- sum(dunif(log_param[1:(4+withQD)], hyperparams$lowerbound[1:(4+withQD)], hyperparams$upperbound[1:(4+withQD)],log = TRUE))
#   }
#   return(d)
# }


#=====================================================================
#Logit et invlogit , from_log_param_to_param and inverse
#======================================================================

#=============
logit = function(x){
  
  y <- x
  w <- which(x==-Inf)
  if(length(w)>0){
    y[w] = 0;
    y[-w] = log(x[-w]/(1-x[-w]))
  }else{
    y = log(x/(1-x))
  }
  return(y)
}
  
#=============
invlogit = function(x){1/(1+exp(-x))}

#=============
from_param_to_log_param <- function(param){
  
  log_param <- param
  if(is.vector(param)){
    log_param[c(1,3,5,6,8)] <- log(param[c(1,3,5,6,8)])
    log_param[c(2,4,7,9)]<- logit(param[c(2,4,7,9)])
  }
  if(is.matrix(param)){
    log_param[,c(1,3,5,6,8)] <- log(param[,c(1,3,5,6,8)])
    log_param[,c(2,4,7,9)]<- logit(param[,c(2,4,7,9)])
  }
  names(log_param) <- paste0(c('log_','logit_','log_','logit_','log_','log_','ogit_','log_','logit_'),names(param))
  return(log_param)
}
#=============
from_logparam_to_param <- function(log_param){
  
  if(is.vector(log_param)){
    param <- log_param
    param[c(1,3,5,6,8)] <- exp(log_param[c(1,3,5,6,8)])
    param[c(2,4,7,9)]<- invlogit(log_param[c(2,4,7,9)])
  }
  if(is.matrix(log_param)){
    param <- log_param
    param[,c(1,3,5,6,8)] <- exp(log_param[,c(1,3,5,6,8)])
    param[,c(2,4,7,9)]<- invlogit(log_param[,c(2,4,7,9)])
  }
  names(param) <- c('lambda_ND_V','piTrunc_ND_V','lambda_ND_R','piTrunc_ND_R','lambda_c','lambda_e','piTrunc_Read', 'lambda_QD','pi_QD')
  
  return(param)
}
#=============



transfo_Gamma_Param <- function(mu,sigma){
  #mu : mean
  # sigma : sd
  # output : moment matching alpha and beta parameters
  beta <- mu/sigma^2
  alpha <- mu*beta
  return(c(alpha,beta))
  
}

#========================================================================
#----------- Estimators of Gamma parameters
#========================================================================
estim_param_Gamma <- function(X){
  EX <- mean(X)
  VX <-var(X)
  alpha_hat <- EX^2/VX
  beta_hat <- alpha_hat/EX
  return(c(alpha_hat,beta_hat))
}


#========================================================================
#----------- Estimators of EMG distribution
#========================================================================
estim_param_emg <- function(X){
  
  m <- mean(X)
  s <- sd(X)
  gamma1 <- skewness(X)
  rho <- (gamma1/2)^(1/3)
  muhat <- m-s*rho
  sigmahat <- s*sqrt(1-rho^2)
  lambdahat <- 1/(s*rho)
  return(c(muhat,sigmahat,lambdahat))
}



#========================================================================
# Expectation of Censored exponential dist
#========================================================================
espExpCensored <-function(lambda,Tmax){
  1/lambda-Tmax*(1-pexp(Tmax,lambda))/pexp(Tmax,lambda)
}
#========================================================================
init_param <- function(data){
#========================================================================
  param <- rep(0,9)
  #------- data Contr V
  L_V <- sort(1/seq(1,500,by=0.01))
  D_V <-  espExpCensored(L_V,Tmax = data$Tmax_Contr_V)
  param[1] <- L_V[which.min(abs(D_V-mean(data$T_Contr_V[data$T_Contr_V<data$Tmax_Contr_V])))]
  param[2] <- mean(data$T_Contr_V == data$Tmax_Contr_V) #pi_trunc_ND_V
  
  #------- data Contr R
  L_R <- sort(1/seq(1,500,by=0.01))
  D_R <-  espExpCensored(L_R,Tmax =data$Tmax_Contr_R)
  param[3] <- L_R[which.min(abs(D_R-mean(data$T_Contr_R[data$T_Contr_R<data$Tmax_Contr_R])))]
  param[4] <- mean(data$T_Contr_R == data$Tmax_Contr_R) #pi_trunc_ND_R
  
  #--------------------------- About lambda_e, lambde_c, pitrunc
  
  FNV <- ecdf(data$T_Exp_V)
  plot(FNV)
  UV <- function(x){
    1-(1-FNV(x))/(1-pExpCensored(x,lambda=param[1],Tmax =data$Tmax_Contr_V,piTrunc = param[2]))
  }
  plot(UV,0,100)
  abs <- 1:(data$Tmax_Contr_V-1)
  DV <- c(0,diff(UV(abs)))
  DV <- DV*(DV>0)
  DV <- DV/sum(DV)  
  MV <- sum(DV*abs)


 
  FNR <- ecdf(data$T_Exp_R)
  UR <- function(x){
    1-(1-FNR(x))/(1-pExpCensored(x,lambda=param[3],Tmax =data$Tmax_Contr_R,piTrunc = param[4]))
  }
  abs <- 1:(data$Tmax_Contr_R-1)
  DR <- c(0,diff(UR(abs)))
  DR <- DR*(DR>0)
  DR <- DR/sum(DR)
  MR <- sum(DR*abs)
  
  
  param[c(7)] <-  0.5*(mean(data$T_Exp_V == data$Tmax_Exp_V)/param[2] + mean(data$T_Exp_R == data$Tmax_Exp_R)/param[4])
  param[6] <- data$kprime/(MR-MV)
  param[5] <- 1/(MV-data$k/param[6]) 
  
  param[c(8)] <- 1
  param[c(9)] <- 0
  names(param) <- c('lambda_ND_V','piTrunc_ND_V','lambda_ND_R','piTrunc_ND_R','lambda_c','lambda_e','piTrunc_Read', 'lambda_QD','pi_QD')
  return(param)
}
