transfo_Gamma_Param <- function(mu,sigma){
  #mu : mean
  # sigma : sd
  # output : moment matching alpha and beta parameters
  beta <- mu/sigma^2
  alpha <- mu*beta
  return(c(alpha,beta))
  
}


#----------------------------------------------------------------
#----------- min(Exponentielle, Tmax)  : simulation
#----------------------------------------------------------------
rExpCensored = function(n,lambda,Tmax = Inf){
  e = rexp(n,lambda)
  e[e>=Tmax] = Tmax
  return(e)
}

#----------------------------------------------------------------
#----------- min(Exponentielle, Tmax)  : simulation
#----------------------------------------------------------------
pExpCensored = function(x,lambda,Tmax = Inf){
  e = pexp(x,lambda)
  e[x>=Tmax] = 1
  return(e)
}

#----------------------------------------------------------------
#----------- min(Exponentielle, Tmax)  : densité
#----------------------------------------------------------------
dExpCensored = function(x,lambda,Tmax = Inf,log = FALSE){
  
  if(length(x)==0){
    return(0)
    }else{
    if(Tmax==Inf){
      d = dexp(x,lambda,log)
      }else{
      d = 0*x
      d[x==Tmax] = 1-pexp(Tmax,lambda)
      d[x<Tmax] = dexp(x[x<Tmax],lambda)
      if(log){d <- log(d)}
  }
  return(d)
  }
}
 
#----------------------------------------------------------------
#--------------------------  densité min(Gamma,Exponentially Modified Gaussian) 
#-----------------------------------------------------------------

# TV_LD ~ min( UV,VV) with UV ~ Exp(lamda_ND), VV ~ Exp(lambda) + N(mu,sigma)   
dminExpExpplusGaussian <- function(x,lambda,mu,sigma,lambda_ND,Tmax= Inf,log = FALSE){
 
  # x: vector of observations
  # lambda: positive real number. Param of the exponential distribution
  # mu:  real number. Mean of the gaussian
  # sigma:  positive real number. SD of the gaussian
  # theta:  vector of 2 positive real number. Param of the gamma  distribution
  # output : vector of length x
  
  
  if(length(x) == 0){
    return(0)
  }else{
    
    F1z <- pemg(x,lambda=lambda,mu=mu,sigma=sigma)
    f1z <- demg(x,lambda=lambda,mu=mu,sigma=sigma)
    F2z <- pexp(x,lambda_ND)
    f2z <- dexp(x,lambda_ND)
    y <- f1z*(1-F2z) + f2z*(1-F1z)
    if(Tmax<Inf){
      FTmax <- pminExpExpplusGaussian(Tmax,lambda,mu,sigma,lambda_ND)
      y[y==Tmax] = 1-FTmax
      y[y>Tmax] = 0
    }
    if(log){y <- log(y)}
    return(y)
  }
}



#----------------------------------------------------------------
#--------------------------  Probability distribution min(Gamma,EMG) 
#----------------------------------------------------------------

# TV  ~ min( UV,VV) with UV ~ Gamma(theta[1], theta[2]), VV ~ Exp(lambda) + N(mu,sigma)   
pminExpExpplusGaussian <- function(x,lambda,mu,sigma,lambda_ND){
  # x: vector of real number
  # lambda: positive real number. Param of the exponential distribution
  # mu:  real number. Mean of the gaussian
  # sigma:  positive real number. SD of the gaussian
  # theta:  vector of 2 positive real number. Param of the gamma  distribution
  # output : vector of length x
  F1z <- pemg(x,lambda = lambda,mu = mu, sigma = sigma)
  F2z <- pexp(x,lambda_ND)
  1-(1-F1z)*(1-F2z)
  
  
  
}

#----------------------------------------------------------------
#--------------------------  Simulation of a sample from  min(Gamma,EMGamma + delta) 
#-----------------------------------------------------------------

# TV_LD ~ min( UV,VV) with UV ~ Exp(lambda_ND), VV ~ Exp(lambda_c) + Normal(k/lambda_e,sqrt(k)/lambda)   
rminExpExpplusGaussian <- function(n,param,k,kprime,color='red',Tmax = Inf){
  
  # n: size of sample
  lambda_ND <- ifelse(color=='red',param[2],param[1])
  lambda_c <- param[3]
  lambda_e <- param[4]
  nbcodons <- ifelse(color=='red',k+kprime,k)
  
  U <- rExpCensored(n,lambda_ND,Tmax)    # mort naturelle lente
  V <- rexp(n,lambda_c) + rnorm(n,nbcodons /lambda_e,sqrt(nbcodons)/lambda_e) # extinction by reading
  Y <- apply(cbind(U,V),1,min)
  Y[Y>Tmax] = Tmax
  return(Y)
  
}


#---------------------------------------------------------------------------------------------------
#----------Simulation of a sample from  pi_QD Gamma() + (1-pi_QD)*min(Gamma,EMGamma)  
#----------------------------------------------------------------------------------------------------

rOurModel <- function(n,param,k,kprime,color='red',Tmax = Inf){
  
  # n: number of simulation 
  # lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k:  Param of the Gamma. Number of codon
  # lambda_e:  positive real number. Param of the Gamma. Inverse of the mean length of lecture per Codon (élongation)
  # theta_ND:  vector of 2 positive real number. Param of the gamma  distribution representing the natural death
  # pi_QD : proportion of particles that will die very quickly
  # lambda_QD: 1/mean of time of quick death phenomenon 
  # output : vector of length n
  
  lambda_c <- param[3]
  lambda_e <- param[4]
  lambda_QD <- param[5]
  pi_QD <- param[6]
  lambda_ND <- ifelse(color=='red',param[2],param[1])
  nbcodons <- ifelse(color=='red',k+kprime,k)
  
 
  Z_QD <- sample(c(1,0),n,replace=TRUE,prob = c(pi_QD, 1-pi_QD)) 
  T_QD <- rExpCensored(n,lambda_QD)
  T_LD <- rminExpExpplusGaussian(n,param,k,kprime,color,Tmax)
  Y <- Z_QD * T_QD  + (1-Z_QD)*T_LD
  Y[Y>Tmax] <- Tmax 
  return(list(Y = Y,Z  = Z_QD))
}

#---------------------------------------------------------------------------------------------------
#----------Density of  pi_QD Gamma() + (1-pi_QD)*min(Gamma,EMGamma)  
#----------------------------------------------------------------------------------------------------
dOurModel <- function(x,param,k,kprime,color='red',Tmax = Inf,log = FALSE){
  
  # param = (lambda_ND,lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k,kprime:  Param of the Gamma. Number of codon
  # Tmax : truncature of data
  # output : vector of same length as x
  
  lambda_c <- param[3]
  lambda_e <- param[4]
  lambda_QD <- param[5]
  pi_QD <- param[6]
  
  lambda_ND <- ifelse(color=='red',param[2],param[1])
  nbcodons <- ifelse(color=='red',k+kprime,k)
  
  f <- rep(0,length(x))
  dQD <- dexp(x[x<Tmax],lambda_QD)
  dNQD <- dminExpExpplusGaussian(x[x<Tmax],lambda_c,mu= nbcodons/lambda_e ,sigma=sqrt(nbcodons)/lambda_e,lambda_ND)
  f[x<Tmax] <- dQD*pi_QD + (1-pi_QD)*dNQD
  
  if(Tmax < Inf){
    FTmax <- pi_QD*pexp(Tmax,lambda_QD) + (1-pi_QD)*pminExpExpplusGaussian(Tmax,lambda_c,mu= nbcodons/lambda_e ,sigma=sqrt(nbcodons)/lambda_e,lambda_ND)
    f[x==Tmax]  = 1- FTmax
  }
  if(log){return(log(f))}else{return(f)}
} 

#========================================================================
#----------------- Likelihood (TV, TR,ZR,ZV)
#========================================================================
log_lik <- function(log_param,data,Z){

  # log_param: log(lambda_ND_V),log(lambda_ND_R), log(lambda_c),log(lambda_e) [log(lambda_QD) pi_QD]

  Tmax <- data$Tmax
  k <- data$k
  kprime <- data$kprime
  

  d = 0
  #--------------------------------------
  lambda_ND_V <- exp(log_param[1])
  lambda_ND_R <- exp(log_param[2])
  lambda_c <- exp(log_param[3])
  lambda_e <- exp(log_param[4])      
  lambda_QD <- exp(log_param[5])
  #--------------------------------------
  d1 <- dExpCensored(data$T_Contr_V, lambda_ND_V,Tmax,log  = TRUE)
  d  <- sum(d1)

  #--------------------------------------
  d2 <- dExpCensored(data$T_Contr_R, lambda_ND_R,Tmax,log  = TRUE)
  d  <- d + sum(d2)
  #--------------------------------------
  ZV <- Z$ZV
  mu_V <- k/lambda_e
  sigma_V <- sqrt(k)/lambda_e
  U3_NQD <- dminExpExpplusGaussian(data$T_Exp_V[ZV==0], lambda = lambda_c,mu = mu_V,sigma = sigma_V,lambda_ND_V,Tmax,log = TRUE)
  d <- d + sum(U3_NQD)
  U3_QD <- dExpCensored(data$T_Exp_V[ZV==1],lambda_QD,Tmax,log=TRUE)
  d <- d  + sum(U3_QD)
  #--------------------------------------
  
  ZR <- Z$ZR
  mu_R <- (k+kprime)/lambda_e
  sigma_R <- sqrt(k+kprime)/lambda_e
  U4_NQD <-  dminExpExpplusGaussian(data$T_Exp_R[ZR==0], lambda = lambda_c,mu = mu_R,sigma = sigma_R,lambda_ND_R,Tmax,log = TRUE)
  d <- d + sum(U4_NQD) 
  U4_QD <- dExpCensored(data$T_Exp_R[ZR==1],lambda_QD,Tmax,log=TRUE)
  d <- d + sum(U4_QD)
  
  return(d)
}

#========================================================================
log_prior = function(log_param,hyperparams,withQD){
#========================================================================
  if('mean' %in% names(hyperparams)){
    d <- sum(dnorm(log_param[1:(4+withQD)], hyperparams$mean[1:(4+withQD)], hyperparams$sd[1:(4+withQD)],log = TRUE))
  }
  if('upperbound' %in% names(hyperparams)){
  d <- sum(dunif(log_param[1:(4+withQD)], hyperparams$lowerbound[1:(4+withQD)], hyperparams$upperbound[1:(4+withQD)],log = TRUE))
  }
  return(d)
}



#========================================================================
# plot Curves
#========================================================================

computationToPlotCurves <- function(x,param,k,kprime,Tmax = Inf){
  
  lambda_ND_V <- param[1]
  lambda_ND_R <- param[2]
  lambda_c <- param[3]
  lambda_e <- param[4]
  lambda_QD <- param[5]
  pi_QD <- param[6]
  
  
  if (pi_QD== 0){nbCurves <- 3*2}else{nbCurves <- 5*2}
  
  P <- as.data.frame(rep(x,nbCurves)); names(P)='time'
  P$Color_Part <- rep(c('Green','Red'),each=nrow(P)/2)
  
  myDensity_V <- matrix(0,length(x),nbCurves/2)
  myDensity_V[,1] <- dExpCensored(x,lambda_ND_V,Tmax)
  myDensity_V[,2] <- demg(x,lambda_c,mu = k/lambda_e,sigma = sqrt(k)/lambda_e)
  myDensity_V[,3] <- dminExpExpplusGaussian(x,lambda_c,mu = k/lambda_e,sigma = sqrt(k)/lambda_e,lambda_ND_V,Tmax)
  

  myDensity_R <- matrix(0,length(x),nbCurves/2)
  myDensity_R[,1] <- dExpCensored(x,lambda_ND_R,Tmax)
  myDensity_R[,2] <- demg(x,lambda_c,mu = (k+kprime)/lambda_e,sigma = sqrt(k+kprime)/lambda_e)
  myDensity_R[,3] <- dminExpExpplusGaussian(x,lambda_c,mu = (k+kprime)/lambda_e,sigma = sqrt(k+kprime)/lambda_e,lambda_ND_R,Tmax)

  
  myCurves = rep(c('1.Natural Death',
                   '2. Sum Arrival + reading',
                   '3. min(ND,Reading)'),each=length(x))
                 
  if(pi_QD>0){
    myDensity_V[,4] <-  myDensity_R[,4]  <- dexp(x,lambda_QD)
    myDensity_V[,5] <- dOurModel(x,param,k,kprime,Tmax,color='green')
    myDensity_R[,5] <- dOurModel(x,param,k,kprime,Tmax,color='red')
    myCurves <- c(myCurves,rep(c('4. Quick Death','5. Final model'),each=length(x)))
  }
  
  
  P$density <- c(c(myDensity_V),c(myDensity_R))
  P$Curves <- as.factor(rep(myCurves,2))
  P$Color_Part <- as.factor(P$Color_Part)
  return(P)
              
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

