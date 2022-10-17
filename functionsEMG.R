transfo_Gamma_Param <- function(mu,sigma){
  #mu : mean
  # sigma : sd
  # output : moment matching alpha and beta parameters
  beta <- mu/sigma^2
  alpha <- mu*beta
  return(c(alpha,beta))
  
}

#----------------------------------------------------------------
#--------------------------  densité min(Gamma,Exponentially Modified Gaussian) 
#-----------------------------------------------------------------

# TV_LD ~ min( UV,VV) with UV ~ Gamma(theta[1], theta[2]), VV ~ Exp(lambda) + N(mu,sigma)   
dminGammaEMGaussian <- function(x,lambda,mu,sigma,theta){
 
  if(length(theta)==1){theta=c(1,theta)}
  
  # x: vector of observations
  # lambda: positive real number. Param of the exponential distribution
  # mu:  real number. Mean of the gaussian
  # sigma:  positive real number. SD of the gaussian
  # theta:  vector of 2 positive real number. Param of the gamma  distribution
  # output : vector of length x
  
  
  F1z <- pemg(x,lambda=lambda,mu=mu,sigma=sigma)
  f1z <- demg(x,lambda=lambda,mu=mu,sigma=sigma)
  
  F2z <- pgamma(x,theta[1],theta[2])
  f2z <- dgamma(x,theta[1],theta[2])
  y <- f1z*(1-F2z) + f2z*(1-F1z)
  return(y)
}



#----------------------------------------------------------------
#--------------------------  Probability distribution min(Gamma,EMG) 
#----------------------------------------------------------------

# TV  ~ min( UV,VV) with UV ~ Gamma(theta[1], theta[2]), VV ~ Exp(lambda) + N(mu,sigma)   
pminGammaEMGaussian <- function(x,lambda,mu,sigma,theta){
  
  if(length(theta)==1){theta=c(1,theta)}
  # x: vector of real number
  # lambda: positive real number. Param of the exponential distribution
  # mu:  real number. Mean of the gaussian
  # sigma:  positive real number. SD of the gaussian
  # theta:  vector of 2 positive real number. Param of the gamma  distribution
  # output : vector of length x

  F1z <- pemg(x,lambda = lambda,mu = mu, sigma = sigma)
  F2z <- pgamma(x,theta[1],theta[2])
  1-(1-F1z)*(1-F2z)
  
  
  
}

#----------------------------------------------------------------
#--------------------------  Simulation of a sample from  min(Gamma,EMGamma + delta) 
#-----------------------------------------------------------------

# TV_LD ~ min( UV,VV) with UV ~ Gamma(theta[1], theta[2]), VV ~ delta + Exp(lambda_c) + Gamma(k,lambda_e)   
rminGammaEMGamma <- function(n,lambda_c,k,lambda_e,theta_ND,delta = 0){
  
  # n: size of sample
  # lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k:  Param of the Gamma. Number of codon
  # lambda_e:  positive real number. Param of the Gamma. Inverse of the mean length of lecture per Codon (élongation)
  # theta_ND:  vector of 2 positive real number. Param of the gamma  distribution representing the natural death
  # delta : parameter of lag 
  # output : vector of length n
  
  if(length(theta_ND)==1){theta_ND=c(1,theta_ND)}
  U <- rgamma(n,shape = theta_ND[1],rate = theta_ND[2])    # mort naturelle lente
  V <- delta + rexp(n,lambda_c) + rnorm(n,k/lambda_e,sqrt(k)/lambda_e) # extinction by reading
  Y <- apply(cbind(U,V),1,min)
  return(Y)
  
}


#---------------------------------------------------------------------------------------------------
#----------Simulation of a sample from  pi_QD Gamma() + (1-pi_QD)*min(Gamma,EMGamma)  
#----------------------------------------------------------------------------------------------------

rOurModel <- function(n,lambda_c,k,lambda_e,theta_ND,delta=0,pi_QD=0,theta_QD=c(1,1)){
  
  
  # x: vector of real number
  # lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k:  Param of the Gamma. Number of codon
  # lambda_e:  positive real number. Param of the Gamma. Inverse of the mean length of lecture per Codon (élongation)
  # theta_ND:  vector of 2 positive real number. Param of the gamma  distribution representing the natural death
  # delta : parameter of lag 
  # pi_QD : proportion of particles that will die very quickly
  # lambda_QD: 1/mean of time of quick death phenomenon 
  # output : vector of length n
  
  if(length(theta_QD)==1){theta_QD=c(1,theta_QD)}
  if(length(theta_ND)==1){theta_ND=c(1,theta_ND)}
  
  Z_QD <- sample(c(1,0),n,replace=TRUE,prob = c(pi_QD, 1-pi_QD)) 
  T_QD <- rgamma(n,shape = theta_QD[1], rate = theta_QD[2])
  T_LD <- rminGammaEMGamma(n,lambda_c,k,lambda_e,theta_ND,delta)
  Y <- Z_QD * T_QD  + (1-Z_QD)*T_LD
  return(list(Y = Y,Z  = Z_QD))
}

#---------------------------------------------------------------------------------------------------
#----------Density of  pi_QD Gamma() + (1-pi_QD)*min(Gamma,EMGamma)  
#----------------------------------------------------------------------------------------------------
dOurModel <- function(x,lambda_c,k,lambda_e,theta_ND,delta=0,pi_QD,theta_QD){
  
  # lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k:  Param of the Gamma. Number of codon
  # lambda_e:  positive real number. Param of the Gamma. Inverse of the mean length of lecture per Codon (élongation)
  # theta_ND:  vector of 2 positive real numbers. Param of the gamma  distribution representing the natural death
  # delta : parameter of lag 
  # pi_QD : proportion of particles that will die very quickly
  # theta_QD: vector of 2 positive real numbers. Param of the gamma  distribution representing the quick death
  # output : vector of same length as x
  
  if (length(theta_ND)==1){theta_ND = c(1,theta_ND)}
  if (length(theta_QD)==1){theta_QD = c(1,theta_QD)}
  d1 <- dgamma(x,shape = theta_QD[1], rate = theta_QD[2])
  #d1 <- dexp(x,theta_QD[2])
  d2 <- dminGammaEMGaussian(x,lambda_c,mu=k/lambda_e + delta,sigma=sqrt(k)/lambda_e,theta = theta_ND)
  return(d1*pi_QD + (1-pi_QD)*d2)
} 

# 
# 


#----------------------------------------------------------------
#----------- Estimators of Gamma parameters
#----------------------------------------------------------------
estim_param_Gamma <- function(X){
  EX <- mean(X)
  VX <-var(X)
  alpha_hat <- EX^2/VX
  beta_hat <- alpha_hat/EX
  return(c(alpha_hat,beta_hat))
}


 

#----------------- Likelihood (TV, TR)
log_lik_withoutQD <- function(log_param,data=list()){
  
  # log_param: log(lambda_ND_V),log(lambda_ND_R), log(lambda_c),log(lambda_e)
  
  lambda_ND_V <- exp(log_param[1])
  lambda_ND_R <- exp(log_param[2])
  lambda_c <- exp(log_param[3])
  lambda_e <- exp(log_param[4])
  
  k <- data$k
  kprime <- data$kprime
  
  d <- 0 
  
  if (!is.null(data$T_Contr_V)){
    d1 <- sum(dexp(data$T_Contr_V, lambda_ND_V,log  = TRUE))
    d<- d + d1
  }
  
  
  if (!is.null(data$T_Contr_R)){
    d2 <- sum(dexp(data$T_Contr_R, lambda_ND_R,log  = TRUE))
    d<- d + d2
  }
  
  
  if (!is.null(data$T_Exp_V)){
    mu_V <- k/lambda_e
    sigma_V <- sqrt(k)/lambda_e
    U3 <- dminGammaEMGaussian(data$T_Exp_V, lambda = lambda_c,mu = mu_V,sigma = sigma_V,c(1,lambda_ND_V))
    d <- d + sum(log(U3)) }
  
  if (!is.null(data$T_Exp_R)){
    mu_R <- (k+kprime)/lambda_e
    sigma_R <- sqrt(k+kprime)/lambda_e
    U4 <- dminGammaEMGaussian(data$T_Exp_R, lambda = lambda_c,mu = mu_R,sigma = sigma_R,c(1,lambda_ND_R))
    d <- d + sum(log(U4)) 
  }
  
  return(d)
  }



#----------------- Likelihood (TV, TR)
log_lik_withQD <- function(log_param,data=list(),ZV,ZR){
  
  # log_param: log(lambda_ND_V),log(lambda_ND_R), log(lambda_c),log(lambda_e)
  
  lambda_ND_V <- exp(log_param[1])
  lambda_ND_R <- exp(log_param[2])
  lambda_c <- exp(log_param[3])
  lambda_e <- exp(log_param[4])
  lambda_QD <- exp(log_param[5])
  pi_QD <- log_param[6]
  
  k <- data$k
  kprime <- data$kprime
  
 
  
  
  d <- 0 
  
  if (!is.null(data$T_Contr_V)){
    d1 <- sum(dexp(data$T_Contr_V, lambda_ND_V,log  = TRUE))
    d<- d + d1
  }
  
  if (!is.null(data$T_Contr_R)){
    d2 <- sum(dexp(data$T_Contr_R, lambda_ND_R,log  = TRUE))
    d<- d + d2
  }
  
  if (!is.null(data$T_Exp_V)){
    mu_V <- k/lambda_e
    sigma_V <- sqrt(k)/lambda_e
    U3_NQD <- log(dminGammaEMGaussian(data$T_Exp_V, lambda = lambda_c,mu = mu_V,sigma = sigma_V,c(1,lambda_ND_V)))
    U3_QD <- dexp(data$T_Exp_V,lambda_QD,log=TRUE)
    d <- d + sum(ZV*U3_QD + (1-ZV)*U3_NQD)
  }
  
  if (!is.null(data$T_Exp_R)){
    mu_R <- (k+kprime)/lambda_e
    sigma_R <- sqrt(k+kprime)/lambda_e
    U4_NQD <- dminGammaEMGaussian(data$T_Exp_R, lambda = lambda_c,mu = mu_R,sigma = sigma_R,c(1,lambda_ND_R))
    U4_QD <- dexp(data$T_Exp_R,lambda_QD,log=TRUE)
    d <- d + sum(ZR*U4_QD + (1-ZR)*U4_NQD)
  } 
  return(d)
}



#----------------------------------------------------------------------
# 
# log_likelihood_ourModel_onecolor <- function(param,theta_ND,k,Times){
#   
#   # param: log(lambda_c),log(lambda_e), qnorm(pi_QD) , log(theta_QD) (of length 2) 
#   
#   lambda_c <- exp(param[1])
#   lambda_e <- exp(param[2])
#   pi_QD <- pnorm(param[3])
#   theta_QD <- exp(param[4:5])
#   
#   LT <- dOurModel(Times,lambda_c,k,lambda_e,theta_ND,delta=0,pi_QD,theta_QD)
#   return(- sum(LT))
# }

#----------------------------------------------------------------------

# log_likelihood_ourModel_allcolors <- function(log_param,k,kprime,data){
#   
#   # param: log(lambda_c),log(lambda_e), qnorm(pi_QD) , log(theta_QD) (of length 2) 
#   # theta_ND  : params of natural death for each color (G then R)
#   lambda_c <- exp(param[3]) 
#   lambda_e <- exp(param[4])
# #  pi_QD <- pnorm(param[3])
# #  theta_QD <- exp(param[4:5])
#   
#   theta_ND_V <- c(1,log_param[1])
#   LTV <- dOurModel(data$TV,lambda_c,k         ,lambda_e,c(1,log_param[1]),delta=0,pi_QD = 0,theta_QD)
#   LTR <- dOurModel(data$TR,lambda_c,k + kprime,lambda_e,theta_ND[3:4],delta=0,pi_QD,theta_QD)
#   
#   return(- sum(LTV) - sum(LTR))
# }

#------------------------------------------------------------------------- 
# plot Curves

computationToPlotCurves <- function(x,lambda_ND,lambda_c,lambda_e,nbCodons,pi_QD = NULL,theta_QD= NULL){
  
  theta_ND <- c(1,lambda_ND)
  if (is.null(theta_QD)){nbCurves <- 5}else{nbCurves <- 7}
  P <- as.data.frame(rep(x,nbCurves)); names(P)='time'
  
  myDensity<- c(dgamma(x,theta_ND[1],theta_ND[2]),
              dexp(x,lambda_c),
              dgamma(x,nbCodons,lambda_e), 
              demg(x,lambda_c,mu = nbCodons/lambda_e + delta,sigma = sqrt(nbCodons)/lambda_e),
              dminGammaEMGaussian(x,lambda_c,mu = nbCodons/lambda_e + delta,sigma = sqrt(nbCodons)/lambda_e,theta=c(1,lambda_ND)))
  
  myCurves = rep(c('0.Natural Death',
                   '1. Arrival Time',
                   '2. Reading k codons',
                   '3. Sum Arrival + reading',
                   '4. min(ND,Reading)'),each=length(x))
                 
  if(!is.null(theta_QD)){
    mydensity  <- c(mydensity ,dgamma(x,theta_QD[1],theta_QD[2]),
                    dOurModel(x,lambda_c,nbCodons,lambda_e,theta_ND,delta,pi_QD,theta_QD))
    myCurves <- c(myCurves,rep(c('5. Quick Death','6. Final model'),each=length(x)))
    
  }
  
  P$density <- myDensity
  P$Curves <- myCurves
  return(P)
              
}



