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
  f1z*(1-F2z) + f2z*(1-F1z)
}



#----------------------------------------------------------------
#--------------------------  Probability distribution min(Gamma,EMG) 
#----------------------------------------------------------------

# TV  ~ min( UV,VV) with UV ~ Gamma(theta[1], theta[2]), VV ~ Exp(lambda) + N(mu,sigma)   
pminGammaEMGaussian <- function(x,lambda,mu,sigma,theta){
  
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
  
  
  U <- rgamma(n,shape = theta_ND[1],rate = theta_ND[2])    # mort naturelle lente
  V <- delta + rexp(n,lambda_c) + rgamma(n,k,lambda_e) # extinction by reading
  Y <- apply(cbind(U,V),1,min)
  return(Y)
  
}


#----------------------------------------------------------------
#--------------------------  Simulation of a sample from  min(Gamma,EMGamma) 
#-----------------------------------------------------------------

rOurModel <- function(n,lambda_c,k,lambda_e,theta_ND,delta,pi_QD,theta_QD){
  
  
  # x: vector of real number
  # lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k:  Param of the Gamma. Number of codon
  # lambda_e:  positive real number. Param of the Gamma. Inverse of the mean length of lecture per Codon (élongation)
  # theta_ND:  vector of 2 positive real number. Param of the gamma  distribution representing the natural death
  # delta : parameter of lag 
  # pi_QD : proportion of particles that will die very quickly
  # lambda_QD: 1/mean of time of quick death phenomenon 
  # output : vector of length n
  
    Z_QD <- sample(c(1,0),n,replace=TRUE,prob = c(pi_QD, 1-pi_QD)) 
    T_QD <- rgamma(n,shape = theta_QD[1], rate = theta_QD[2])
    T_LD <- rminGammaEMGamma(n,lambda_c,k,lambda_e,theta_ND,delta)
    Y <- Z_QD * T_QD  + (1-Z_QD)*T_LD
    return(list(Y = Y,Z  = Z_QD))
}


dOurModel <- function(x,lambda_c,k,lambda_e,theta_ND,delta,pi_QD,theta_QD){
  
  # lambda_c: Param of the exponential distribution for the time to go to the ADN
  # k:  Param of the Gamma. Number of codon
  # lambda_e:  positive real number. Param of the Gamma. Inverse of the mean length of lecture per Codon (élongation)
  # theta_ND:  vector of 2 positive real numbers. Param of the gamma  distribution representing the natural death
  # delta : parameter of lag 
  # pi_QD : proportion of particles that will die very quickly
  # theta_QD: vector of 2 positive real numbers. Param of the gamma  distribution representing the quick death
  # output : vector of same length as x
  
  d1 <- dgamma(x,shape = theta_QD[1], rate = theta_QD[2])
  #d1 <- dexp(x,theta_QD[2])
  d2 <- dminGammaEMGaussian(x,lambda_c,mu=k/lambda_e + delta,sigma=sqrt(k/lambda_e^2),theta = theta_ND)
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




log_lik_TV_TR <- function(log_param,theta_N,k,kprime,data=list(TV=TV, TR = TR)){
  
  # log_param: log(lambda_c),log(lambda_e)
  
  lambda_c <- exp(log_param[1])
  lambda_e <- exp(log_param[2])
  mu_V <- k/lambda_e
  mu_R <- (k+kprime)/lambda_e
  sigma_V <- sqrt(mu_V/lambda_e)
  sigma_R <- sqrt(mu_R/lambda_e)
 
  LogLTV <- dens_minGammaEMG(data$TV, lambda = lambda_c,mu = mu_V,sigma = sigma_V,theta_N)
  LogLTR <- dens_minGammaEMG(data$TR, lambda = lambda_c,mu = mu_R,sigma = sigma_R,theta_N)
  - sum(log(LogLTV )) - sum(log(LogLTR))
}


