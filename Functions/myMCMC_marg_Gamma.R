my_mcmc_marg_Gamma_onechain  = function(data,log_param_init,
                             hyperparams = list(),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:8), withQD  = FALSE)){

  #param_init : log(alpha_ND_V, beta_ND_V, alpha_ND_R,beta_ND_R, lambda_c,lambda_e,[log lambda_QD,pi_DQ])
  
  if(is.null(data$Tmax_R)){data$Tmax_R = Inf}
  if(is.null(data$Tmax_V)){data$Tmax_V = Inf}
  
  
   
  #--------------------------------------
  withQD <- paramsChains$withQD
  if(is.null(withQD)){stop('With Quick Death or not?')}
  #----------- checks
  if(!withQD){# default values 
    log_param_init[8] <- 0
    log_param_init[7] <- 0
  }
  #---------------------------- Params to sample
  if(is.null(paramsChains$paramsToSample)){paramsChains$paramsToSample <- 1:ifelse(withQD,8,6)}
  whereRW <- whereToUpdate_marg(paramsChains$paramsToSample)
  
  #--------------------- INIT and Stockage
  myPostSample <- matrix(0,paramsChains$nMCMC-paramsChains$nBurnin,8) 
  log_param <- log_param_init
   
  LL <- log_lik_marg_Gamma(log_param,data)
  logprior <- log_prior_marg_Gamma(log_param,hyperparams,withQD) 
  
  
  for (i in 1:paramsChains$nMCMC){
    
    if(i%%1000==0){print(paste0('Iteration ', i))}
    
    #-------------------------------------------------------------
    # log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD ,logit
    #-------------------------------------------------------------
    if(length(whereRW)>0){
      log_param_c <- log_param 
      w = sample(whereRW,1) 
      s <- sample(c(0.01,0.1,1),1)*paramsChains$rho
      log_param_c[w] <- log_param[w] + rnorm(1,0,s[w])
      logprior_c <- log_prior_marg_Gamma(log_param_c,hyperparams,withQD)
      if(logprior_c!=-Inf){
        LL_c <- log_lik_marg_Gamma(log_param_c,data)
        alpha <- LL_c + logprior_c - LL - logprior
        if (log(runif(1))<alpha){
          log_param <- log_param_c
          LL <- LL_c
          logprior <-  logprior_c 
        }
      }
    }
    
    
    if(i>paramsChains$nBurnin){myPostSample[i-paramsChains$nBurnin,] = log_param}
  }
  

  #============================ output
   
  myLogPostSample <- as.data.frame(myPostSample)
  names_myLogPostSample <- c('log_alpha_ND_V', 'log_beta_ND_V', 'log_alpha_ND_R', 'log_beta_ND_R', 'log_lambda_c','log_lambda_e','log_lambda_QD','inv_logit_pi_QD')
  names(myLogPostSample) <- names_myLogPostSample
  
  myPostSample <- myLogPostSample
  myPostSample[,1:7] <- exp(myLogPostSample[,1:7])
  myPostSample[,8] <- 1/(1+exp(-myLogPostSample[,8]))
  
  names_myPostSample <- c('alpha_ND_V', 'beta_ND_V', 'alpha_ND_R', 'beta_ND_R', 'lambda_c','lambda_e','lambda_QD','pi_QD')
  names(myPostSample) <- names_myPostSample
  

  return(list(myPostSample  = myPostSample,myLogPostSample  = myLogPostSample))
}


#================================================================
whereToUpdate_marg = function(vect){
  
  myvect = vect
  oneParam <- (length(myvect)==1)
  if (oneParam){
    myvect = c(myvect,myvect)
  }else{
    
  }
  return(myvect)
  
}

#========================================================================
#----------------- Likelihood (TV, TR)
#========================================================================
log_lik_marg_Gamma <- function(log_param,data){
  
  #log_param: log(alpha_ND_V),log(beta_ND_V), log(alpha_ND_R),log(beta_ND_R), log(lambda_c),log(lambda_e) [log(lambda_QD) pi_QD]
  
  Tmax_R <- data$Tmax_R
  Tmax_V <- data$Tmax_V
  k <- data$k
  kprime <- data$kprime
  
  
  d = 0
  #--------------------------------------
  alpha_ND_V <- exp(log_param[1])
  beta_ND_V <- exp(log_param[2])
  alpha_ND_R <- exp(log_param[3])
  beta_ND_R <- exp(log_param[4])
  lambda_c <- exp(log_param[5])
  lambda_e <- exp(log_param[6])      
  lambda_QD <- exp(log_param[7])
  pi_QD <- 1/(1+ exp(-log_param[8]))
  #--------------------------------------
  d1 <- dGammaCensored(data$T_Contr_V, alpha_ND_V,beta_ND_V,Tmax_V,log  = TRUE)
  d  <- sum(d1)
  
  #--------------------------------------
  d2 <- dGammaCensored(data$T_Contr_R, alpha_ND_R,beta_ND_R,Tmax_R,log  = TRUE)
  d  <- d + sum(d2)
  #--------------------------------------
  mu_V <- k/lambda_e
  sigma_V <- sqrt(k)/lambda_e
  U3_NQD <- dminGammaExpplusGaussian(data$T_Exp_V, lambda = lambda_c,mu = mu_V,sigma = sigma_V,alpha_ND_V,beta_ND_V,Tmax_V)
  U3_QD <- dExpCensored(data$T_Exp_V,lambda_QD,Tmax_V)
  U3 <- pi_QD*U3_QD + (1-pi_QD)*U3_NQD
  if(length(U3)==1){
    if(U3==0){
      U3 = 1
    }
  }
  
  
  d <- d  + sum(log(U3))
  #--------------------------------------

  mu_R <- (k+kprime)/lambda_e
  sigma_R <- sqrt(k+kprime)/lambda_e
  U4_NQD <-  dminGammaExpplusGaussian(data$T_Exp_R, lambda = lambda_c,mu = mu_R,sigma = sigma_R,alpha_ND_R,beta_ND_R,Tmax_R)
  U4_QD <- dExpCensored(data$T_Exp_R,lambda_QD,Tmax_R)
  U4 <- pi_QD*U4_QD + (1-pi_QD)*U4_NQD
  if(length(U4)==1){
    if(U4==0){
      U4 = 1
    }
    }
  

  
  
  
  return(d)
}

#========================================================================
log_prior_marg_Gamma = function(log_param,hyperparams,withQD){
  #========================================================================
  r <- 6+2*withQD
  if('mean' %in% names(hyperparams)){
    r <- 6+2*withQD
    d <- sum(dnorm(log_param[1:r], hyperparams$mean[1:r], hyperparams$sd[1:r],log = TRUE))
  }
  if('upperbound' %in% names(hyperparams)){
    d <- sum(dunif(log_param[1:r], hyperparams$lowerbound[1:r], hyperparams$upperbound[1:r],log = TRUE))
  }
  return(d)
}



