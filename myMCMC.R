my_mcmc_onechain  = function(data,log_param_init,
                             hyperparams = list(lowerbound=-10,upperbound=0,a=1,b=1),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:6),withQD  = FALSE)){

  #param_init : log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,[log lambda_QD,pi_DQ])
  
  #--------------------------------------
  withQD <- paramsChains$withQD
  if(is.null(withQD)){withQD  = FALSE}
  q <- length(log_param_init)
  n_Exp_V <- length(data$T_Exp_V)
  n_Exp_R <- length(data$T_Exp_R)
  n_Exp <- n_Exp_V+n_Exp_R
  #----------- checks
  if(!withQD){log_param_init <- log_param_init[1:4]}
  if(withQD & (q!=6)){stop('With QD requires 6 parameters')}
  
  checkLB <- length(hyperparams$lowerbound)== 4 + withQD
  checkUB <- length(hyperparams$lowerbound)== 4 + withQD
  if((!checkLB) | (!checkUB)){
    stop('Specify upper bounds and lower bounds for all the parameters')}

  #---------------------------- Params to sample
  if(is.null(paramsChains$paramsToSample)){paramsChains$paramsToSample <- 1:q}
  whereRW <- whereToUpdate(paramsChains$paramsToSample)
  
  #--------------------- INIT and Stockage
  myPostSample <- matrix(0,paramsChains$nMCMC-paramsChains$nBurnin,q) 
  log_param <- log_param_init
  
  if(withQD){
    pi_QD <- log_param_init[6]
    Z <- list()
    Z$ZR <- dbinom(n_Exp_R,1,pi_QD)
    Z$ZV <- dbinom(n_Exp_V,1,pi_QD)
  }else{Z = NULL}
  LL <- log_lik(log_param,data=mydata,Z)
  
  logprior <- sum(dunif(log_param[1:(4+withQD)], hyperparams$lowerbound[1:(4+withQD)], hyperparams$upperbound[1:(4+withQD)],log = TRUE))
  if(withQD){
    logprior <- logprior + dbeta(log_param[6],hyperparams$a,hyperparams$b,log = TRUE)
  }
  
  
  for (i in 1:paramsChains$nMCMC){
    
    if(i%%1000==0){print(paste0('Iteration ', i))}
    
    #-------------------------------------------------------------
    # log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD 
    #-------------------------------------------------------------
    log_param_c <- log_param 
    w = sample(whereRW,1) 
    s <- sample(c(0.01,0.1,1),1)*paramsChains$rho
    log_param_c[w] <- log_param[w] + rnorm(1,0,s[w])
    
    
    logprior_c <- sum(dunif(log_param_c[1:(4+withQD)], hyperparams$lowerbound[1:(4+withQD)], hyperparams$upperbound[1:(4+withQD)],log = TRUE))
    if(withQD){
      logprior_c <- logprior_c + dbeta(log_param_c[6],hyperparams$a,hyperparams$b,log = TRUE)
    }
    
    if(logprior_c!=-Inf){
      LL_c <- log_lik(log_param_c,data=mydata,Z)
      alpha <- LL_c + logprior_c - LL - logprior
      if(is.na(LL_c)){browser()}
      #if(LL_c==Inf){LL_c  = -Inf}
      if (log(runif(1))<alpha){
        log_param <- log_param_c
        LL <- LL_c
        logprior <-  logprior_c 
      }
    }
    
    
    #-------------------------------------------------------------
    # Z  and pi_QD
    #-------------------------------------------------------------
    if(withQD){
      
      Z <- sample_Z_QD(log_param,data,pi_QD)
      S <- sum(sapply(Z, sum))
      log_param[6] = rbeta(1,hyperparams$a+S, hyperparams$b+ n_Exp - S)
    }
    
    if(i>paramsChains$nBurnin){myPostSample[i-paramsChains$nBurnin,] = log_param}
  }
  

  #============================ output
  myPostSample <- as.data.frame(myPostSample)
  names_myPostSample <- c('log_lambda_ND_V', 'log_lambda_ND_R', 'log_lambda_c','log_lambda_e')
  if(withQD){names_myPostSample = c(names_myPostSample,'log_lambda_QD','pi_QD')}
  names(myPostSample) <- names_myPostSample

  return(myPostSample)
}

#=================================================================================================
my_mcmc  = function(data,list_log_param_init,
                    hyperparams = list(mu=0,sigma=10,a=1,b=1),
                    paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,nChains= 3)){
  
  
  myListPostSamplelapply(1:paramsChains$nChains,function(m){
    my_mcmc_onechain(data,
                     list_log_param_init[[i]],
                     hyperparams,
                     paramsChains)})
  
  return(myListPostSample)
}

#================================================================
whereToUpdate = function(vect){
  
  if(max(vect)==6){myvect = vect[-which(vect==6)]}else{myvect = vect}
  oneParam <- (length(myvect)==1)
  if (oneParam){
    myvect = c(myvect,myvect)
  }else{
    
  }
  return(myvect)
  
}

#================================================================
sample_Z_QD <- function(log_param,data,pi_QD){
  
  lambda_ND_V <- exp(log_param[1])
  lambda_ND_R <- exp(log_param[2])
  lambda_c <- exp(log_param[3])
  lambda_e <- exp(log_param[4])
  lambda_QD <- exp(log_param[5])
  pi_QD <- log_param[6]

  # the green
  pZV_QD <- pi_QD*dexp(data$T_Exp_V,lambda_QD)
  pZV_NQD <- (1-pi_QD)*dminGammaEMGaussian(data$T_Exp_V,
                                           lambda = lambda_c,
                                           mu = data$k/lambda_e,
                                           sigma=sqrt(data$k)/lambda_e,
                                           theta =c(1,lambda_ND_V) )
  prob_ZV <-  pZV_QD / (pZV_QD+pZV_NQD)
  ZV <- rbinom(length(data$T_Exp_V),1,prob_ZV)
  
  # the red
  pZR_QD <- pi_QD*dexp(data$T_Exp_R,lambda_QD)
  pZR_NQD <- (1-pi_QD)*dminGammaEMGaussian(data$T_Exp_R,
                                           lambda = lambda_c,
                                           mu = (data$k+ data$kprime)/lambda_e,
                                           sigma=sqrt(data$k+data$kprime)/lambda_e,
                                           theta =c(1,lambda_ND_R) )
  prob_ZR <-  pZR_QD / (pZR_QD+pZR_NQD)
  ZR <- rbinom(length(data$T_Exp_R),1,prob_ZR)
  return(list(ZV = ZV,ZR = ZR))
}
