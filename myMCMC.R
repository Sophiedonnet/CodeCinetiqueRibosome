my_mcmc_onechain  = function(data,log_param_init,
                             hyperparams = list(lowerbound=-10,upperbound=0,a=1,b=1),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:6),withQD  = FALSE)){

  #param_init : log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,[log lambda_QD,pi_DQ])
  
  if(is.null(data$Tmax)){data$Tmax = Inf}
   
  #--------------------------------------
  withQD <- paramsChains$withQD
  if(is.null(withQD)){withQD  = FALSE}
  n_Exp_V <- length(data$T_Exp_V)
  n_Exp_R <- length(data$T_Exp_R)
  n_Exp <- n_Exp_V+n_Exp_R
  if(n_Exp==0){withQD = FALSE; withExp = FALSE}
  #----------- checks
  if(!withQD){
    log_param_init[6] <- 0
    log_param_init[5] <- 0
    hyperparams$lowerbound[5] <- -10
    hyperparams$upperbound[5] <- 0
    hyperparams$a = 1
    hyperparams$b = 1
  }

  #---------------------------- Params to sample
  if(is.null(paramsChains$paramsToSample)){paramsChains$paramsToSample <- 1:ifelse(withQD,6,4)}
  whereRW <- whereToUpdate(paramsChains$paramsToSample)
  
  #--------------------- INIT and Stockage
  myPostSample <- matrix(0,paramsChains$nMCMC-paramsChains$nBurnin,6) 
  log_param <- log_param_init
  Z <- list()
  if(withQD){
    Z$ZR <- rbinom(n_Exp_R,1,log_param[6])
    Z$ZV <- rbinom(n_Exp_V,1,log_param[6])
  }else{
    Z$ZR = rep(0,n_Exp_R)
    Z$ZV = rep(0,n_Exp_V)
  }
  LL <- log_lik(log_param,data,Z,withQD)
  logprior <- log_prior(log_param,hyperparams,withQD) 
  
  
  for (i in 1:paramsChains$nMCMC){
    
    if(i%%1000==0){print(paste0('Iteration ', i))}
    
    #-------------------------------------------------------------
    # log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD 
    #-------------------------------------------------------------
    if(length(whereRW>0)){
      log_param_c <- log_param 
      w = sample(whereRW,1) 
      s <- sample(c(0.01,0.1,1),1)*paramsChains$rho
      log_param_c[w] <- log_param[w] + rnorm(1,0,s[w])
      logprior_c <- log_prior(log_param_c,hyperparams,withQD)
      if(logprior_c!=-Inf){
        LL_c <- log_lik(log_param_c,data=mydata,Z,withQD)
        alpha <- LL_c + logprior_c - LL - logprior
        if(is.na(LL_c)){browser()}
      #if(LL_c==Inf){LL_c  = -Inf}
        if (log(runif(1))<alpha){
          if (w==5){print(c(i,'accept'))}
          log_param <- log_param_c
          LL <- LL_c
          logprior <-  logprior_c 
        }
      }
    }
    
    
    #-------------------------------------------------------------
    # Z  and pi_QD
    #-------------------------------------------------------------
    if(withQD){
      Z <- sample_Z_QD(log_param,data)
      S <- sum(sapply(Z, sum))
      if(6 %in% paramsChains$paramsToSample){
        log_param[6] = rbeta(1,hyperparams$a+S, hyperparams$b+ n_Exp - S)
      }
    }
    if(i>paramsChains$nBurnin){myPostSample[i-paramsChains$nBurnin,] = log_param}
  }
  

  #============================ output
   
  myLogPostSample <- as.data.frame(myPostSample)
  names_myLogPostSample <- c('log_lambda_ND_V', 'log_lambda_ND_R', 'log_lambda_c','log_lambda_e','log_lambda_QD','pi_QD')
  names(myLogPostSample) <- names_myLogPostSample
  
  myPostSample <- myLogPostSample
  myPostSample[,1:5] <- exp(myLogPostSample[,1:5])
  names_myPostSample <- c('lambda_ND_V', 'lambda_ND_R', 'lambda_c','lambda_e','lambda_QD','pi_QD')
  names(myPostSample) <- names_myPostSample
  

  return(list(myPostSample  = myPostSample,myLogPostSample  = myLogPostSample))
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
sample_Z_QD <- function(log_param,data){
  
  n_Exp_V <- data$n_Exp_V
  n_Exp_R <- data$n_Exp_R
  
  pi_QD <- log_param[6]
  if(pi_QD==1){
    Z <- list(ZV =rep(1,n_Exp_V) ,ZR = rep(1,n_Exp_R))
  }
  if(pi_QD==0){
    Z <- list(ZV =rep(0,n_Exp_V) ,ZR = rep(0,n_Exp_R))
  }
  if ((pi_QD< 1) & (pi_QD>0)){
    lambda_ND_V <- exp(log_param[1])
    lambda_ND_R <- exp(log_param[2])
    lambda_c <- exp(log_param[3])
    lambda_e <- exp(log_param[4])
    lambda_QD <- exp(log_param[5])
    # the green
    pZV_QD <- pi_QD*dexp(data$T_Exp_V,lambda_QD)
    pZV_NQD <- (1-pi_QD)*dminGammaEMGaussian(data$T_Exp_V,
                                           lambda = lambda_c,
                                           mu = data$k/lambda_e,
                                           sigma=sqrt(data$k)/lambda_e,
                                           theta =c(1,lambda_ND_V) )
    prob_ZV <-  pZV_QD / (pZV_QD+pZV_NQD)
    ZV <- rbinom(n_Exp_V,1,prob_ZV)
    # the red
    pZR_QD <- pi_QD*dexp(data$T_Exp_R,lambda_QD)
    pZR_NQD <- (1-pi_QD)*dminGammaEMGaussian(data$T_Exp_R,
                                           lambda = lambda_c,
                                           mu = (data$k+ data$kprime)/lambda_e,
                                           sigma=sqrt(data$k+data$kprime)/lambda_e,
                                           theta =c(1,lambda_ND_R) )
    prob_ZR <-  pZR_QD / (pZR_QD+pZR_NQD)
    ZR <- rbinom(length(data$T_Exp_R),1,prob_ZR)
    Z <- list(ZV = ZV, ZR= ZR)
  }
    return(Z)
}
