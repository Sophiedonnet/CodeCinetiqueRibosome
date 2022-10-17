my_mcmc_onechain  = function(data,log_param_init,
                             hyperparams = list(lowerbound=-10,upperbound=0,a=1,b=1),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:6),withQD  = FALSE)){
  
  
  withQD <- paramsChains$withQD
  if(is.null(withQD)){withQD  = FALSE}
  
  
  if(is.null(paramsChains$paramsToSample)){
    if(withQD){
      paramsChains$paramsToSample = 1:5
    }else{
      paramsChains$paramsToSample = 1:4}
  }
  
  #param_init : log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD,pi_DQ
  myPostSample <- matrix(0,paramsChains$nMCMC-paramsChains$nBurnin,length(log_param_init)) 
  q <- length(log_param_init)
  log_param <- log_param_init
  
  if(withQD){
    pi_QD <- log_param_init[6]
    ZR <- dbinom(n_Exp_R,1,pi_QD)
    ZV <- dbinom(n_Exp_V,1,pi_QD)
    LL <- log_lik_withQD(log_param,data=mydata,ZV,ZR)
  }else{
    LL <- log_lik_withoutQD(log_param,data=mydata)  
  }
  
  logprior <- sum(dunif(log_param[1:4], hyperparams$lowerbound[1:4], hyperparams$upperbound[1:4],log = TRUE))
  if(withQD){
    logprior <- logprior + dunif(log_param[5], hyperparams$lowerbound[5], hyperparams$upperbound[5],log = TRUE) 
    logprior <- logprior + dbeta(log_param[6],a,b,log = TRUE)
  }
  
  
  
  for (i in 1:paramsChains$nMCMC){
    
    if(i%%1000==0){print(paste0('Iteration ', i))}
    
    #-------------------------------------------------------------
    # log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD 
    #-------------------------------------------------------------
    log_param_c <- log_param 
    s <- sample(c(0.01,0.1,1),1)*paramsChains$rho
    if(length(paramsChains$paramsToSample)==1){
      w = paramsChains$paramsToSample}else{
        w <- sample(paramsChains$paramsToSample,1)
      }
    log_param_c[w] <- log_param[w] + rnorm(1,0,s)
    
    
    logprior_c <- sum(dunif(log_param_c[1:4], hyperparams$lowerbound[1:4], hyperparams$upperbound[1:4],log = TRUE))
    if(withQD){
      logprior_c <- logprior + dunif(log_param_c[5], hyperparams$lowerbound[5], hyperparams$upperbound[5],log = TRUE) 
    }
    
    if(logprior_c!=-Inf){
      
      if(withQD){
        LL <- log_lik_withQD(log_param_c,data=mydata,ZV,ZR)
      }else{
        LL_c <- log_lik_withoutQD(log_param_c,data=mydata)
      }
      alpha <- LL_c + logprior_c - LL - logprior
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
      lambda_QD <- exp(log_param[5])
      lambda_c <- exp(log_param[3])
      lambda_e <- exp(log_param[3])
      lambda_ND_R <- exp(log_param[2])
      lambda_ND_V <- exp(log_param[1])
      pi_QD <- log_param[6]
      
      # the green
      pZV_QD <- pi_QD*dexp(data$T_Exp_V,lambda_QD)
      pZV_NQD <- (1-pi_QD)*dminGammaEMGaussian(data$T_Exp_V,lambda = lambda_c,mu = data$k/lambda_e,sigma=data$k/lambda_e,theta =c(1,lambda_ND_V) )
      prob_ZV <-  pZV_QD / (pZV_QD+pZV_NQD)
      ZV <- rbinom(n_Exp_V,1,prob_ZV)
      
      # the red
      pZR_QD <- pi_QD*dexp(data$T_Exp_R,lambda_QD)
      pZR_NQD <- (1-pi_QD)*dminGammaEMGaussian(data$T_Exp_R,lambda = lambda_c,mu = (data$k+ data$kprime)/lambda_e,sigma=(data$k+data$kprime)/lambda_e,theta =c(1,lambda_ND_R) )
      prob_ZR <-  pZR_QD / (pZR_QD+pZR_NQD)
      ZR <- rbinom(n_Exp_R,1,prob_ZR)
      
      SR <- sum(ZR)
      SV <- sum(ZV)
      log_param[6] = rbeta(1,a+SR+SV, b+ n_Exp_V+n_Exp_R-SR-SV)
    }
    if(i>paramsChains$nBurnin){
      myPostSample[i-paramsChains$nBurnin,] = log_param
    }
  }
  
  myPostSample <- as.data.frame(myPostSample)
  if(withQD){
    names(myPostSample) <- c('log_lambda_ND_V', 'log_lambda_ND_R', 'log_lambda_c','log_lambda_e','log_lambda_QD','pi_QD')
  }else{
    names(myPostSample) <- c('log_lambda_ND_V', 'log_lambda_ND_R', 'log_lambda_c','log_lambda_e')
  }
  
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
