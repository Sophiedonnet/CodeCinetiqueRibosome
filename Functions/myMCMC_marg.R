my_mcmc_marg  = function(data,log_param_init,
                             hyperparams = list(),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:9))){
 
  
  # log_param_init : log(lambda_ND_UP) ,logit(piTrunc_ND_UP),log(lambda_ND_DN),logit(piTrunc_ND_DN), log(lambda_c), log(lambda_e), logit(piTrunc_Read), log(lambda_QD),logit(pi_QD))
  
   

  #---------------------------- Params to sample
  if(is.null(paramsChains$paramsToSample)){paramsChains$paramsToSample <- 1:9}
  whereRW <- whereToUpdate_marg(paramsChains$paramsToSample)
  #----------------------------- Param fr kernel 
  if( is.null(paramsChains$rho) ){paramsChains$rho <- rep(1,9)}
  if( length(paramsChains$rho)==1 ){paramsChains$rho <- rep(paramsChains$rho,9)}
  
  
  #--------------------- INIT and Stockage
  myPostSample <- matrix(0,paramsChains$nMCMC-paramsChains$nBurnin,9) 
  log_param <- log_param_init
  logprior <- log_prior_marg(log_param,hyperparams,whereRW) 
  LL <- log_lik_marg(log_param,data)
  
  
  for (iter in 1:paramsChains$nMCMC){
    
    if(iter%%1000==0){print(paste0('Iteration ', iter))}
    
    #-------------------------------------------------------------
    # log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD ,logit
    #-------------------------------------------------------------
    if(length(whereRW)>0){
      log_param_c <- log_param 
      w <- sample(whereRW,1) 
      sw <- sample(c(0.01,0.1,1),1)*paramsChains$rho[w]
      log_param_c[w] <- log_param[w] + rnorm(1,0,sw)
      logprior_c <- log_prior_marg(log_param_c,hyperparams,whereRW)
      if(logprior_c!=-Inf){
        LL_c <- log_lik_marg(log_param_c,data)
        alpha <- LL_c + logprior_c - LL - logprior
        if (log(runif(1))<alpha){
          log_param <- log_param_c
          LL <- LL_c
          logprior <-  logprior_c 
          #print('accept')
        }
      }
    }
    
    
    if(iter>paramsChains$nBurnin){myPostSample[iter-paramsChains$nBurnin,] = log_param}
  }
  

  #============================ output
   
  myLogPostSample <- as.data.frame(myPostSample)
  names(myLogPostSample) <- names(log_param_init)
  
  myPostSample <- myLogPostSample
  myPostSample[,c(1,3,5,6,8)] <- exp(myLogPostSample[,c(1,3,5,6,8)])
  myPostSample[,c(2,4,7,9)] <- invlogit(myLogPostSample[,c(2,4,7,9)])
  
  return(list(myPostSample  = myPostSample,myLogPostSample  = myLogPostSample))
}


#================================================================
whereToUpdate_marg = function(vect){
  
  myvect <- vect
  if (length(myvect)==1){
    myvect = c(myvect,myvect)
  }
  return(myvect)
  
}

#========================================================================
#----------------- Likelihood (TV, TR)
#========================================================================
log_lik_marg <- function(log_param,data){
  
  # log_param : log(lambda_ND_UP) ,logit(piTrunc_ND_UP),log(lambda_ND_DN),logit(piTrunc_ND_DN), log(lambda_c), log(lambda_e), logit(piTrunc_Read), log(lambda_QD),logit(pi_QD))

  param <- from_logparam_to_param(log_param)
  
  #-----------------------------Data controle  UP ---------
  d1 <- dExpCensored(data$Tctr_UP, lambda = param[1], Tmax = data$Tmax_Tctr_UP,piTrunc = param[2],log  = TRUE)
 
  #-----------------------------Data controle DN ---------
  d2 <- dExpCensored(data$Tctr_DN, lambda = param[3],Tmax = data$Tmax_Tctr_DN,piTrunc = param[4],log  = TRUE)
  #------------------------------Data Exp UP ---------
  d3 <- dOurModelExp(data$Texp_UP,param,data$k,data$kprime,UPDN ='UP',Tmax = data$Tmax_Texp_UP,log = TRUE) 
  #------------------------------Data Exp DN ---------
  d4 <- dOurModelExp(data$Texp_DN,param,data$k,data$kprime,UPDN='DN',Tmax = data$Tmax_Texp_DN,log = TRUE) 
 
  res <- sum(d1) + sum(d2) + sum(d3) + sum(d4) 
  return(res)
}

#========================================================================
log_prior_marg = function(log_param,hyperparams,whereRW){
  #========================================================================
 
  if('mean' %in% names(hyperparams)){
    d <- sum(dnorm(log_param[whereRW] , hyperparams$mean[whereRW], hyperparams$sd[whereRW],log = TRUE))
  }
  if('upperbound' %in% names(hyperparams)){
    d <- sum(dunif(log_param[whereRW], hyperparams$lowerbound[whereRW], hyperparams$upperbound[whereRW],log = TRUE))
  }
  return(d)
}



