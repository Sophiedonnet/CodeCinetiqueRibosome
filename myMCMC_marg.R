
my_mcmc_marg_onechain  = function(data,log_param_init,
                             hyperparams = list(),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:9))){

  #param_init : log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,[log lambda_QD,pi_DQ])
  

   
  
  
   

  #---------------------------- Params to sample
  if(is.null(paramsChains$paramsToSample)){paramsChains$paramsToSample <- 1:9}
  whereRW <- whereToUpdate_marg(paramsChains$paramsToSample)
  
  #--------------------- INIT and Stockage
  myPostSample <- matrix(0,paramsChains$nMCMC-paramsChains$nBurnin,9) 
  log_param <- log_param_init
   
  LL <- log_lik_marg(log_param,data)
  logprior <- log_prior_marg(log_param,hyperparams,whereRW) 
  
  
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
      logprior_c <- log_prior_marg(log_param_c,hyperparams,whereRW)
      #if(logprior_c!=-Inf){
        LL_c <- log_lik_marg(log_param_c,data)
        alpha <- LL_c + logprior_c - LL - logprior
        if (log(runif(1))<alpha){
          log_param <- log_param_c
          LL <- LL_c
          logprior <-  logprior_c 
        }
      #}
    }
    
    
    if(i>paramsChains$nBurnin){myPostSample[i-paramsChains$nBurnin,] = log_param}
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
log_lik_marg <- function(log_param,data){
  
  # log_param: log(lambda_ND_V),log(lambda_ND_R), log(lambda_c),log(lambda_e) [log(lambda_QD) pi_QD]
  

  param <- from_logparam_to_param(log_param)
  
  d = 0
  
  
  #-----------------------------Data controle Vertes---------
  d1 <- dExpCensored(data$T_Contr_V, lambda = param[1], Tmax = data$Tmax_Contr_V,piTrunc = param[2],log  = TRUE)
  d  <- sum(d1)
  
  #-----------------------------Data controle Rouges---------
  d2 <- dExpCensored(data$T_Contr_R, lambda = param[3],Tmax = data$Tmax_Contr_R,piTrunc = param[4],log  = TRUE)
  d  <- d + sum(d2)
  #------------------------------Data Exp Green---------
  d3 <- dOurModelExp(data$T_Exp_V,param,data$k,data$kprime,color='green',Tmax = data$Tmax_Exp_V,log = TRUE) 
  d <- d  + sum(d3)
  #------------------------------Data Exp Red ---------
  d4 <- dOurModelExp(data$T_Exp_R,param,data$k,data$kprime,color='red',Tmax = data$Tmax_Exp_R,log = TRUE) 
  d <- d  + sum(d4)
  
  return(d)
}

#========================================================================
log_prior_marg = function(log_param,hyperparams,w){
  #========================================================================
 
  if('mean' %in% names(hyperparams)){
    d <- sum(dnorm(log_param[w] , hyperparams$mean[w], hyperparams$sd[w],log = TRUE))
  }
  if('upperbound' %in% names(hyperparams)){
    d <- sum(dunif(log_param[w], hyperparams$lowerbound[w], hyperparams$upperbound[w],log = TRUE))
  }
  return(d)
}



