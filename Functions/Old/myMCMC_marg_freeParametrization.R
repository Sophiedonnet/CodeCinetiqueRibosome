my_mcmc_marg_freeParametrization  = function(theData,log_param_init,
                             hyperparams = list(),
                             paramsChains = list(nMCMC=50000,rho=1,nBurnin=1000,paramsToSample=c(1:nbParam))){
 
  
  # param 
  #   1 - "lambda_ND_UP"    
  #   2 - "piTrunc_ND_UP"   
  #   3 - "lambda_ND_DN"
  #   4 - "piTrunc_ND_DN"
  #   5-  "lambda_c"
  #   6 - "mu_e_UP" 
  #   7 - "sigma_e_UP"
  #   8 - "piTrunc_Read_UP"
  #   9 - "mu_e_DN" 
  #   10 - "sigma_e_DN"
  #   11 - "piTrunc_Read_DN"
  
   
  nbParam <- length(log_param_init)

  #---------------------------- Params to sample
  if(is.null(paramsChains$paramsToSample)){paramsChains$paramsToSample <- 1:nbParam}
  whereRW <- whereRWFunction(paramsChains$paramsToSample)
  #----------------------------- Param fr kernel 
  if( is.null(paramsChains$rho) ){paramsChains$rho <- rep(1,nbParam)}
  if( length(paramsChains$rho)==1 ){paramsChains$rho <- rep(paramsChains$rho,nbParam)}
  
  
  #--------------------- INIT and Stockage
  myPostSample <- matrix(0,paramsChains$nMCMC,nbParam) 
  myPostSample[1,] <- log_param <- log_param_init
  logprior <- log_prior_marg(log_param,hyperparams) 
  LL <- log_lik_marg_freeParametrization(log_param,theData)
  
  
  for (iter in 2:paramsChains$nMCMC){
    
    if(iter%%100==0){print(paste0('Iteration ', iter))}
    
    #------------------ param[2] 
    NbTmax_ctr_UP <- sum(theData$Tctr_UP==theData$Tmax_Tctr_UP)
    log_param[2] <- rbeta(1,hyperparams$param1[2]+NbTmax_ctr_UP,hyperparams$param2[2] + length(theData$Tctr_UP)-NbTmax_ctr_UP)
    
    #------------------ param[4] 
    NbTmax_ctr_DN <- sum(theData$Tctr_DN==theData$Tmax_Tctr_DN)
    log_param[4] <- rbeta(1,hyperparams$param1[4]+NbTmax_ctr_DN,hyperparams$param2[4] + length(theData$Tctr_DN)-NbTmax_ctr_DN)
    
    #------------------ param[8] 
    NbTmax_exp_UP <- sum(theData$Texp_UP==theData$Tmax_Texp_UP)
    piTrunc_Exp_UP <- rbeta(1,hyperparams$param1[8]+NbTmax_exp_UP,hyperparams$param2[8] + length(theData$Texp_UP)-NbTmax_exp_UP)
    log_param[8] <- piTrunc_Exp_UP/log_param[2]
    
    
    #------------------ param[11] 
    NbTmax_exp_DN <- sum(theData$Texp_DN==theData$Tmax_Texp_DN)
    piTrunc_Exp_DN <- rbeta(1,hyperparams$param1[11]+NbTmax_exp_DN,hyperparams$param2[11] + length(theData$Texp_DN)-NbTmax_exp_DN)
    log_param[11] <- piTrunc_Exp_DN/log_param[4]
    
    
    
    #-------------------------------------------------------------
    # log(lambda_ND_R, lambda_ND_V, lambda_c,lambda_e,log lambda_QD ,logit
    #-------------------------------------------------------------
    ord_w <- sample(whereRW,length(whereRW) ,replace=FALSE) 
    for (w in ord_w){ 
        log_param_c <- log_param 
        sw <- sample(c(0.01,0.1,1),1)*paramsChains$rho[w]
        log_param_c[w] <- log_param[w] + rnorm(1,0,sw)
        logprior_c <- log_prior_marg(log_param_c,hyperparams)
        LL_c <- log_lik_marg_freeParametrization(log_param_c,theData)
        alpha <- LL_c + logprior_c - LL - logprior
        if (log(runif(1))<alpha){
          log_param <- log_param_c
          LL <- LL_c
          logprior <-  logprior_c 
          #print('accept')
        }# fin accept
    }# fin sampling over parameters
    myPostSample[iter,] = log_param
  }
  
  

  #============================ output
   
  myLogPostSample <- as.data.frame(myPostSample)
  names(myLogPostSample) <- names(log_param_init)
  
  #-------------------------- 
  myPostSample <- myLogPostSample
  myPostSample[,c(1,3,5,7,10)] <- exp(myLogPostSample[,c(1,3,5,7,10)])
  myPostSample[,c(2,4,8,11)] <- myLogPostSample[,c(2,4,8,11)]
  
  return(list(myPostSample  = myPostSample,myLogPostSample  = myLogPostSample))
}


#================================================================
whereRWFunction = function(vect){
  
  myvect <- vect
  w <- which(vect %in% c(2,4,8,11))
  myvect <- myvect[-w]
  
  if (length(myvect)==1){
    myvect = c(myvect,myvect)
  }
  return(myvect)
  
}

#========================================================================
#----------------- Likelihood (TV, TR)
#========================================================================
log_lik_marg_freeParametrization <- function(log_param,theData){
  
  # log_param_init : 
  #   - log(lambda_ND_UP) 
  #   - logit(piTrunc_ND_UP),
  #   - log(lambda_ND_DN) 
  #   - logit(piTrunc_ND_DN), 
  #   - log(lambda_c)
  #   - mu_c_UP 
  #   - log(sigma_c_UP),
  #   - mu_c_DN 
  #   - log(sigma_c_DN),
  
  
  
  param <- from_logparam_to_param_freeParametrization(log_param)
  
  #-----------------------------Data controle  UP ---------
  d1 <- dExpCensored(theData$Tctr_UP, lambda = param[1], piTrunc = param[2],Tmax = theData$Tmax_Tctr_UP, log  = TRUE)
 
  #-----------------------------Data controle DN ---------
  d2 <- dExpCensored(theData$Tctr_DN, lambda = param[3], piTrunc = param[4],Tmax = theData$Tmax_Tctr_DN, log  = TRUE)
  #------------------------------Data Exp UP ---------
  d3 <- dOurModelExp_freeParametrization(theData$Texp_UP,param,UPDN ='UP',Tmax = theData$Tmax_Texp_UP,log = TRUE) 
  #------------------------------Data Exp DN ---------
  d4 <- dOurModelExp_freeParametrization(theData$Texp_DN,param,UPDN='DN',Tmax = theData$Tmax_Texp_DN,log = TRUE) 
 
  res <- sum(d1) + sum(d2) + sum(d3) + sum(d4) 
  return(res)
}


###############################################################" 
#========================================================================
log_prior_marg = function(log_param,hyperparams){
#========================================================================
  
  whichParam <- 1:length(log_param)
  w <- c(2,4,8,11)
  wbarre <- c(1,3,5,6,7,10,11)
  d + sum(dbeta(log_param[w], hyperparams$param1[w], hyperparams$param2[w],log = TRUE))
  d <- d + sum(dnorm(log_param[wbarre], hyperparams$param1[wbarre], hyperparams$param2[wbarre],log = TRUE))
  return(d)
}


