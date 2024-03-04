

#=============================================================
from_param_to_log_param_freeParametrization <- function(param){
#=============================================================  
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
  
  
  log_param <- param
  if(is.vector(param)){
    log_param[c(1,3,5,7,10)] <- log(param[c(1,3,5,7,10)])
      }
  if(is.matrix(param)){
    log_param[,c(1,3,5,7,10)] <- log(param[,c(1,3,5,7,10)])
      }
  names(log_param) <- paste0(c('log_','','log_','','log_','','log_','', '','log_',''),names(param))
  return(log_param)
}
#=============
from_logparam_to_param_freeParametrization <- function(log_param,constMu = TRUE){
  
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
  #   12 - "lambda_e if constMu
  
  
  
  if(is.vector(log_param)){
    param <- log_param
    param[c(1,3,5,7,10)] <- exp(log_param[c(1,3,5,7,10)])
  }
  if(is.matrix(log_param)){
    param <- log_param
    param[,c(1,3,5,7,10)] <- exp(log_param[,c(1,3,5,7,10)])
  }
  names(param) <- c('lambda_ND_UP','piTrunc_ND_UP','lambda_ND_DN','piTrunc_ND_DN','lambda_c','mu_e_UP', 'sigma_e_UP','piTrunc_Read_UP','mu_e_DN', 'sigma_e_DN','piTrunc_Read_DN')
  
  return(param)
}
#=================================================================================




#################################################################### 
dOurModelExp_freeParametrization <- function(x,param,UPDN='DN',Tmax = Inf,log = FALSE){
  ####################################################################  

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
  
  # Tmax : truncature of data
  # output : vector of same length as x
  
  
  lambda_ND <- ifelse(UPDN =='UP',param[1],param[3])
  piTrunc_ND <- ifelse(UPDN =='UP',param[2],param[4])
  mu_e <- ifelse(UPDN =='UP',param[6],param[9])
  sigma_e <- ifelse(UPDN =='UP',param[7],param[10])
  piTrunc_Read <- ifelse(UPDN =='UP',param[8],param[11])
  lambda_c <- param[5]
  
  f <- dminExpExpplusGaussian(x,mu = mu_e,sigma = sigma_e,lambda = lambda_c,piTrunc = piTrunc_Read,lambda_ND = lambda_ND,piTrunc_ND = piTrunc_ND,Tmax,log) 
  
  
  return(f)
} 
#################################################################### 
rOurModelExp_freeParametrization <- function(n,param,UPDN='DN',Tmax = Inf){
  ####################################################################  
  
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
  
  # Tmax : truncature of data
  # output : vector of same length as x


  
  lambda_ND <- ifelse(UPDN =='UP',param[1],param[3])
  piTrunc_ND <- ifelse(UPDN =='UP',param[2],param[4])
  mu_e <- ifelse(UPDN =='UP',param[6],param[9])
  sigma_e <- ifelse(UPDN =='UP',param[7],param[10])
  piTrunc_Read <- ifelse(UPDN =='UP',param[8],param[11])
  lambda_c <- param[5]
  
  E <- rminExpExpplusGaussian(n,mu = mu_e,sigma = sigma_e,lambda = lambda_c,piTrunc = piTrunc_Read,lambda_ND = lambda_ND,piTrunc_ND = piTrunc_ND,Tmax) 
  return(E)
}


