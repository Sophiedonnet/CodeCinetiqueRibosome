rm(list=ls())
library(e1071)
##################################################
################"" Lois exponentielles
###################################################

n <- n_C  <-1000

alpha_N <- 100
beta_N  <- 2
theta_N <-c(alpha_N,beta_N) 
lambda_e <- 1/0.5
lambda_c <- 1/15
k <- 16
kprime <- 30
param <- c(lambda_c,lambda_e)

###données de contrôle
TC =   rgamma(n_C,alpha_N,beta_N)
hist(TC)


###données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e)
UV  <-  rgamma(n,alpha_N,beta_N)
VV <-   rexp(n,lambda_c) + rgamma(n,k,lambda_e)
TV <- apply(cbind(UV,VV),1,min)

UR  <- rgamma(n,alpha_N,beta_N)
VR <- VV  + rgamma(n,kprime,lambda_e)
TR <-  apply(cbind(UR,VR),1,min)


#####  verif distributions 

#------- TV 

dens_TV <- function(x,lambda_c,lambda_e,k,theta_N){
  
  R <- sum(c(1,k)/c(lambda_c,lambda_e))
  Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
  k_V<- R^2/Q
  lambda_V<- k_V/R 
  theta_V <- c(k_V,lambda_V)
  F1z <- pgamma(x,theta_V[1],theta_V[2])
  f1z <- dgamma(x,theta_V[1],theta_V[2])
  F2z <- pgamma(x,theta_N[1],theta_N[2])
  f2z <- dgamma(x,theta_N[1],theta_N[2])
  f1z*(1-F2z) + f2z*(1-F1z)
}

hist(TV,nclass = n/10,freq = FALSE)
curve(dens_TV(x,lambda_c,lambda_e,k,theta_N),add=TRUE,col='red')

prob_TV <- function(x,lambda_c,lambda_e,k,theta_N){
  R <- sum(c(1,k)/c(lambda_c,lambda_e))
  Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
  k_V<- R^2/Q
  lambda_V<- k_V/R 
  theta_V <- c(k_V,lambda_V)
  
  F1z <- pgamma(x,theta_V[1],theta_V[2])
  F2z <- pgamma(x,theta_N[1],theta_N[2])
  1-(1-F1z)*(1-F2z)
  
}

plot(ecdf(TV))
curve(prob_TV(x,lambda_c,lambda_e,k,theta_N),add=TRUE,col='red')


#------- TR 

dens_VR <- function(x,lambda_c,lambda_e,k,kprime){
  R <- sum(c(1,k+kprime)/c(lambda_c,lambda_e))
  Q <- sum(c(1,k+kprime)/c(lambda_c,lambda_e)^2)
  k_R<- R^2/Q
  lambda_R<- k_R/R 
  theta_R <- c(k_R,lambda_R)
  dgamma(x,k_R,lambda_R)
}

plot(density(VR))
curve(dens_VR(x,lambda_c,lambda_e,k,kprime),add=TRUE,col='red')


prob_VR <- function(x,lambda_c,lambda_e,k,kprime){
  R <- sum(c(1,k+kprime)/c(lambda_c,lambda_e))
  Q <- sum(c(1,k+kprime)/c(lambda_c,lambda_e)^2)
  k_R<- R^2/Q
  lambda_R<- k_R/R 
  theta_R <- c(k_R,lambda_R)
  pgamma(x,k_R,lambda_R)
}

plot(ecdf(VR))
curve(prob_VR(x,lambda_c,lambda_e,k,kprime),add=TRUE,col='red')


prob_TR <- function(x,lambda_c,lambda_e,k,kprime,theta_N){
  1-(1-prob_VR(x,lambda_c,lambda_e,k,kprime))*(1-pgamma(x,theta_N[1],theta_N[2]))
}

plot(ecdf(TR))
curve(prob_TR(x,lambda_c,lambda_e,k,kprime,theta_N),add=TRUE,col='red')



dens_TR <- function(x,lambda_c,lambda_e,k,kprime,theta_N){
  
  F_VR <- prob_VR(x,lambda_c,lambda_e,k,kprime)
  f_VR <- dens_VR(x,lambda_c,lambda_e,k,kprime)
  
  F_UR <- pgamma(x,theta_N[1],theta_N[2])
  f_UR <- dgamma(x,theta_N[1],theta_N[2])
  f_VR*(1-F_UR) + f_UR*(1-F_VR)
  
}


hist(TR,freq=FALSE,nclass=n/10)
lines(density(TR),col='orange')
curve(dens_TR(x,lambda_c,lambda_e,k,kprime,theta_N),add=TRUE,col='red')





################################################""
#------------ Estimations
#################################################


estim_param_Gamma <- function(X){
  EX <- mean(X)
  VX <-var(X)
  alpha_hat <- EX^2/VX
  beta_hat <- alpha_hat/EX
  return(c(alpha_hat,beta_hat))
}


estim_param_SumGamma <- function(X,k,kprime){
  
  param_gamma <- estim_param_Gamma(X)
  A <- k^2+k
  R <-param_gamma[1]/param_gamma[2]
  B <- -2*R*k
  C <- R^2-R/param_gamma[2]
  DEL <- B^2-4*A*C
  SOL <- 1/(2*A)*(-B+sqrt(DEL)*c(1,-1))
  lambda_e_hat <- 1/SOL
  lambda_c_hat <- 1/(R-k/lambda_e_hat)
  Theta_hat <- rbind(lambda_c_hat,lambda_e_hat)
  w <- which(apply(Theta_hat>0,2,sum)==2)
  return(Theta_hat[,w])
}

#-----------------------------------  ESTIM theta_N  = (alpha_N, beta_N)---- 
theta_N_hat<- estim_param_Gamma(TC)
cbind(theta_N,theta_N_hat)
#------------------------------------- ESTIM l,lambda_e, ----- 


log_lik_TV <- function(log_param,theta_N,k,TV){
  lambda_c <- exp(log_param[1])
  lambda_e <- exp(log_param[2])
  UTV <- dens_TV(TV, lambda_c,lambda_e,k,theta_N)
  - sum(log(UTV)) 
}


log_lik_TR <- function(log_param,theta_N,k,kprime,TR){
  lambda_c <- exp(log_param[1])
  lambda_e <- exp(log_param[2])
  UTR <- dens_TR(TR, lambda_c,lambda_e,k,kprime,theta_N)
  - sum(log(UTR)) 
}


log_lik_TV_TR <- function(log_param,theta_N,k,kprime,TV,TR){
  log_lik_TR(log_param,theta_N,k,kprime,TR)+ log_lik_TV(log_param,theta_N,k,TV)
}



param_init_from_TV <- estim_param_SumGamma(TV,k,kprime=0)
param_init_from_TR <- estim_param_SumGamma(TR,k,kprime)




log_param_hat_from_TV <- optim(log(param_init_from_TV), log_lik_TV,theta_N = theta_N_hat,k=k,TV = TV)$par 
log_param_hat_from_TR <- optim(log(param_init_from_TR), log_lik_TR,theta_N = theta_N_hat,k=k,kprime = kprime,TR = TR)$par 
log_param_hat_from_TV_TR <- optim(log(param_init_from_TV), log_lik_TV_TR,theta_N = theta_N_hat,k=k,kprime = kprime,TV = TV,TR= TR)$par 

param_hat_from_TV <- exp(log_param_hat_from_TV)
param_hat_from_TR <- exp(log_param_hat_from_TR)
param_hat_from_TV_TR <- exp(log_param_hat_from_TV_TR)

hist(TV,freq=FALSE,nclass=min(n/10,100))
lines(density(TV),col='orange')
curve(dens_TV(x,lambda_c,lambda_e,k,theta_N),add=TRUE,col='red')
#curve(dens_TV(x,lambda_c,lambda_e,k,theta_N_hat),add=TRUE,col='blue')
curve(dens_TV(x,param_hat_from_TV[1],param_hat_from_TV[2],k,theta_N_hat),add=TRUE,col='magenta')
curve(dens_TV(x,param_hat_from_TR[1],param_hat_from_TR[2],k,theta_N_hat),add=TRUE,col='green')
curve(dens_TV(x,param_hat_from_TV_TR[1],param_hat_from_TV_TR[2],k,theta_N_hat),add=TRUE,col='blue')


hist(TR,freq=FALSE,nclass=min(n/10,100))
lines(density(TR),col='orange')
curve(dens_TR(x,lambda_c,lambda_e,k,kprime,theta_N),add=TRUE,col='red')
#curve(dens_TV(x,lambda_c,lambda_e,k,theta_N_hat),add=TRUE,col='blue')
curve(dens_TR(x,param_hat_from_TV[1],param_hat_from_TV[2],k,kprime,theta_N_hat),add=TRUE,col='magenta')
curve(dens_TR(x,param_hat_from_TR[1],param_hat_from_TR[2],k,kprime,theta_N_hat),add=TRUE,col='green')
curve(dens_TR(x,param_hat_from_TV_TR[1],param_hat_from_TV_TR[2],k,kprime,theta_N_hat),add=TRUE,col='blue')


cbind(param,param_hat_from_TV,param_hat_from_TR,param_hat_from_TV_TR)



names(param_hat)<-c('lambda_c','lambda_e')


#--------------------- Estim from VV et VR (do we manage to remove the effect of natural death? )


param_hat_fromVV <- estim_param_SumGamma(VV,k,kprime=0)
param_hat_fromVR <- estim_param_SumGamma(VR,k,kprime)



#---------------- ESTIM from 



cbind(param,param_hat_fromVV,param_hat)


