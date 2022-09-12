rm(list=ls())
#library(e1071)
library(ggplot2)
##################################################
################ Illustration de l'approx de somme de Gamma
###################################################

n <- n_C  <-1000

alpha_N <- 100
beta_N  <- 2
theta_N <-c(alpha_N,beta_N)
lambda_e <- 1/0.5
lambda_c <- 1/15
k <- 16
kprime <- 30

alpha_N_R <- 1000


### données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e)
VV <-   rgamma(n,1,lambda_c) + rgamma(n,k,lambda_e)
UV  <- rgamma(n,alpha_N,beta_N)
TV <- apply(cbind(UV,VV),1,min)
summary(TV)

UR  <- rgamma(n,alpha_N_R,beta_N)
VR <- VV  + rgamma(n,kprime,lambda_e)
TR <-  apply(cbind(UR,VR),1,min)
summary(UR)
summary(TR)
summary(VR)

### approx 2;
R <- sum(c(1,k)/c(lambda_c,lambda_e))
Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
k_sum <- R^2/Q
lambda_sum <- k_sum/R

plot(density(VV))
#curve(dgamma(x,k_sum, lambda_sum),add=TRUE,col='green')
m <- emg.mle(VV)
curve(demg(x,m@coef[1],m@coef[2],m@coef[3]),add=TRUE,col='cyan')




R <- sum(c(1,k+kprime)/c(lambda_c,lambda_e))
Q <- sum(c(1,k+prime)/c(lambda_c,lambda_e)^2)
k_sum <- R^2/Q
lambda_sum <- k_sum/R
m <- emg.mle(VR)

plot(density(VR))
curve(dgamma(x,k_sum, lambda_sum),add=TRUE,col='red')
curve(demg(x,m@coef[1],m@coef[2],m@coef[3]),add=TRUE,col='cyan')




#VV <-   rgamma(n,l,lambda_e) + rgamma(n,k,lambda_e)




#####  verif distributions

#------- TV

# dens_TV <- function(x,lambda_c,lambda_e,k,theta_N){
# 
# #  theta_L <- c(k+l,lambda_e)
#   R <- sum(c(1,k)/c(lambda_c,lambda_e))
#   Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
#   k_sum <- R^2/Q
#   lambda_sum <- k_sum/R
#   theta_L <- c(k_sum,lambda_sum)
#   F1z <- pgamma(x,theta_L[1],theta_L[2])
#   f1z <- dgamma(x,theta_L[1],theta_L[2])
#   

#   F2z <- pgamma(x,theta_N[1],theta_N[2])
#   f2z <- dgamma(x,theta_N[1],theta_N[2])
#   f1z*(1-F2z) + f2z*(1-F1z)
# }

dens_TV <- function(x,lambda_c,mu,sigma,theta_N){
  
  #  theta_L <- c(k+l,lambda_e)
  #R <- sum(c(1,k)/c(lambda_c,lambda_e))
  #Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
  #k_sum <- R^2/Q
  #lambda_sum <- k_sum/R
  #theta_L <- c(k_sum,lambda_sum)
  #F1z <- pgamma(x,theta_L[1],theta_L[2])
  #f1z <- dgamma(x,theta_L[1],theta_L[2])
  
  F1z <- pemg(x,lambda=lambda_c,mu=mu,sigma=sigma)
  f1z <- demg(x,lambda=lambda_c,mu=mu,sigma=sigma)
  
  F2z <- pgamma(x,theta_N[1],theta_N[2])
  f2z <- dgamma(x,theta_N[1],theta_N[2])
  f1z*(1-F2z) + f2z*(1-F1z)
}



hist(TV,nclass = n/10,freq = FALSE)
lines(density(TV),col='orange')
m <- emg.mle(VV)
lambda_c_est <- m@coef[3]
mu <- m@coef[1]
sigma <- m@coef[2]
curve(dens_TV(x,lambda_c_est,mu,sigma,theta_N),add=TRUE,col='green')


# prob_TV <- function(x,lambda_c,lambda_e,k,theta_N){
#   #theta_L <- c(k+l,lambda_e)
#   R <- sum(c(1,k)/c(lambda_c,lambda_e))
#   Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
#   k_sum <- R^2/Q
#   lambda_sum <- k_sum/R
#   theta_L <- c(k_sum,lambda_sum)
# 
#   F1z <- pgamma(x,theta_L[1],theta_L[2])
#   F2z <- pgamma(x,theta_N[1],theta_N[2])
#   1-(1-F1z)*(1-F2z)
# 
# }
prob_TV <- function(x,lambda_c,lambda_e,mu_sigma,theta_N){
  #theta_L <- c(k+l,lambda_e)
  #R <- sum(c(1,k)/c(lambda_c,lambda_e))
  #Q <- sum(c(1,k)/c(lambda_c,lambda_e)^2)
  #k_sum <- R^2/Q
  #lambda_sum <- k_sum/R
  #theta_L <- c(k_sum,lambda_sum)

  F1z <- pemg(x,lambda_e,mu,sigma)
  F2z <- pgamma(x,theta_N[1],theta_N[2])
  1-(1-F1z)*(1-F2z)

}



plot(ecdf(TV))
curve(prob_TV(x,lambda_c,lambda_e_est,mu,sigma,theta_N),add=TRUE,col='green')


#------- TR

#dens_VR <- function(x,lambda_e,l,k,kprime){
#  dgamma(x,k+kprime+l,lambda_e)
#}

#plot(density(VR))
#curve(dens_VR(x,lambda_e,l,k,kprime),add=TRUE,col='red')


#prob_VR <- function(x,lambda_e,l,k,kprime){
#  pgamma(x,k+kprime+l,lambda_e)
#}


dens_TR <- function(x,lambda_c,lambda_e,k,kprime,theta_N){
  dens_TV(x,lambda_c,lambda_e,k+kprime,theta_N)
}
  
prob_TR <- function(x,lambda_c,lambda_e,k,kprime,theta_N){
  prob_TV(x,lambda_c,lambda_e,k+kprime,theta_N)
}




plot(ecdf(TR))
theta_N_R <- c(alpha_N_R,beta_N)
curve(prob_TR(x,lambda_c,lambda_e,k,kprime,theta_N_R),add=TRUE,col='red')



#dens_TR <- function(x,lambda_e,l,k,kprime,theta_N){
#
#  F_VR <- prob_VR(x,lambda_e,l,k,kprime)
#  f_VR <- dens_VR(x,lambda_e,l,k,kprime)
#  F_UR <- pgamma(x,theta_N[1],theta_N[2])
#  f_UR <- dgamma(x,theta_N[1],theta_N[2])
#  f_VR*(1-F_UR) + f_UR*(1-F_VR)
#}



hist(TR,freq=FALSE,nclass=n/10)
lines(density(TR),col='orange')
curve(dens_TR(x,lambda_c,lambda_e,k,kprime,theta_N_R),add=TRUE,col='red')


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


#-----------------------------------  ESTIM theta_N  = (alpha_N, beta_N)----
TC <-  rgamma(n_C,alpha_N,beta_N)
theta_N_hat<- estim_param_Gamma(TC)
cbind(theta_N,theta_N_hat)
#------------------------------------- ESTIM l,lambda_e, -----


log_lik_TV_TR <- function(log_param,theta_N,k,kprime,TV,TR){
  l <- exp(log_param[1])
  lambda_e <- exp(log_param[2])

  UTV <- dens_TV(TV, lambda_e,l,k,theta_N)
  UTR <- dens_TR(TR,lambda_e,l,k,kprime,theta_N)

  - sum(log(UTV)) - sum(log(UTR))
}


lambda_e_init <- 0.5*(estim_param_Gamma(TV)[2] + estim_param_Gamma(TR)[2])
l_init <- 22#% 0.5*(estim_param_Gamma(TV)[1]-k + estim_param_Gamma(TR)[1]-k-kprime)
log_param_init <- c(log(l_init),log(lambda_e_init))
log_param<- log(c(l,lambda_e))
log_lik_TV_TR(log_param,theta_N,k,kprime,TV,TR)
log_lik_TV_TR(log_param_init,theta_N_hat,k,kprime,TV,TR)


log_param_hat <- optim(log_param_init, log_lik_TV_TR,theta_N = theta_N_hat,k=k,kprime=kprime,TV = TV,TR = TR)$par
param_hat <- exp(log_param_hat)
hist(TR,freq=FALSE,nclass=min(n/10,100))
lines(density(TR),col='orange')
curve(dens_TR(x,lambda_e,l,k,kprime,theta_N),add=TRUE,col='red')
curve(dens_TR(x,exp(log_param_hat[2]),exp(log_param_hat[1]),k,kprime,theta_N_hat),add=TRUE,col='magenta')

names(param_hat)<-c('l','lamdda_e')

all_param<- c(theta_N,param,param[2]/param[1])
names(all_param) = c('alpha_N','beta_N','l','lamdda_e','lambda_c')

all_param_hat <- c(theta_N_hat,param_hat,param_hat[2]/param_hat[1]);
cbind(all_param,all_param_hat)

