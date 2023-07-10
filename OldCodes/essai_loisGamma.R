library(e1071)
##################################################
################"" Lois exponentielles
###################################################

n <- n_C  <- 2000

Delta_C <- 0
Delta <- 0

alpha_N <- 80
beta_N  <- 2
theta_N <-c(alpha_N,beta_N) 
alpha_L <- 40
beta_L <- 2
theta_L <-c(alpha_L,beta_L) 
nu_L <- 1/10

###données de contrôle
TC = Delta_C + rgamma(n_C,alpha_N,beta_N)

###données expérimentales 
UV  <- Delta + rgamma(n,alpha_N,beta_N)
VV <- Delta + rgamma(n,alpha_L,beta_L)
TV <- apply(cbind(UV,VV),1,min)

UR  <- Delta + rgamma(n,alpha_N,beta_N)
VR <- VV  + rgamma(n,nu_L,beta_L)
TR <-  apply(cbind(UR,VR),1,min)


#####  verif distributions 

#------- TV 
dens_TV <- function(x,theta_L,theta_N,Delta=0){
  z <- x- Delta  
  F1z <- pgamma(z,theta_L[1],theta_L[2])
  f1z <- dgamma(z,theta_L[1],theta_L[2])
  F2z <- pgamma(z,theta_N[1],theta_N[2])
  f2z <- dgamma(z,theta_N[1],theta_N[2])
  f1z*(1-F2z) + f2z*(1-F1z)
}

hist(TV,nclass = n/10,freq = FALSE)
curve(dens_TV(x,theta_L,theta_N),add=TRUE,col='red')

prob_TV <- function(x,theta_L,theta_N,Delta=0){
  z <- x- Delta
  F1z <- pgamma(z,theta_L[1],theta_L[2])
  F2z <- pgamma(z,theta_N[1],theta_N[2])
  1-(1-F1z)*(1-F2z)

}

plot(ecdf(TV))
curve(prob_TV(x,theta_L,theta_N),add=TRUE,col='red')


#------- TR 

dens_VR <- function(x,theta_L,nu_L,Delta=0){
  z <- x- Delta
  dgamma(z,theta_L[1] + nu_L,theta_L[2])
}

plot(density(VR))
curve(dens_VR(x,theta_L,nu_L,Delta=0),add=TRUE,col='red')


prob_VR <- function(x,theta_L,nu_L,Delta=0){
  z <- x- Delta
  pgamma(z,theta_L[1] + nu_L,theta_L[2])
}

plot(ecdf(VR))
curve(prob_VR(x,theta_L,nu_L),add=TRUE,col='red')


prob_TR <- function(x,theta_L,nu_L,theta_N,Delta=0){
  1-(1-prob_VR(x,theta_L,nu_L,Delta))*(1-pgamma(x-Delta,theta_N[1],theta_N[2]))
}

plot(ecdf(TR))
curve(prob_TR(x,theta_L,nu_L,theta_N),add=TRUE,col='red')



dens_TR <- function(x,theta_L,nu_L,theta_N,Delta=0){

  F_VR <- prob_VR(x,theta_L,nu_L,Delta)
  f_VR <- dens_VR(x,theta_L,nu_L,Delta)
  
  F_UR <- pgamma(x-Delta,theta_N[1],theta_N[2])
  f_UR <- dgamma(x-Delta,theta_N[1],theta_N[2])
  f_VR*(1-F_UR) + f_UR*(1-F_VR)

}


hist(TR,freq=FALSE,nclass=n/10)
lines(density(TR),col='orange')
curve(dens_TR(x,theta_L,nu_L,theta_N),add=TRUE,col='red')


################################################""
#------------ Estimations
#################################################



#-----------------------------------  ESTIM alpha_N, beta_N---- 
theta_N_hat<- estim_param_Gamma(TC)

#------------------------------------- ESTIM alpha_L, beta_Lt----- 


log_lik_TV <- function(theta_L,theta_N,TV){
  
  U <- dens_TV(TV,theta_L,theta_N)
  U <- U[U!=0]
  - sum(log(U))
}


theta_L_estim_init <- estim_param_Gamma(TV)


log_lik_TV(theta_L,theta_N,TV)
log_lik_TV(theta_L,theta_N_hat,TV)


theta_L_hat <- optim(theta_L_estim_init, log_lik_TV,theta_N = theta_N_hat,TV = TV)$par 

#------------------------------------- ESTIM nu_L----- 


log_lik_TR <- function(grid_nu_L,theta_L,theta_N,TR){
  m <- length(grid_nu_L)
  LL <- rep(0,m)
  for (i in 1:m){
    LL[i] <-   sum(log(dens_TR(TR,theta_L,grid_nu_L[i],theta_N)))
  }
  return(-LL)
}

grid_nu_L <- seq(0.0001,0.7,length=100)
LL <-log_lik_TR(grid_nu_L,theta_L,theta_N,TR)
plot(grid_nu_L,LL,type='l')
abline(v = nu_L,col='red')
nu_L_hat <- optimize(log_lik_TR,c(0.0000001, 0.7),tol = 0.0001,theta_L = theta_L, theta_N = theta_N,TR = TR) 
abline(v = nu_L_hat,col='orange')

LL_hat <- log_lik_TR(grid_nu_L,theta_L_hat,theta_N_hat,TR)
nu_L_hat_hat <- optimize(log_lik_TR,c(0.0000001, 0.7),tol = 0.0001,theta_L = theta_L_hat,theta_N = theta_N_hat,TR = TR) 
plot(grid_nu_L,LL,type='l',col='red',ylim = range(c(LL,LL_hat)))
lines(grid_nu_L,LL_hat,col='blue')
abline(v = nu_L_hat,col='red')
abline(v = nu_L_hat_hat,col='blue')
abline(v = nu_L,col='green')


1/nu_L_hat_hat$minimum
1/nu_L






