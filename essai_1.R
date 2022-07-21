##################################################
################"" Lois exponentielles
###################################################


Delta = 5 
theta = 2
n = 1000
Y = Delta + rexp(n,theta) 
 
hist(Y)

lik_expDelta <- function(seq_Delta,theta,Y){
  n <- length(Y)
  L <- (n*theta - sum(Y)*log(theta) + n*seq_Delta* theta) * (min(Y)>seq_Delta) 
  return(L)
}

theta_trial = 0.5
seq_Delta = seq(1,10, by= 0.1)
plot(seq_Delta,lik_expDelta(seq_Delta,theta_trial,Y), type='l')



######################" Trial my model, check density 
theta_N <- 1/40
theta_L <- 1/20
nu_L <- 1/10
Delta_C <- 2
Delta <- 3
n = n_C  = 200

###données de contrôle
TC = Delta_C + rexp(n_C,theta_N)

###données expérimentales 
UV  <- Delta + rexp(n,theta_N)
VV <- Delta + rexp(n,theta_L)
TV <- apply(cbind(UV,VV),1,min)

UR  <- Delta + rexp(n,theta_N)
VR <- VV  + rexp(n,nu_L)
TR <-  apply(cbind(UR,VR),1,min)


#####  verif distributions 

#------- TV 
dens_TV <- function(x,theta_L,theta_N,Delta){
   z <- x- Delta  
   dexp(z,theta_L + theta_N)
}

hist(TV,nclass = n/10,freq = FALSE)
curve(dens_TV(x,theta_L,theta_N,Delta),add=TRUE,col='red')

prob_TV <- function(x,theta_L,theta_N,Delta){
  z <- x- Delta  
  pexp(z,theta_L + theta_N)
}

plot(ecdf(TV))
curve(prob_TV(x,theta_L,theta_N,Delta),add=TRUE,col='red')


#------- TR 

dens_VR <- function(x,theta_L,nu_L,Delta){
  z <- x- Delta
  -  (z>0)*theta_L * nu_L/(theta_L-nu_L) * (exp(-theta_L*z) - exp(-nu_L*z))
}

plot(density(VR))
curve(dens_VR(x,theta_L,nu_L,Delta),add=TRUE,col='red')


prob_VR <- function(x,theta_L,nu_L,Delta){
  z <- x- Delta
  (z>=0)*(1-1/(nu_L-theta_L) * (nu_L*exp(-theta_L*z) - theta_L*exp(-nu_L*z)))
}

plot(ecdf(VR))
curve(prob_VR(x,theta_L,nu_L,Delta),add=TRUE,col='red')


prob_TR <- function(x,theta_L,nu_L,theta_N,Delta){
  1-(1-prob_VR(x,theta_L,nu_L,Delta))*(1-pexp(x-Delta,theta_N))
}

plot(ecdf(TR))
curve(prob_TR(x,theta_L,nu_L,theta_N,Delta),add=TRUE,col='red')



dens_TR <- function(x,theta_L,nu_L,theta_N,Delta){
  z <- x- Delta
  R1 <- theta_L * (theta_N + nu_L)*exp(-(theta_N+nu_L)*z)
  R2 <-nu_L * (theta_L + theta_N)*exp(-(theta_N+theta_L)*z)
  (z>=0)/(theta_L-nu_L) * (R1-R2)
}


hist(TR,freq=FALSE,nclass=n/10)
lines(density(TR),col='orange')
curve(dens_TR(x,theta_L,nu_L,theta_N,Delta),add=TRUE,col='red')


################################################""
#------------ Estimations
#################################################

Delta_C_hat <- min(TC)
theta_N_hat <- 1/mean(TC-Delta_C_hat)
print(c(1/theta_N_hat, 1/theta_N))



Delta_hat <- min(c(TR,TV))
theta_L_hat <-1/mean(TV-Delta_hat) - theta_N_hat 
print(c(1/theta_L_hat, 1/theta_L))


log_lik_TR <- function(grid_nu_L,theta_L,theta_N,Delta,Y){

  m <- length(grid_nu_L)
  LL <- rep(0,m)
  for (i in 1:m){
    LL[i] <-   sum(log(dens_TR(Y,theta_L,grid_nu_L[i],theta_N,Delta)))
  }
  return(-LL)
}

grid_nu_L <- seq(0.0001,0.2,length=100)
LL <-log_lik_TR(grid_nu_L,theta_L,theta_N,Delta,TR)
plot(grid_nu_L,LL,type='l')
abline(v = nu_L,col='red')
nu_L_hat <- optimize(log_lik_TR,c(0.0000001, 0.2),tol = 0.0001,theta_L = theta_L,theta_N = theta_N,Delta= Delta,Y = TR) 
abline(v = nu_L_hat,col='orange')

LL_hat <- log_lik_TR(grid_nu_L,theta_L_hat,theta_N_hat,Delta_hat,TR)
nu_L_hat_hat <- optimize(log_lik_TR,c(0.0000001, 0.2),tol = 0.0001,theta_L = theta_L_hat,theta_N = theta_N_hat,Delta= Delta_hat,Y = TR) 
plot(grid_nu_L,LL,type='l',col='red')
lines(grid_nu_L,LL_hat,col='blue')
abline(v = nu_L_hat,col='red')
abline(v = nu_L_hat_hat,col='blue')
abline(v = nu_L,col='green')


1/nu_L_hat_hat$minimum
1/nu_L



#rhat <- n*sum(X*log(X)) - sum(log(X))*sum(X)
#alpha_hat <- n*sum(X) /rhat
#beta_hat <- n^2 / rhat

#plot(density(X))
#curve(dgamma(x,shape = alpha_hat,rate = beta_hat),add=TRUE,col='red')

#Y <- rgamma(n,shape= alpha_hat, rate = beta_hat)
#plot(density(Y))

### approx Gamma


