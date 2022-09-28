rm(list=ls())
library(ggplot2)
library(emg)
source('functionsEMG.R')

##################################################
################ Params  for simulation
###################################################


#------------------- Params of natural death 
mu_ND_V <- 70 # natural death 
sd_ND_V <-400
#mu_ND_R <- 85
#sd_ND_R <-10
theta_ND <- transfo_Gamma_Param(mu_ND_V,sd_ND_V)
theta_ND <- c(1,1/mu_ND_V)

#------------------- Params of ADN reading
lambda_e <- 1/3
lambda_c <- 1/3
k <- 16
kprime <- 30
delta = 0

#------------------- Params of Quick death
pi_QD <-0.1 # probability of quick death

#lambda_QD = 1/5; theta_QD <-  c(1,lambda_QD)
theta_QD <- transfo_Gamma_Param(10,10)

#------------------- Nb d'observations
n_V <- n_R <- 10000 # expériences
n_C  <-1000 # controle



########################################################################### 
#--------------------- Compute densities of all phenomenon
########################################################################### 
nbCurves <- 7
len_abs <- 1000
x <- seq(0,200,length = len_abs )
P <- as.data.frame(rep(x,nbCurves)); names(P)='time'
P$density<- c(dgamma(x,theta_ND[1],theta_ND[2]),
              dexp(x,lambda_c),
              dgamma(x,k,lambda_e), 
              demg(x,lambda_c,mu = k/lambda_e + delta,sigma = sqrt(k/lambda_e^2)),
              dminGammaEMGaussian(x,lambda_c,mu = k/lambda_e + delta,sigma = sqrt(k/lambda_e^2),theta=theta_ND),
              dgamma(x,theta_QD[1],theta_QD[2]),
              dOurModel(x,lambda_c,k,lambda_e,theta_ND,delta,pi_QD,theta_QD))

P$Curves = rep(c('0.Natural Death',
                 '1. Arrival Time',
                 '2. Reading k codons',
                 '3. Sum Arrival + reading',
                 '4. min(ND,Reading)',
                 '5. Quick Death',
                 '6. Final model'
                 ),each=len_abs)



########################################################################## 
#--------------------- Plots various component of the model
########################################################################### 
library(dplyr)

Ptemp <- P %>% filter(Curves %in% c('1. Arrival Time','2. Reading k codons','3. Sum Arrival + reading'))
ggplot(data = Ptemp, aes(x=time,y=density,colour=Curves))+ geom_line(size=1.1)              


Ptemp <- P %>% filter(Curves %in% c('0.Natural Death','3. Sum Arrival + reading','4. min(ND,Reading)')) 
ggplot(data = Ptemp, aes(x=time,y=density,colour=Curves))+ geom_line(size=1.1)              


Ptemp <- P %>% filter(Curves %in% c('4. min(ND,Reading)',
                                    '6. Final model')) 
ggplot(data = Ptemp, aes(x=time,y=density,colour=Curves))+ geom_line(size=1.1)+ggtitle('With posible Quicke Death')              




#------------------ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) without Quick death

delta  <- 0
T_V <- rminGammaEMGamma(n_V,lambda_c,k,lambda_e,theta_ND , delta)

par(mfrow=c(2,1))
hist(T_V,freq = FALSE,nclass = n_V/100,main="Without Quick Death")  
curve(dminGammaEMGaussian(x,lambda_c,mu=k/lambda_e + delta,sigma=sqrt(k/lambda_e^2),theta = theta_ND),add=TRUE,col='orange',lwd=2)

#------------------ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) with  additional Quick death


T_V_QD <- rOurModel(n_V,lambda_c,k,lambda_e,theta_ND,delta,pi_QD,theta_QD)
hist(T_V_QD$Y,freq = FALSE,nclass = n_V/100,main="With Quick Death")  



# 
# 
# #################################################################################""
# ################ Ajustements par modele EMG sur VV et VR (version facile du probleme : pas de mort naturelle)
# ##############################################################################################
# m <- emg.mle(VV)
# theta_hat <- estim_param_Gamma(UV)
# 
# par(mfrow= c(2,1))
# hist(VV,nclass = n/10,freq = FALSE)
# lines(density(VV),col='orange',lwd = 2)
# curve(demg(x,lambda=m@coef[3],mu = m@coef[1],sigma =  m@coef[2]),add=TRUE,col='green',lwd = 2)
# 
# hist(TV,nclass = n/10,freq = FALSE)
# lines(density(TV),col='orange',lwd = 2)
# curve(dens_minGammaEMG(x,lambda=m@coef[3],mu = m@coef[1],sigma =  m@coef[2], theta = theta_hat),add=TRUE,col='green',lwd = 2)
# 
# 
# 
# ################################################""
# #------------ Estimations
# #################################################
# 
# 
# #-------------------------------- Estim de theta_N
# theta_ND_hat<- estim_param_Gamma(TC)
# cbind(theta_ND_TRUE,theta_ND_hat)
# 
# #------------------------------------- ESTIM lambda_c,lambda_e, -----
# log_param_0 <- c()
# m_TV <- emg.mle(TV)
# m_TR <- emg.mle(TR) 
# log_param_0[1] <-log(m_TV@coef[3])
# log_param_0[2] <-log(k/m_TV@coef[1])
# 
# 
# #log_param_0[1] <-log((m_TV@coef[3]+ m_TR@coef[3])/2)
# #log_param_0[2] <- log(0.5) + log(k/m_TV@coef[1] + (k+kprime)/m_TR@coef[1] )
# 
# cbind(log(c(lambda_c,lambda_e)),log_param_0)
# 
# 
# log_param_hat <- optim(log_param_0, log_lik_TV_TR,theta_N = theta_ND_hat,k=k,kprime=kprime,data = list(TV = TV,TR = TR)) $par
# 
# 
# param_hat <- exp(log_param_hat)
# hist(TR,freq=FALSE,nclass=min(n/10,100))
# lines(density(TR),col='orange')
# curve(dens_TR(x,lambda_e,l,k,kprime,theta_N),add=TRUE,col='red')
# curve(dens_TR(x,exp(log_param_hat[2]),exp(log_param_hat[1]),k,kprime,theta_ND_hat),add=TRUE,col='magenta')
# 
# names(param_hat)<-c('l','lamdda_e')
# 
# all_param<- c(theta_N,param,param[2]/param[1])
# names(all_param) = c('alpha_N','beta_N','l','lamdda_e','lambda_c')
# 
# all_param_hat <- c(theta_ND_hat,param_hat,param_hat[2]/param_hat[1]);
# cbind(all_param,all_param_hat)
# 
