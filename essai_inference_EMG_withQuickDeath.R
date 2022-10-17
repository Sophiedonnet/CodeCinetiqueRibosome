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
theta_ND_true <- c(1,1/mu_ND_V)
theta_ND_TRUE  <- theta_ND

#------------------- Params of ADN reading
lambda_e <- 1/3
lambda_c <- 1/3
k <- 16
kprime <- 30
delta <- 0

#------------------- Params of Quick death
pi_QD <-0.0 # probability of quick death

#lambda_QD = 1/5; theta_QD <-  c(1,lambda_QD)
theta_QD <- transfo_Gamma_Param(10,10)

#------------------- Nb d'observations
n_V <- n_R <- 200 # expériences
n_C  <-100 # controle



#------------------ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) with  additional Quick death

T_R_C <- rgamma(n_C,theta_ND[1],theta_ND[2])
T_V_C <- rgamma(n_C,theta_ND[1],theta_ND[2])


T_V_QD <- rOurModel(n_V,lambda_c,k,lambda_e,theta_ND,delta,pi_QD,theta_QD)
T_R_QD <- rOurModel(n_V,lambda_c,k+kprime,lambda_e,theta_ND,delta,pi_QD,theta_QD)
hist(T_V_QD$Y,freq = FALSE,nclass = n_V/10,main="With Quick Death")  
hist(T_R_QD$Y,freq = FALSE,nclass = n_V/10,main="With Quick Death")  

#####################################################################################
########"" Parameters inference
#####################################################################################

T_R <- T_R_QD$Y
T_V <- T_V_QD$Y

# #-------------------------------- Estim de theta_ND

theta_ND_R_hat<- estim_param_Gamma(T_R_C)
theta_ND_V_hat<- estim_param_Gamma(T_V_C)



#############################################################
###########"  ESTIMATION for the model without Quick death
#############################################################
my_data <- list(TV= T_V, TR=T_R)

c(lambda_c,lambda_e)

log_lik_withoutQD(log_param = log(c(lambda_c,lambda_e)),theta_ND  = c(theta_ND_V_hat,theta_ND_R_hat),k,kprime,data=my_data)







# #################################################################################""
# ################ Ajustements par modele EMG sur VV et VR (version facile du probleme : pas de mort naturelle)
# ##############################################################################################


theta_QD_init_V = estim_param_Gamma(T_V[T_V < quantile(T_V,0.1)])
theta_QD_init_R = estim_param_Gamma(T_R[T_R < quantile(T_R,0.1)])
theta_QD_init <- (theta_QD_init_R + theta_QD_init_V)*0.5





cbind(theta_ND_TRUE,theta_ND_V_hat,theta_ND_R_hat)
# #------------------------------------- ESTIM lambda_c,lambda_e, -----
param_init <- c()
m_TV <- emg.mle(T_V[T_V > quantile(T_V,0.1)])
#m_TR <- emg.mle(T_R[T_R > quantile(T_R,0.1)])
param_init[1] <-log(m_TV@coef[3])  # log lambda_c 
param_init[2] <-log(k/m_TV@coef[1]) #  log lambda_e 
param_init[3] <- qnorm(0.1)
param_init[4:5] <- log(theta_QD_init_V)

param_hat <- optim(param_init, log_likelihood_ourModel_allcolors,theta_ND = c(theta_ND_V_hat, theta_ND_R_hat),k=k,kprime=kprime,data = list(TV = T_V,TR = T_R)) $par
pi_QD_hat <- pnorm(param_hat[3])
theta_QD_hat <- exp(param_hat[4:5])
pi_QD_hat <- pnorm(param_hat[3])



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
