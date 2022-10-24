rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC.R')


##################################################
################ Params  for simulation
###################################################


#------------------- Params of natural death 
mu_ND_V <- 70 # natural death 
lambda_ND_V <- 1/mu_ND_V
mu_ND_R <- 85
lambda_ND_R <- 1/mu_ND_R

#------------------- Params of ADN reading
lambda_c <- 1/10
lambda_e <- 1/3
k <- 16
kprime <- 30
delta <- 0

param_true = c(lambda_ND_V, lambda_ND_R,lambda_c,lambda_e)
log_param_true <- log(param_true)

#------------------- Params of Quick death
pi_QD <-0.1# probability of quick death
lambda_QD = 1/3;

if(pi_QD>0){
  param_true <- c(param_true,lambda_QD,pi_QD)
  log_param_true <- c(log_param_true,log(lambda_QD),pi_QD)
}

#------------------- Nb d'observations
n_V <- n_R <- 500 # expériences
n_C  <-500 # controle

#------------------ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) with  additional Quick death
#------- Natural death
T_R_C <- rexp(n_C,1/mu_ND_R)
T_V_C <- rexp(n_C,1/mu_ND_V)
hist(T_R_C,freq = FALSE,nclass = n_V/10,main="Natural Death")  
hist(T_V_C,freq = FALSE,nclass = n_V/10,main="Natural Death")  

#------- Exp with quick death 
T_Exp_V <- rOurModel(n_V,lambda_c,k,lambda_e,lambda_ND_V,pi_QD,lambda_QD)$Y
T_Exp_R <- rOurModel(n_R,lambda_c,k+kprime,lambda_e,lambda_ND_R,pi_QD,lambda_QD)$Y

par(mfrow=c(1,2))
begin_title <- ifelse(pi_QD > 0,'With','Without')
hist(T_Exp_V,freq=FALSE,nclass =n_V/10, main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(T_Exp_V))
curve(dminGammaEMGaussian(x,param_true[3],mu = k/param_true[4],sigma = sqrt(k/param_true[4]^2),theta=c(1,param_true[1])),add=TRUE,col='green',lwd=2)
curve(dOurModel(x,param_true[3],k,param_true[4],theta_ND=c(1,param_true[1]),pi_QD,lambda_QD),add=TRUE,col='green',lwd=2,lty = 2)


hist(T_Exp_R,freq=FALSE,nclass =n_R/10, main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(T_Exp_R))
curve(dminGammaEMGaussian(x,param_true[3],mu = (k+kprime)/param_true[4],sigma = sqrt((k+kprime)/param_true[4]^2),theta=c(1,param_true[2])),add=TRUE,col='red',lwd=2)
curve(dOurModel(x,param_true[3],k+kprime,param_true[4],theta_ND=c(1,param_true[2]),pi_QD,lambda_QD),add=TRUE,col='red',lwd=2,lty = 2)



################################# ESTIM ND  
mydata <- list()
mydata$T_Contr_R = T_R_C
mydata$T_Contr_V = T_V_C
#mydata$n_Contr_R <- length(mydata$T_Contr_R)
#mydata$n_Contr_V <- length(mydata$T_Contr_V)

mydata$T_Exp_R <- T_Exp_R
mydata$T_Exp_V <- T_Exp_V
mydata$n_Exp_R <- length(mydata$T_Exp_R)
mydata$n_Exp_V <- length(mydata$T_Exp_V)
mydata$k = k
mydata$kprime = kprime


#param_init = log(c(lambda_ND_V, lambda_ND_R,lambda_c,lambda_e,lambda_QD,piQD))
paramsChains = list(nMCMC=5000,rho=rep(10,6))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$rho[4] <- 1
paramsChains$paramsToSample =c(6) # 1:4 #
log_param_init = c(-1,-1,-2,-1)
for (i in 1:6){
  if (!(i %in% paramsChains$paramsToSample)){
  log_param_init[i] <- log(param_true[i])
  }
}


paramsChains$withQD  = TRUE
if(paramsChains$withQD){
  log_param_init[5] = -log(3)
  log_param_init[6] = 0.5
}




hyperparams= list(lowerbound=c(-10,-10,-10,-10,-10),upperbound=rep(0,5),a=1,b=1)
hyperparams$lowerbound[3] <- -log(20)  #lambda_c
hyperparams$upperbound[3] <- - log(5)

hyperparams$lowerbound[4] <- -log(5) #lambda_e
hyperparams$upperbound[4] <- -log(1) 

hyperparams$lowerbound[5] <- -log(10) #lambda_QD
hyperparams$upperbound[5] <- -log(1) 

####################################################
myPostSample <- my_mcmc_onechain(mydata,log_param_init,hyperparams,paramsChains)
#######################################################""""

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(myPostSample[,m],type='l',xlab='',main=names(myPostSample)[m],ylab=''); 
  abline(h=log_param_true[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(density(myPostSample[,m]),main=names(myPostSample)[m],xlim=c(hyperparams$lowerbound[m],hyperparams$upperbound[m]))
  abline(v= log_param_true[m],col='green',lwd=2)
  curve(dunif(x,hyperparams$lowerbound[m],hyperparams$upperbound[m]),col='orange',add=TRUE)
  
}

if((3 %in% paramsChains$paramsToSample) & (4 %in% paramsChains$paramsToSample)){
  par(mfrow=c(1,1))
  plot(myPostSample[,4],myPostSample[,3],xlab=names(myPostSample)[4],ylab= names(myPostSample)[3])
  points(log_param_true[4],log_param_true[3],col='red',pch=4)
}

 


############## THE GREEN
abs <- seq(0,max(c(mydata$T_Exp_V,mydata$T_Exp_R)),len= 1000)
P <- computationToPlotCurves(abs, exp(log_param_true[1]),exp(log_param_true[3]),exp(log_param_true[4]),nbCodons=k)
P$comp <- "True"
Pestim <- computationToPlotCurves(abs, exp(mean(myPostSample[,1])),exp(mean(myPostSample[,3])),exp(mean(myPostSample[,4])),nbCodons=k)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)
Ptemp <- allP %>% filter(Curves %in% c('0.Natural Death','3. Sum Arrival + reading','4. min(ND,Reading)')) 
ggplot(data = Ptemp, aes(x=time,y=density,colour=Curves))+ geom_line(aes(linetype=comp))  


############## THE RED
abs <- seq(0,max(mydata$T_Exp_R),len= 1000)
P <- computationToPlotCurves(abs, exp(log_param_true[2]),exp(log_param_true[3]),exp(log_param_true[4]),nbCodons=k+kprime)
P$comp <- "True"
Pestim <- computationToPlotCurves(abs, exp(mean(myPostSample[,2])),exp(mean(myPostSample[,3])),exp(mean(myPostSample[,4])),nbCodons=k+kprime)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)
Ptemp <- allP %>% filter(Curves %in% c('0.Natural Death','3. Sum Arrival + reading','4. min(ND,Reading)')) 
ggplot(data = Ptemp, aes(x=time,y=density,colour=Curves))+ geom_line(aes(linetype=comp))  


par(mfrow=c(1,2))
param_estim <- apply(myPostSample,2,mean)
param_estim[-5] <- exp(param_estim[-5])


hist(T_Exp_V,freq=FALSE,nclass =100, main=paste(begin_title," Quick Death Green",sep=' '))  
lines(density(T_Exp_V))
#curve(dminGammaEMGaussian(x,param_estim[3],mu = k/param_estim[4],sigma = sqrt(k)/param_estim[4],theta=c(1,param_estim[1])),add=TRUE,col='green',lwd=2,lty=2,)
#curve(dminGammaEMGaussian(x,param_true[3],mu = k/param_true[4],sigma = sqrt(k)/param_true[4],theta=c(1,param_true[1])),add=TRUE,col='green',lwd=2)
curve(dOurModel(x,param_true[3],k,param_true[4],theta_ND=c(1,param_true[1]),param_true[6],param_true[5]),add=TRUE,col='green',lwd=2,lty = 2)
curve(dOurModel(x,param_estim[3],k,param_estim[4],theta_ND=c(1,param_estim[1]),param_estim[6],param_estim[5]),add=TRUE,col='green',lwd=2,lty = 2)




hist(T_Exp_R,freq=FALSE,nclass =100, main=paste(begin_title," Quick Death Red",sep=' '))
lines(density(T_Exp_R))
curve(dminGammaEMGaussian(x,param_estim[3],mu = (k+kprime)/param_estim[4],sigma = sqrt(k+kprime)/param_estim[4],theta=c(1,param_estim[2])),add=TRUE,col='red',lty=2,lwd=2)
curve(dminGammaEMGaussian(x,param_true[3],mu = (k+kprime)/param_true[4],sigma = sqrt(k+kprime)/param_true[4],theta=c(1,param_true[2])),add=TRUE,col='red',lwd=2)

