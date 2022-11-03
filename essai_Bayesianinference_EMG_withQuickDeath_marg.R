rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC_marg.R')

##################################################
################ Params  for simulation
###################################################

#-  Params of natural death 
lambda_ND_V <- 1/70
lambda_ND_R <- 1/85

#-  Params of ADN reading
lambda_c <- 1/10
lambda_e <- 1/0.5
k <- 16
kprime <- 30

#-  Params of Quick death
pi_QD <- 0.3# probability of quick death
lambda_QD <- 1;

#- all Param
param_true <-  c(lambda_ND_V, lambda_ND_R,lambda_c,lambda_e,lambda_QD,pi_QD)
log_param_true <- param_true
log_param_true[1:5] <- log(param_true[1:5])
log_param_true[6] <- -log(1-param_true[6]) + log(param_true[6])
mu_true <-param_true
mu_true[1:5] <- 1/param_true[1:5]


#------------------- Nb d'observations
Tmax <- 100
n_V <- n_R <- 1000 # expériences
n_C  <- 1000 # controle

##################################################################################################
################ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) with  additional Quick death
##################################################################################################

#------- Natural death
T_R_C <- rExpCensored(n_C,lambda_ND_R,Tmax)
T_V_C <- rExpCensored(n_C,lambda_ND_V,Tmax)

par(mfrow=c(2,2))
hist(T_V_C,freq = FALSE,nclass = n_C/10,main="Natural Death")  
curve(dExpCensored(x,lambda_ND_V,Tmax),add=TRUE,col='green')
hist(T_R_C,freq = FALSE,nclass = n_C/10,main="Natural Death")  
curve(dExpCensored(x,lambda_ND_R,Tmax),add=TRUE,col='red')

plot(ecdf(T_V_C))
curve(pExpCensored(x,lambda_ND_V,Tmax),add=TRUE,col='green')
plot(ecdf(T_R_C))
curve(pExpCensored(x,lambda_ND_R,Tmax),add=TRUE,col='red')


#------- Exp with quick death 
data_Exp_V <- rOurModel(n_V,param_true,k,kprime,color='green',Tmax)
T_Exp_V <- data_Exp_V$Y
data_Exp_R <- rOurModel(n_R,param_true,k,kprime,color='red',Tmax)
T_Exp_R <- data_Exp_R$Y


par(mfrow=c(1,2))
begin_title <- ifelse(pi_QD > 0,'With','Without')
hist(T_Exp_V,freq=FALSE,nclass =min(n_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
curve(dOurModel(x,param_true,k,kprime,color='green',Tmax),add=TRUE,col='green',lwd=2,lty = 1)

hist(T_Exp_R,freq=FALSE,nclass= min(n_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
curve(dOurModel(x,param_true,k,kprime,color='red',Tmax),add=TRUE,col='red',lwd=2,lty = 1)


mydata  <- list()
mydata$T_Contr_R <- T_R_C
mydata$T_Contr_V <- T_V_C
mydata$k <- k
mydata$kprime <- kprime
mydata$Tmax <- Tmax
mydata$T_Exp_R <- T_Exp_R
mydata$T_Exp_V <- T_Exp_V
mydata$n_Exp_R <- length(mydata$T_Exp_R)
mydata$n_Exp_V <- length(mydata$T_Exp_V)

save(mydata, file='essais_datasimu.Rdata') 
save(param_true,log_param_true,file='essais_datasimu_param.Rdata')
########################################################################
################################# ESTIM ND  
##########################################################################



rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC_marg.R')

load('essais_datasimu.Rdata')
paramsChains = list(nMCMC=20000,rho=rep(1,6))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$withQD  = TRUE

paramsChains$paramsToSample = c(1:6) #

log_param_init = c(-1,-1,-2,-1,1,0)
hyperparams= list(mean=rep(0,6) ,sd = rep(3,6))
for (i in 1:6){
  if (!(i %in% paramsChains$paramsToSample)){
  log_param_init[i] <- log_param_true[i]
  }
}
 




####################################################
 
mySample <- my_mcmc_marg_onechain(mydata,log_param_init,hyperparams,paramsChains)

load('essais_datasimu_param.Rdata')

 


################################# Log param


par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(mySample$myLogPostSample[,m],type='l',xlab='',main=names(mySample$myPostSample)[m],ylab=''); 
  abline(h=log_param_true[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
    plot(density(mySample$myLogPostSample[,m]),main=names(mySample$myLogPostSample)[m])
    abline(v= log_param_true[m],col='green',lwd=2)
    curve(dnorm(x,hyperparams$mean[m],hyperparams$sd[m]),add=TRUE,col='orange',lwd=2)
}



######################## 
log_param_estim <- apply(mySample$myLogPostSample,2,mean)
abs <- seq(0,mydata$Tmax,len= 1000)
param_estim <- log_param_estim
param_estim[1:5] <- exp(log_param_estim[1:5])
param_estim[6] <- 1/(1+exp(-log_param_estim[6]))

P <- computationToPlotCurves(abs, param_true,mydata$k,mydata$kprime)
P$comp <- "True"
Pestim <- computationToPlotCurves(abs, param_estim,mydata$k,mydata$kprime)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

levels(P$Curves)
Ptemp <- allP %>% filter(Curves %in% c('2. Sum Arrival + reading', '3. min(ND,Reading)'))

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  

#-------------------- On data 


par(mfrow=c(1,2))
begin_title <- ifelse(param_true[6] > 0,'With','Without')
hist(mydata$T_Exp_V,freq=FALSE,nclass =min(mydata$n_Exp_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(mydata$T_Exp_V))
curve(dOurModel(x,param_true,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModel(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Exp_R,freq=FALSE,nclass= min(mydata$n_Exp_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(mydata$T_Exp_R))
curve(dOurModel(x,param_true,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModel(x,param_estim,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax),add=TRUE,col='orange',lwd=2,lty = 2)



