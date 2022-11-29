rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC_marg.R')

##################################################
################ Load data
###################################################



T_R_C <- read.table("DataKaren/TRC.txt", quote="\"", comment.char="")$V1
T_V_C <- read.table("DataKaren/TVC.txt", quote="\"", comment.char="")$V1
T_Exp_R <- read.table("DataKaren/TRExp.txt", quote="\"", comment.char="")$V1
T_Exp_V <- read.table("DataKaren/TVExp.txt", quote="\"", comment.char="")$V1






par(mfrow=c(2,2))
hist(T_V_C,freq = FALSE,nclass =50,main="Natural Death. GREEN")  
hist(T_R_C,freq = FALSE,nclass =50,main="Natural Death. RED")  
hist(T_Exp_V,freq=FALSE,nclass =50, main=paste("Experi Data. GREEN",sep=' '))  
hist(T_Exp_R,freq=FALSE,nclass =50, main=paste("Experi Data. RED",sep=' '))  


k <- 16
kprime <- 30



mydata  <- list()
mydata$T_Contr_R <- T_R_C*2
mydata$T_Contr_V <- T_V_C
mydata$k <- k
mydata$kprime <- kprime
mydata$T_Exp_R <- T_Exp_R*2
mydata$T_Exp_V <- T_Exp_V
mydata$Tmax_Contr_R <- max(T_R_C*2)
mydata$Tmax_Contr_V <- max(T_V_C)
mydata$Tmax_Exp_R <- max(T_Exp_R)
mydata$Tmax_Exp_V <- max(T_Exp_V)


########################################################################
################################# ESTIM ND  
##########################################################################

paramsChains = list(nMCMC=20000,rho=rep(1,9))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$paramsToSample = c(1:7) #

param_init <- init_param(mydata)

param_init[9] = 0
log_param_init <- from_param_to_log_param(param_init)
log_param_init
hyperparams= list(mean=rep(0,length(param_init)) ,sd = rep(3,length(param_init)))





####################################################
mySample <- my_mcmc_marg_onechain(mydata,log_param_init,hyperparams,paramsChains)





################################# Log param


par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(mySample$myLogPostSample[,m],type='l',xlab='',main=names(mySample$myPostSample)[m],ylab=''); 
#  abline(h=log_param_true[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(density(mySample$myLogPostSample[,m]),main=names(mySample$myLogPostSample)[m])
  #abline(v= log_param_true[m],col='green',lwd=2)
  curve(dnorm(x,hyperparams$mean[m],hyperparams$sd[m]),add=TRUE,col='orange',lwd=2)
}



######################## 
log_param_estim <- apply(mySample$myLogPostSample,2,mean)
param_estim <- from_logparam_to_param(log_param_estim)
abs <- seq(0,max(c(mydata$T_Contr_R,mydata$T_Contr_V)),len= 1000)


P <- computationToPlotCurves(abs, param_init,mydata)
P$comp <- "Init"
Pestim <- computationToPlotCurves(abs, param_estim,mydata)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

levels(P$Curves)
Ptemp <- allP %>% filter(Curves %in% c("1.Natural Death"))

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  

#-------------------- On data  Controle


par(mfrow=c(1,2))
hist(mydata$T_Contr_V,freq = FALSE,nclass =50,main="Natural Death. GREEN")  
lines(density(mydata$T_Contr_V))
curve(dExpCensored(x,lambda = param_init[1],piTrunc = param_init[2],Tmax = mydata$Tmax_Contr_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dExpCensored(x,lambda = param_estim[1],piTrunc = param_estim[2],Tmax = mydata$Tmax_Contr_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Contr_R,freq = FALSE,nclass =50,main="Natural Death. Red")  
lines(density(mydata$T_Contr_R))
curve(dExpCensored(x,lambda = param_init[3],piTrunc = param_init[4],Tmax = mydata$Tmax_Contr_R),add=TRUE,col='red',lwd=2,lty = 1)
curve(dExpCensored(x,lambda = param_estim[3],piTrunc = param_estim[4],Tmax = mydata$Tmax_Contr_R),add=TRUE,col='orange',lwd=2,lty = 2)

#-------------------- On data  Exp
par(mfrow=c(1,2))
begin_title <- ifelse(param_init[9] > 0,'With','Without')
hist(mydata$T_Exp_V,freq=FALSE,nclass =min(mydata$n_Exp_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(mydata$T_Exp_V))
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax_Exp_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_Exp_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Exp_R,freq=FALSE,nclass= min(mydata$n_Exp_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(mydata$T_Exp_R))
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='red',Tmax=mydata$Tmax_Exp_R),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax_Exp_R),add=TRUE,col='orange',lwd=2,lty = 2)





hist(mydata$T_Exp_V,freq=FALSE,nclass =min(mydata$n_Exp_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(mydata$T_Exp_V))
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax_Exp_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_Exp_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Exp_R,freq=FALSE,nclass= min(mydata$n_Exp_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(mydata$T_Exp_R))
curve(dOurModelExp(x,init,mydata$k,mydata$kprime,color='red',Tmax=mydata$Tmax_Exp_R),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax_Exp_R),add=TRUE,col='orange',lwd=2,lty = 2)

