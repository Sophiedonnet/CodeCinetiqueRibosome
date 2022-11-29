rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC_marg.R')

##################################################
################ Data
###################################################
T_R_C <- read.table("DataKaren/TRC.txt", quote="\"", comment.char="")$V1
T_V_C <- read.table("DataKaren/TVC.txt", quote="\"", comment.char="")$V1
T_Exp_R <- read.table("DataKaren/TRExp.txt", quote="\"", comment.char="")$V1
T_Exp_V <- read.table("DataKaren/TVExp.txt", quote="\"", comment.char="")$V1

#---------------- param expérimentaux
k <- 16
kprime <- 32

##################################################
################ Etude data de controle
######################################################

mydata_Contr  <- list()
mydata_Contr$T_Contr_R <- T_R_C*2.8
mydata_Contr$T_Contr_V <- T_V_C*1.4
mydata_Contr$k <- k
mydata_Contr$kprime <- kprime
mydata_Contr$Tmax_R <-max(mydata_Contr$T_Contr_R)
mydata_Contr$Tmax_V <-max(mydata_Contr$T_Contr_V)

par(mfrow=c(1,2))
hist(mydata_Contr$T_Contr_V,nclass=50,freq = FALSE,main="Natural Death. GREEN")  
hist(mydata_Contr$T_Contr_R,nclass=50,freq = FALSE,main="Natural Death. RED")  



#-------------- tuning param  
 
paramsChains = list(nMCMC=20000,rho=rep(1,6))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$withQD  = TRUE
paramsChains$paramsToSample = c(1:2) #
log_param_init = c(-1,-1,-2,-1,1,0)
log_param_init[1] <- log(1/mean(mydata_Contr$T_Contr_V[mydata_Contr$T_Contr_V<mydata_Contr$Tmax_V]))
log_param_init[2] <- log(1/mean(mydata_Contr$T_Contr_R[mydata_Contr$T_Contr_R<mydata_Contr$Tmax_R]))
hyperparams= list(mean=rep(0,6) ,sd = rep(3,6))



#-------------- MCMC
mySample_Contr <- my_mcmc_marg_onechain(mydata_Contr,log_param_init,hyperparams,paramsChains)

save(mySample_Contr, file='DataKaren/results_estim_data_Contr.Rdata')

#---------------- plot post sample

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(mySample_Contr$myLogPostSample[,m],type='l',xlab='',main=names(mySample_Contr$myPostSample)[m],ylab=''); 
  abline(h=log_param_init[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(density(mySample_Contr$myLogPostSample[,m]),main=names(mySample_Contr$myLogPostSample)[m])
  abline(v= log_param_init[m],col='green',lwd=2)
  curve(dnorm(x,hyperparams$mean[m],hyperparams$sd[m]),add=TRUE,col='orange',lwd=2)
}

#----------------- Fitting
log_param_estim <- apply(mySample_Contr$myLogPostSample,2,mean)
abs <- sort(c(seq(0,max(mydata_Contr$Tmax_R,mydata_Contr$Tmax_V+10),len= 1000),mydata_Contr$Tmax_R,mydata_Contr$Tmax_V))
param_estim <- log_param_estim
param_estim[1:5] <- exp(log_param_estim[1:5])
param_estim[6] <- 1/(1+exp(-log_param_estim[6]))


param_init <- log_param_estim
param_init[1:5] <- exp(log_param_init[1:5])
param_init[6] <- 1/(1+exp(-log_param_init[6]))

P <- computationToPlotCurves(abs, param_init,mydata_Contr$k,mydata_Contr$kprime,Tmax=list(R=mydata_Contr$Tmax_R,V=mydata_Contr$Tmax_V))
P$comp <- "Init"
Pestim <- computationToPlotCurves(abs, param_estim,mydata_Contr$k,mydata_Contr$kprime,Tmax=list(R=mydata_Contr$Tmax_R,V=mydata_Contr$Tmax_V))
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

Ptemp <- allP %>% filter(Curves %in% c('1.Natural Death'))

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  

#-------------------- On data 


par(mfrow=c(1,2))

hist(mydata_Contr$T_Contr_V,freq=FALSE,nclass =min(length(mydata_Contr$T_Contr_V)/10,50), main=paste("Quick Death. Green",sep=' '),xlim=c(0,178),ylim=c(0,0.02))  
lines(density(mydata_Contr$T_Contr_V))
#curve(dOurModel(x,param_true,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax),add=TRUE,col='green',lwd=2,lty = 1)
curve(dExpCensored(x,param_estim[1],Tmax = mydata_Contr$Tmax_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dExpCensored(x,param_init[1],Tmax = mydata_Contr$Tmax_V),add=TRUE,col='green',lwd=2,lty = 3)


hist(mydata_Contr$T_Contr_R,freq=FALSE,nclass= min(length(mydata_Contr$T_Contr_R)/10,50), main=paste("Quick Death. Red",sep=' '),xlim=c(0,178),ylim=c(0,0.02))  
lines(density(mydata_Contr$T_Contr_R))
curve(dExpCensored(x,param_estim[2]),add=TRUE,col='red',lwd=2,lty = 1)
curve(dExpCensored(x,param_init[2],Tmax = mydata_Contr$Tmax_R),add=TRUE,col='red',lwd=2,lty = 3)






########################################################" 
hist(mydata$T_Exp_V,freq=FALSE,nclass =min(mydata_Contr$/10,100), main=paste("Quick Death. Green",sep=' '))  
lines(density(mydata$T_Exp_V))
#curve(dOurModel(x,param_true,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModel(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Exp_R,freq=FALSE,nclass= min(mydata$n_Exp_R/10,100), main=paste("Quick Death. Red",sep=' '))  
lines(density(mydata$T_Exp_R))
curve(dOurModel(x,param_estim,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax_R),add=TRUE,col='orange',lwd=2,lty = 2)



#mydata$T_Exp_R <- T_Exp_R*2.8
#mydata$T_Exp_V <- T_Exp_V*1.4
#mydata$n_Exp_R <- length(mydata$T_Exp_R)
#mydata$n_Exp_V <- length(mydata$T_Exp_V)
hist(mydata$T_Exp_V,nclass=50,freq=FALSE, main=paste("Experi Data. GREEN",sep=' '))  
hist(mydata$T_Exp_R,nclass=50,freq=FALSE, main=paste("Experi Data. RED",sep=' '))  


