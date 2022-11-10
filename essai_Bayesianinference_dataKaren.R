rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC_marg.R')

##################################################
################ Load data
###################################################

load('/data/TCR.txt')

T_R_C <- read.table("DataKaren/TRC.txt", quote="\"", comment.char="")$V1
T_V_C <- read.table("DataKaren/TVC.txt", quote="\"", comment.char="")$V1
T_Exp_R <- read.table("DataKaren/TRExp.txt", quote="\"", comment.char="")$V1
T_Exp_V <- read.table("DataKaren/TVExp.txt", quote="\"", comment.char="")$V1


k <- 16
kprime <- 30





par(mfrow=c(2,2))
hist(T_V_C,freq = FALSE,main="Natural Death. GREEN")  
hist(T_R_C,freq = FALSE,main="Natural Death. RED")  
hist(T_Exp_V,freq=FALSE,nclass =min(n_V/10,100), main=paste("Experi Data. GREEN",sep=' '))  
hist(T_Exp_R,freq=FALSE,nclass =min(n_V/10,100), main=paste("Experi Data. RED",sep=' '))  




mydata  <- list()
mydata$T_Contr_R <- T_R_C*2
mydata$T_Contr_V <- T_V_C
mydata$k <- k
mydata$kprime <- kprime
mydata$Tmax <- 140
mydata$T_Exp_R <- T_Exp_R*2
mydata$T_Exp_V <- T_Exp_V
mydata$n_Exp_R <- length(mydata$T_Exp_R)
mydata$n_Exp_V <- length(mydata$T_Exp_V)

########################################################################
################################# ESTIM ND  
##########################################################################

paramsChains = list(nMCMC=20000,rho=rep(1,6))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$withQD  = TRUE
paramsChains$paramsToSample = c(1:6) #
log_param_init = c(-1,-1,-2,-1,1,0)
hyperparams= list(mean=rep(0,6) ,sd = rep(3,6))





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



