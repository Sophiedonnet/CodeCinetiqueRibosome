rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC_marg.R')
source('functionsPlotPostInference.R')

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


#------------ format data
mydata  <- list()
mydata$T_Contr_R <- T_R_C*2
mydata$T_Contr_V <- T_V_C
mydata$k <-16
mydata$kprime <- 30
mydata$Tmax_Contr_R <- max(T_R_C*2)
mydata$Tmax_Contr_V <- max(T_V_C)
mydata$Tmax_Exp_R <- max(T_Exp_R*2)
mydata$Tmax_Exp_V <- max(T_Exp_V)

myCompleteData <- mydata
myCompleteData$T_Exp_R <- T_Exp_R*2
myCompleteData$T_Exp_V <- T_Exp_V


########################################################################
################################# ESTIM Contrôle Data 
##########################################################################

paramsChains = list(nMCMC=20000,rho=rep(1,9))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$paramsToSample = c(1:4) #


param_init <- init_param(myCompleteData)
log_param_init <- from_param_to_log_param(param_init)

hyperparams= list(mean=rep(0,length(param_init)) ,sd = rep(3,length(param_init)))


#---------------- MCMC
Estim_Control_Data  <- my_mcmc_marg_onechain(mydata,log_param_init,hyperparams,paramsChains)
#--------------- Plot 
plotPostInf(Estim_Control_Data, paramChains,log_param_ref = log_param_init)

#----------------- Adjust to data 
log_param_estim <- apply(Estim_Control_Data$myLogPostSample,2,mean)
param_estim <- from_logparam_to_param(log_param_estim)
par(mfrow=c(1,2))
X <- mydata$T_Contr_V[mydata$T_Contr_V<max(mydata$T_Contr_V)]
hist(X,freq = FALSE,nclass =50,main="Natural Death. GREEN")  
lines(density(X))
curve(dExpCensored(x,lambda = param_init[1],piTrunc = 0,Tmax = mydata$Tmax_Contr_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dExpCensored(x,lambda = param_estim[1],piTrunc =0,Tmax = mydata$Tmax_Contr_V),add=TRUE,col='orange',lwd=2,lty = 2)

X <- mydata$T_Contr_R[mydata$T_Contr_R<max(mydata$T_Contr_R)]
hist(X,freq = FALSE,nclass =50,main="Natural Death. Red")  
lines(density(X))
curve(dExpCensored(x,lambda = param_init[3],piTrunc = 0,Tmax = mydata$Tmax_Contr_R),add=TRUE,col='red',lwd=2,lty = 1)
curve(dExpCensored(x,lambda = param_estim[3],piTrunc = 0,Tmax = mydata$Tmax_Contr_R),add=TRUE,col='orange',lwd=2,lty = 2)


########################################################################
################################# ESTIM Exp Data 
##########################################################################


#---------------- MCMC
paramsChains$paramsToSample = c(5:7) #
log_param_init[1:4]<-log_param_estim[1:4]
log_param_init[9] <- -Inf
Estim_Complete_Data  <- my_mcmc_marg_onechain(myCompleteData,log_param_init,hyperparams,paramsChains)
#--------------- Plot 
plotPostInf(Estim_Complete_Data, paramChains,log_param_ref = log_param_init)

################################# Log param



abs <- seq(0,max(c(myCompleteData$T_Exp_R,myCompleteData$T_Exp_V)),len= 1000)


P <- computationToPlotCurves(abs, param_init,mydata)
P$comp <- "Init"
Pestim <- computationToPlotCurves(abs, param_estim,mydata)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

levels(P$Curves)
Ptemp <- allP #%>% filter(Curves %in% c("1.Natural Death"))

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  
#-------------------- On data  Exp
par(mfrow=c(1,2))
begin_title <- ifelse(param_init[9] > 0,'With','Without')
hist(myCompleteData$T_Exp_V,freq=FALSE,nclass =50, main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(myCompleteData$T_Exp_V))
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax_Exp_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_Exp_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(myCompleteData$T_Exp_R,freq=FALSE,nclass= 50, main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(myCompleteData$T_Exp_R))
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

