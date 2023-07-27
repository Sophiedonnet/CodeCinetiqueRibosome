rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC.R')


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

########################################################################
################################# ESTIM ND  
##########################################################################
mydata_Contr <- list()
mydata_Contr$T_Contr_R <- T_R_C
mydata_Contr$T_Contr_V <- T_V_C
mydata_Contr$k <- k
mydata_Contr$kprime <- kprime
mydata_Contr$Tmax <- Tmax
mydata <- mydata_Contr
mydata$T_Exp_R <- T_Exp_R
mydata$T_Exp_V <- T_Exp_V
mydata$n_Exp_R <- length(mydata$T_Exp_R)
mydata$n_Exp_V <- length(mydata$T_Exp_V)

 


paramsChains = list(nMCMC=50000,rho=rep(10,6))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$rho[4] <- 1
paramsChains$rho[5] <- 1
paramsChains$paramsToSample = c(1:6) #
log_param_init = c(-1,-1,-2,-1,1,0.5)
paramsChains$withQD  = (pi_QD>0)

for (i in 1:6){
  if (!(i %in% paramsChains$paramsToSample)){
  log_param_init[i] <- log_param_true[i]
  }
}
 



hyperparams= list(mean=log_param_true[1:5],sd = rep(0.5,5),a=1,b=1)

#hyperparams= list(lowerbound=rep(-10,5),upperbound=rep(0,5),a=1,b=1)

#hyperparams$lowerbound[1:2] <- -log(200)
#hyperparams$upper[1:2] <- -log(30)
#hyperparams$lowerbound[3] <- -log(20)  #lambda_c
#hyperparams$upperbound[3] <- log(5)
#hyperparams$lowerbound[4] <- -log(20) #lambda_e
#hyperparams$upperbound[4] <- log(2)
#hyperparams$lowerbound[5] <- log(1/20)  #lambda_QD
#hyperparams$upperbound[5] <- log(2) 


#-------------------------

hyperparams$b <- 10
mode_piQD <- 0.3
hyperparams$a <- (1+(hyperparams$b-2)*mode_piQD)/(1-mode_piQD)
curve(dbeta(x,hyperparams$a,hyperparams$b))

U <- rbeta(1000,hyperparams$a, hyperparams$b)
hyperparams$mean[6] <- mean(U)
hyperparams$sd[6] <- sd(U)


####################################################
mySample_Contr <- my_mcmc_onechain(mydata_Contr,log_param_init,hyperparams,paramsChains)
hyperparams$mean[1] <- mean(mySample_Contr$myLogPostSample$log_lambda_ND_V)
hyperparams$mean[2] <- mean(mySample_Contr$myLogPostSample$log_lambda_ND_R)
hyperparams$sd[1] <- sd(mySample_Contr$myLogPostSample$log_lambda_ND_V)
hyperparams$sd[2] <- sd(mySample_Contr$myLogPostSample$log_lambda_ND_R)

mySample <- my_mcmc_marg_onechain(mydata,log_param_init,hyperparams,paramsChains)




myPostSample <- mySample$myPostSample
myLogPostSample <- mySample$myLogPostSample


#######################################################""""


par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(myLogPostSample[,m],type='l',xlab='',main=names(myPostSample)[m],ylab=''); 
  abline(h=log_param_true[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  
  if(m %in% c(1:5)){
    sprior <- 1/exp(rnorm(10000,hyperparams$mean[m],hyperparams$sd[m]))
    plot(density(1/myPostSample[,m]),main=names(myPostSample)[m],xlim=range(sprior))
    abline(v= 1/param_true[m],col='green',lwd=2)
    lines(density(sprior),col='orange',lwd = 2)
  }else{
    plot(density(myPostSample[,m]),main=names(myPostSample)[m],xlim=c(0,1))
    abline(v= param_true[m],col='green',lwd=2)
    curve(dbeta(x,hyperparams$a,hyperparams$b),col='orange',lwd = 2,add = TRUE)
  }
}

par(mfrow=c(1,2))
curve(dnorm(x,hyperparams$mean[1],hyperparams$sd[1]),-4.5,-4)
lines(density(mySample_Contr$myLogPostSample$log_lambda_ND_V),col='magenta')
lines(density(mySample$myLogPostSample$log_lambda_ND_V),col='red')
abline(v=log_param_true[1],col='green')
curve(dnorm(x,hyperparams$mean[2],hyperparams$sd[2]),-5,-4)
lines(density(mySample_Contr$myLogPostSample$log_lambda_ND_R),col='magenta')
lines(density(mySample$myLogPostSample$log_lambda_ND_R),col='red')
abline(v=log_param_true[2],col='green')




if((3 %in% paramsChains$paramsToSample) & (4 %in% paramsChains$paramsToSample)){
  par(mfrow=c(1,1))
  plot(myPostSample[,4],myPostSample[,3],xlab=names(myPostSample)[4],ylab= names(myPostSample)[3])
  points(param_true[4],param_true[3],col='red',pch=4)
}



log_param_estim <- apply(myLogPostSample,2,mean)
R_true <- sample_Z_QD(log_param_true,mydata)  
R_estim <- sample_Z_QD(log_param_estim,mydata)  
plot(R_true$prob_ZR,R_estim$prob_ZR)
abline(a=0,b=1,col='red')

plot(R_true$prob_ZV,R_estim$prob_ZV)
abline(a=0,b=1,col='red')



############## 
abs <- seq(0,Tmax,len= 1000)











param_estim <- log_param_estim
param_estim[1:5] <- exp(log_param_estim[1:5])


P <- computationToPlotCurves(abs, param_true,k,kprime)
P$comp <- "True"
Pestim <- computationToPlotCurves(abs, param_estim,k,kprime)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

Ptemp <- allP #%>% filter(Curves %in% c('1.Natural Death','')) 

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  

#-------------------- On data 


par(mfrow=c(1,2))
begin_title <- ifelse(pi_QD > 0,'With','Without')
hist(T_Exp_V,freq=FALSE,nclass =min(n_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(T_Exp_V))
curve(dOurModel(x,param_true,k,kprime,color='green'),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModel(x,param_estim,k,kprime,color='green'),add=TRUE,col='green',lwd=2,lty = 2)

hist(T_Exp_R,freq=FALSE,nclass= min(n_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(T_Exp_R))
curve(dOurModel(x,param_true,k,kprime,color='red'),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModel(x,param_estim,k,kprime,color='red'),add=TRUE,col='red',lwd=2,lty = 2)



