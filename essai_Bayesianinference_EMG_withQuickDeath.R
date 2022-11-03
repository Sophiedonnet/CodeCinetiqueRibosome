rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('myMCMC.R')


##################################################
################ Params  for simulation
###################################################

Tmax <- 100
#------------------- Params of natural death 
lambda_ND_V <- 1/70
lambda_ND_R <- 1/85

#------------------- Params of ADN reading
lambda_c <- 1/10
lambda_e <- 1/0.5
k <- 16
kprime <- 30

param_true = c(lambda_ND_V, lambda_ND_R,lambda_c,lambda_e)
log_param_true <- log(param_true)
mu_true <- 1/param_true

#------------------- Params of Quick death
pi_QD <-0.3# probability of quick death
lambda_QD = 1/1;

param_true <- c(param_true,lambda_QD,pi_QD)
log_param_true <- c(log_param_true,log(lambda_QD),pi_QD)
mu_true <- 1/param_true
mu_true[6] <- pi_QD
  

#------------------- Nb d'observations
n_V <- n_R <- 1000 # expériences
n_C  <- 5000 # controle


#------------------ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) with  additional Quick death
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


################################# ESTIM ND  

mydata <- list()
mydata$T_Contr_R = T_R_C
mydata$T_Contr_V = T_V_C
mydata$k = k
mydata$kprime = kprime
mydata$T_Exp_R <- T_Exp_R
mydata$T_Exp_V <- T_Exp_V
mydata$n_Exp_R <- length(mydata$T_Exp_R)
mydata$n_Exp_V <- length(mydata$T_Exp_V)
mydata$Tmax <- Tmax
mydata$Z <- list(ZV = data_Exp_V$Z, ZR = data_Exp_R$Z) 
######################################  Post Z  | Y param
Z <- sample_Z_QD(log_param_true,mydata)

par(mfrow=c(1,1))
plot(data_Exp_V$Y,Z$prob_ZV, col=data_Exp_V$Z+1)

table(Z$ZV,data_Exp_V$Z)/1000




paramsChains = list(nMCMC=20000,rho=rep(10,6))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$rho[4] <- 1
paramsChains$rho[5] <- 0.01
paramsChains$paramsToSample = c(5) #
log_param_init = c(-1,-1,-2,-1,1,0)
paramsChains$withQD  = (pi_QD>0)
if(paramsChains$withQD){
  log_param_init[5] = log(1)
  log_param_init[6] = 0.5
}

for (i in 1:6){
  if (!(i %in% paramsChains$paramsToSample)){
  log_param_init[i] <- log_param_true[i]
  }
}
 



hyperparams= list(lowerbound=rep(-10,5),upperbound=rep(0,5),a=1,b=1)
hyperparams$lowerbound[1:2] <- -log(200)
hyperparams$upper[1:2] <- -log(50)
hyperparams$lowerbound[3] <- -log(20)  #lambda_c
hyperparams$upperbound[3] <- - log(5)
hyperparams$lowerbound[4] <- -log(5) #lambda_e
hyperparams$upperbound[4] <- 2 
hyperparams$lowerbound[5] <- -log(5) #lambda_QD
hyperparams$upperbound[5] <- 2 

####################################################
mySample <- my_mcmc_onechain(mydata,log_param_init,hyperparams,paramsChains)
myPostSample <- mySample$myPostSample
myLogPostSample <- mySample$myLogPostSample
#######################################################""""


par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  plot(myPostSample[,m],type='l',xlab='',main=names(myPostSample)[m],ylab=''); 
  abline(h=param_true[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
for(m in paramsChains$paramsToSample){
  
  if(m %in% c(1:5)){
    sprior <- 1/exp(runif(10000,hyperparams$lowerbound[m],hyperparams$upperbound[m]))
  }else{
    sprior <- rbeta(10000,hyperparams$a,hyperparams$b)
  }
  plot(density(1/myPostSample[,m]),main=names(myPostSample)[m],xlim=range(sprior))
  abline(v= 1/param_true[m],col='green',lwd=2)
  lines(density(sprior),col='orange',lwd = 2)
  
}


if((3 %in% paramsChains$paramsToSample) & (4 %in% paramsChains$paramsToSample)){
  par(mfrow=c(1,1))
  plot(myPostSample[,4],myPostSample[,3],xlab=names(myPostSample)[4],ylab= names(myPostSample)[3])
  points(param_true[4],param_true[3],col='red',pch=4)
}

 


############## 
abs <- seq(0,200,len= 1000)
param_estim <- apply(myPostSample,2,mean)
P <- computationToPlotCurves(abs, param_true,k,kprime)
P$comp <- "True"
Pestim <- computationToPlotCurves(abs, param_estim,k,kprime)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

Ptemp <- allP %>% filter(Curves %in% c('1.Natural Death')) 

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



