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
piTrunc_ND_V <- 0.3

lambda_ND_R <- 1/85
piTrunc_ND_R <- 0.2

#-  Params of ADN reading
lambda_c <- 1/10
lambda_e <- 1/0.5
k <- 16
kprime <- 30
piTrunc_Read <- 0.2

#-  Params of Quick death
pi_QD <- 0# probability of quick death
lambda_QD <- 1;

#- all Param
param_true <-  c(lambda_ND_V, piTrunc_ND_V, lambda_ND_R,piTrunc_ND_R,lambda_c,lambda_e,piTrunc_Read,lambda_QD,pi_QD)
names(param_true) <- c('lambda_ND_V','piTrunc_ND_V','lambda_ND_R','piTrunc_ND_R','lambda_c','lambda_e','piTrunc_Read', 'lambda_QD','pi_QD')
log_param_true <- from_param_to_log_param(param_true)

#------------------- Nb d'observations
Tmax <- 100
n_V <- n_R <- 1000 # expériences
n_C  <- 1000 # controle

##################################################################################################
################ données expérimentales  SIMULATION  exp(lambda_c) + Gamma(k, lambda_e) with  additional Quick death
##################################################################################################

#------- Natural death
T_V_C <- rExpCensored(n_C,lambda_ND_V,Tmax,piTrunc = piTrunc_ND_V)
T_R_C <- rExpCensored(n_C,lambda_ND_R,Tmax,piTrunc = piTrunc_ND_R)

par(mfrow=c(2,2))
hist(T_V_C,freq = FALSE,nclass = n_C/10,main="Natural Death")  
curve(dExpCensored(x,lambda_ND_V,Tmax,piTrunc_ND_V),add=TRUE,col='green')
hist(T_R_C,freq = FALSE,nclass = n_C/10,main="Natural Death")  
curve(dExpCensored(x,lambda_ND_R,Tmax,piTrunc_ND_R),add=TRUE,col='red')

plot(ecdf(T_V_C))
curve(pExpCensored(x,lambda_ND_V,Tmax,piTrunc_ND_V),add=TRUE,col='green')
plot(ecdf(T_R_C))
curve(pExpCensored(x,lambda_ND_R,Tmax,piTrunc_ND_R),add=TRUE,col='red')


#------- Exp with quick death 
data_Exp_V <- rOurModelExp(n_V,param_true,k,kprime,color='green',Tmax)
T_Exp_V <- data_Exp_V$Y
data_Exp_R <- rOurModelExp(n_R,param_true,k,kprime,color='red',Tmax)
T_Exp_R <- data_Exp_R$Y


par(mfrow=c(1,2))
begin_title <- ifelse(pi_QD > 0,'With','Without')
hist(T_Exp_V,freq=FALSE,nclass =min(n_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
curve(dOurModelExp(x,param_true,k,kprime,color='green',Tmax),add=TRUE,col='green',lwd=2,lty = 1)

hist(T_Exp_R,freq=FALSE,nclass= min(n_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
curve(dOurModelExp(x,param_true,k,kprime,color='red',Tmax),add=TRUE,col='red',lwd=2,lty = 1)


mydata  <- list()
mydata$T_Contr_R <- T_R_C
mydata$T_Contr_V <- T_V_C
mydata$Tmax_Contr_R <- max(T_R_C)
mydata$Tmax_Contr_V <- max(T_V_C)
mydata$T_Exp_R <- T_Exp_R
mydata$T_Exp_V <- T_Exp_V
mydata$Tmax_Exp_R <- max(T_Exp_R)
mydata$Tmax_Exp_V <- max(T_Exp_V)

mydata$k <- k
mydata$kprime <- kprime

mycompletedata <- mydata
mycompletedata$T_Exp_R <- T_Exp_R
mycompletedata$T_Contr_R <- T_R_C
mycompletedata$T_Contr_V <- T_V_C

# save(mydata, file='essais_datasimu.Rdata') 
# save(param_true,log_param_true,file='essais_datasimu_param.Rdata')
########################################################################
################################# ESTIM ND  
##########################################################################

data <- mycompletedata
FNV <- ecdf(data$T_Exp_V)
U <- function(x){
  1-(1-FNV(x))/(1-pExpCensored(x,lambda=param_true[1],Tmax =data$Tmax_Contr_V,piTrunc = param_true[2]))
}
V  <- function(x){
    pemgCensored(x,mu=data$k/param_true[6],sigma = sqrt(data$k)/param_true[6],lambda = param_true[5],Tmax=data$Tmax_Exp_V,piTrunc=param_true[7])
}

V2  <- function(x){
  demgCensored(x,mu=data$k/param_true[6],sigma = sqrt(data$k)/param_true[6],lambda = param_true[5],Tmax=data$Tmax_Exp_V,piTrunc= 0)
}

abs <- 1:(data$Tmax_Contr_V-1)
D <- diff(U(abs))
abs<- abs[-1]
D <- D*(D>0)
D <- D*(abs<60)
D <- D/sum(D)

M <- sum(D*abs)
V <- sum(D*abs^2)-M^2
  
EchaEMG <- remgCensored(10000,mu=data$k/param_true[6],sigma = sqrt(data$k)/param_true[6],lambda = param_true[5],Tmax=data$Tmax_Exp_V,piTrunc= 0)
mean(EchaEMG)

par(mfrow=c(2,1))
plot(abs,D,type='l')
hist(Echan,freq=FALSE,add=TRUE,col='red',nclass=50)
curve(V2,0,100,col='blue',add=TRUE)

 


plot(D*(D>0),type='l')
  
log_lik_marg(log_param_true,mydata)
log_lik_marg(log_param_true,mycompletedata)





paramsChains = list(nMCMC=20000,rho=rep(1,9))
paramsChains$nBurnin <- paramsChains$nMCMC/10
paramsChains$withQD  = TRUE

paramsChains$paramsToSample = c(1:7) #





param_init <- init_param(mycompletedata)
log_param_init <- from_param_to_log_param(param_init)
hyperparams= list(mean=rep(0,length(param_init)) ,sd = rep(3,length(param_init)))
for (i in 1:length(param_init)){
  if (!(i %in% paramsChains$paramsToSample)){
  log_param_init[i] <- log_param_true[i]
  }
}
 




####################################################
 
mySample <- my_mcmc_marg_onechain(mydata,log_param_init,hyperparams,paramsChains)


 


################################# Log param


par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/3),3))
for(m in paramsChains$paramsToSample){
  plot(mySample$myLogPostSample[,m],type='l',xlab='',main=names(mySample$myPostSample)[m],ylab=''); 
  abline(h=log_param_true[m],col='green',lwd=2)
}

par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/3),3))
for(m in paramsChains$paramsToSample){
    plot(density(mySample$myLogPostSample[,m]),main=names(mySample$myLogPostSample)[m])
    abline(v= log_param_true[m],col='green',lwd=2)
    curve(dnorm(x,hyperparams$mean[m],hyperparams$sd[m]),add=TRUE,col='orange',lwd=2)
}



######################## 
log_param_estim <- apply(mySample$myLogPostSample,2,mean)
abs <- seq(0,mydata$Tmax_Contr_R,len= 1000)
param_estim <- from_logparam_to_param(log_param_estim)

P <- computationToPlotCurves(abs, param_true,mycompletedata)
P$comp <- "True"
Pestim <- computationToPlotCurves(abs, param_estim,mycompletedata)
Pestim$comp <- "Estim"
allP <- rbind(P,Pestim)

levels(P$Curves)
Ptemp <- allP %>% filter(Curves %in% c("1.Natural Death",'2. Sum Arrival + reading', '3. min(ND,Reading)'))

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  

#-------------------- On data 


par(mfrow=c(1,2))
begin_title <- ifelse(param_true[9] > 0,'With','Without')
hist(mydata$T_Exp_V,freq=FALSE,nclass =min(mydata$n_Exp_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(mydata$T_Exp_V))
curve(dOurModelExp(x,param_true,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax_Exp_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_Exp_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Exp_R,freq=FALSE,nclass= min(mydata$n_Exp_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(mydata$T_Exp_R))
curve(dOurModelExp(x,param_true,mydata$k,mydata$kprime,color='red',Tmax=mydata$Tmax_Exp_R),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax_Exp_R),add=TRUE,col='orange',lwd=2,lty = 2)



