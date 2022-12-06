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

mydata$T_Exp_R <- T_Exp_R*2
mydata$T_Exp_V <- T_Exp_V


###########################################################"" 
data <- mydata
FNV <- ecdf(data$T_Exp_V)
param <- rep(0,9)
#------- data Contr V
L_V <- sort(1/seq(1,500,by=0.01))
D_V <-  espExpCensored(L_V,Tmax = data$Tmax_Contr_V)
param[1] <- L_V[which.min(abs(D_V-mean(data$T_Contr_V[data$T_Contr_V<data$Tmax_Contr_V])))]
param[2] <- mean(data$T_Contr_V == data$Tmax_Contr_V) #pi_trunc_ND_V

#------- data Contr R
L_R <- sort(1/seq(1,500,by=0.01))
D_R <-  espExpCensored(L_R,Tmax =data$Tmax_Contr_R)
param[3] <- L_R[which.min(abs(D_R-mean(data$T_Contr_R[data$T_Contr_R<data$Tmax_Contr_R])))]
param[4] <- mean(data$T_Contr_R == data$Tmax_Contr_R) #pi_trunc_ND_R

#----------- data Exp V
FNV <- ecdf(data$T_Exp_V)
UV <- function(x){
  1-(1-FNV(x))/(1-pExpCensored(x,lambda=param[1],Tmax =data$Tmax_Contr_V,piTrunc = param[2]))
}
plot(UV,0,100)
abs <- 1:(data$Tmax_Contr_V-1)
DV <- c(0,diff(UV(abs)))
DV <- DV*(DV>0)
DV <- DV/sum(DV)
EchanV <- sample(abs,10000,prob = DV,replace=TRUE)
TUP = 60
theta_hat_V <- estim_param_emg(EchanV[EchanV<TUP])
plot(density(EchanV[EchanV<data$Tmax_Exp_V]),main='Exp green corrected by ND')
curve(demg(x,theta_hat_V[1],theta_hat_V[2],theta_hat_V[3]),add=TRUE,col='red')


#----------- data Exp R
FNR <- ecdf(data$T_Exp_R)
UR <- function(x){
  1-(1-FNR(x))/(1-pExpCensored(x,lambda=param[3],Tmax =data$Tmax_Contr_R,piTrunc = param[4]))
}
abs <- 1:(data$Tmax_Contr_R-1)
DR <- c(0,diff(UR(abs)))
DR <- DR*(DR>0)
DR <- DR/sum(DR)
EchanR <- sample(abs,100000,prob = DR,replace=TRUE)
theta_hat_R <- estim_param_emg(EchanR[EchanR<TUP])
thetagamma_hat_R <- estim_param_Gamma(EchanR)

par(mfrow=c(1,1))
plot(density(EchanR),lwd=2)
curve(demg(x,lambda=theta_hat_V[3],mu=theta_hat_V[1]*(data$k+data$kprime)/data$k,sigma = theta_hat_V[2]*sqrt((data$k+data$kprime)/data$k)),add=TRUE,col='green')
curve(demg(x,theta_hat_R[1],theta_hat_R[2],theta_hat_R[3]),add=TRUE,col='red')
curve(dgamma(x,thetagamma_hat_R[1],thetagamma_hat_R[2]),add=TRUE,col='magenta')
curve(dexp(x,rate = 1/mean(EchanR)),add=TRUE,col='blue')
curve(dlnorm(x,mean(log(EchanR)),sd(log(EchanR))),add=TRUE,col='orange')



WR <- function(x){
  1-(1-FNR(x))/(1-pemg(x,lambda=theta_hat_V[3],mu=theta_hat_V[1]*(data$k+data$kprime)/data$k,sigma = theta_hat_V[2]*sqrt((data$k+data$kprime)/data$k)))
}
plot(1:100,WR(1:100),type='l')
DR <- c(0,diff(WR(abs)))
DR <- DR*(DR>0)
DR <- DR/sum(DR)
EchanR_ND <- sample(abs,100000,prob = DR,replace=TRUE)
par(mfrow=c(1,1))
plot(density(EchanR_ND))





¨

TUP2 = 100
FNR2 <- ecdf(EchanR[EchanR<TUP2])
UR2 <- function(x){
  1-(1-FNR2(x))/(1-)
}
abs <- 1:(TUP2-1)

plot(abs,UR2(abs),type='l')
DR2 <- c(0,diff(UR2(abs)))
DR2 <- DR*(DR2>0)
DR2 <- DR2/sum(DR2)
EchanR_SD <- sample(abs,100000,prob = DR2,replace=TRUE)
par(mfrow=c(1,1))
plot(density(EchanR_SD))

lambda_c <- 0.5*(theta_hat_R[3]+theta_hat_V[3])
1/lambda_c



data$k/theta_hat_V[1]
(data$k+data$kprime)/theta_hat_R[1]



param_init <- init_param(data)
log_param_init <- from_param_to_log_param(param_init)

hyperparams= list(mean=rep(0,length(param_init)) ,sd = rep(3,length(param_init)))



par(mfrow=c(1,2))
X <- mydata$T_Contr_V[mydata$T_Contr_V<max(mydata$T_Contr_V)]
hist(X,freq = FALSE,nclass =50,main="Natural Death. GREEN")  
lines(density(X))
curve(dExpCensored(x,lambda = param_init[1],piTrunc = 0,Tmax = mydata$Tmax_Contr_V),add=TRUE,col='green',lwd=2,lty = 1)

X <- mydata$T_Contr_R[mydata$T_Contr_R<max(mydata$T_Contr_R)]
hist(X,freq = FALSE,nclass =50,main="Natural Death. Red")  
lines(density(X))
curve(dExpCensored(x,lambda = param_init[3],piTrunc = 0,Tmax = mydata$Tmax_Contr_R),add=TRUE,col='red',lwd=2,lty = 1)


########################################################################
################################# ESTIM Exp Data 
##########################################################################




abs <- seq(0,max(c(myCompleteData$T_Exp_R,myCompleteData$T_Exp_V)),len= 1000)


P <- computationToPlotCurves(abs, param_init,mydata)
P$comp <- "Init"
levels(P$Curves)
Ptemp <- P #%>% filter(Curves %in% c("1.Natural Death"))

g <- ggplot(data = Ptemp, aes(x=time,y=density,colour = Curves))+ geom_line(aes(colour=Curves, linetype=comp))#+  geom_point(aes(shape=comp)) 
g + facet_wrap(~Color_Part)  
#-------------------- On data  Exp
par(mfrow=c(1,2))
begin_title <- ifelse(param_init[9] > 0,'With','Without')
hist(myCompleteData$T_Exp_V,freq=FALSE,nclass =50, main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(myCompleteData$T_Exp_V))
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax_Exp_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dminExpExpplusGaussian(x,theta_hat_V[1],theta_hat_V[2],theta_hat_V[3],param_init[7],param_init[1],param_init[2],Tmax= Inf,log = FALSE),add=TRUE,col='green',lwd=2,lty = 1)

#curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_Exp_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(myCompleteData$T_Exp_R,freq=FALSE,nclass= 50, main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(myCompleteData$T_Exp_R))
curve(dminExpExpplusGaussian(x,theta_hat_R[1],theta_hat_R[2],theta_hat_R[3],param_init[7],param_init[3],param_init[4],Tmax= Inf,log = FALSE),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='red',Tmax=mydata$Tmax_Exp_R),add=TRUE,col='red',lwd=2,lty = 1)





hist(mydata$T_Exp_V,freq=FALSE,nclass =min(mydata$n_Exp_V/10,100), main=paste(begin_title,"Quick Death. Green",sep=' '))  
lines(density(mydata$T_Exp_V))
curve(dOurModelExp(x,param_init,mydata$k,mydata$kprime,color='green',Tmax=mydata$Tmax_Exp_V),add=TRUE,col='green',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='green',Tmax = mydata$Tmax_Exp_V),add=TRUE,col='orange',lwd=2,lty = 2)

hist(mydata$T_Exp_R,freq=FALSE,nclass= min(mydata$n_Exp_R/10,100), main=paste(begin_title,"Quick Death. Red",sep=' '))  
lines(density(mydata$T_Exp_R))
curve(dOurModelExp(x,init,mydata$k,mydata$kprime,color='red',Tmax=mydata$Tmax_Exp_R),add=TRUE,col='red',lwd=2,lty = 1)
curve(dOurModelExp(x,param_estim,mydata$k,mydata$kprime,color='red',Tmax = mydata$Tmax_Exp_R),add=TRUE,col='orange',lwd=2,lty = 2)

