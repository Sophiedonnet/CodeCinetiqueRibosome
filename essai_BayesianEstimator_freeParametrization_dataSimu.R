rm(list=ls()) ### efface toutes les variables en mémoire


##########################################################################
########################### installation des packages nécéssaires
#####################################################################"#### 
#install.packages("ggplot2") # pour les graphes
#install.packages("dplyr") #  manipulation de tableau de données
#install.packages("emg") #


##########################################################################
######################################### Loader les packages
###########################################################################
library(ggplot2)
library(dplyr)
library(emg)
source('Functions/functionsEMG.R')
source('Functions/functionsMomentEstimations.R')
source('Functions/functionsPlotPostInference.R')
#source('Functions/myMCMC_marg_freeParametrization.R')
#source('Functions/functionsLikelihood_freeParametrization.R')
source('Functions/myMCMC_marg_constMu.R')
##################################################
################ Data simulation
###################################################

# 
where_data <- c('DataKarenComplete/FormattedData/') # où sont les données envoyées par Karen 
names_data <- list.files(where_data) # liste les fichiers de Rdata
Nbfiles <- length(names_data)
all_directories <- paste0(where_data,names_data)



load('DataKarenComplete/resEstim/res_EstimMoment_allData.Rdata')


d = 1; 
load(paste0(where_data,names_data[d]))
data.d <- data.i
data.sim <- data.d
  

param_sim <- rep(0,12)
names(param_sim) <- c("lambda_ND_UP", "piTrunc_ND_UP", 
                      "lambda_ND_DN" , "piTrunc_ND_DN", 
                      "lambda_c", 
                      "mu_e_UP" , "sigma_e_UP", "piTrunc_Read_UP", 
                      "mu_e_DN" , "sigma_e_DN", "piTrunc_Read_DN",
                      "lambda_e")

param_sim[c(1,2)] <- param_estim_UP_all[d,c(1,2)] # lambda_ND_UP ,piTrunc_ND_UP
param_sim[c(3,4)] <- param_estim_DN_all[d,c(1,2)] # # lambda_ND_DN ,piTrunc_ND_DN
param_sim[5] <- 0.5*(param_estim_UP_all[d,3] + param_estim_DN_all[d,3]) # lambda_c
param_sim[c(6,7,8)] <- param_estim_UP_all[d,c(4,5,6)]  # mu_e_UP, sigma_e_UP, piTrunc_Read_UP
param_sim[c(9,10,11)] <- param_estim_DN_all[d,c(4,5,6)]  # mu_e_DN, sigma_e_DN, piTrunc_Read_DN
param_sim[12] = 0.5*(data.d$k/param_sim[6]  + (data.d$k + data.d$kprime)/param_sim[9])
param_sim[6] <-data.d$k/param_sim[12]
param_sim[9] <- (data.d$k + data.d$kprime)/param_sim[12]



fact.ndata  = 1
data.sim$Texp_UP <- rOurModelExp_constMu(n = length(data.d$Texp_UP)*fact.ndata,param= param_sim,UPDN='UP',Tmax = data.d$Tmax_Texp_UP)
data.sim$Texp_DN <- rOurModelExp_constMu(n = length(data.d$Texp_DN)*fact.ndata,param= param_sim,UPDN='DN',Tmax = data.d$Tmax_Texp_DN)
data.sim$Tctr_UP <- rExpCensored(n = length(data.d$Tctr_UP)*fact.ndata,lambda = param_sim[1],Tmax = data.d$Tmax_Tctr_UP, piTrunc=param_sim[2])
data.sim$Tctr_DN <- rExpCensored(n = length(data.d$Tctr_DN)*fact.ndata,lambda = param_sim[3],Tmax = data.d$Tmax_Tctr_DN, piTrunc=param_sim[4])

log_param_sim <- from_param_to_log_param_constMu(param_sim)
############################################# 
resEstimUP <- estim_param_moment(data.sim,'UP')
resEstimDN <- estim_param_moment(data.sim,'DN')

resEstimUP <- estim_param_maxlik(data.sim,'UP')
resEstimDN <- estim_param_maxlik(data.sim,'DN')


#####################
#### 
mean(resEstimUP$echan_exp_corr)
mean(remgCensored(1000,param_sim[6],param_sim[7],param_sim[5],Tmax= 98, param_sim[8]))

param_estim_moment <- rep(0,12)
param_estim_moment[c(1,2)] <- resEstimUP$param_estim[c(1,2)]
param_estim_moment[c(3,4)] <- resEstimDN$param_estim[c(1,2)]
param_estim_moment[5] <-  0.5*(resEstimUP$param_estim[3] + resEstimDN$param_estim[3])
param_estim_moment[c(6,7,8)] <- resEstimUP$param_estim[c(4,5,6)]
param_estim_moment[c(9,10,11)] <- resEstimDN$param_estim[c(4,5,6)]
param_estim_moment[12] = 0.5*(data.d$k/param_estim_moment[6]  + (data.d$k + data.d$kprime)/param_estim_moment[9])
names(param_estim_moment) <- names(param_sim)

rbind(param_sim,param_estim_moment)

 

##############################################"
###############" MCMC Estimation starting from moment estimator
############################################ 
#--------------- run MCMC


log_param_init <- log_param_estim_moment <- from_param_to_log_param_constMu(param_estim_moment)
param_init <- param_estim_moment
nbParam <- length(log_param_init)
paramsChains <- list(nMCMC=1000,rho=rep(1,nbParam),nBurnin=1,paramsToSample=c(1:12)[-c(2,4,8,11)])
paramsChains$rho[c(12)] <- 1
paramsChains$rho[c(1,3)] <- 0.1


log_param_init[-paramsChains$paramsToSample] <- from_param_to_log_param_constMu(param_sim)[-paramsChains$paramsToSample]
param_init[-paramsChains$paramsToSample] <- param_sim[-paramsChains$paramsToSample]


hyperparams <- list(param1=rep(0,nbParam) ,param2 = rep(3,nbParam))
hyperparams$param1[5] = -3;
hyperparams$param2[5] = 1;

hyperparams$param1[12] = 2; 
hyperparams$param2[12] = 1; 

hyperparams$param1[c(2,4,8,11)] = 1; 
hyperparams$param2[c(2,4,8,11)] = 1; 

resMCMC <- my_mcmc_marg_constMu(data.sim,log_param_init,
                                               hyperparams = hyperparams,
                                               paramsChains = paramsChains)


thinning <- 1
burnin = 0
extr <- seq(burnin+1,paramsChains$nMCMC,by=thinning)
par(mfrow=c(4,3))
for (p in c(6,7,8,9,10,11,1,2,3,4,5,12)){
  U <- resMCMC$myLogPostSample[extr,p]
  plot(resMCMC$myLogPostSample[extr,p],type='l',main=names(log_param_sim)[p],ylab = '',xlab = 'iter',ylim=range(c(U,log_param_sim[p],log_param_estim_moment[p]))); 
  abline(h=log_param_sim[p],col='red',lwd=2)  
  abline(h=log_param_estim_moment[p],col='green',lwd=2,lty=2)  
  #abline(h=log_param_init[p],col='magenta',lwd=2,lty=3)  
  
}

par(mfrow=c(4,3))
for (p in c(6,7,8,9,10,11,1,2,3,4,5,12)){
  U <- resMCMC$myPostSample[extr,p]
  plot(resMCMC$myPostSample[extr,p],type='l',main=names(param_sim)[p],ylab = '',xlab = 'iter',ylim=range(c(U,param_sim[p],param_estim_moment[p]))); 
  abline(h=param_sim[p],col='red',lwd=2)  
  abline(h=param_estim_moment[p],col='green',lwd=2,lty=2)  
}



  
##############################################"
###############" PLOT fit Exp distri to  data Ctr
#############################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  
  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/dataCtr/PDF/",gsub("\\..*","",names_data[i]),"_plotFitCtr.png")
  load(paste0(where_data,names_data[i]))
  g <- plot_Fit_Ctr(data.i,param_estim_UP_all[i,],param_estim_DN_all[i,],which.curve = 'pdf')
  g
  ggsave(where_plot_fit_Ctr)
  
  
  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/dataCtr/Density/",gsub("\\..*","",names_data[i]),"_plotFitCtr.png")
  load(paste0(where_data,names_data[i]))
  g <- plot_Fit_Ctr(data.i,param_estim_UP_all[i,],param_estim_DN_all[i,],which.curve = 'density')
  g
  ggsave(where_plot_fit_Ctr)
  
  
}


#####################################################################################"
###############" PLOT fit  emg data Exp corrected
#######################################################################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  where_plot_fit_emg <- paste0("DataKarenComplete/plotFit/dataExpCorrected/Density/",gsub("\\..*","",names_data[i]),"_plotFitEmg")
  load(paste0(where_data,names_data[i]))
  
  png(file=paste0(where_plot_fit_emg,'_dens.png'))
  par(mfrow=c(2,1))
  if(!is.null(data.i$Texp_UP)){
    hist(echan_exp_corr_UP[[i]],nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. UP'))
    lines(density(echan_exp_corr_UP[[i]]),col = data.i$colorUP,lwd = 2)
    curve(demg(x,mu=param_estim_UP_all[i,4],sigma=param_estim_UP_all[i,5],lambda=param_estim_UP_all[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
    legend('topright', legend=c("Density estim", "EMG fit"),col=c(data.i$colorUP, data.i$colorUP), lty=1:2, cex=0.8)

      }
  if(!is.null(data.i$Texp_DN)){
    hist(echan_exp_corr_DN[[i]],nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
    lines(density(echan_exp_corr_DN[[i]]),col = data.i$colorDN,lwd = 2)
    curve(demg(x,mu=param_estim_DN_all[i,4],sigma=param_estim_DN_all[i,5],lambda=param_estim_DN_all[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
    legend('topright', legend=c("Density estim", "EMG fit"),col=c(data.i$colorDN, data.i$colorDN), lty=1:2, cex=0.8)
    
  }
  dev.off()

  where_plot_fit_emg <- paste0("DataKarenComplete/plotFit/dataExpCorrected/PDF/",gsub("\\..*","",names_data[i]),"_plotFitEmg")
  
  png(file=paste0(where_plot_fit_emg,'_pdf.png'))
  par(mfrow=c(2,1))
  if(!is.null(data.i$Texp_UP)){
    plot(ecdf(echan_exp_corr_UP[[i]]),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. UP'))
    curve(pemg(x,mu=param_estim_UP_all[i,4],sigma=param_estim_UP_all[i,5],lambda=param_estim_UP_all[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
    legend('bottomright', legend=c("Empir pdf", "EMG fit"),col=c('black', data.i$colorUP), lty=1:2, cex=0.8)
    
  }
  if(!is.null(data.i$Texp_DN)){
    plot(ecdf(echan_exp_corr_DN[[i]]),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
    curve(pemg(x,mu=param_estim_DN_all[i,4],sigma=param_estim_DN_all[i,5],lambda=param_estim_DN_all[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
    legend('bottomright', legend=c("Empir pdf", "EMG fit"),col=c('black', data.i$colorDN), lty=1:2, cex=0.8)
    
  }
  dev.off()

}


#################################################################################""
###############" PLOT fit  to Experi Data
####################################################################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  load(paste0(where_data,names_data[i]))
  
  where_plot_fit_exp<- paste0("DataKarenComplete/plotFit/dataExp/PDF/",gsub("\\..*","",names_data[i]),"_plotExpData")
  
  ##############         PDF 
  D <- c()
  if(!is.null(data.i$Texp_UP)){
    FnTExp <- ecdf(data.i$Texp_UP)
    FnTCtr <- ecdf(data.i$Tctr_UP)
    x <- seq(1,data.i$Tmax_Texp_UP)
    A <-  mean(data.i$Texp_UP==data.i$Tmax_Texp_UP)
    DUP <- as.data.frame(rep(x,3)); names(DUP) <- 't'
    YExp <- 1-pemg(x,mu=param_estim_UP_all[i,4],sigma=param_estim_UP_all[i,5],lambda=param_estim_UP_all[i,3])
    Y <- 1-(1-FnTCtr(x))*YExp
    Y <- Y*(1-A)
    Y2 <-  1-(1-pExpCensored(x,lambda=param_estim_UP_all[i,1],piTrunc=param_estim_UP_all[i,2],Tmax = data.i$Tmax_Tctr_UP))*YExp
    Y2 <- Y2*(1-A)
    DUP$pdf <- c(FnTExp(x),Y,Y2)
    DUP$curve <- as.factor(rep(c('Data','Fit Emg + empir','Fit Emg + Exp'),each = length(x)))
    DUP$color <- as.factor(data.i$colorUP)
    DUP$UPDN <- as.factor('UP')
    D <- rbind(D,DUP)
  }
  
  if(!is.null(data.i$Texp_DN)){
    FnTExp <- ecdf(data.i$Texp_DN)
    FnTCtr <- ecdf(data.i$Tctr_DN)
    x <- seq(1,data.i$Tmax_Texp_DN)
    A <- mean(data.i$Texp_DN==data.i$Tmax_Texp_DN)
    DDN <- as.data.frame(rep(x,3))
    names(DDN) <- 't'
    YExp <- 1-pemg(x,mu=param_estim_DN_all[i,4],sigma=param_estim_DN_all[i,5],lambda=param_estim_DN_all[i,3])
    Y <- 1-(1-FnTCtr(x))*YExp
    Y <- Y*(1-A)
    Y2 <-  1-(1-pExpCensored(x,lambda=param_estim_DN_all[i,1],piTrunc=param_estim_DN_all[i,2],Tmax = data.i$Tmax_Tctr_DN))*YExp
    Y2 <- Y2*(1-A)
    DDN$pdf <- c(FnTExp(x),Y,Y2)
    DDN$curve <- as.factor(rep(c('Data','Fit Emg + empir','Fit Emg + Exp'),each = length(x)))
    DDN$color <- as.factor(data.i$colorDN)
    DDN$UPDN <- as.factor('DN')
    D <- rbind(D,DDN)
  }
  
  p <- ggplot(D,aes(x=t,y=pdf,group=curve,col=color))+geom_line(aes(linetype=curve))+facet_grid(rows = vars(UPDN) )+ggtitle(data.i$name_data)
  p
  ggsave(paste0(where_plot_fit_exp,'.png'))
  
  ##############        hist 
  where_plot_fit_exp<- paste0("DataKarenComplete/plotFit/dataExp/Density/",gsub("\\..*","",names_data[i]),"_plotExpData")
  png(file=paste0(where_plot_fit_exp,'.png'))
  par(mfrow=c(2,1))
  
  if(!is.null(data.i$Texp_UP)){
    hist(data.i$Texp_UP,nclass=50,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg+Ctr corr to Exp data. UP'),freq = FALSE)
    x <- seq(1,data.i$Tmax_Texp_UP)
    Y2 <- dminExpExpplusGaussian(x,mu=param_estim_UP_all[i,4],sigma=param_estim_UP_all[i,5],lambda = param_estim_UP_all[i,3],piTrunc = param_estim_UP_all[i,6],lambda_ND = param_estim_UP_all[i,1],piTrunc_ND = param_estim_UP_all[i,2],Tmax=  data.i$Tmax_Tctr_UP,log = FALSE)
    lines(x,Y2,col=data.i$colorUP,lty = 3,lwd=2)
    legend('topright', legend=c("EMG fit + Exp Ctr corr"),col=c(data.i$colorUP), lty=2, cex=0.8)
    
  }
  if(!is.null(data.i$Texp_DN)){
    hist(data.i$Texp_DN,nclass=50,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg+Ctr corr to Exp data. DN'),freq = FALSE)
    x <- seq(1,data.i$Tmax_Texp_DN)
    Y2 <- dminExpExpplusGaussian(x,mu=param_estim_DN_all[i,4],sigma=param_estim_DN_all[i,5],lambda = param_estim_DN_all[i,3],piTrunc = param_estim_DN_all[i,6],lambda_ND = param_estim_DN_all[i,1],piTrunc_ND = param_estim_DN_all[i,2],Tmax=  data.i$Tmax_Tctr_DN,log = FALSE)
    lines(x,Y2,col=data.i$colorDN,lty = 3,lwd=2)
    legend('topright', legend=c("EMG fit + Exp Ctr corr"),col=c(data.i$colorDN), lty=2, cex=0.8)
    
    
    
  }
  dev.off()

    
  # 
  # par(mfrow=c(2,1))
  # 
  # if(!is.null(data.i$Texp_UP)){
  #   hist(data.it$Texp_UP,nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. UP'))
  #   lines(density(echan_exp_corr_UP[[i]]),col = data.i$colorUP,lwd = 2)
  #   curve(demg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
  #   legend('topright', legend=c("Density estim", "EMG fit"),col=c(data.i$colorUP, data.i$colorUP), lty=1:2, cex=0.8)
  #   
  # }
  # if(!is.null(data.i$Texp_DN)){
  #   hist(echan_exp_corr_DN[[i]],nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
  #   lines(density(echan_exp_corr_DN[[i]]),col = data.i$colorDN,lwd = 2)
  #   curve(demg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
  #   legend('topright', legend=c("Density estim", "EMG fit"),col=c(data.i$colorDN, data.i$colorDN), lty=1:2, cex=0.8)
  #   
  # }
  
 
  
}


############################### transfo param_estim

param_estim <- cbind(param_estim_UP_all[,-c(2,4,5,6)], param_estim_DN_all[,-c(2,4,5,6)])
param_estim <- 1/param_estim
colnames(param_estim) <-  c('mean_ND_UP','mean_c_UP','mean_e_UP','mean_ND_DN','mean_c_DN','mean_e_DN')
param_estim <- as.data.frame(param_estim)
write.csv(param_estim,file='DataKarenComplete/param_estim.csv',row.names  = TRUE)


