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
source('Functions/functionsLikelihood_LogNormOnCtrData.R')
source('Functions/functionsMomentEstimations.R')
source('Functions/functionsPlotPostInference.R')
source('Functions/functionsEMG.R')


##################################################
################ Load and format data
###################################################

# 
#where_data <- c('DataKarenComplete/RawData/') # où sont les données envoyées par Karen 

#toCorrect <- c("DataKarenComplete/RawData/UP049RDN205G_CtrA0","DataKarenComplete/RawData/UP180RDN205G_CtrA0") # 2 fichiers que j'ai repérés pour lesquels les temps interimages n'étaient pas égaux? 


# names_data <- list.files(where_data) # liste les fichiers de Rdata
# Nbfiles <- length(names_data)
# all_directories <- paste0(where_data,names_data)
# for (i in 1:Nbfiles){
#   names.i <- list.files(all_directories[i])
#   print(names.i)
#   print(all_directories[i])
#   names.i.1 <- paste0(all_directories[i],'/',names.i)
#   data.i <- lapply(names.i.1,function(f){read.table(f, quote="\"", comment.char="")$V1})
#   names(data.i)<- sapply(names.i,function(s){substr(s,1,7)})
# #  if(all_directories[i] %in% toCorrect){
# #    data.i$Tctr_DN  = data.i$Tctr_DN*2
# #    data.i$Texp_DN  = data.i$Texp_DN*2
# #  }
#   
#   Tmax <- lapply(data.i,max)
#   names(Tmax) <- paste0(rep('Tmax_',length(names(data.i))),names(data.i))
#   data.i <- c(data.i,Tmax)
# 
# 
#   data.i$colorUP <-ifelse(substr(names_data[i],6,6)=='G','green','red')
#   data.i$colorDN <-ifelse(substr(names_data[i],12,12)=='G','green','red')
#   data.i$k <- as.numeric(substr(names_data[i],3,5) )
#   data.i$kprime <- as.numeric(substr(names_data[i],9,11) ) - data.i$k
#   data.i$name_data <- gsub("\\..*","",names_data[i])
#   save(data.i,file = paste0('DataKarenComplete/FormattedData/',names_data[i],'.Rdata' ))
# 
# 
################################################"
############ Plot all data
##################################################
where_data <- c('DataKarenComplete/FormattedData/')
names_data <- list.files(where_data)
nbData <- 1 # length(names_data)
for (i in 1:nbData){
# 
   print(i)
   w <- paste0(where_data,names_data[i])
#   where_plot <- paste0("DataKarenComplete/plotData/",gsub("\\..*","",names_data[i]),"_data.png")
   load(w)
   p <- plot_FN_all_data(data.i) ### plot_FN_all_data est une fonction adhoc pour tracer les données. p est un objet ggplot2 . 
   p
   #   ggsave(where_plot)
#  
}

################################################"
############ Fit Mixture on Ctr
##################################################
where_data <- c('DataKarenComplete/FormattedData/')
names_data <- list.files(where_data)

paramCtrDN <- paramCtrUP <- vector("list",nbData)


for (i in 1:nbData){
  # 
  print(i)
  w <- paste0(where_data,names_data[i])
  load(w)
  
  #------------------- 
  
  Tctr_UP <- data.i$Tctr_UP
  res_LogNormTrun_UP <-  fitLogNormTruncCtr(Tctr_UP,Tmax =data.i$Tmax_Tctr_UP)
  paramCtrUP[[i]] <- res_LogNormTrun_UP$paramCtr
  res_ExpoTrunc_UP <-  fitExpoTruncCtr(Tctr_UP,Tmax =data.i$Tmax_Tctr_UP)
  lambda_hat_UP <- res_ExpoTrunc_UP$lambda_hat
  logLikExp <- res_ExpoTrunc_UP$loglik
  logLikLogNormal <-  res_LogNormTrun_UP$loglik
  if(logLikExp>logLikLogNormal){bestmodelUP = 'exp'}else{bestmodelUP='lognormal'}
  print(bestmodelUP)
  #------------------- 
  Tctr_DN <- data.i$Tctr_DN
  res_LogNormTrun_DN <-  fitLogNormTruncCtr(Tctr_DN,Tmax =data.i$Tmax_Tctr_DN)
  paramCtrDN[[i]] <- res_LogNormTrun_DN$paramCtr
  res_ExpoTrunc_DN <-  fitExpoTruncCtr(Tctr_DN,Tmax =data.i$Tmax_Tctr_DN)
  lambda_hat_DN <- res_ExpoTrunc_DN$lambda_hat
  logLikExp <- res_ExpoTrunc_DN$loglik
  logLikLogNormal <-  res_LogNormTrun_DN$loglik
  if(logLikExp>logLikLogNormal){bestmodelDN = 'exp'}else{bestmodelDN='lognormal'}
  print(bestmodelUP)
  
  
  abs = seq(0,data.i$Tmax_Tctr_UP,len = 1000)
  pi_ctrUP <- (1-mean(data.i$Tctr_UP==data.i$Tmax_Tctr_UP))
  pi_ctrDN <- (1-mean(data.i$Tctr_DN==data.i$Tmax_Tctr_DN))
  
  par(mfrow= c(2,2))
  plot(ecdf(data.i$Tctr_UP), main = paste(i,'UP ', bestmodelUP,sep=' . '))
  lines(abs,pCtrData(abs,paramCtr=paramCtrUP[[i]],Tmax = data.i$Tmax_Tctr_UP),col='red',lwd = 2)
  lines(abs,pi_ctrUP*pexpTrunc(abs,lambda_hat_UP,Tmax = data.i$Tmax_Tctr_UP),col='purple')
  
  
  
  hist(Tctr_UP,freq= FALSE,nclass=50,main = 'UP')
  lines(abs,dCtrData(abs,paramCtr=paramCtrUP[[i]],Tmax = data.i$Tmax_Tctr_UP),col='red',lwd = 2)
  lines(abs,pi_ctrUP*dexpTrunc(abs,lambda_hat_UP,Tmax = data.i$Tmax_Tctr_UP),col='purple',lwd = 2)

  
  plot(ecdf(data.i$Tctr_DN), main = paste(i,'DN ', bestmodelDN,sep=' . '))
  lines(abs,pCtrData(abs,paramCtr=paramCtrDN[[i]],Tmax = data.i$Tmax_Tctr_DN),col='red',lwd = 2)
  lines(abs,pi_ctrDN*pexpTrunc(abs,lambda_hat_DN,Tmax = data.i$Tmax_Tctr_UP),col='purple',lwd = 2)
  
  hist(Tctr_DN,freq= FALSE,nclass=50,main = 'DN')
  lines(abs,dCtrData(abs,paramCtr=paramCtrDN[[i]],Tmax = data.i$Tmax_Tctr_DN),col='red',lwd = 2)
  lines(abs,pi_ctrDN*dexpTrunc(abs,lambda_hat_DN,Tmax = data.i$Tmax_Tctr_DN),col='purple',lwd = 2)

  }






##############################################"
###############" Estim moment
#############################################


for (i in 1:ndData){

  print(c(i,names_data[i]))

  load(paste0(where_data,names_data[i]))
  
  # pi_ctrUP <- (1-mean(data.i$Tctr_UP==data.i$Tmax_Tctr_UP))
  # pi_ctrDN <- (1-mean(data.i$Tctr_DN==data.i$Tmax_Tctr_DN))
  # pi_expUP <- (1-mean(data.i$Texp_UP==data.i$Tmax_Texp_UP))
  # pi_expDN <- (1-mean(data.i$Texp_DN==data.i$Tmax_Texp_DN))
  # 
  # par(mfrow=c(2,2))
  # hist(data.i$Texp_UP,freq = FALSE,nclass= 50); abline(v = mean(data.i$Texp_UP))
  # hist(data.i$Tctr_UP,freq = FALSE,nclass= 50,add = TRUE,col='green'); abline(v = mean(data.i$Tctr_UP),col='red')
  # plot(ecdf(data.i$Texp_UP)); abline(v = mean(data.i$Texp_UP))
  # plot(ecdf(data.i$Tctr_UP),add = TRUE,col='green'); abline(v = mean(data.i$Tctr_UP),col='green')
  # hist(data.i$Texp_DN,freq = FALSE,nclass= 50); abline(v = mean(data.i$Texp_UP))
  # hist(data.i$Tctr_DN,freq = FALSE,nclass= 50,add = TRUE,col='green'); abline(v = mean(data.i$Tctr_UP),col='red')
  # plot(ecdf(data.i$Texp_DN)); abline(v = mean(data.i$Texp_UP))
  # plot(ecdf(data.i$Tctr_DN),add = TRUE,col='green'); abline(v = mean(data.i$Tctr_UP),col='green')
  # 
  
  if(!is.null(data.i$Texp_UP)){
    Texp  <- data.i$Texp_UP
    Tctr  <- data.i$Tctr_UP
    paramCtr <- paramCtrUP[[i]]
    Tmax <- data.i$Tmax_Texp_UP
    resEstim_UP <- estimProc_Up_or_Dn(Texp,Tctr,Tmax, paramCtr)
  }

  if(!is.null(data.i$Texp_DN)){
    Texp_DN <- data.i$Texp_DN#[data.i$Texp_DN< data.i$Tmax_Texp_DN]
    Tctr_DN  <- data.i$Tctr_DN
    paramCtr_DN <- paramCtrDN[[i]]
    Tmax_DN <- data.i$Tmax_Texp_DN
    resEstim_DN <- estimProc_Up_or_Dn(Texp_DN,Tctr_DN,Tmax_DN, paramCtr_DN)
  }
}

  par(mfrow = c(2,2))
  plot(ecdf(Texp_UP[Texp_UP<Tmax_UP]),main = 'UP')
  x <- seq(0,98,by=0.1)
  lines(x, pExpData(x,resEstim_UP$moment,paramCtr_UP,Tmax_UP),col='blue',lwd= 2)
  lines(x, pExpData(x,resEstim_UP$loglik,paramCtr_UP,Tmax_UP),col='red',lwd= 2,lty = 2)
  lines(x, pExpData(x,resEstim_UP$fitFn,paramCtr_UP,Tmax_UP),col='orange',lwd= 2,lty = 3)
  plot(ecdf(Tctr_UP),add = TRUE,col='green')
  curve(pi_ctrUP*pCtrData(x,paramCtr_UP,Tmax_UP),col='green',lwd= 2,lty = 3,add = TRUE)

  curve(dRead(x,resEstim_UP$fitFn),0,Tmax,col='orange',main='Reading dist',lwd = 2,ylab='')
  curve(dRead(x,resEstim_UP$loglik),0,Tmax,col='red',lwd= 2,add = TRUE)
  #curve(dRead(x,resEstim_UP$moment),0,Tmax,col='blue',lwd= 2,add = TRUE)
  abline(v = (resEstim_UP$loglik[1] + 1/resEstim_UP$loglik[3]),col='red')
  abline(v = (resEstim_UP$fitFn[1] + 1/resEstim_UP$fitFn[3]),col='orange')
  #abline(v = (resEstim_UP$moment[1] + 1/resEstim_UP$moment[3]),col='blue')
  
  plot(ecdf(Texp_DN))
  #curve(pi_expUP*pminCtrEMG(x,resEstim_DN$moment,paramCtr_DN,Tmax_DN),0,Tmax,col='blue',lwd= 2,add = TRUE)
  curve(pi_expDN*pminCtrEMG(x,resEstim_DN$loglik,paramCtr_DN,Tmax_DN),0,Tmax,col='red',lwd= 2,lty = 2,add = TRUE)
  curve(pi_expDN*pminCtrEMG(x,resEstim_DN$fitFn,paramCtr_DN,Tmax_DN),0,Tmax,col='orange',lwd= 2,lty = 3,add = TRUE)
  plot(ecdf(Tctr_DN),add = TRUE,col='green')
  curve(pi_ctrDN*pCtrData(x,paramCtr_DN,Tmax_DN),col='green',lwd= 2,lty = 3,add = TRUE)
  
  curve(dRead(x,resEstim_DN$fitFn),0,Tmax,col='orange',main='Reading dist',lwd = 2,ylab='')
  curve(dRead(x,resEstim_DN$loglik),0,Tmax,col='red',lwd= 2,add = TRUE)
  #curve(dRead(x,resEstim_DN$moment),0,Tmax,col='blue',lwd= 2,add = TRUE)
  abline(v = (resEstim_DN$loglik[1] + 1/resEstim_DN$loglik[3]),col='red')
  abline(v = (resEstim_DN$fitFn[1] + 1/resEstim_DN$fitFn[3]),col='orange')
  #abline(v = (resEstim_DN$moment[1] + 1/resEstim_DN$moment[3]),col='blue')
  
}

save(param_estim_DN_all, param_estim_UP_all,  echan_exp_corr_UP,  echan_exp_corr_DN,file='DataKarenComplete/res_EstimMoment_allData.Rdata')
  
load(file='DataKarenComplete/res_EstimMoment_allData.Rdata') 
##############################################"
###############" PLOT fit Exp distri to  data Ctr
#############################################
for (i in 1:length(names_data)){

  print(names_data[i])

  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/MomentEstimation/dataCtr/",gsub("\\..*","",names_data[i]),"_Ctr_PDF.png")
  load(paste0(where_data,names_data[i]))
  g <- plot_Fit_Ctr(data.i,param_estim_UP_all[i,],param_estim_DN_all[i,],which.curve = 'pdf')
  g
  ggsave(where_plot_fit_Ctr)


  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/MomentEstimation/dataCtr/",gsub("\\..*","",names_data[i]),"_Ctr_density.png")
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
  where_plot_fit_emg <- paste0("DataKarenComplete/plotFit/MomentEstimation/dataExpCorrected/",gsub("\\..*","",names_data[i]),"_FitEmg_density")
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

  where_plot_fit_emg <- paste0("DataKarenComplete/plotFit/MomentEstimation/dataExpCorrected/",gsub("\\..*","",names_data[i]),"_FitEmg_pdf")
  
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
###############" PLOT fit  to Experimental Data
####################################################################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  load(paste0(where_data,names_data[i]))
  
  where_plot_fit_exp<- paste0("DataKarenComplete/plotFit/MomentEstimation/dataExp/",gsub("\\..*","",names_data[i]),"_ExpData_pdf")
  
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
  where_plot_fit_exp<- paste0("DataKarenComplete/plotFit/MomentEstimation/dataExp/",gsub("\\..*","",names_data[i]),"_ExpData_density")
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
write.csv(param_estim,file='DataKarenComplete/res_EstimMoment.csv',row.names  = TRUE)


