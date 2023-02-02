rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('functionsMomentEstimations.R')
##################################################
################ Load and format data
###################################################

# 
where_data <- c('DataKarenComplete/RawData/')

toCorrect <- c("DataKarenComplete/RawData/UP049RDN205G_CtrA0","DataKarenComplete/RawData/UP180RDN205G_CtrA0")


names_data <- list.files(where_data)
Nbfiles <- length(names_data)
all_directories <- paste0(where_data,names_data)
for (i in 1:Nbfiles){
  names.i <- list.files(all_directories[i])
  print(names.i)
  print(all_directories[i])
  names.i.1 <- paste0(all_directories[i],'/',names.i)
  data.i <- lapply(names.i.1,function(f){read.table(f, quote="\"", comment.char="")$V1})
  names(data.i)<- sapply(names.i,function(s){substr(s,1,7)})
  if(all_directories[i] %in% toCorrect){
    data.i$Tctr_DN  = data.i$Tctr_DN*2
    data.i$Texp_DN  = data.i$Texp_DN*2
  }
  
  Tmax <- lapply(data.i,max)
  names(Tmax) <- paste0(rep('Tmax_',length(names(data.i))),names(data.i))
  data.i <- c(data.i,Tmax)


  data.i$colorUP <-ifelse(substr(names_data[i],6,6)=='G','green','red')
  data.i$colorDN <-ifelse(substr(names_data[i],12,12)=='G','green','red')
  data.i$k <- as.numeric(substr(names_data[i],3,5) )
  data.i$kprime <- as.numeric(substr(names_data[i],9,11) ) - data.i$k
  data.i$name_data <- gsub("\\..*","",names_data[i])
  save(data.i,file = paste0('DataKarenComplete/FormattedData/',names_data[i],'.Rdata' ))
}
# 
# ################################################"
# ############ Plot all data
# #################################################
where_data <- c('DataKarenComplete/FormattedData/')
names_data <- list.files(where_data)
for (i in 1:length(names_data)){

  print(i)
  w <- paste0(where_data,names_data[i])
  where_plot <- paste0("DataKarenComplete/plotData/",gsub("\\..*","",names_data[i]),"_plotdata.png")
  load(w)
  data.i
  p <- plot_FN_all_data(data.i)
  ggsave(where_plot)

}

##############################################"
###############" Estim moment
#############################################
where_data <- c('DataKarenComplete/FormattedData/')
names_data <- list.files(where_data)


#--------------------------------------------
param_estim_UP <- matrix(NA,length(names_data),7)
echan_exp_corr_UP <- vector("list", length(names_data))
param_estim_DN <- matrix(NA,length(names_data),7)
echan_exp_corr_DN <-vector("list", length(names_data))
rownames(param_estim_UP) <- rownames(param_estim_DN) <- gsub("\\..*","",names_data)
names(echan_exp_corr_DN) <-names(echan_exp_corr_UP)<- gsub("\\..*","",names_data)
colnames(param_estim_UP) <- paste0(c('lambda_ND','piTrunc_ND','lambda_c','mu_emg','sigma_emg','piTrunc_Read','lambda_e'),'_UP')
colnames(param_estim_DN) <- paste0(c('lambda_ND','piTrunc_ND','lambda_c','mu_emg','sigma_emg','piTrunc_Read','lambda_e'),'_DN')
for (i in 1:length(names_data)){

  print(names_data[i])
  load(paste0(where_data,names_data[i]))
  if(!is.null(data.i$Texp_UP)){
    resEstimUP <- estim_param_moment(data.i,'UP')
    param_estim_UP[i,]<- resEstimUP$param_estim
    echan_exp_corr_UP[[i]] <- resEstimUP$echan_exp_corr
  }

  if(!is.null(data.i$Texp_DN)){
    resEstimDN <- estim_param_moment(data.i,'DN')
    param_estim_DN[i,]<- resEstimDN$param_estim
    echan_exp_corr_DN[[i]] <- resEstimDN$echan_exp_corr
  }
}

save(param_estim_DN, param_estim_UP,  echan_exp_corr_UP,  echan_exp_corr_DN,file='DataKarenComplete/res_EstimMoment_allData.Rdata')
  
load(file='DataKarenComplete/res_EstimMoment_allData.Rdata') 
##############################################"
###############" PLOT fit Exp distri to  data Ctr
#############################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  
  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/dataCtr/PDF/",gsub("\\..*","",names_data[i]),"_plotFitCtr.png")
  load(paste0(where_data,names_data[i]))
  g <- plot_Fit_Ctr(data.i,param_estim_UP[i,],param_estim_DN[i,],which.curve = 'pdf')
  g
  ggsave(where_plot_fit_Ctr)
  
  
  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/dataCtr/Density/",gsub("\\..*","",names_data[i]),"_plotFitCtr.png")
  load(paste0(where_data,names_data[i]))
  g <- plot_Fit_Ctr(data.i,param_estim_UP[i,],param_estim_DN[i,],which.curve = 'density')
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
    curve(demg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
    legend('topright', legend=c("Density estim", "EMG fit"),col=c(data.i$colorUP, data.i$colorUP), lty=1:2, cex=0.8)

      }
  if(!is.null(data.i$Texp_DN)){
    hist(echan_exp_corr_DN[[i]],nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
    lines(density(echan_exp_corr_DN[[i]]),col = data.i$colorDN,lwd = 2)
    curve(demg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
    legend('topright', legend=c("Density estim", "EMG fit"),col=c(data.i$colorDN, data.i$colorDN), lty=1:2, cex=0.8)
    
  }
  dev.off()

  where_plot_fit_emg <- paste0("DataKarenComplete/plotFit/dataExpCorrected/PDF/",gsub("\\..*","",names_data[i]),"_plotFitEmg")
  
  png(file=paste0(where_plot_fit_emg,'_pdf.png'))
  par(mfrow=c(2,1))
  if(!is.null(data.i$Texp_UP)){
    plot(ecdf(echan_exp_corr_UP[[i]]),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. UP'))
    curve(pemg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
    legend('bottomright', legend=c("Empir pdf", "EMG fit"),col=c('black', data.i$colorUP), lty=1:2, cex=0.8)
    
  }
  if(!is.null(data.i$Texp_DN)){
    plot(ecdf(echan_exp_corr_DN[[i]]),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
    curve(pemg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
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
  png(file=paste0(where_plot_fit_exp,'.png'))
  par(mfrow=c(2,1))
  
  if(!is.null(data.i$Texp_UP)){
    plot(ecdf(data.i$Texp_UP),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg+Ctr corr to Exp data. UP'))
    x <- seq(1,data.i$Tmax_Texp_UP)
    FnTCtr <- ecdf(data.i$Tctr_UP) 
    Y <- 1-(1-FnTCtr(x))*(1-pemg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3]))
    Y <- Y*(1-mean(data.i$Texp_UP==data.i$Tmax_Texp_UP))
    lines(x,Y,col=data.i$colorUP,lty = 2,lwd=2)
    Y2 <-  1-(1-pExpCensored(x,lambda=param_estim_UP[i,1],piTrunc=param_estim_UP[i,2],Tmax = data.i$Tmax_Tctr_UP))*(1-pemg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3]))
    Y2 <- Y2*(1-mean(data.i$Texp_UP==data.i$Tmax_Texp_UP))
   lines(x,Y2,col=data.i$colorUP,lty = 3,lwd=2)
    legend('bottomright', legend=c("Empir pdf", "EMG fit + Fn Ctr corr","EMG fit + Exp Ctr corr"),col=c('black', data.i$colorUP,data.i$colorUP), lty=1:3, cex=0.8)
    
  }
  if(!is.null(data.i$Texp_DN)){
    plot(ecdf(data.i$Texp_DN),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg+Ctr corr to Exp data. DN'))
    x <- seq(1,data.i$Tmax_Texp_DN)
    FnTCtr <- ecdf(data.i$Tctr_DN) 
    Y <- 1-(1-FnTCtr(x))*(1-pemg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3]))
    Y <- Y*(1-mean(data.i$Texp_DN==data.i$Tmax_Texp_DN))
    lines(x,Y,col=data.i$colorDN,lty = 2,lwd=2)
    Y2 <-  1-(1-pExpCensored(x,lambda=param_estim_DN[i,1],piTrunc=param_estim_DN[i,2],Tmax = data.i$Tmax_Tctr_DN))*(1-pemg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3]))
    Y2 <- Y2*(1-mean(data.i$Texp_DN==data.i$Tmax_Texp_DN))
    lines(x,Y2,col=data.i$colorDN,lty = 3,lwd=2)
    legend('bottomright', legend=c("Empir pdf", "EMG fit + Fn Ctr corr","EMG fit + Exp Ctr corr"),col=c('black', data.i$colorDN,data.i$colorDN), lty=1:3, cex=0.8)
    
    
  }
  dev.off()

  ##############        hist 
  where_plot_fit_exp<- paste0("DataKarenComplete/plotFit/dataExp/Density/",gsub("\\..*","",names_data[i]),"_plotExpData")
  png(file=paste0(where_plot_fit_exp,'.png'))
  par(mfrow=c(2,1))
  
  if(!is.null(data.i$Texp_UP)){
    hist(data.i$Texp_UP,nclass=50,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg+Ctr corr to Exp data. UP'),freq = FALSE)
    x <- seq(1,data.i$Tmax_Texp_UP)
    Y2 <- dminExpExpplusGaussian(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda = param_estim_UP[i,3],piTrunc = param_estim_UP[i,6],lambda_ND = param_estim_UP[i,1],piTrunc_ND = param_estim_UP[i,2],Tmax=  data.i$Tmax_Tctr_UP,log = FALSE)
    lines(x,Y2,col=data.i$colorUP,lty = 3,lwd=2)
    legend('topright', legend=c("EMG fit + Exp Ctr corr"),col=c(data.i$colorUP), lty=2, cex=0.8)
    
  }
  if(!is.null(data.i$Texp_DN)){
    hist(data.i$Texp_DN,nclass=50,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg+Ctr corr to Exp data. DN'),freq = FALSE)
    x <- seq(1,data.i$Tmax_Texp_DN)
    Y2 <- dminExpExpplusGaussian(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda = param_estim_DN[i,3],piTrunc = param_estim_DN[i,6],lambda_ND = param_estim_DN[i,1],piTrunc_ND = param_estim_DN[i,2],Tmax=  data.i$Tmax_Tctr_DN,log = FALSE)
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

param_estim <- cbind(param_estim_UP[,-c(2,4,5,6)], param_estim_DN[,-c(2,4,5,6)])
param_estim <- 1/param_estim
colnames(param_estim) <-  c('mean_ND_UP','mean_c_UP','mean_e_UP','mean_ND_DN','mean_c_DN','mean_e_DN')
param_estim <- as.data.frame(param_estim)
write.csv(param_estim,file='DataKarenComplete/param_estim.csv',row.names  = TRUE)


