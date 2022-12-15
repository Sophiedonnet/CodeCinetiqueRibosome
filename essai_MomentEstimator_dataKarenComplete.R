rm(list=ls())
library(ggplot2)
library(dplyr)
library(emg)
source('functionsEMG.R')
source('functionsMomentEstimations.R')
##################################################
################ Load and format data
###################################################


where_data <- c('DataKarenComplete/RawData/')
names_data <- list.files(where_data)
Nbfiles <- length(names_data)
all_directories <- paste0(where_data,names_data)
for (i in 1:Nbfiles){
  
  
  names.i <- list.files(all_directories[i])
  print(all_directories[i])
  names.i.1 <- paste0(all_directories[i],'/',names.i)
  data.i <- lapply(names.i.1,function(f){read.table(f, quote="\"", comment.char="")$V1})
  names(data.i)<- sapply(names.i,function(s){substr(s,1,7)})
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

################################################"
############ Plot all data
#################################################
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
param_estim_UP <- matrix(NA,length(names_data),6)
echan_exp_corr_UP <- vector("list", length(names_data))
param_estim_DN <- matrix(NA,length(names_data),6)
echan_exp_corr_DN <-vector("list", length(names_data))
rownames(param_estim_UP) <- rownames(param_estim_DN) <- gsub("\\..*","",names_data)
names(echan_exp_corr_DN) <-names(echan_exp_corr_UP)<- gsub("\\..*","",names_data)

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

  
  
##############################################"
###############" PLOT fit  data Ctr
#############################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  
  where_plot_fit_Ctr <- paste0("DataKarenComplete/plotFit/",gsub("\\..*","",names_data[i]),"_plotFitCtr.png")
  load(paste0(where_data,names_data[i]))
  g <- plot_Fit_Ctr(data.i,param_estim_UP[i,],param_estim_DN[i,])
  g
  ggsave(where_plot_fit_Ctr)
  
 #g <- plot_Fit_Exp(data.i,gsub("\\..*","",names_data[i]),param_estim[i,])
  #ggsave(where_plot_fit)
  
}

##############################################"
###############" PLOT fit  emg data Exp corrected
#############################################
for (i in 1:length(names_data)){
  
  print(names_data[i])
  where_plot_fit_emg <- paste0("DataKarenComplete/plotFit/",gsub("\\..*","",names_data[i]),"_plotFitEmg.png")
  load(paste0(where_data,names_data[i]))
  
  par(mfrow=c(2,1))
  
  if(!is.null(data.i$Texp_UP)){
    hist(echan_exp_corr_UP[[i]],nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. UP'))
    lines(density(echan_exp_corr_UP[[i]]),col = data.i$colorUP,lwd = 2)
    curve(demg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
  }
  if(!is.null(data.i$Texp_DN)){
    hist(echan_exp_corr_DN[[i]],nclass = 100,freq = FALSE,main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
    lines(density(echan_exp_corr_DN[[i]]),col = data.i$colorDN,lwd = 2)
    curve(demg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
  }
  
  par(mfrow=c(2,1))
  
  if(!is.null(data.i$Texp_UP)){
    plot(ecdf(echan_exp_corr_UP[[i]]),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. UP'))
    curve(pemg(x,mu=param_estim_UP[i,4],sigma=param_estim_UP[i,5],lambda=param_estim_UP[i,3] ),add=TRUE,col=data.i$colorUP,lty = 2,lwd=2)
  }
  if(!is.null(data.i$Texp_DN)){
    plot(ecdf(echan_exp_corr_DN[[i]]),main=paste0(gsub("\\..*","",names_data[i]),'. Fit emg to data Exp corrected. DN'))
    curve(pemg(x,mu=param_estim_DN[i,4],sigma=param_estim_DN[i,5],lambda=param_estim_DN[i,3] ),add=TRUE,col=data.i$colorDN,lty = 2,lwd=2)
  }

}

############################### transfo param_estim
param_estim = param_estim[,-c(8,9)]
param_estim[,c(1,3,5,6)] <- 1/param_estim[,c(1,3,5,6)]
colnames_param_estim = rep(0,7)
colnames_param_estim[1] = 'Mean Time ND UP'
colnames_param_estim[3] = 'Mean Time ND DN'
colnames_param_estim[2] = 'pi TRUNC UP'
colnames_param_estim[4] = 'pi TRUNC DN'
colnames_param_estim[5] = '1/lambda_c'
colnames_param_estim[6] = '1/lambda_e'
colnames_param_estim[7] = 'pi TRUNC Reading'

colnames(param_estim) <- colnames_param_estim




