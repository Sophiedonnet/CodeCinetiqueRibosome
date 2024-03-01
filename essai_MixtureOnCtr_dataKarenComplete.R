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
source('Functions/functionsLikelihood_MixtureOnCtrData.R')
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
for (i in 1:length(names_data)){
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
Kmax = 1
paramCtrUP <- vector("list",length(names_data))
paramCtrDN <- vector("list",length(names_data))

for (i in 1:length(names_data)){
  # 
  print(i)
  w <- paste0(where_data,names_data[i])
  load(w)
  
  #-------------------
  resFitMixture <- vector("list",2)
  names(resFitMixture) <- c('UP','DN')
  #------------------- 
  
  Tctr_UP <- data.i$Tctr_UP[data.i$Tctr_UP < data.i$Tmax_Tctr_UP & data.i$Tctr_UP>0]
  paramCtrUP[[i]] <-  fitMixtureCtr(Tctr_UP,Tmax =data.i$Tmax_Tctr_UP ,Kmax = Kmax, select = TRUE)
  
  logLikExp <- sum(dexp(Tctr_UP,1/mean(Tctr_UP),log = TRUE))
  logLikLogNormal <-  paramCtrUP[[i]]$loglik
  if(logLikExp>logLikLogNormal){bestmodelUP = 'exp'}else{bestmodelUP='lognormal'}
  
  #------------------- 
  Tctr_DN <- data.i$Tctr_DN[data.i$Tctr_DN < data.i$Tmax_Tctr_DN & data.i$Tctr_DN>0]
  paramCtrDN[[i]] <-  fitMixtureCtr(Tctr_DN,Tmax =data.i$Tmax_Tctr_DN,Kmax = Kmax,select   = TRUE)
  
  logLikExp <- sum(dexp(Tctr_DN,1/mean(Tctr_DN),log = TRUE))
  logLikLogNormal <-  paramCtrDN[[i]]$loglik
  if(logLikExp>logLikLogNormal){bestmodelDN = 'exp'}else{bestmodelDN='lognormal'}
  
  
  abs = seq(0,data.i$Tmax_Tctr_UP,len = 1000)
  par(mfrow= c(2,2))
  
  plot(ecdf(Tctr_UP), main = paste(i,'UP ', bestmodelUP,sep=' . '))
  lines(abs,pCtrData(abs,param=paramCtrUP[[i]],Tmax = data.i$Tmax_Tctr_UP) ,col='red')
  lines(abs,pexp(abs,1/mean(Tctr_UP)),col='green')
  
  hist(Tctr_UP,freq= FALSE,nclass=50,main = 'UP')
  lines(abs,dCtrData(abs,param=paramCtrUP[[i]],Tmax = data.i$Tmax_Tctr_UP),col='red',add = TRUE)
  lines(abs,dexp(abs,1/mean(Tctr_UP)),col='green')

  plot(ecdf(Tctr_DN), main = paste(i,'DN ', bestmodelDN,sep=' . '))
  lines(abs,pCtrData(abs,param=paramCtrDN[[i]],Tmax = data.i$Tmax_Tctr_DN),col='red')
  lines(abs,pexp(abs,1/mean(Tctr_DN)),col='green')
  
  hist(Tctr_DN,freq= FALSE,nclass=50,main = 'DN')
  lines(abs,dCtrData(abs,param=paramCtrDN[[i]],Tmax = data.i$Tmax_Tctr_DN),col='red')
  lines(abs,dexp(abs,1/mean(Tctr_DN)),col='green')
}






##############################################"
###############" Estim moment
#############################################


#-------------------------------------------- Stoackage param estimés
param_estim_all <- matrix(NA,length(names_data),3)
rownames(param_estim_all) <- rownames(param_estim_all ) <- gsub("\\..*","",names_data)
colnames(param_estim_all) <- paste0(c('lambda_c','lambda_e','sigma'))
echan_exp_corr <- vector("list", length(names_data))
echan_exp_corr <-vector("list", length(names_data))
names(echan_exp_corr) 

for (i in 1:length(names_data)){

  print(names_data[i])

  load(paste0(where_data,names_data[i]))
  if(!is.null(data.i$Texp_UP)){
    
    Texp <- data.i$Texp_UP[data.i$Texp_UP< data.i$Tmax_Texp_UP]
    Tctr <- data.i$Tctr_UP[data.i$Tctr_UP< data.i$Tmax_Tctr_UP]
    paramCtr <- paramCtrUP[[i]]
    Tmax <- data.i$Tmax_Texp_UP
    hist(Texp,freq= FALSE,nclass = 100)
    hist(Tctr,freq= FALSE,nclass = 100,add = TRUE, col='red')
    
    
    resEstimUP <- estim_param_moment(Texp, Tmax,paramCtr)
    paramExp <- resEstimUP$paramExpMoment
    
    hist(Texp,freq = FALSE,nclass= 100)
    curve(dminCtrEMG(x,paramExp,paramCtr,Tmax),0,Tmax,col='red',add = TRUE)
    curve(demg(x,mu=paramExp[1],sigma = paramExp[2],lambda = paramExp[3]),0,Tmax,col='green',add = TRUE)
    
    logparamExp <- paramExp
    logparamExp[c(2,3)] <-log(paramExp[c(2,3)])
    resOptim <- optim(par = logparamExp, fn = loglik_Texp_mixture, paramCtr= paramCtr, Texp = Texp,Tmax = Tmax)
    paramExp_estim <-resOptim$par
    paramExp_estim[c(2,3)] <-exp(resOptim$par[c(2,3)])
    
    hist(Texp,freq = FALSE,nclass= 100)
    curve(dminCtrEMG(x,paramExp,paramCtr,Tmax),0,Tmax,col='red',add = TRUE)
    curve(demg(x,mu=paramExp[1],sigma = paramExp[2],lambda = paramExp[3]),0,Tmax,col='green',add = TRUE)
    curve(dminCtrEMG(x,paramExp_estim,paramCtr,Tmax),0,Tmax,col='red',lty = 2,add = TRUE)
    curve(demg(x,mu=paramExp[1],sigma = paramExp[2],lambda = paramExp[3]),0,Tmax,col='green',add = TRUE)
    curve(demg(x,mu=paramExp_estim[1],sigma = paramExp_estim[2],lambda = paramExp_estim[3]),lty = 2,Tmax,col='darkgreen',add = TRUE)
    
    
  }

  if(!is.null(data.i$Texp_DN)){
    resEstimDN <- estim_param_moment(data.i,'DN')
    #resEstimDN <- estim_param_maxlik(data.i,'DN')
    param_estim_DN_all[i,]<- resEstimDN$param_estim
    echan_exp_corr_DN[[i]] <- resEstimDN$echan_exp_corr
  }
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


