rm(list=ls()) ### efface toutes les variables en mémoire


##########################################################################
########################### installation des packages nécéssaires
#####################################################################"#### 
#install.packages("ggplot2") # pour les graphes
#install.packages("dplyr") #  manipulation de tableau de données
#install.packages("emg") #

library(ggplot2)
library(dplyr)
library(emg)

##########################################################################
######################################### Loader nos functions 
##########################################################################
source('Functions/functionsSimulationRead.R')
source('Functions/functionsLikelihood_LogNormOnCtrData.R')
source('Functions/functionsMomentEstimations.R')
#source('Functions/functionsPlotPostInference.R')
source('Functions/functionsEMG.R')


##################################################
################ Load and format data
###################################################

# 

list_Ctr <-list.files('DataKarenComplete_WGE/RawData')
nbFiles <- 0
for (i in 1:length(list_Ctr)){
  which_SubFile <- paste0(list_Ctr[i],'/')
  where_data <- paste('DataKarenComplete_WGE/RawData',which_SubFile,sep='/') # où sont les données envoyées par Karen 
  nbFiles  <- nbFiles +length(list.files(where_data))
}
all_estim <- data.frame(file= rep(0,nbFiles),
                        ctr=rep(0,nbFiles), 
                        tIres=rep(0,nbFiles), 
                        scanningRate=rep(0,nbFiles)) 

for( j in 1:length(list_Ctr)){
  
  which_SubFile <- paste0(list_Ctr[j],'/')
  where_Rawdata <- paste('DataKarenComplete_WGE/RawData',which_SubFile,sep='/') # où sont les données envoyées par Karen
  names_data <- list.files(where_Rawdata) # liste les fichiers de Rdata
  nbData <- length(names_data)
  
  all_directories <- paste0(where_Rawdata,names_data)
  for (i in 1:nbData){
    names.i <- list.files(all_directories[i])
    print(names.i)
    print(all_directories[i])
    names.i.1 <- paste0(all_directories[i],'/',names.i)
    data.i <- lapply(names.i.1,function(f){read.table(f, quote="\"", comment.char="")$V1})
    names(data.i)<- sapply(names.i,function(s){substr(s,1,7)})
    # #  if(all_directories[i] %in% toCorrect){
    # #    data.i$Tctr_DN  = data.i$Tctr_DN*2
    # #    data.i$Texp_DN  = data.i$Texp_DN*2
    # #  }
    #
    Tmax <- lapply(data.i,max)
    names(Tmax) <- paste0(rep('Tmax_',length(names(data.i))),names(data.i))
    data.i <- c(data.i,Tmax)
    #
    #
    data.i$colorUP <-ifelse(substr(names_data[i],6,6)=='G','green','red')
    data.i$colorDN <-ifelse(substr(names_data[i],12,12)=='G','green','red')
    data.i$k <- as.numeric(substr(names_data[i],3,5) )
    data.i$kprime <- as.numeric(substr(names_data[i],9,11) ) - data.i$k
    data.i$name_data <- gsub("\\..*","",names_data[i])
    #dir.create(file.path(paste0('DataKarenComplete_WGE/FormattedData/',which_SubFile)))
               
    save(data.i,file = paste0('DataKarenComplete_WGE/FormattedData/',which_SubFile,names_data[i],'.Rdata' ))
  }
}

################################################"
############ Plot all data
##################################################
for( j in 1:length(list_Ctr)){
  
  which_SubFile <- paste0(list_Ctr[j],'/')
  where_data <- paste('DataKarenComplete_WGE/FormattedData',which_SubFile,sep='/') # où sont les données envoyées par Karen
  names_data <- list.files(where_data) # liste les fichiers de Rdata
  nbData <- length(names_data)
  
  for (i in 1:nbData){
    #
    print(paste0('Plot_data ',i,'--', which_SubFile))
    w <- paste0(where_data,names_data[i])
    where_plot <- paste0("DataKarenComplete_WGE/plotData/",which_SubFile,gsub("\\..*","",names_data[i]),"_data.png")
    load(w)
    p <- plot_FN_all_data(data.i) ### plot_FN_all_data est une fonction adhoc pour tracer les données. p est un objet ggplot2 .
    p
    ggsave(where_plot)
    #
  }
}

  # 
################################################"
############ Fit LogNormale on Ctr
##################################################
for( j in 1:length(list_Ctr)){
  
  which_SubFile <- paste0(list_Ctr[j],'/')
  where_data <- paste('DataKarenComplete_WGE/FormattedData',which_SubFile,sep='/') # où sont les données envoyées par Karen
  names_data <- list.files(where_data) # liste les fichiers de Rdata
  nbData <- length(names_data)
  
  paramCtrDN <- paramCtrUP <- vector("list",nbData)
  for (i in 1:nbData){
    #
    print(paste0('Fit Ctr ',i,'--', which_SubFile))

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

    where_plot <- paste0("DataKarenComplete_WGE/plotFit/",which_SubFile,gsub("\\..*","",names_data[i]),"_fitCtr.png")
    png(where_plot, width = 2000,height=2000)

    par(mfrow= c(2,2))
    plot(ecdf(data.i$Tctr_UP), main = paste(i,'UP ', bestmodelUP,sep=' . '))
    lines(abs,pCtrData(abs,paramCtr=paramCtrUP[[i]],Tmax = data.i$Tmax_Tctr_UP),col='red',lwd = 2)
    lines(abs,pi_ctrUP*pexpTrunc(abs,lambda_hat_UP,Tmax = data.i$Tmax_Tctr_UP),col='purple')

    hist(Tctr_UP[Tctr_UP<data.i$Tmax_Tctr_UP],freq= FALSE,nclass=50,main = 'UP')
    lines(abs,1/pi_ctrUP*dCtrData(abs,paramCtr=paramCtrUP[[i]],Tmax = data.i$Tmax_Tctr_UP-1),col='red',lwd = 2)
    lines(abs,dexpTrunc(abs,lambda_hat_UP,Tmax = data.i$Tmax_Tctr_UP),col='purple',lwd = 2)


    plot(ecdf(data.i$Tctr_DN), main = paste(i,'DN ', bestmodelDN,sep=' . '))
    lines(abs,pCtrData(abs,paramCtr=paramCtrDN[[i]],Tmax = data.i$Tmax_Tctr_DN),col='red',lwd = 2)
    lines(abs,pi_ctrDN*pexpTrunc(abs,lambda_hat_DN,Tmax = data.i$Tmax_Tctr_UP),col='purple',lwd = 2)

    hist(Tctr_DN[Tctr_DN<data.i$Tmax_Tctr_DN],freq= FALSE,nclass=50,main = 'DN')
    lines(abs,1/pi_ctrDN*dCtrData(abs,paramCtr=paramCtrDN[[i]],Tmax = data.i$Tmax_Tctr_DN-1),col='red',lwd = 2)
    lines(abs,dexpTrunc(abs,lambda_hat_DN,Tmax = data.i$Tmax_Tctr_DN),col='purple',lwd = 2)

    dev.off()

  }

  whereRes <- paste0('DataKarenComplete_WGE/resEstim/estimParamCtr_',gsub('/','.',which_SubFile),'Rdata')
  resCtr = list(paramCtrUp=paramCtrUP,paramCtrDN=paramCtrDN)
  save(resCtr,file = whereRes)
  # 
  #load(whereRes)
  #paramCtrUP <- resCtr$paramCtrUp
  #paramCtrDN <- resCtr$paramCtrDN
}

  
  
  
##############################################"
###############" Estim moment
#############################################
iter = 0 
for( j in 1:length(list_Ctr)){
  
  which_SubFile <- paste0(list_Ctr[j],'/')
  where_data <- paste('DataKarenComplete_WGE/FormattedData',which_SubFile,sep='/') # où sont les données envoyées par Karen
  names_data <- list.files(where_data) # liste les fichiers de Rdata
  nbData <- length(names_data)
  
  whereRes <- paste0('DataKarenComplete_WGE/resEstim/estimParamCtr_',gsub('/','.',which_SubFile),'Rdata')
  load(whereRes)
  paramCtrUP <- resCtr$paramCtrUp
  paramCtrDN <- resCtr$paramCtrDN
  
  
  #-------------------------------------------- Stoackage param estimés
  for (i in 1:length(names_data)){
    
    iter = iter + 1 
    print(iter)
    all_estim$file[iter] = names_data[i]
    all_estim$ctr[iter] = which_SubFile
    
    print(paste0('Fit Read ',i,'--', which_SubFile))
    
    load(paste0(where_data,names_data[i]))
    
    if(!is.null(data.i$Texp_UP)){
      Texp_UP  <- data.i$Texp_UP
      Tctr_UP  <- data.i$Tctr_UP
      paramCtr_UP <- paramCtrUP[[i]]
      FnCtr_UP <-   ecdf(Tctr_UP)
      Tmax_UP <- data.i$Tmax_Texp_UP
      echanRead_UP <- simu_readData(Texp_UP, Tctr_UP, Tmax_UP,paramCtr_UP,MC = 10000)
      
      #--------------- check on figures 
      #------------------ Ctr
      
      
      FestimExp <- function(t,paramCtr_UP,echanRead_UP){
        Fn <- ecdf(echanRead_UP)
        #V <- FRead(t,Tmax,piCtr,paramCtr)
        y <-  1-(1-pCtrData(t,paramCtr = paramCtr_UP,Tmax = Tmax_UP))*(1-Fn(t))
        
        return(y)
      }
    }
    
    if(!is.null(data.i$Texp_DN)){
      Texp_DN  <- data.i$Texp_DN
      Tctr_DN  <- data.i$Tctr_DN
      paramCtr_DN <- paramCtrDN[[i]]
      FnCtr_DN <-   ecdf(Tctr_DN)
      Tmax_DN <- data.i$Tmax_Texp_DN
      echanRead_DN <- simu_readData(Texp_DN, Tctr_DN, Tmax_DN,paramCtr_DN,MC = 10000)
    }
    #--------------------------------------------
    
    mean_readUP <- mean(echanRead_UP[echanRead_UP<Tmax_UP])
    mean_readDN <- mean(echanRead_DN[echanRead_DN<Tmax_DN])
    
    lambda_e <- data.i$kprime/(mean_readDN - mean_readUP)
    lambda_c <- 1/(mean_readUP - data.i$k/lambda_e)
    
    all_estim$tIres[iter] = 1/lambda_c
    all_estim$scanningRate[iter] = lambda_e
    
    
    #----------------------------------------
    where_plot <- paste0("DataKarenComplete_WGE/plotFit/",which_SubFile,gsub("\\..*","",names_data[i]),"_fitRead.png")
    png(where_plot, width = 2000,height=2000)
    
    par(mfrow=c(1,3))
    curve(FnCtr_UP,0,Tmax_UP,col=data.i$colorUP,lwd=2,main=paste('Ctr',names_data[i]))
    curve(pCtrData(x,paramCtr = paramCtr_UP,Tmax = Tmax_UP),add = TRUE,lty = 2,col=data.i$colorUP)
    curve(FnCtr_DN,0,Tmax_DN,col=data.i$colorDN,lwd =2 ,add = TRUE)
    curve(pCtrData(x,paramCtr = paramCtr_DN,Tmax = Tmax_DN),add = TRUE,lty = 2,col=data.i$colorDN)
    legend('topleft', legend=c("UP", "DN"),
           col=c(data.i$colorUP, data.i$colorDN),lty=1,lwd=2)
    
    plot(ecdf(Texp_UP),main='Exp',col=data.i$colorUP)
    curve(FestimExp(x,paramCtr_UP,echanRead_UP),add=TRUE,col=data.i$colorUP)
    plot(ecdf(Texp_DN),add = TRUE,col=data.i$colorDN)
    curve(FestimExp(x,paramCtr_DN,echanRead_DN),add=TRUE,col=data.i$colorDN)
    plot(ecdf(echanRead_UP),main=paste('Read. mu_c = ',round(1/lambda_c,1),'lambda_e = ', round(lambda_e,2)),col=data.i$colorUP)
    plot(ecdf(echanRead_DN),col=data.i$colorDN,add = TRUE)
    dev.off()
    #est_LogNormTrun_Read_UP <-  fitLogNormTruncCtr(echanRead_UP,Tmax =data.i$Tmax_Tctr_UP)
    
    #plot(density(echanRead_UP[echanRead_UP<data.i$Tmax_Texp_UP],bw = 10),main='Read',col=data.i$colorUP)
    #lines(density(echanRead_DN[echanRead_DN<data.i$Tmax_Texp_DN],bw = 10),col=data.i$colorDN)
    
    
  }
}

save(all_estim,file='DataKarenComplete_WGE/res_EstimMoment_allData.Rdata')

