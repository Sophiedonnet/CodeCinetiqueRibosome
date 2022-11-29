plotPostInf <- function(mySample, paramChains,log_param_ref = NULL){
  
  par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
  for(m in paramsChains$paramsToSample){
    plot(mySample$myLogPostSample[,m],type='l',xlab='',main=names(mySample$myPostSample)[m],ylab=''); 
    if(!is.null(log_param_ref)){
      abline(h=log_param_ref[m],col='green',lwd=2)
    }
  }

  par(mfrow=c(ceiling(length(paramsChains$paramsToSample)/2),2))
  for(m in paramsChains$paramsToSample){
    plot(density(mySample$myLogPostSample[,m]),main=names(mySample$myLogPostSample)[m])
    if(!is.null(log_param_ref)){
      abline(v=log_param_ref[m],col='green',lwd=2)
    }
    curve(dnorm(x,hyperparams$mean[m],hyperparams$sd[m]),add=TRUE,col='orange',lwd=2) 
  }
}

#========================================================================
# plot Curves
#========================================================================

computationToPlotCurves <- function(x,param,data){
  
  
  param[is.na(param)]=0
  nbCurves <- 3*2
  
  P <- as.data.frame(rep(x,nbCurves)); names(P)='time'
  P$Color_Part <- rep(c('Green','Red'),each=nrow(P)/2)
  
  myDensity_V <- matrix(0,length(x),nbCurves/2)
  myDensity_V[,1] <- dExpCensored(x,lambda=param[1],Tmax = data$Tmax_Contr_V,piTrunc =param[2] )
  myDensity_V[,2] <- demgCensored(x,lambda=param[5],mu = data$k/param[6],sigma = sqrt(data$k)/param[6],piTrunc = param[7],Tmax =data$Tmax_Exp_V)
  myDensity_V[,3] <- dOurModelExp(x,param,k = data$k, kprime = data$kprime,color='green',Tmax = data$Tmax_Exp_V)
  
  
  myDensity_R <- matrix(0,length(x),nbCurves/2)
  myDensity_R[,1] <-  dExpCensored(x,lambda=param[3],Tmax = data$Tmax_Contr_R,piTrunc =param[4])
  myDensity_R[,2] <- demgCensored(x,param[5],mu = (data$k+data$kprime)/param[6],sigma = sqrt(data$k+data$kprime)/param[6],piTrunc = param[7],Tmax =data$Tmax_Exp_V)
  myDensity_R[,3] <- dOurModelExp(x,param,k = data$k, kprime = data$kprime,color='red',Tmax = data$Tmax_Exp_R)
  
  
  myCurves = rep(c('1.Natural Death',
                   '2. Sum Arrival + reading',
                   '3. min(ND,Reading)'),each=length(x))
  
  # myDensity_V[,4] <-  myDensity_R[,4]  <- dexp(x,param[8])
  # myDensity_V[,5] <- dOurModel(x,param,k,kprime,Tmax,color='green')
  # myDensity_R[,5] <- dOurModel(x,param,k,kprime,Tmax,color='red')
  # 
  # myCurves <- c(myCurves,rep(c('4. Quick Death','5. Final model'),each=length(x)))
  #}
  
  
  P$density <- c(c(myDensity_V),c(myDensity_R))
  P$Curves <- as.factor(rep(myCurves,2))
  P$Color_Part <- as.factor(P$Color_Part)
  return(P)
  
}
