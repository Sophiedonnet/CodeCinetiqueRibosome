
#========================================================================
# plot Chains
#========================================================================
plotChains <- function(resMCMC, paramsChains, log_param_estim_moment, log_param_sim = NULL){
  

  extr <- seq(paramsChains$nBurnin+1,paramsChains$nMCMC,by= paramsChains$thinning)
  
  param_estim_moment <- from_logparam_to_param_constMu(log_param_estim_moment)
  if(!is.null(log_param_sim)){param_sim <- from_logparam_to_param_constMu(log_param_sim)}
  
  
  par(mfrow=c(3,2))
  for (p in c(1,3,5,12,7,10)){
    U <- resMCMC$myLogPostSample[extr,p]
    plot(U,type='l',main=names(log_param_estim_moment)[p],ylab = '',xlab = 'iter'); 
    abline(h=log_param_estim_moment[p],col='green',lwd=2,lty=2)  
    if(!is.null(log_param_sim)){
      abline(h=log_param_sim[p],col='red',lwd=2)
    }
  }
  
  par(mfrow=c(3,2))
  for (p in c(1,3,5,12,7,10)){
    U <- resMCMC$myPostSample[extr,p]
    plot(U,type='l',main=names(param_estim_moment)[p],ylab = '',xlab = 'iter',ylim=range(c(U,param_estim_moment[p]))); 
    abline(h=param_estim_moment[p],col='green',lwd=2,lty=2)  
    if(!is.null(log_param_sim)){
      abline(h=param_sim[p],col='red',lwd=2)
    }
    
  }
}

#========================================================================
# plot post distrib
#========================================================================
plotPosteriorDistr <- function(resMCMC, paramsChains, hyperparams, log_param_estim_moment, log_param_sim = NULL){
  
  extr <- seq(paramsChains$nBurnin+1,paramsChains$nMCMC,by= paramsChains$thinning)
  
  par(mfrow=c(3,2))
  for (p in c(1,3,5,12,7,10)){
    U <- resMCMC$myLogPostSample[extr,p]
    plot(density(U),ylab = '',main=names(log_param_estim_moment)[p]); 
    curve(dnorm(x,hyperparams$param1[p],hyperparams$param2[p]),add=TRUE,col='blue',lty=2)
   
    abline(v=log_param_estim_moment[p],col='green',lwd=2,lty=2) 
    if(!is.null(log_param_sim)){
      abline(v=log_param_sim[p],col='red',lwd=2)  
    }
  }
  par(mfrow=c(3,2))
  for (p in c(1,3,5,12,7,10)){
    U <- resMCMC$myPostSample[extr,p]
    plot(density(U),ylab = '',main=names(param_estim_moment)[p]); 
  }
  
}

#========================================================================
# plot post curves with Confidence bandwith
#========================================================================
plotCurves <- function(resMCMC, theData,paramsChains, hyperparams, param_estim_moment){
  
  extr <- seq(paramsChains$nBurnin+1,paramsChains$nMCMC,by= paramsChains$thinning)
  
  
  D_CTR_UP  <- density(theData$Tctr_UP[theData$Tctr_UP< theData$Tmax_Tctr_UP])
  m_CTR_UP <- length(D_CTR_UP$x)
  D_CTR_DN  <- density(theData$Tctr_DN[theData$Tctr_DN< theData$Tmax_Tctr_DN])
  m_CTR_DN <- length(D_CTR_DN$x)
  D_EXP_UP  <- density(theData$Texp_UP[theData$Texp_UP< theData$Tmax_Texp_UP])
  m_EXP_UP <- length(D_EXP_UP$x)
  D_EXP_DN  <- density(theData$Texp_DN[theData$Texp_DN< theData$Tmax_Texp_DN])
  m_EXP_DN <- length(D_EXP_DN$x)
  
  CURVES_DATA<- as.data.frame(c(D_CTR_UP$x,D_CTR_DN$x,D_EXP_UP$x,D_EXP_DN$x))
  names(CURVES_DATA) <- c('t')  
  CURVES_DATA$UPDN <- as.factor(c(rep('UP',m_CTR_UP),rep('DN',m_CTR_UP),rep('UP',m_EXP_UP),rep('DN',m_EXP_DN)))
  CURVES_DATA$EXP_CTR <- as.factor(c(rep('CTR',m_CTR_UP+ m_CTR_DN),rep('EXP',m_EXP_UP+ m_EXP_DN)))
  CURVES_DATA$TYPE <- 'DATA'
  CURVES_DATA$DENSITY  <- c(D_CTR_UP$y,D_CTR_DN$y,D_EXP_UP$y,D_EXP_DN$y)
  
  #### 
  CURVES_MOMENT_ESTIM<- CURVES_DATA
  CURVES_MOMENT_ESTIM$TYPE <- 'MOMENT'
  dCTR_UP <- dExpCensored(D_CTR_UP$x,lambda=param_estim_moment[1],Tmax = Inf,piTrunc= param_estim_moment[2],log = FALSE)
  dCTR_DN <- dExpCensored(D_CTR_DN$x,lambda=param_estim_moment[3],Tmax = Inf,piTrunc= param_estim_moment[4],log = FALSE)
  dEXP_UP <- dminExpExpplusGaussian(D_EXP_UP$x,mu =param_estim_moment[6],sigma = param_estim_moment[7],lambda =  param_estim_moment[5],piTrunc = param_estim_moment[8],lambda_ND = param_estim_moment[1],piTrunc_ND = param_estim_moment[2], Inf,log = FALSE) 
  dEXP_DN <- dminExpExpplusGaussian(D_EXP_DN$x,mu =param_estim_moment[9],sigma = param_estim_moment[10],lambda =  param_estim_moment[5],piTrunc = param_estim_moment[11],lambda_ND = param_estim_moment[3],piTrunc_ND = param_estim_moment[4], Inf,log = FALSE) 
  CURVES_MOMENT_ESTIM$DENSITY  <- c(dCTR_UP,dCTR_DN,dEXP_UP,dEXP_DN)
  
  
  CURVES_ESTIM <- as.data.frame(rep(CURVES_DATA$t,3)); 
  names(CURVES_ESTIM) <- c('t') 
  
  pl <- 0 
  CURVES <- list()
  CURVES$CTR_UP <- matrix(0,length(extr),m_CTR_UP)
  CURVES$CTR_DN <- matrix(0,length(extr),m_CTR_DN)
  CURVES$EXP_UP <- matrix(0,length(extr),m_EXP_UP)
  CURVES$EXP_DN <- matrix(0,length(extr),m_EXP_DN)

  for (i in extr){
    pl = pl+1
    param <- as.numeric(resMCMC$myPostSample[i,])
    lambda_c <- param[5]
    
    # UPDN = 'UP'
    # lambda_ND <-     ifelse(UPDN =='UP',param[1],param[3])
    # piTrunc_ND <-    ifelse(UPDN =='UP',param[2],param[4])
    # mu_e <-          ifelse(UPDN =='UP',param[6],param[9])
    # sigma_e <-       ifelse(UPDN =='UP',param[7],param[10])
    # piTrunc_Read <-  ifelse(UPDN =='UP',param[8],param[11])
    
    CURVES$CTR_UP[pl,] <- dExpCensored(D_CTR_UP$x,lambda=param[1],Tmax = Inf,piTrunc= param[2],log = FALSE)
    CURVES$CTR_DN[pl,] <- dExpCensored(D_CTR_DN$x,lambda=param[3],Tmax = Inf,piTrunc= param[4],log = FALSE)
    CURVES$EXP_UP[pl,] <- dminExpExpplusGaussian(D_EXP_UP$x,mu =param[6],sigma = param[7],lambda =  param[5],piTrunc = param[8],lambda_ND = param[1],piTrunc_ND = param[2], Inf,log = FALSE) 
    CURVES$EXP_DN[pl,] <- dminExpExpplusGaussian(D_EXP_DN$x,mu =param[9],sigma = param[10],lambda =  param[5],piTrunc = param[11],lambda_ND = param[3],piTrunc_ND = param[4], Inf,log = FALSE) 
  }
  
  CURVES_MEAN <- lapply(CURVES,function(M){apply(M,2,mean)})
  CURVES_Q95<- lapply(CURVES,function(M){apply(M, MARGIN = 2, FUN = quantile, probs = c(0.95))})
  CURVES_Q05<- lapply(CURVES,function(M){apply(M, MARGIN = 2, FUN = quantile, probs = c(0.05))})
  
  CURVES_ESTIM$UPDN <- rep(CURVES_DATA$UPDN,3)
  CURVES_ESTIM$EXP_CTR <- rep(CURVES_DATA$EXP_CTR,3)
  CURVES_ESTIM$TYPE <- rep(c('MEAN','Q95','Q05'),each =nrow(CURVES_DATA))
  CURVES_ESTIM$DENSITY <- c(unlist(CURVES_MEAN),unlist(CURVES_Q95),unlist(CURVES_Q05)) 
  
  CURVES_DF <- rbind(CURVES_DATA,CURVES_ESTIM,CURVES_MOMENT_ESTIM)
  g <- ggplot(CURVES_DF,aes(x=t,y=DENSITY,linetype=TYPE)) + geom_line(aes(color=UPDN)) + facet_wrap( UPDN  ~ EXP_CTR)
  
  
  par(mfrow=c(1,1))
  hist(theData$Tctr_UP[theData$Tctr_UP<theData$Tmax_Tctr_UP],freq = FALSE,nclass = 100,main='CTR_UP',xlab='')
  lines(D_CTR_UP$x,CURVES_MEAN$CTR_UP,col='red',lty=1)
  lines(D_CTR_UP$x,CURVES_Q05$CTR_UP,col='red',lty=2)
  lines(D_CTR_UP$x,CURVES_Q95$CTR_UP,col='red',lty=2)
  lines(D_CTR_UP$x,dCTR_UP,col='magenta',lty=2)
  
  hist(theData$Tctr_DN[theData$Tctr_DN<theData$Tmax_Tctr_DN],freq = FALSE,nclass = 100,main='CTR_DN',xlab='')
  lines(D_CTR_DN$x,CURVES_MEAN$CTR_DN,col='green',lty=1)
  lines(D_CTR_DN$x,CURVES_Q05$CTR_DN,col='green',lty=2)
  lines(D_CTR_DN$x,CURVES_Q95$CTR_DN,col='green',lty=2)
  lines(D_CTR_DN$x,dCTR_DN,col='magenta',lty=2)
  
  
  hist(theData$Texp_UP[theData$Texp_UP<theData$Tmax_Texp_UP],freq = FALSE,nclass = 100,main='EXP_UP',xlab='')
  lines(D_EXP_UP$x,CURVES_MEAN$EXP_UP,col='red',lty=1)
  lines(D_EXP_UP$x,CURVES_Q05$EXP_UP,col='red',lty=2)
  lines(D_EXP_UP$x,CURVES_Q95$EXP_UP,col='red',lty=2)
  lines(D_EXP_UP$x,dEXP_UP,col='magenta',lty=2)
  
  
  hist(theData$Texp_DN[theData$Texp_DN<theData$Tmax_Texp_DN],freq = FALSE,nclass = 100,main='EXP_DN',xlab='')
  lines(D_EXP_DN$x,CURVES_MEAN$EXP_DN,col='green',lty=1)
  lines(D_EXP_DN$x,CURVES_Q05$EXP_DN,col='green',lty=2)
  lines(D_EXP_DN$x,CURVES_Q95$EXP_DN,col='green',lty=2)
  lines(D_EXP_DN$x,dEXP_DN,col='magenta',lty=2)
  
  return(g)
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
