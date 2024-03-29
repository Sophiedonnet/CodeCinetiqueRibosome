#========================================================================
estim_param_moment <- function(Texp, Tctr, Tmax,paramCtr,plot = FALSE){
#========================================================================
  
 
  paramRead <- rep(NaN,3)
  names(paramRead) <- c('mu','sigma','lambda')
  
  
  
  
  #--------------------------- data Exp: About mu_e, lambde_c
  par(mfrow=c(1,1))
  FN_exp <- ecdf(Texp)
  abs <- c(sort(unique(Texp)),100)
  plot(abs,FN_exp(abs),type='s')
  plot(ecdf(Tctr),add = TRUE,col='red')
  lines(abs,pCtrData(abs,paramCtr = paramCtr,Tmax = Tmax),col='red')
  
  FnRead <- function(x,Tmax,paramCtr){
    y <- 1-(1-FN_exp(x))/(1-pCtrData(x,paramCtr = paramCtr,Tmax = Tmax))
    y[x>=Tmax] = 1
    return(y)
  }
  # vérif 
  
  abs <- 0:Tmax 
  U <- FRead(abs,Tmax,paramCtr)
  plot(abs,U,type='l')
  plot(diff(U)/diff(abs))
  P <- c(0,diff(U))
  w <- which(P<0)
  P[w] <- 0
  
  Z <- rbinom(10000,1,P[Tmax+1])
  Echan <- sample(0:(Tmax-1),10000,prob = P[-(Tmax+1)],replace = TRUE)
  EchanRead <- Z*Tmax + (1-Z)*Echan
  curve(FRead(x,Tmax,paramCtr),0,Tmax)
  plot(ecdf(EchanRead),add  =TRUE)
  
  FnRead2  <- ecdf(EchanRead)
  plot(abs,FnRead2(abs),type='l')
  
  FnExp2 <- function(x,Tmax,paramCtr){
    y <- 1-(1-FnRead2(x))*(1-pCtrData(x,paramCtr = paramCtr,Tmax = Tmax))
    y[x>=Tmax] = 1
    return(y)
  }
  plot(abs,FN_exp(abs),type='s')
  lines(abs,FnExp2(abs,Tmax,paramCtr),col='red')
  estFexp <- function(t,paramCtr,piCtr,EchanRead){
    return(y)
  }
  
  plot(ecdf(Texp_UP))
  curve(estFexp(x,paramCtr,piCtr,EchanRead),col='red',add = TRUE)
  
  TUpperBound <- Tmax
  while(is.na(paramRead[2])){
    paramRead <- estim_param_emg(Echan[Echan<TUpperBound])
    TUpperBound <- TUpperBound-10
  }
  curve(pemg(x,paramRead[1],paramRead[2],paramRead[3]),add  = TRUE,col='red')
  
  
  
  myF <-function(logparamRead,Tmax){
    t <- seq(0,Tmax,by = 0.1)
    paramRead <- exp(logparamRead)
    yEmp <- FRead(t,Tmax,piUp)
    yP <- pemg(t,mu,sigma,lambda) - pemg(0,mu,sigma,lambda)
    C <- pemg(Tmax,mu,sigma,lambda) - pemg(0,mu,sigma,lambda)
    
    LL <- sum((yEmp-yMod)^2)
    return(LL)
  }
  
  
  
  
  plot(ecdf(Echan))
  curve((1-sum(Echan==Tmax))*pemg(x,mu=paramRead[1],sigma = paramRead[2],lambda = paramRead[3]),add = TRUE)
  curve(plnorm(x,meanlog = mean(log(Echan)),sdlog =  sd(log(Echan))),add = TRUE,col='green')
  
  return(list(paramRead = paramRead, echan_exp_corr = Echan))
}  


#========================================================================
estim_param_maxlik <- function(data,UP.or.DN = 'UP'){
  #========================================================================
  
  
  
  
  if(UP.or.DN == 'UP'){Tctr <- data$Tctr_UP}else{Tctr <- data$Tctr_DN}
  if(UP.or.DN == 'UP'){Texp <- data$Texp_UP}else{Texp <- data$Texp_DN}
  Tmax_Tctr <- ifelse(UP.or.DN == 'UP',data$Tmax_Tctr_UP,data$Tmax_Tctr_DN)
  Tmax_Texp <- ifelse(UP.or.DN == 'UP',data$Tmax_Texp_UP,data$Tmax_Texp_DN)
  
  
  
  param <- rep(0,7)
  names(param) <- c('lambda_ND','piTrunc_ND','lambda_c','mu_emg','sigma_emg','piTrunc_Read','lambda_e')
  
  
  
  
  #------- data Contr
  L <- sort(1/seq(1,500,by=0.01))
  D <-  espExpCensored(L,Tmax = Tmax_Tctr)
  lambda_ND <- param[1] <-  L[which.min(abs(D-mean(Tctr[Tctr<Tmax_Tctr])))]
  pi_trunc_ND <-param[2]  <- mean(Tctr == Tmax_Tctr) 
  
  pi_trunc_Read <- param[6] <- mean(Texp == Tmax_Texp)/param[2]
  
  
  
  #--------------------------- Init theta 
  FN_exp <- ecdf(Texp)
  FN_ctr <-  ecdf(Tctr)
  abs <- seq(0,Tmax_Texp-1,len=1000)
  U <- function(x){
    1-(1-FN_exp(x))/(1-FN_ctr(x))
  }
  D <- c(0,diff(U(abs)))
  D <- D*(D>0)/sum(D*(D>0)) 
  Echan <- sample(abs,10000,prob = D,replace=TRUE)
  theta_hat <- rep(NaN,3)
  TUpperBound <- Tmax_Tctr
  while(is.na(theta_hat[2])){
    theta_hat <- estim_param_emg(Echan[Echan<TUpperBound])
    TUpperBound <- TUpperBound-10
  }
  param[5] <- theta_hat[2] # sigma
  param[3] <- theta_hat[3] # lambda_e
  ######################################## 
  
  LL <- function(log_mu,param){
    - sum(dminExpExpplusGaussian(Texp,
                               mu  = exp(log_mu),
                               sigma = param[5], 
                               lambda = param[3], 
                               piTrunc = param[6], 
                               lambda_ND = param[1], 
                               piTrunc_ND = param[2],
                               Tmax = Tmax_Texp,
                               log = TRUE ))
  }
  init_log_mu <- log(theta_hat[c(1)])
  res_optim <- optim(init_log_mu,LL,param=param,method='Brent',lower =init_log_mu-10,upper = init_log_mu +10  )
  # lambda_e 
  param[4] <- exp(res_optim$par) # mu
  
  param[7] <- ifelse(UP.or.DN == 'UP', data$k,data$kprime + data$k)/param[4]
  
  return(list(param_estim = param, echan_exp_corr = Echan))
}  


# DD <- density(data$Tctr_UP,bw=1)
# x <- DD$x
# F1z <- pemgCensored(x,mu=theta_hat_UP[1],sigma=theta_hat_UP[2],lambda=theta_hat_UP[3],Tmax=data$Tmax_Texp_UP,piTrunc=param[7])
# f1z <- demgCensored(x,mu=theta_hat_UP[1],sigma=theta_hat_UP[2],lambda=theta_hat_UP[3],Tmax=data$Tmax_Texp_UP,piTrunc=param[7])
# F2z <- FN_Ctr_UP(x)
# f2z <- DD$y
# y <- f1z*(1-F2z) + f2z*(1-F1z)
# 
# par(mfrow=c(1,2))
# plot(density(Echan_UP[Echan_UP<data$Tmax_Texp_UP]),main=paste0('Exper UP corrected by Contr UP. ',name_data))
# curve(demg(x,theta_hat_UP[1],theta_hat_UP[2],theta_hat_UP[3]),add=TRUE,col=data.i$colorUP)
# 
# hist(data$Texp_UP,freq =FALSE,nclass = 50)
# lines(x,y,col=data.i$colorUP,lwd=2)



#========================================================================
#---------- plot fonctionde repartition de toutes les données
#========================================================================
plot_FN_all_data <- function(data){ 
  
  
  
  Times <- data[which(sapply(data,function(s){length(s)})>1)]
  
  DTT <- lapply(1:length(Times),function(k){
    t <- 1:max(Times[[k]])
    Fn <- ecdf(Times[[k]]);  # fonction de répartition. 
    DT <- as.data.frame(t)
    DT$Fn <- Fn(t)
    DT$data <- names(Times)[k]
    DT$exp <- substr(names(Times)[k],2,4)
    DT$UPDN <- substr(names(Times)[k],6,7)
    DT$color <- ifelse(substr(names(Times[k]),6,7)=='UP',data$colorUP,data$colorDN)
    
    return(DT)
  }
  )
  DTT <- do.call("rbind", DTT)
  p <- ggplot(DTT,aes(x=t,y=Fn,group=data,col=UPDN))+geom_line(aes(linetype=exp))+ggtitle(data$name_data)
  return(p)
}


#========================================================================
#---------- plot Fit des fonction de repartition de toutes les données
#========================================================================
plot_Fit_Ctr <- function(data,param_estim_UP,param_estim_DN,which.curve='pdf'){ 
  
  
  
    DT_UP <- c()
    DT_DN <- c()
    ################# UP
    if(!is.null(data$Tctr_UP)){
      if(which.curve=='pdf'){ 
        t <- 1:max(data$Tmax_Tctr_UP)
        Fn <- ecdf(data$Tctr_UP)
        Y_UP <- Fn(t)
        Fit_UP <- pExpCensored(t,lambda=param_estim_UP[1],piTrunc=param_estim_UP[2],Tmax = data$Tmax_Tctr_UP)
      }
      if(which.curve=='density'){
        w <- which(data$Tctr_UP< data$Tmax_Tctr_UP)
        pr <- length(w)/length(data$Tctr_UP)
        d_UP <- density(data$Tctr_UP[w])
        t <- c(d_UP$x)
        Y_UP <- d_UP$y*pr
        Fit_UP <- dExpCensored(t,lambda=param_estim_UP[1],piTrunc=param_estim_UP[2],Tmax = data$Tmax_Tctr_UP)
      }
      
    
      DT_UP <- as.data.frame(rep(t,2))
      names(DT_UP) = 't'
      DT_UP$Prob <- c(Y_UP,Fit_UP)
      DT_UP$Curve <- paste('UP.', rep(c('Data','Fit Exp'),each=length(t)) )
      DT_UP$data <- 'UP'
      DT_UP$type <- rep(c('Data','Fit Exp'),each=length(t))
    }
  
    ################ DN
    if(!is.null(data$Tctr_DN)){
    
      if(which.curve=='pdf'){ 
        t <- 1:max(data$Tmax_Tctr_DN)
        Fn <- ecdf(data$Tctr_DN)
        Y_DN <- Fn(t)
        Fit_DN <- pExpCensored(t,lambda=param_estim_DN[1],piTrunc=param_estim_DN[2],Tmax = data$Tmax_Tctr_DN)
      }
      if(which.curve=='density'){
        w <- which(data$Tctr_DN < data$Tmax_Tctr_DN)
        pr <- length(w)/length(data$Tctr_DN)
        d_DN <- density(data$Tctr_DN[w])
        t <- c(d_DN$x)
        Y_DN <- d_DN$y*pr
        Fit_DN <- dExpCensored(t,lambda=param_estim_DN[1],piTrunc=param_estim_DN[2],Tmax = data$Tmax_Tctr_DN)
      }
      
      DT_DN <- as.data.frame(rep(t,2))
      names(DT_DN) = 't'
      DT_DN$Prob <- c(Y_DN,Fit_DN)
      DT_DN$Curve <- paste('DN.', rep(c('Data','Fit Exp'),each=length(t)) ) 
      DT_DN$data <- 'DN'
      DT_DN$type <- rep(c('Data','Fit Exp'),each=length(t))
    }
    DT <- rbind(DT_UP,DT_DN)
    
    
    q <- ggplot(DT,aes(x=t,y=Prob,group=Curve,color=data)) + geom_line(aes(linetype=type)) + ggtitle(paste0(data$name_data,'. Exponential Fit on Ctr Data'))
    return(q)
}

#========================================================================
#---------- plot Fit des fonction de repartition de toutes les données
#========================================================================
plot_Fit_Exp <- function(data,name_data,param_estim){ 
  
  
  
  ################# UP
  t <- 1:max(data$Tmax_Tctr_UP)
  Fn <- ecdf(data$Tctr_UP); 
  DT_UP <- as.data.frame(rep(t,2))
  names(DT_UP) = 't'
  DT_UP$Prob <- c(Fn(t),pExpCensored(t,param_estim[1],param_estim[2],Tmax = data$Tmax_Tctr_UP))
  DT_UP$Curve <- paste('UP.', rep(c('Data','Fit'),each=length(t)) )
  DT_UP$data <- 'UP'
  DT_UP$type <- rep(c('Data','Fit'),each=length(t))
  DT <- DT_UP
  ################ DN
  if(!is.null(data$Tctr_DN)){
    
    t <- 1:max(data$Tmax_Tctr_DN)
    Fn <- ecdf(data$Tctr_DN); 
    DT_DN <- as.data.frame(rep(t,2))
    names(DT_DN) = 't'
    DT_DN$Prob <- c(Fn(t),pExpCensored(t,param_estim[3],param_estim[4],Tmax = data$Tmax_Tctr_DN))
    DT_DN$Curve <- paste('DN.', rep(c('Data','Fit'),each=length(t)) ) 
    DT_DN$data <- 'DN'
    DT_DN$type <- rep(c('Data','Fit'),each=length(t))
    DT <- rbind(DT,DT_DN)
  }
  
  
  q <- ggplot(DT,aes(x=t,y=Prob,group=Curve,color=data)) + geom_line(aes(linetype=type)) + ggtitle(paste0(name_data,' .Fit Data Ctr'))
  return(q)
}
