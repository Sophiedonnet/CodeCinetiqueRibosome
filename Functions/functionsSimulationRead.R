#========================================================================
simu_readData <- function(Texp, Tctr, Tmax,paramCtr,MC = 10000){
#========================================================================
  

  FN_exp <- ecdf(Texp)
  FN_ctr <- ecdf(Tctr)
  
  absc = seq(1,Tmax-1)
  plot(absc,FN_exp(absc),type='s')
  lines(absc,pCtrData(absc,paramCtr = paramCtr,Tmax = Tmax),col = 'red')
  FRead <- function(x,Tmax,paramCtr){
    y <- 1-(1-FN_exp(x))/(1-pCtrData(x,paramCtr = paramCtr,Tmax = Tmax))
    return(y)
  }
  FRead <- function(x,Tmax,paramCtr){
    y <- 1-(1-FN_exp(x))/(1-FN_ctr(x))
    return(y)
  }
  U <- c(FRead(absc,Tmax,paramCtr), 1)
  P <- c(0,diff(U))
  w <- which(P<0)
  P[w] <- 0
  P = P/sum(P)
  EchanRead <- sample(c(absc,Tmax),MC,prob = P,replace = TRUE)
  return(EchanRead)
} 

