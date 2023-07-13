source('Functions/functionsEMG.R')
n = 10000
mu = 10
sigma=5
lambda= 1/3

Tmax = 20
piTrunc = 1- pemg(Tmax,mu,sigma,lambda) 
piTrunc  = 0.2
E <- remgCensored(n,mu, sigma,lambda,Tmax, piTrunc)

Esp <- Tmax-(1-piTrunc)/myF(Tmax)*integrate(myF,-1000,Tmax)$value
mean(E==Tmax)
piTrunc
mean(E)
hist(E)

hist(E,freq = FALSE)
curve(demgCensored(x,mu,sigma,lambda,Tmax,piTrunc),add=TRUE,col='red')

mean(E[E<Tmax])


myF <- function(x){pemg(x,mu,sigma,lambda)}


#myFLim <- function(x){x*pemg(x,mu,sigma,lambda)}
mean(E)
integrate(pemg())
