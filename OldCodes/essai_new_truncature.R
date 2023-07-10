Tmax <- 80
lambda <- 1/70
n   <- 100000



Y <- rExpCensored(n,lambda,Tmax)

abs <- seq(0,Tmax,by=1)
dens <- dExpCensored(abs,lambda,Tmax)

hist(Y,freq = FALSE,nclass = 100)
lines(abs,dens,col='red')

####################################
piTrunc = 0.1
Y <- rExpCensored_2(n,lambda,Tmax,piTrunc)


print(c(mean(Y[Y<Tmax]),espExpCensored_2(lambda,Tmax)))

hist(Y,freq = FALSE,nclass = 100)
dens2 <- dExpCensored_2(abs,lambda,Tmax,piTrunc)
lines(abs,dens2,col='red')

plot(ecdf(Y))
lines(abs,pExpCensored_2(abs,lambda,Tmax,piTrunc),col='red')
 

mean(Y)
mean(Y[Y<Tmax])
1/lambda-Tmax*(1-pexp(Tmax,lambda))/pexp(Tmax,lambda)

curve(dbeta(x,sum(Y==Tmax),sum(Y<Tmax)))


################################ROUGES 

Y <- read.table("DataKaren/TRC.txt", quote="\"", comment.char="")$V1*2; 
Y <- read.table("DataKaren/TVC.txt", quote="\"", comment.char="")$V1; 

Tmax <- max(Y)
n <- length(Tmax)

sum(Y == Tmax)
piTrunc <- mean(Y == Tmax)

hist(Y,freq = FALSE,nclass = 50)
abs <- seq(0,Tmax,by=1)
dens2 <- dExpCensored_2(abs,lambda=1/mean(Y[Y<Tmax]),Tmax,piTrunc=mean(Y<Tmax))
lines(abs,dens2,col='red')



L = sort(1/seq(50,500,by=0.01))
espExpCensored <-function(lambda,Tmax){
  1/lambda-Tmax*(1-pexp(Tmax,lambda))/pexp(Tmax,lambda)}

D = espExpCensored(L,Tmax = Tmax)
lambda.hat <- L[which.min(abs(D-mean(Y[Y<Tmax])))]



hist(Y,freq = FALSE,nclass = 50)
abs <- seq(0,Tmax,by=1)
dens2 <- dExpCensored(abs,lambda=1/mean(Y[Y<Tmax]),Tmax,piTrunc)
dens3 <- dExpCensored(abs,lambda=lambda.hat,Tmax,piTrunc)
lines(abs,dens2,col='red')
lines(abs,dens3,col='green')

plot(ecdf(Y))
lines(abs,pExpCensored(abs,lambda=1/mean(Y[Y<Tmax]),Tmax,piTrunc),col='red')
lines(abs,pExpCensored(abs,lambda=lambda.hat,Tmax,piTrunc),col='green')

N = 100000
rExpCensored(N,lambda.hat,Tmax,piTrunc)



#################################"
n = 100000
mu = 10
sigma = 1
lambda = 1/4
lambda_ND <- 1/7
piTrunc_ND <- 0.2
Tmax = 20
Y <- rminExpExpplusGaussian(n,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax)
hist(Y,freq = FALSE,nclass=50)  
abs <- seq(0,20,len=200)
d <- dminExpExpplusGaussian(abs,mu,sigma,lambda,piTrunc,lambda_ND,piTrunc_ND,Tmax)
lines(abs,d,col="red")

mean(Y==Tmax)
piTrunc*piTrunc_ND
