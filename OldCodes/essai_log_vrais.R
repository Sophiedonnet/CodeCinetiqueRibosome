param = rep(10,4)

lambda_ND_V <- exp(log_param_true[1])
lambda_ND_R <- exp(log_param_true[2])
lambda_c <- exp(param[3])
lambda_e <- exp(param[4])

mu_V <- mydata$k/lambda_e
sigma_V <- sqrt(mu_V/lambda_e)

x <- mydata$T_Exp_V
d3 <- sum(log(dminGammaEMGaussian(x, lambda = lambda_c,mu = mu_V,sigma = sigma_V,c(1,lambda_ND_V))))

print(d3)


curve(dmy)