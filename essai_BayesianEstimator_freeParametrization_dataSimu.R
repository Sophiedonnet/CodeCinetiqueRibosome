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
source('Functions/functionsEMG.R')
source('Functions/functionsMomentEstimations.R')
source('Functions/functionsPlotPostInference.R')
source('Functions/myMCMC_marg_constMu.R')
##################################################
################ Data simulation
###################################################

# 
where_data <- c('DataKarenComplete/FormattedData/') # où sont les données envoyées par Karen 
names_data <- list.files(where_data) # liste les fichiers de Rdata
Nbfiles <- length(names_data)
all_directories <- paste0(where_data,names_data)
load('DataKarenComplete/resEstim/res_EstimMoment_allData.Rdata')


d = 1; 
load(paste0(where_data,names_data[d]))
data.d <- data.i
# data.sim <- data.d
#   
# 
# param_sim <- rep(0,12)
# names(param_sim) <- c("lambda_ND_UP", "piTrunc_ND_UP", 
#                       "lambda_ND_DN" , "piTrunc_ND_DN", 
#                       "lambda_c", 
#                       "mu_e_UP" , "sigma_e_UP", "piTrunc_Read_UP", 
#                       "mu_e_DN" , "sigma_e_DN", "piTrunc_Read_DN",
#                       "lambda_e")
# 
# param_sim[c(1,2)] <- param_estim_UP_all[d,c(1,2)] # lambda_ND_UP ,piTrunc_ND_UP
# param_sim[c(3,4)] <- param_estim_DN_all[d,c(1,2)] # # lambda_ND_DN ,piTrunc_ND_DN
# param_sim[5] <- 0.5*(param_estim_UP_all[d,3] + param_estim_DN_all[d,3]) # lambda_c
# param_sim[c(6,7,8)] <- param_estim_UP_all[d,c(4,5,6)]  # mu_e_UP, sigma_e_UP, piTrunc_Read_UP
# param_sim[c(9,10,11)] <- param_estim_DN_all[d,c(4,5,6)]  # mu_e_DN, sigma_e_DN, piTrunc_Read_DN
# param_sim[12] = 0.5*(data.d$k/param_sim[6]  + (data.d$k + data.d$kprime)/param_sim[9])
# param_sim[6] <-data.d$k/param_sim[12]
# param_sim[9] <- (data.d$k + data.d$kprime)/param_sim[12]
# 
# 
# 
# fact.ndata  = 1
# data.sim$Texp_UP <- rOurModelExp_constMu(n = length(data.d$Texp_UP)*fact.ndata,param= param_sim,UPDN='UP',Tmax = data.d$Tmax_Texp_UP)
# data.sim$Texp_DN <- rOurModelExp_constMu(n = length(data.d$Texp_DN)*fact.ndata,param= param_sim,UPDN='DN',Tmax = data.d$Tmax_Texp_DN)
# data.sim$Tctr_UP <- rExpCensored(n = length(data.d$Tctr_UP)*fact.ndata,lambda = param_sim[1],Tmax = data.d$Tmax_Tctr_UP, piTrunc=param_sim[2])
# data.sim$Tctr_DN <- rExpCensored(n = length(data.d$Tctr_DN)*fact.ndata,lambda = param_sim[3],Tmax = data.d$Tmax_Tctr_DN, piTrunc=param_sim[4])
# 
# log_param_sim <- from_param_to_log_param_constMu(param_sim)

theData <- data.d


##################################################
################ ESTIMATION MOMENT
###################################################

resEstimUP <- estim_param_maxlik(theData,'UP')
resEstimDN <- estim_param_maxlik(theData,'DN')



param_estim_moment <- rep(0,12)
param_estim_moment[c(1,2)] <- resEstimUP$param_estim[c(1,2)]
param_estim_moment[c(3,4)] <- resEstimDN$param_estim[c(1,2)]
param_estim_moment[5] <-  0.5*(resEstimUP$param_estim[3] + resEstimDN$param_estim[3])
param_estim_moment[c(6,7,8)] <- resEstimUP$param_estim[c(4,5,6)]
param_estim_moment[c(9,10,11)] <- resEstimDN$param_estim[c(4,5,6)]
param_estim_moment[12] = 0.5*(theData$k/param_estim_moment[6]  + (theData$k + theData$kprime)/param_estim_moment[9])
names(param_estim_moment) <-  c("lambda_ND_UP", "piTrunc_ND_UP", 
                               "lambda_ND_DN" , "piTrunc_ND_DN", 
                                 "lambda_c", 
                              "mu_e_UP" , "sigma_e_UP", "piTrunc_Read_UP",
                                "mu_e_DN" , "sigma_e_DN", "piTrunc_Read_DN",
                                "lambda_e")

rbind(param_estim_moment)

 

##############################################"
###############" MCMC Estimation starting from moment estimator
############################################ 
#--------------- run MCMC


log_param_init <- log_param_estim_moment <- from_param_to_log_param_constMu(param_estim_moment)
param_init <- param_estim_moment
nbParam <- length(log_param_init)
paramsChains <- list(nMCMC=5000,rho=rep(1,nbParam),nBurnin=1,paramsToSample=c(1:12)[-c(2,4,8,11)])
paramsChains$rho[c(12)] <- 0.7
paramsChains$rho[c(5)] <- 0.7
paramsChains$rho[c(1,3)] <- 0.1
paramsChains$rho[c(7,10)] <- 0.1


#-------------- pi trunc non estime
#log_param_init[-paramsChains$paramsToSample] <- from_param_to_log_param_constMu(param_sim)[-paramsChains$paramsToSample]
#param_init[-paramsChains$paramsToSample] <- param_sim[-paramsChains$paramsToSample]
log_param_init[c(2,4,8,11)] <- from_param_to_log_param_constMu(param_estim_moment)[c(2,4,8,11)]
param_init[c(2,4,8,11)] <- param_estim_moment[c(2,4,8,11)]

#------------------ Hyperparams de la prior
hyperparams <- list(param1=rep(0,nbParam) ,param2 = rep(3,nbParam))
hyperparams$param1[c(1,3)] = -3;
hyperparams$param2[c(1,3)] =2;
hyperparams$param1[5] = -3;
hyperparams$param2[5] =0.5;
hyperparams$param1[12] = 2; 
hyperparams$param2[12] = 0.5; 
hyperparams$param1[7] = 2; 
hyperparams$param2[7] = 0.5; 
hyperparams$param1[10] = 2; 
hyperparams$param2[10] = 0.5; 

#------------------Run MCMC
resMCMC <- my_mcmc_marg_constMu(theData,log_param_init,hyperparams = hyperparams,paramsChains = paramsChains)
 


##############################################"
#------------------------- PLOT RESULTAS
##############################################"
paramsChains$nBurnin = 1001
paramsChains$thinning = 1

plotChains(resMCMC, paramsChains, log_param_estim_moment, log_param_sim = NULL)

plotPosteriorDistr(resMCMC, paramsChains, hyperparams, log_param_estim_moment, log_param_sim = NULL)

g<- plotCurves(resMCMC, theData,paramsChains, hyperparams, param_estim_moment)
  


