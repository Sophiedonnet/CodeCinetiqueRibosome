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
library(mixtools)
library()
source('Functions/functionsEMG.R')
source('Functions/functionsMomentEstimations.R')
source('Functions/functionsPlotPostInference.R')





where_data <- c('DataKaren20230602/TvieB_Controles_concat')
names_data <- list.files(where_data)

for (i in 1:length(names_data)){
  print(names_data[i])
  data_C <- read.table(paste0(where_data,'/',names_data[i]))
  untruncated <- data_C$V1[data_C$V1< max(data_C$V1) & data_C$V1>0]
  hist(untruncated,freq= FALSE,nclass=50)
  
  res <- lapply(2:4, function(K){ 
    mixtools::normalmixEM(log(untruncated),maxit = 10000,k=K)
    })
  loglik <- lapply(res,function(r){r$loglik})
  res <- res[[which.max(loglik)-1]]
  
  abs = seq(min(untruncated),max(untruncated),len = 1000)
  fit <- apply(sapply(1:which.max(loglik),function(k){res$lambda[k]*plnorm(abs,res$mu[k],res$sigma[k])}),1,sum)
  plot(ecdf(untruncated))
  lines(abs,fit,col='red',add = TRUE)
  }
