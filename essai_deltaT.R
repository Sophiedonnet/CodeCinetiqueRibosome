

#---------------------------------------------------------
wExp = which((data.i$Texp_UP<data.i$Tmax_Texp_UP) & (data.i$Texp_DN<data.i$Tmax_Texp_DN))
length(wExp)/length(data.i$Texp_UP)



par(mfrow=c(1,1))

plot(data.i$Texp_DN,data.i$Texp_UP)
abline(a=0,b=1)


par(mfrow = c(2,1))
hist(data.i$Texp_DN[wExp]-data.i$Texp_UP[wExp])
abline(v =  mean(data.i$Texp_DN[wExp]-data.i$Texp_UP[wExp]),col='red',lwd = 3)
abline(v =  0,col='blue',lwd = 3)


wCtr = which((data.i$Tctr_UP<data.i$Tmax_Tctr_UP) & (data.i$Tctr_DN<data.i$Tmax_Tctr_DN))
length(wCtr)/length(data.i$Tctr_UP)
       
hist(data.i$Tctr_DN[wCtr]-data.i$Tctr_UP[wCtr])
abline(v =  mean(data.i$Tctr_DN[wExp]-data.i$Tctr_UP[wExp]),col='red',lwd= 3)
abline(v =  0,col='blue',lwd= 3)

       
cor(data.i$Texp_UP,data.i$Texp_DN)

cor(data.i$Tctr_UP,data.i$Tctr_DN)
